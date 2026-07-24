"""Tests for the `display` section of the merged spec document.

The section moves twelve hand-written frontend `case` bodies into data. Three
things have to hold:

  * it is *consistent* -- every `from`/`compose` field reference resolves to a
    real parse plugin and a field that plugin actually produces, checked at load
    time like the config<->parsing half (this is what stops the labels drifting
    from the parsers);
  * scopes are *derived*, never authored -- the display rows name a plugin only,
    and the allele-vs-transcript answer comes from `parsing`;
  * the load side is *defensive* -- a spec pinned before this section existed has
    no `display` key, must still load, and its job must still render (falling
    back to the current genome's display spec).
"""

import json

import pytest
from pydantic import FilePath, ValidationError

from app.vep.models.display_spec_model import (
    DisplayGroupBlock,
    DisplayListBlock,
    DisplayRowsBlock,
)
from app.vep.models.merged_spec_model import MergedSpec
from app.vep.utils.spec_loader import (
    SPEC_SIDECAR_FILE,
    load_merged_spec,
    write_spec_sidecar,
)
from app.vep.utils.vcf_results import (
    _load_pinned_merged_spec,
    _resolve_display_payload,
)

SPEC = load_merged_spec("human_grch38")

# The options moved off the frontend switch in this change.
SPEC_DRIVEN_OPTIONS = {
    "hgvs", "hgvsg", "spdi", "alphamissense", "revel", "clinpred", "cadd", "spliceai",
    "loeuf", "dosage_sensitivity", "utrannotator", "nmd", "riboseqorfs", "eve",
    # `list`-block options (repeat + truncate, migrated off frontend overrides)
    "phenotypes", "go", "mavedb", "nearest_gene", "nearest_exon_jb",
    # sub-option rows (Show-all enumeration)
    "mutfunc",
    # multi-cell list items under an option heading (GWAS + QTL groups)
    "opentargets",
    # conditional (`when`) + group + list-as-rows breakdown
    "clinvar",
    # view gating (default vs Show all) + row/item link builder + count
    "protvar",
    # ENSP parse plugin + a "Protein ID" row with an app_popup builder link
    "protein",
    # default-vs-Show-all views + count + sub-option count rows (no new operator)
    "intact",
    # a fields-less gff-overlap custom (new Regulatory panel)
    "gencode_promoters",
}


def _vcf_path(tmp_path) -> FilePath:
    path = tmp_path / "results.vcf"
    path.write_text("##fileformat=VCFv4.2\n")
    return FilePath(path)


def _doc(display, plugins=None):
    """A minimal merged document carrying `display`, with one parse plugin whose
    fields the display rows can reference."""
    return {
        "genome": {"assembly": "GRCh38"},
        "config": {"entries": []},
        "parsing": {
            "plugins": plugins
            if plugins is not None
            else [
                {
                    "plugin": "revel",
                    "scope": "transcript",
                    "output": "revel",
                    "csq_fields": ["REVEL"],
                    "targets": [
                        # a float, as REVEL scores are — so a `num` format over it
                        # is type-compatible (see the format<->type checks below)
                        {"field": "score", "from": "REVEL", "transform": "scalar", "type": "float"}
                    ],
                }
            ]
        },
        "display": display,
    }


def _display(*rows, **block):
    block.setdefault("kind", "rows")
    block["rows"] = list(rows)
    return {"options": [{"option_id": "revel", "blocks": [block]}]}


# --- the shipped document ---------------------------------------------------


def test_bundled_spec_has_a_display_section_for_the_moved_options():
    assert SPEC.display is not None
    assert {o.option_id for o in SPEC.display.options} == SPEC_DRIVEN_OPTIONS


def test_bundled_display_references_resolve():
    """Belt and braces: load_merged_spec already runs the check, but state the
    invariant — every display ref resolves: a fixed row's `plugin.field`, a
    block's `when` field, a list block's `plugin.listField`, and each list
    element's item fields (label + cells). Groups are flattened by iter_blocks."""
    targets = {
        plugin.plugin: {t.field: t for t in plugin.targets}
        for plugin in SPEC.parse_plugins()
    }

    def resolves(ref: str) -> bool:
        plugin, _, field = ref.partition(".")
        return plugin in targets and field in targets[plugin]

    for option in SPEC.display.options:
        for block in option.iter_blocks():
            if block.when:
                assert resolves(block.when.field_ref), (
                    option.option_id, block.when
                )
            if isinstance(block, DisplayGroupBlock):
                continue
            if isinstance(block, DisplayListBlock):
                plugin, list_field = block.list_ref()
                assert list_field in targets[plugin], (
                    option.option_id, plugin, list_field
                )
                item_fields = set(targets[plugin][list_field].item_fields or [])
                for item_field in block.item.item_field_refs():
                    assert item_field in item_fields, (
                        option.option_id, plugin, list_field, item_field
                    )
            elif isinstance(block, DisplayRowsBlock):
                for row in block.rows:
                    for ref in row.field_refs():
                        assert resolves(ref), (option.option_id, ref)


# --- the static reference check ---------------------------------------------


def test_unknown_plugin_reference_raises():
    doc = _doc(_display({"label": "REVEL", "from": "not_a_plugin.score"}))
    with pytest.raises(ValidationError, match="unknown parse plugin 'not_a_plugin'"):
        MergedSpec.model_validate(doc)


def test_unknown_field_reference_raises():
    """The deliberate-typo case: right plugin, field it does not produce."""
    doc = _doc(_display({"label": "REVEL", "from": "revel.scores"}))
    with pytest.raises(
        ValidationError,
        match="references field 'scores' that parse plugin 'revel' does not produce",
    ):
        MergedSpec.model_validate(doc)


def test_compose_field_references_are_checked_too():
    doc = _doc(
        _display(
            {
                "label": "REVEL",
                "compose": {
                    "format": "with_score",
                    "classification": "revel.score",
                    "score": "revel.nope",
                },
            }
        )
    )
    with pytest.raises(ValidationError, match="'nope'"):
        MergedSpec.model_validate(doc)


def test_unknown_requires_plugin_raises():
    doc = _doc(_display({"label": "REVEL", "from": "revel.score"}, requires="ghost"))
    with pytest.raises(ValidationError, match="requires unknown parse plugin 'ghost'"):
        MergedSpec.model_validate(doc)


def test_when_field_reference_is_checked():
    """A block's `when` reads a field, so a bad ref there fails like a row's."""
    doc = _doc(
        _display(
            {"label": "REVEL", "from": "revel.score"},
            when={"present": "revel.nope"},
        )
    )
    with pytest.raises(ValidationError, match="'nope'"):
        MergedSpec.model_validate(doc)


def test_group_subblock_references_are_checked():
    """The check flattens groups, so a bad ref inside a group's block fails."""
    doc = _doc(
        {
            "options": [
                {
                    "option_id": "revel",
                    "blocks": [
                        {
                            "kind": "group",
                            "heading": "Group",
                            "blocks": [
                                {
                                    "kind": "rows",
                                    "rows": [{"label": "R", "from": "revel.nope"}],
                                }
                            ],
                        }
                    ],
                }
            ]
        }
    )
    with pytest.raises(ValidationError, match="'nope'"):
        MergedSpec.model_validate(doc)


def test_when_needs_exactly_one_of_present_or_empty():
    doc = _doc(
        _display(
            {"label": "REVEL", "from": "revel.score"},
            when={"present": "revel.score", "empty": "revel.score"},
        )
    )
    with pytest.raises(
        ValidationError, match="exactly one of `present` or `empty`"
    ):
        MergedSpec.model_validate(doc)


def test_item_label_needs_exactly_one_of_from_or_template():
    doc = _doc(
        {
            "options": [
                {
                    "option_id": "revel",
                    "blocks": [
                        {
                            "kind": "list",
                            "from": "revel.score",
                            "item": {
                                "label": {"from": "x", "template": "y"},
                                "cells": [{"from": "x"}],
                            },
                        }
                    ],
                }
            ]
        }
    )
    with pytest.raises(
        ValidationError, match="exactly one of `from` or `template`"
    ):
        MergedSpec.model_validate(doc)


def test_row_needs_exactly_one_of_from_or_compose():
    with pytest.raises(ValidationError, match="exactly one of `from` or `compose`"):
        MergedSpec.model_validate(_doc(_display({"label": "REVEL"})))


def test_unknown_row_key_is_rejected():
    """extra="forbid", like the rest of the spec document: a key we do not
    understand is a spec we do not understand."""
    with pytest.raises(ValidationError):
        MergedSpec.model_validate(
            _doc(_display({"label": "REVEL", "from": "revel.score", "italic": True}))
        )


def test_valid_display_loads():
    spec = MergedSpec.model_validate(
        _doc(_display({"label": "REVEL", "from": "revel.score", "format": "num"}))
    )
    assert spec.display.options[0].blocks[0].rows[0].source == "revel.score"


# --- the static format <-> type compatibility check -------------------------
#
# A `format` assumes a value shape; applying it to the wrong parsing type crashes
# the renderer (`num` -> `.toPrecision`, `join`/`humanize_join` -> `.join`/`.map`).
# One plugin with a field of every shape the checks distinguish.

_TYPED_PLUGIN = [
    {
        "plugin": "p",
        "scope": "transcript",
        "output": "p",
        "csq_fields": ["C"],
        "targets": [
            {"field": "score", "from": "C", "transform": "scalar", "type": "float"},
            {"field": "name", "from": "C", "transform": "scalar", "type": "string"},
            {"field": "terms", "from": "C", "transform": "list", "type": "string"},
            {
                "field": "assays",
                "from": "C",
                "transform": "chunk",
                "size": 2,
                "as": [
                    {"field": "urn", "type": "string"},
                    {"field": "sc", "type": "float"},
                ],
                "item_fields": ["urn", "sc"],
            },
        ],
    }
]


def _typed(*rows):
    return _doc(_display(*rows), plugins=_TYPED_PLUGIN)


def _typed_list(item):
    return _doc(
        {"options": [{"option_id": "p", "blocks": [
            {"kind": "list", "from": "p.assays", "item": item}
        ]}]},
        plugins=_TYPED_PLUGIN,
    )


def test_num_format_over_a_string_field_raises():
    """The motivating case: `num` calls `.toPrecision` and throws on a string."""
    with pytest.raises(
        ValidationError,
        match=r"formats 'p.name' as 'num'.*needs a numeric field",
    ):
        MergedSpec.model_validate(
            _typed({"label": "N", "from": "p.name", "format": "num"})
        )


def test_num_format_over_a_list_field_raises():
    with pytest.raises(ValidationError, match=r"formats 'p.terms' as 'num'"):
        MergedSpec.model_validate(
            _typed({"label": "N", "from": "p.terms", "format": "num"})
        )


def test_list_format_over_a_scalar_raises():
    """`humanize_join` maps `.replace` over the elements, so it needs a list."""
    with pytest.raises(
        ValidationError,
        match=r"formats 'p.name' as 'humanize_join'.*needs a list of strings",
    ):
        MergedSpec.model_validate(
            _typed({"label": "N", "from": "p.name", "format": "humanize_join"})
        )


def test_count_over_a_numeric_field_raises():
    """`count` is for a list or a delimited string; a number always drops."""
    with pytest.raises(ValidationError, match=r"formats 'p.score' as 'count'"):
        MergedSpec.model_validate(
            _typed({"label": "N", "from": "p.score", "format": "count"})
        )


def test_list_item_cell_format_is_type_checked():
    """A cell's `format` is checked against the element field's declared `as`
    type — `urn` is a string, so `num` over it would crash."""
    with pytest.raises(ValidationError, match=r"formats 'p.assays.urn' as 'num'"):
        MergedSpec.model_validate(
            _typed_list({"cells": [{"from": "urn", "format": "num"}]})
        )


def test_with_score_requires_a_numeric_score():
    """The `with_score` compose renders `num(score)`, so a string score is a
    crash the same way a bad row format is."""
    with pytest.raises(
        ValidationError,
        match=r"'p.name' in a with_score value, which needs a numeric field",
    ):
        MergedSpec.model_validate(
            _typed(
                {
                    "label": "N",
                    "compose": {
                        "format": "with_score",
                        "classification": "p.name",
                        "score": "p.name",
                    },
                }
            )
        )


def test_compatible_formats_over_matching_types_load():
    """Every pairing the real spec uses: `num` over a float, a list format over a
    list, `count` over a delimited string, and a cell `num` over a float item
    field — all load."""
    doc = _doc(
        {"options": [{"option_id": "p", "blocks": [
            {"kind": "rows", "rows": [
                {"label": "S", "from": "p.score", "format": "num"},
                {"label": "T", "from": "p.terms", "format": "humanize_join"},
                {"label": "C", "from": "p.name", "format": "count"},
            ]},
            {"kind": "list", "from": "p.assays", "item": {"cells": [
                {"from": "urn", "format": "text"},
                {"from": "sc", "format": "num"},
            ]}},
        ]}]},
        plugins=_TYPED_PLUGIN,
    )
    spec = MergedSpec.model_validate(doc)
    assert spec.display is not None


# --- scopes are derived, not authored ---------------------------------------


def test_plugin_scopes_come_from_the_parsing_plugins():
    scopes = SPEC.plugin_scopes()
    # The four allele-scoped plugins the moved options read.
    assert scopes["cadd"] == "allele"
    assert scopes["spdi"] == "allele"
    assert scopes["hgvsg"] == "allele"
    assert scopes["hgvs"] == "transcript"
    assert set(scopes) == {p.plugin for p in SPEC.parse_plugins()}


def test_display_payload_carries_options_and_scopes():
    payload = SPEC.display_payload()
    assert {o.option_id for o in payload.options} == SPEC_DRIVEN_OPTIONS
    assert payload.plugin_scopes == SPEC.plugin_scopes()


def test_display_section_serialises_with_the_authored_key_names():
    """The wire format uses `from`, as authored -- not the Python field name."""
    dumped = SPEC.display_payload().model_dump(mode="json", by_alias=True)
    hgvs = next(o for o in dumped["options"] if o["option_id"] == "hgvs")
    assert hgvs["blocks"][0]["rows"][0]["from"] == "hgvs.transcript"


# --- the legacy fallback ----------------------------------------------------


def _legacy_document() -> dict:
    """The bundled spec as it was written before this change: same document with
    the `display` key genuinely absent."""
    payload = json.loads(
        (SPEC.model_dump_json(exclude={"spec_version": True}))
    )
    payload.pop("display")
    assert "display" not in payload
    return payload


def test_a_spec_without_a_display_key_still_loads():
    spec = MergedSpec.model_validate(_legacy_document())
    assert spec.display is None
    assert spec.display_payload() is None


def test_legacy_pinned_sidecar_loads_and_falls_back_to_the_current_display(tmp_path):
    (tmp_path / SPEC_SIDECAR_FILE).write_text(json.dumps(_legacy_document()))
    pinned = _load_pinned_merged_spec(_vcf_path(tmp_path))
    assert pinned is not None and pinned.display is None

    payload = _resolve_display_payload(pinned)
    assert payload is not None
    assert {o.option_id for o in payload.options} == SPEC_DRIVEN_OPTIONS
    # Scopes still describe the *pinned* parsers -- only the layout is current.
    assert payload.plugin_scopes == pinned.plugin_scopes()


def test_a_current_pinned_sidecar_uses_its_own_display(tmp_path):
    write_spec_sidecar(tmp_path, SPEC)
    pinned = _load_pinned_merged_spec(_vcf_path(tmp_path))
    payload = _resolve_display_payload(pinned)
    assert payload is not None
    assert {o.option_id for o in payload.options} == SPEC_DRIVEN_OPTIONS


def test_no_pinned_spec_means_no_display_payload():
    assert _resolve_display_payload(None) is None


def test_unknown_assembly_falls_back_to_nothing_rather_than_raising():
    legacy = MergedSpec.model_validate(
        {**_legacy_document(), "genome": {"assembly": "Nothing_v1"}}
    )
    assert _resolve_display_payload(legacy) is None
