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
    CellSpec,
    DEFAULT_TRUNCATE_VISIBLE_COUNT,
    DisplayGroupBlock,
    DisplayListBlock,
    DisplayRowsBlock,
    DisplayTableBlock,
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
    "hgvs", "spdi", "alphamissense", "revel", "clinpred", "cadd", "avi", "spliceai",
    "loeuf", "pli", "dosage_sensitivity", "utrannotator", "nmd", "riboseqorfs",
        "tss_distance", "eve", "gerp",
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
    # a count linked through a sibling field (`link_from` on a row)
    "geno2mp",
    # `map_rows`: rows discovered from the job's own population vocabulary
    # rather than written into the spec. These were the last options drawn by
    # hand-written frontend components instead of the display spec.
    "gnomad_exomes", "gnomad_genomes", "allofus", "gnomad_sv", "gnomad_cnv",
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


def test_af_map_rows_blocks_declare_the_plugin_they_read():
    """`plugin_refs` decides which options a genome's assembled spec offers — an
    option is included only when every plugin it reads is present. A `map_rows`
    block that contributed no refs would put gnomAD SV (GRCh38-only) into a spec
    for an assembly that has no such plugin, and the option would render empty
    rather than being absent."""
    by_id = {o.option_id: o for o in SPEC.display.options}

    for option_id, plugin in (
        ("gnomad_exomes", "gnomad_exomes"),
        ("gnomad_genomes", "gnomad_genomes"),
        ("allofus", "all_of_us"),
        ("gnomad_sv", "gnomad_sv"),
        ("gnomad_cnv", "gnomad_cnv"),
    ):
        assert plugin in by_id[option_id].plugin_refs(), option_id


def test_af_map_rows_reads_the_overall_beside_the_population_dict():
    """The two halves of a source's frequencies live in different targets — the
    all-ancestry figure is a scalar, the rest a dict — so the block names both.
    Without `overall_from` the vocabulary's "All" row has nowhere to read."""
    option = next(
        o for o in SPEC.display.options if o.option_id == "gnomad_genomes"
    )
    block = next(b for b in option.iter_blocks() if b.kind == "map_rows")

    assert block.source == "gnomad_genomes.populations"
    assert block.overall_from == "gnomad_genomes.overall"
    assert block.vocabulary == "af_populations"
    assert block.scope == "gnomad_genomes"


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


def _clinvar_label_doc(label):
    return _doc(
        {
            "options": [
                {
                    "option_id": "revel",
                    "blocks": [
                        {
                            "kind": "list",
                            "from": "revel.score",
                            "item": {"label": label, "cells": [{"from": "count"}]},
                        }
                    ],
                }
            ]
        }
    )


def test_item_label_wrap_needs_a_from():
    """`wrap` surrounds a formatted `from` value — there is nothing to wrap
    around a free-text `template`."""
    doc = _clinvar_label_doc({"template": "y", "wrap": "Submitters reporting \"{}\""})
    with pytest.raises(ValidationError, match="`wrap` needs `from`"):
        MergedSpec.model_validate(doc)


def test_item_label_wrap_needs_a_placeholder():
    """`wrap` without a `{}` slot would drop the value entirely — reject it."""
    doc = _clinvar_label_doc({"from": "significance", "wrap": "Submitters reporting"})
    with pytest.raises(ValidationError, match="`wrap` needs a `\\{\\}` placeholder"):
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


# --- the `table` block ------------------------------------------------------


def _typed_table(*columns):
    return _doc(
        {"options": [{"option_id": "p", "blocks": [
            {"kind": "table", "from": "p.assays", "columns": list(columns)}
        ]}]},
        plugins=_TYPED_PLUGIN,
    )


def test_valid_table_loads():
    spec = MergedSpec.model_validate(
        _typed_table(
            {"label": "URN", "from": "urn"},
            {"label": "Score", "from": "sc", "format": "num"},
        )
    )
    block = spec.display.options[0].blocks[0]
    assert block.kind == "table"
    assert [c.label for c in block.columns] == ["URN", "Score"]


def test_a_column_may_state_its_alignment():
    """The house rule derives alignment from the format — a numeric column reads
    right without saying so. `align` exists for the case a format cannot
    express: a number the source publishes pre-formatted as a string, like
    OpenTargets' p-value, where `format: num` would be a lie the load-time type
    check rightly rejects."""
    spec = MergedSpec.model_validate(
        _typed_table(
            {"label": "URN", "from": "urn"},
            {"label": "p-value", "from": "urn", "align": "right"},
        )
    )
    columns = spec.display.options[0].blocks[0].columns
    assert [c.align for c in columns] == [None, "right"]


def test_alignment_is_restricted_to_left_or_right():
    """Centre-aligned figures do not line up on their digits, which is the whole
    point of the rule."""
    with pytest.raises(ValidationError):
        MergedSpec.model_validate(
            _typed_table({"label": "x", "from": "urn", "align": "centre"})
        )


def test_table_column_ref_must_be_an_item_field():
    """A column reads a field of the list element, checked like a list cell."""
    with pytest.raises(
        ValidationError, match=r"item field 'nope' not in p.assays"
    ):
        MergedSpec.model_validate(_typed_table({"label": "X", "from": "nope"}))


def test_table_column_format_is_type_checked():
    """A column's `format` is checked against the element field's `as` type —
    `urn` is a string, so `num` over it would crash."""
    with pytest.raises(ValidationError, match=r"formats 'p.assays.urn' as 'num'"):
        MergedSpec.model_validate(
            _typed_table({"label": "URN", "from": "urn", "format": "num"})
        )


def _grouped_table(group_by):
    return _doc(
        {"options": [{"option_id": "p", "blocks": [
            {
                "kind": "table",
                "from": "p.assays",
                "group_by": group_by,
                "columns": [{"label": "Score", "from": "sc", "format": "num"}],
            }
        ]}]},
        plugins=_TYPED_PLUGIN,
    )


def test_group_by_labels_rename_headings():
    """`labels` renames the headings the data supplies, for the values it names —
    phenotypes group on `type` but read as "Variant associated", not
    "Variation"."""
    spec = MergedSpec.model_validate(
        _grouped_table({"field": "urn", "labels": {"a": "Ay"}})
    )
    assert spec.display.options[0].blocks[0].group_by.labels == {"a": "Ay"}


def test_group_by_labels_are_optional():
    spec = MergedSpec.model_validate(_grouped_table({"field": "urn"}))
    assert spec.display.options[0].blocks[0].group_by.labels is None


def test_group_by_field_must_be_an_item_field():
    """The field the rows group on is checked like a column's."""
    with pytest.raises(
        ValidationError, match=r"item field 'nope' not in p.assays"
    ):
        MergedSpec.model_validate(_grouped_table({"field": "nope"}))


# --- the `table` block, fixed / matrix mode (SpliceAI) ----------------------


def _typed_matrix(*rows, columns=None):
    columns = columns or [{"label": "Metric"}, {"label": "Value", "format": "num"}]
    return _doc(
        {"options": [{"option_id": "p", "blocks": [
            {"kind": "table", "columns": columns, "rows": list(rows)}
        ]}]},
        plugins=_TYPED_PLUGIN,
    )


def test_valid_fixed_table_loads():
    spec = MergedSpec.model_validate(
        _typed_matrix({"label": "Score", "values": ["p.score"]})
    )
    block = spec.display.options[0].blocks[0]
    assert block.kind == "table"
    assert block.rows[0].label == "Score"
    assert block.rows[0].values == ["p.score"]


def test_table_needs_from_or_rows():
    """Exactly one of `from` (list) or `rows` (fixed) — neither is invalid."""
    with pytest.raises(
        ValidationError, match=r"table needs exactly one of `from` or `rows`"
    ):
        MergedSpec.model_validate(
            _doc(
                {"options": [{"option_id": "p", "blocks": [
                    {"kind": "table", "columns": [{"label": "X"}]}
                ]}]},
                plugins=_TYPED_PLUGIN,
            )
        )


def test_fixed_table_row_value_count_must_match_columns():
    """Two columns means one value column, so a row with two values is wrong."""
    with pytest.raises(ValidationError, match=r"but there are 1 value column"):
        MergedSpec.model_validate(
            _typed_matrix({"label": "Score", "values": ["p.score", "p.name"]})
        )


def test_fixed_table_columns_take_no_from():
    with pytest.raises(
        ValidationError, match=r"a fixed table's columns take no `from`"
    ):
        MergedSpec.model_validate(
            _typed_matrix(
                {"label": "Score", "values": ["p.score"]},
                columns=[{"label": "Metric"}, {"label": "Value", "from": "x"}],
            )
        )


def test_fixed_table_value_ref_is_checked():
    with pytest.raises(ValidationError, match=r"'nope'"):
        MergedSpec.model_validate(
            _typed_matrix({"label": "X", "values": ["p.nope"]})
        )


def test_fixed_table_value_format_is_type_checked():
    """The value column's `format` (num) is checked against each row value's
    field type — `p.name` is a string, which would crash."""
    with pytest.raises(ValidationError, match=r"formats 'p.name' as 'num'"):
        MergedSpec.model_validate(
            _typed_matrix({"label": "X", "values": ["p.name"]})
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
    # hgvs is now a "HGVS" group with per-param (requires_selected) sub-blocks.
    assert hgvs["blocks"][0]["blocks"][0]["rows"][0]["from"] == "hgvs.transcript"


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
    # And the layout is *whole*: the rating scales its options refer to travel
    # with them, or these jobs would render ClinVar's review status unrated.
    assert "clinvar_submission" in payload.rating_scales


def test_a_current_pinned_sidecar_uses_its_own_display(tmp_path):
    write_spec_sidecar(tmp_path, SPEC)
    pinned = _load_pinned_merged_spec(_vcf_path(tmp_path))
    payload = _resolve_display_payload(pinned)
    assert payload is not None
    assert {o.option_id for o in payload.options} == SPEC_DRIVEN_OPTIONS


def test_no_pinned_spec_means_no_display_payload():
    assert _resolve_display_payload(None) is None


def test_unknown_assembly_falls_back_to_the_base_display():
    """A job pinned before the display section existed, on a genome with no spec
    of its own, now renders the base options rather than nothing — the same
    fallback the submit path uses."""
    legacy = MergedSpec.model_validate(
        {**_legacy_document(), "genome": {"assembly": "Nothing_v1"}}
    )
    payload = _resolve_display_payload(legacy)
    assert payload is not None
    base = load_merged_spec("base")
    assert {o.option_id for o in payload.options} == {
        o.option_id for o in base.display.options
    }


# --- linked table columns (IntAct) ------------------------------------------


def _intact_table():
    """The IntAct interactions table from the shipped spec."""
    spec = load_merged_spec("human_grch38")
    option = next(o for o in spec.display.options if o.option_id == "intact")
    return next(b for b in option.blocks if b.kind == "table")


def test_intact_table_columns_are_in_the_agreed_order():
    assert [c.label for c in _intact_table().columns] == [
        "Interaction AC",
        "Feature Type",
        "Interaction Participants",
        "Feature short label",
        "Affected Protein",
        "PubMed Links",
    ]


def test_only_the_always_emitted_intact_columns_are_ungated():
    """interaction_ac and feature_type come back whatever is selected; the other
    four are sub-option columns, so the table is 2 to 6 columns wide."""
    columns = {c.label: c for c in _intact_table().columns}
    ungated = [label for label, c in columns.items() if c.sub_option is None]
    assert sorted(ungated) == ["Feature Type", "Interaction AC"]
    assert columns["Interaction Participants"].sub_option.id == (
        "intact_interaction_participants"
    )
    assert columns["Affected Protein"].sub_option.id == "intact_ap_ac"


def test_intact_uniprot_columns_split_and_require_the_uniprotkb_prefix():
    """Participants pack several accessions into one value, and a value without
    the `uniprotkb:` prefix is not an accession at all — it must not be linked."""
    columns = {c.label: c for c in _intact_table().columns}

    participants = columns["Interaction Participants"]
    assert participants.split == "_and_"
    assert participants.link_prefix == "uniprotkb:"
    assert participants.link.template == "https://www.uniprot.org/uniprotkb/{value}/entry"

    affected = columns["Affected Protein"]
    assert affected.link_prefix == "uniprotkb:"
    assert affected.split is None  # a single accession


def test_intact_identifier_columns_link_to_their_own_resources():
    columns = {c.label: c for c in _intact_table().columns}
    assert columns["Interaction AC"].link.template == (
        "https://www.ebi.ac.uk/intact/details/interaction/{value}"
    )
    assert columns["PubMed Links"].link.template == (
        "https://europepmc.org/article/MED/{value}"
    )
    # These are never prefixed, so they are linked unconditionally.
    assert columns["Interaction AC"].link_prefix is None
    assert columns["PubMed Links"].link_prefix is None


def test_split_or_prefix_without_a_link_is_rejected_at_load():
    # The rule belongs to a rendered value, so it holds wherever one appears --
    # a column, a cell of a repeated item, a line of a list-valued cell.
    with pytest.raises(ValidationError, match="only apply to a linked value"):
        DisplayTableBlock.model_validate({
            "kind": "table",
            "from": "intact.interactions",
            "columns": [{"label": "Participants", "from": "x", "split": "_and_"}],
        })
    with pytest.raises(ValidationError, match="only apply to a linked value"):
        CellSpec.model_validate({"from": "x", "link_prefix": "uniprotkb:"})


# --- house style: repeating blocks are capped ------------------------------


def _repeating_blocks(spec):
    """Every block that renders one row per data element, with its option id."""
    found = []

    def walk(blocks, option_id):
        for block in blocks:
            if block.kind == "group":
                walk(block.blocks, option_id)
            elif block.kind == "list" or (
                block.kind == "table" and block.source is not None
            ):
                found.append((option_id, block))

    for option in spec.display.options:
        walk(option.blocks, option.option_id)
    return found


def test_every_repeating_block_is_capped():
    """House style: a block that repeats shows three rows and hides the rest
    behind a show-more toggle. Applied as a model default, so a new block gets
    it without having to remember — this asserts none has slipped through."""
    spec = load_merged_spec("human_grch38")
    uncapped = [
        f"{option_id}: {block.heading or block.kind}"
        for option_id, block in _repeating_blocks(spec)
        if block.truncate is None
    ]
    assert not uncapped, "repeating blocks with no cap:\n  " + "\n  ".join(uncapped)


def test_the_default_cap_is_three():
    spec = load_merged_spec("human_grch38")
    counts = {
        block.truncate.visible_count for _, block in _repeating_blocks(spec)
    }
    assert counts == {DEFAULT_TRUNCATE_VISIBLE_COUNT} == {3}


def test_a_fixed_matrix_table_is_not_capped():
    """A fixed table states its own rows, so its height is known and small —
    SpliceAI's four splicing events should not gain a show-more toggle."""
    spec = load_merged_spec("human_grch38")
    spliceai = next(o for o in spec.display.options if o.option_id == "spliceai")
    matrices = [
        b
        for group in spliceai.blocks
        for b in (group.blocks if group.kind == "group" else [group])
        if b.kind == "table" and b.rows is not None
    ]
    assert matrices, "expected SpliceAI to carry a fixed matrix table"
    assert all(b.truncate is None for b in matrices)


def test_every_block_kind_carries_the_shared_gates():
    """The point of the shared base: a block kind cannot arrive missing a gate.

    Each of the four used to redeclare `heading`/`requires_selected`/`when`/
    `view` for itself, and `group_by` had already been added to two of them
    independently. This fails if a fifth kind is written without inheriting."""
    from app.vep.models.display_spec_model import (
        DisplayGroupBlock,
        DisplayListBlock,
        DisplayRowsBlock,
        DisplayTableBlock,
        _GatedBlock,
    )

    gates = {"heading", "requires_selected", "when", "view"}
    for block in (
        DisplayRowsBlock,
        DisplayListBlock,
        DisplayTableBlock,
        DisplayGroupBlock,
    ):
        assert issubclass(block, _GatedBlock), block.__name__
        assert gates <= set(block.model_fields), block.__name__


def test_the_shared_gates_still_reject_an_unknown_field():
    """Inheriting must not loosen `extra=forbid` — a typo'd gate has to fail at
    load rather than being silently ignored."""
    import pytest
    from pydantic import ValidationError

    from app.vep.models.display_spec_model import DisplayRowsBlock

    with pytest.raises(ValidationError):
        DisplayRowsBlock.model_validate(
            {"kind": "rows", "rows": [], "whn": {"present": "x.y"}}
        )


# --- star ratings -----------------------------------------------------------


_SCALE = {
    "confidence": {"out_of": 4, "ratings": {"reviewed by expert panel": 3}}
}


def _rated(display):
    display.setdefault("rating_scales", _SCALE)
    return _doc(display)


def _starred_column(scale):
    """A table column whose items carry a star rating on the named scale."""
    return _doc(
        {"options": [{"option_id": "p", "blocks": [{
            "kind": "table",
            "from": "p.summary",
            "columns": [{"label": "Verdict", "from": "verdicts",
                         "items": {"from": "status", "stars": scale}}],
        }]}]},
        plugins=_LIST_PLUGIN,
    )


def test_stars_must_name_a_known_scale():
    """A typo'd scale would otherwise render no stars — which is exactly what an
    unrecognised term legitimately does, so it would read as data, not a bug."""
    with pytest.raises(ValidationError, match="unknown rating scale"):
        MergedSpec.model_validate(_starred_column("noscale"))


def test_an_item_can_carry_a_known_scale():
    doc = _starred_column("confidence")
    doc["display"]["rating_scales"] = _SCALE
    spec = MergedSpec.model_validate(doc)
    column = spec.display.options[0].blocks[0].columns[0]
    assert column.items.stars == "confidence"


def test_stars_inside_an_expanded_cell_are_checked():
    """The reference can sit two levels down — a column's `items`, and the cells
    of the detail those items expand onto."""
    doc = _doc(
        {
            "options": [{"option_id": "revel", "blocks": [{
                "kind": "table",
                "from": "revel.score",
                "columns": [{
                    "label": "Score",
                    "from": "score",
                    "items": {
                        "from": "score",
                        "expand": {
                            "from": "detail",
                            "cells": [{"from": "status", "stars": "noscale"}],
                        },
                    },
                }],
            }]}]
        }
    )
    with pytest.raises(ValidationError, match="unknown rating scale"):
        MergedSpec.model_validate(doc)


def test_a_rating_cannot_exceed_its_scale():
    from app.vep.models.display_spec_model import RatingScale

    with pytest.raises(ValidationError, match=r"outside 0\.\.4"):
        RatingScale.model_validate({"out_of": 4, "ratings": {"impossible": 5}})


def test_the_payload_carries_the_scales():
    """The frontend renders the stars, so the term -> rating table has to reach
    it — this is the only path it travels."""
    payload = SPEC.display_payload()
    assert payload.rating_scales["clinvar_aggregate"].out_of == 4
    assert (
        payload.rating_scales["clinvar_aggregate"].ratings["practice guideline"] == 4
    )


# --- a row that stacks a list ----------------------------------------------


_LIST_PLUGIN = [
    {
        "plugin": "p",
        "scope": "allele",
        "output": "p",
        "csq_fields": ["C"],
        "targets": [
            {
                "field": "summary",
                "from": "C",
                "transform": "chunk",
                "size": 2,
                "as": [
                    {"field": "kind", "type": "string"},
                    {"field": "verdict", "type": "string"},
                ],
                "item_fields": ["kind", "verdict"],
            },
            {
                "field": "detail",
                "from": "C",
                "transform": "chunk",
                "size": 2,
                "as": [
                    {"field": "kind", "type": "string"},
                    {"field": "status", "type": "string"},
                ],
                "item_fields": ["kind", "status"],
            },
        ],
        # A cell can only hold a list of objects if something put one there, and
        # a join is the one thing that does -- so a fixture for an itemised
        # column needs one, exactly as the real specs do.
        "joins": [
            {
                "into": "summary",
                "from": "detail",
                "left_key": "kind",
                "right_key": "kind",
                "as": "verdicts",
            }
        ],
    }
]


def _stacking_row(*cells):
    return _doc(
        {
            "options": [{"option_id": "p", "blocks": [{
                "kind": "rows",
                "rows": [{
                    "label": "Classification",
                    "from": "p.summary",
                    "item": {"cells": list(cells)},
                }],
            }]}]
        },
        plugins=_LIST_PLUGIN,
    )


def test_a_row_can_stack_a_list_of_items():
    spec = MergedSpec.model_validate(
        _stacking_row({"from": "kind"}, {"from": "verdict"})
    )
    row = spec.display.options[0].blocks[0].rows[0]
    assert row.list_ref() == ("p", "summary")
    assert [c.source for c in row.item.cells] == ["kind", "verdict"]


def test_a_stacking_rows_item_fields_are_checked():
    """The same check a list block's item gets — otherwise a typo'd field just
    renders an empty line."""
    with pytest.raises(ValidationError, match="item field 'verdcit' not in"):
        MergedSpec.model_validate(_stacking_row({"from": "verdcit"}))


def test_stars_from_and_template_fields_are_item_refs_too():
    """Both name *fields of the element*, so both are typo-checked -- including
    the `{field}` placeholders of a cell's template."""
    with pytest.raises(ValidationError, match="item field 'scale' not in"):
        MergedSpec.model_validate(
            _stacking_row({"from": "kind", "stars_from": "scale"})
        )
    with pytest.raises(ValidationError, match="item field 'total' not in"):
        MergedSpec.model_validate(
            _stacking_row({"from": "kind", "template": "{kind} of {total}"})
        )


def test_the_fixture_omits_only_what_the_defaults_put_back():
    """The guard that replaces byte-identity.

    The generated fixture is dumped with `exclude_none`, so it no longer lists
    every field a model could have. That is only safe if what it leaves out is
    exactly what loading it back fills in -- otherwise the frontend renders from
    a document that says less than the spec does.

    Byte-identity used to serve this purpose, but it cannot survive a model
    gaining a field, which is precisely what converging the value types does.

    `exclude_defaults` would be the tighter filter and is the wrong one: it also
    drops non-None defaults, and the frontend cannot put those back. The house
    truncation is one, and dropping it stopped every long list truncating.
    """
    payload = SPEC.display_payload()
    assert payload is not None
    lean = payload.model_dump(mode="json", by_alias=True, exclude_none=True)
    # `type(payload)`, not an import: this suite reaches the models as
    # `app.vep.models...` while the loader builds them as `vep.models...`, and
    # those are two different classes, on which pydantic equality is always
    # False.
    assert type(payload).model_validate(lean) == payload
