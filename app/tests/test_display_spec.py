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

from app.vep.models.display_spec_model import DisplayListBlock
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
    "hgvs", "hgvsg", "spdi", "alphamissense", "revel", "cadd", "spliceai",
    "loeuf", "dosage_sensitivity", "utrannotator", "riboseqorfs", "eve",
    # `list`-block options (repeat + truncate, migrated off frontend overrides)
    "phenotypes", "go", "mavedb",
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
                        {"field": "score", "from": "REVEL", "transform": "scalar"}
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
    invariant — every display ref (fixed-row `plugin.field`, a list block's
    `plugin.listField`, and each list cell's item field) resolves."""
    targets = {
        plugin.plugin: {t.field: t for t in plugin.targets}
        for plugin in SPEC.parse_plugins()
    }
    for option in SPEC.display.options:
        for plugin, field in option.scalar_field_refs():
            assert field in targets[plugin], (option.option_id, plugin, field)
        for block in option.blocks:
            if not isinstance(block, DisplayListBlock):
                continue
            plugin, list_field = block.list_ref()
            assert list_field in targets[plugin], (
                option.option_id, plugin, list_field
            )
            item_fields = set(targets[plugin][list_field].item_fields or [])
            for cell in block.item.cells:
                for item_field in cell.item_field_refs():
                    assert item_field in item_fields, (
                        option.option_id, plugin, list_field, item_field
                    )


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
