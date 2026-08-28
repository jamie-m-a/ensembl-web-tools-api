"""Tests for the config<->parsing consistency check on MergedSpec
(merged_spec_model.py, design section 6.1).

The check runs as a model_validator at load time, so a bad merged document fails
loudly the moment it is loaded rather than producing empty annotations later.
"""

import logging

import pytest
from pydantic import ValidationError

from app.vep.models.merged_spec_model import MergedSpec
from app.vep.models.pipeline_model import ConfigIniParams
from app.vep.utils.spec_loader import load_merged_spec

SPEC = load_merged_spec("human_grch38")


def _expected(**options):
    """expected_csq_columns for a submission in this option state (over the
    bundled spec), through ConfigIniParams so sub-option defaults are realistic.

    Laid over the resolved map rather than replacing it, because a couple of the
    states worth testing are not ones a *client* can send: `hgvsg` has no form
    control and is switched on by ProtVar's `forces_on`, so it never appears in
    a submitted payload but is very much a state the check must handle.
    """
    params = ConfigIniParams(
        genome_id="g", assembly_name="GRCh38", options=options
    )
    return SPEC.expected_csq_columns({**params.options, **options})


def _plugin(plugin_id, csq_fields, *, scope="transcript"):
    return {
        "plugin": plugin_id,
        "scope": scope,
        "output": plugin_id,
        "csq_fields": csq_fields,
        "targets": [],
    }


def _doc(config_entries, parse_plugins):
    return {
        "genome": {"assembly": "GRCh38"},
        "config": {"entries": config_entries},
        "parsing": {"plugins": parse_plugins},
    }


# --- the shipped document ---------------------------------------------------


def test_bundled_merged_spec_is_consistent():
    # load_merged_spec runs the consistency check; a bad spec would raise here.
    spec = load_merged_spec("human_grch38")
    assert len(spec.config_entries()) == 38
    # One more parse plugin than config entries: the Phenotypes option feeds two,
    # splitting gene-associated phenotypes (narrowed to the row's own gene) from
    # variant-associated ones (narrowed to the row's allele).
    assert len(spec.parse_plugins()) == 39


# --- reference integrity ----------------------------------------------------


def test_unknown_parse_plugin_reference_raises():
    doc = _doc(
        [
            {
                "id": "revel",
                "order": 1,
                "parsed_as": ["does_not_exist"],
                "config": {"emit": "plugin", "name": "REVEL", "params": {"file": "x"}},
            }
        ],
        [_plugin("revel", ["REVEL"])],
    )
    with pytest.raises(ValidationError, match="unknown parse plugin 'does_not_exist'"):
        MergedSpec.model_validate(doc)


# --- config-only entries need no parser -------------------------------------


def test_config_only_entry_needs_no_parser():
    doc = _doc(
        [{"id": "spdi", "order": 1, "parsed_as": [], "config": {"emit": "flag", "keyword": "spdi"}}],
        [],
    )
    MergedSpec.model_validate(doc)  # no raise


# --- the non-1:1 relations the sibling-section shape exists to support -------


def test_one_config_to_many_parse_is_valid():
    # eve -> {eve, popeve}
    doc = _doc(
        [
            {
                "id": "eve",
                "order": 1,
                "parsed_as": ["eve", "popeve"],
                "config": {"emit": "plugin", "name": "EVE", "params": {"file": "x"}},
            }
        ],
        [_plugin("eve", ["EVE_CLASS"]), _plugin("popeve", ["popEVE_SCORE"])],
    )
    MergedSpec.model_validate(doc)  # no raise


def test_many_config_to_one_parse_is_valid():
    # {hgvs, hgvsg} -> hgvs
    doc = _doc(
        [
            {"id": "hgvs", "order": 1, "parsed_as": ["hgvs"], "config": {"emit": "flag", "keyword": "hgvs"}},
            {"id": "hgvsg", "order": 2, "parsed_as": ["hgvs"], "config": {"emit": "flag", "keyword": "hgvsg"}},
        ],
        [_plugin("hgvs", ["HGVSg", "HGVSc", "HGVSp"])],
    )
    MergedSpec.model_validate(doc)  # no raise


# --- custom column-level check ----------------------------------------------


def test_custom_literal_column_mismatch_raises():
    doc = _doc(
        [
            {
                "id": "clinvar",
                "order": 1,
                "parsed_as": ["clinvar"],
                "config": {
                    "emit": "custom",
                    "params": {"file": "x", "short_name": "ClinVar", "format": "vcf"},
                    "fields": {"literal": ["CLNSIG", "WRONG"]},
                },
            }
        ],
        [_plugin("clinvar", ["ClinVar_CLNSIG", "ClinVar_CLNSIGCONF"], scope="allele")],
    )
    with pytest.raises(ValidationError, match="ClinVar_WRONG"):
        MergedSpec.model_validate(doc)


def test_custom_literal_columns_that_match_are_valid():
    doc = _doc(
        [
            {
                "id": "clinvar",
                "order": 1,
                "parsed_as": ["clinvar"],
                "config": {
                    "emit": "custom",
                    "params": {"file": "x", "short_name": "ClinVar", "format": "vcf"},
                    "fields": {"literal": ["CLNSIG", "CLNSIGCONF"]},
                },
            }
        ],
        [_plugin("clinvar", ["ClinVar_CLNSIG", "ClinVar_CLNSIGCONF"], scope="allele")],
    )
    MergedSpec.model_validate(doc)  # no raise


def test_custom_builder_short_name_mismatch_raises():
    doc = _doc(
        [
            {
                "id": "gnomad_exomes",
                "order": 1,
                "parsed_as": ["gnomad_exomes"],
                "config": {
                    "emit": "custom",
                    "params": {"file": "x", "short_name": "WRONG_NAME", "format": "vcf"},
                    "fields": {
                        "builder": "gnomad_ancestry_sex",
                        "base": "AF",
                        "ancestries": [{"option": "gnomad_exomes_all", "code": ""}],
                        "sexes": [{"suffix": "both", "code": ""}],
                    },
                },
            }
        ],
        [_plugin("gnomad_exomes", ["gnomAD_exomes_AF"], scope="allele")],
    )
    with pytest.raises(ValidationError, match="WRONG_NAME"):
        MergedSpec.model_validate(doc)


def test_custom_builder_short_name_prefix_match_is_valid():
    doc = _doc(
        [
            {
                "id": "gnomad_exomes",
                "order": 1,
                "parsed_as": ["gnomad_exomes"],
                "config": {
                    "emit": "custom",
                    "params": {"file": "x", "short_name": "gnomAD_exomes", "format": "vcf"},
                    "fields": {
                        "builder": "gnomad_ancestry_sex",
                        "base": "AF",
                        "ancestries": [{"option": "gnomad_exomes_all", "code": ""}],
                        "sexes": [{"suffix": "both", "code": ""}],
                    },
                },
            }
        ],
        [_plugin("gnomad_exomes", ["gnomAD_exomes_AF"], scope="allele")],
    )
    MergedSpec.model_validate(doc)  # no raise


# --- soft: an unreachable parser ---------------------------------------------


def test_parse_plugin_with_no_config_is_a_soft_warning(caplog):
    doc = _doc(
        [
            {
                "id": "revel",
                "order": 1,
                "parsed_as": ["revel"],
                "config": {"emit": "plugin", "name": "REVEL", "params": {"file": "x"}},
            }
        ],
        [_plugin("revel", ["REVEL"]), _plugin("orphan", ["ORPHAN"])],
    )
    with caplog.at_level(logging.WARNING):
        MergedSpec.model_validate(doc)  # no raise — a soft signal, not an error
    assert "orphan" in caplog.text


# --- display `list` blocks: item-field refs are checked -----------------------


def _go_like_doc(cells):
    """A minimal merged doc with a `go`-like list plugin (item_fields id/name)
    and a display `list` block whose item is `cells`."""
    plugins = [
        {
            "plugin": "go",
            "scope": "transcript",
            "output": "go_terms",
            "csq_fields": ["GO"],
            "targets": [
                {
                    "field": "go_terms",
                    "from": "GO",
                    "transform": "list",
                    "item_fields": ["id", "name"],
                }
            ],
        }
    ]
    config = [
        {
            "id": "go",
            "order": 1,
            "parsed_as": ["go"],
            "config": {"emit": "flag", "keyword": "go"},
        }
    ]
    doc = _doc(config, plugins)
    doc["display"] = {
        "options": [
            {
                "option_id": "go",
                "blocks": [
                    {
                        "kind": "list",
                        "heading": "Gene Ontology",
                        "from": "go.go_terms",
                        "item": {"cells": cells},
                    }
                ],
            }
        ]
    }
    return doc


def test_display_list_block_valid_item_refs_load():
    doc = _go_like_doc(
        [
            {"from": "id", "link": {"kind": "external", "template": "x/{id}"}},
            {"from": "name"},
        ]
    )
    MergedSpec.model_validate(doc)  # no raise


def test_display_list_cell_unknown_item_field_raises():
    doc = _go_like_doc([{"from": "bogus"}])
    with pytest.raises(ValidationError, match="item field 'bogus'"):
        MergedSpec.model_validate(doc)


def test_display_list_link_template_unknown_item_field_raises():
    doc = _go_like_doc(
        [{"from": "id", "link": {"kind": "external", "template": "x/{missing}"}}]
    )
    with pytest.raises(ValidationError, match="item field 'missing'"):
        MergedSpec.model_validate(doc)


def test_display_list_unknown_list_field_raises():
    doc = _go_like_doc([{"from": "id"}])
    doc["display"]["options"][0]["blocks"][0]["from"] = "go.not_a_target"
    with pytest.raises(ValidationError, match="not_a_target"):
        MergedSpec.model_validate(doc)


def _set_item(doc, item):
    doc["display"]["options"][0]["blocks"][0]["item"] = item
    return doc


def test_display_list_item_rows_valid_refs_load():
    doc = _set_item(
        _go_like_doc([{"from": "id"}]),
        {"rows": [{"label": "Id", "from": "id"}, {"label": "Name", "from": "name"}]},
    )
    MergedSpec.model_validate(doc)  # no raise


def test_display_list_item_rows_unknown_field_raises():
    doc = _set_item(
        _go_like_doc([{"from": "id"}]), {"rows": [{"label": "X", "from": "bogus"}]}
    )
    with pytest.raises(ValidationError, match="item field 'bogus'"):
        MergedSpec.model_validate(doc)


def test_display_list_item_needs_exactly_one_of_cells_or_rows():
    # both cells and rows
    doc = _set_item(
        _go_like_doc([{"from": "id"}]),
        {"cells": [{"from": "id"}], "rows": [{"label": "Id", "from": "id"}]},
    )
    with pytest.raises(ValidationError, match="exactly one of `cells` or `rows`"):
        MergedSpec.model_validate(doc)
    # neither
    doc = _set_item(_go_like_doc([{"from": "id"}]), {})
    with pytest.raises(ValidationError, match="exactly one of `cells` or `rows`"):
        MergedSpec.model_validate(doc)


def test_bundled_display_has_list_options():
    spec = load_merged_spec("human_grch38")
    ids = {o.option_id for o in spec.display.options}
    assert {"phenotypes", "go"} <= ids


# --- expected_csq_columns (the per-job basis for the missing-field check) -----


def test_defaults_expect_nothing():
    # No annotation option is on by default, so nothing is required.
    assert _expected() == set()


def test_simple_plugin_expects_its_csq_fields():
    assert _expected(revel=True) == {"REVEL"}
    assert _expected(cadd=True) == {"CADD_PHRED", "CADD_RAW"}


CLINVAR_SHORT_COLUMNS = {
    "ClinVar",  # the bare match column a custom always emits
    "ClinVar_CLNDN",
    "ClinVar_CLNDNINCL",
    "ClinVar_CLNDISDB",
    "ClinVar_CLNDISDBINCL",
    "ClinVar_CLNREVSTAT",
    "ClinVar_CLNSIG",
    "ClinVar_CLNSIGINCL",
    "ClinVar_ONCDN",
    "ClinVar_ONCDNINCL",
    "ClinVar_ONCDISDB",
    "ClinVar_ONCDISDBINCL",
    "ClinVar_ONC",
    "ClinVar_ONCINCL",
    "ClinVar_ONCREVSTAT",
    "ClinVar_ONCSCV",
    "ClinVar_ONCCONF",
    "ClinVar_CLNSUBA",
    "ClinVar_CLNPMID",
    "ClinVar_CLNSUBN",
    "ClinVar_CLNRCV",
    "ClinVar_SCI",
    "ClinVar_SCIREVSTAT",
    "ClinVar_SCIDN",
    "ClinVar_SCIDISDB",
    "ClinVar_GENEINFO",
}

CLINVAR_SV_COLUMNS = {"ClinVar_SV", "ClinVar_SV_CLNSIG", "ClinVar_SV_ORIGIN"}


def test_custom_literal_expects_exact_columns():
    # The bare `short_name` match column is always emitted by a custom, so it is
    # expected too, alongside the literal `short_name_<field>` columns. The
    # germline custom rides in on Phenotypes (`forces_on`), which is what makes
    # its columns expected — alongside Phenotypes' own column.
    assert _expected(phenotypes=True) == CLINVAR_SHORT_COLUMNS | {"PHENOTYPES"}


def test_clinvar_short_expects_nothing_without_phenotypes():
    # A stale `clinvar_short` restored by edit/rerun expects nothing on its own —
    # the `requires: ["phenotypes"]` gate holds.
    assert _expected(clinvar_short=True) == set()


def test_clinvar_structural_stands_on_its_own():
    # It used to need a `clinvar` master that, once the germline data moved to
    # Phenotypes, gated nothing else -- so ticking the one control that means
    # something now brings the structural columns by itself.
    assert _expected(clinvar_sv=True) == CLINVAR_SV_COLUMNS
    # The two ClinVar sources stay independent: Phenotypes brings the germline
    # columns, this brings the structural ones, and neither implies the other.
    assert _expected(phenotypes=True) == CLINVAR_SHORT_COLUMNS | {"PHENOTYPES"}
    assert (
        _expected(phenotypes=True, clinvar_sv=True)
        == CLINVAR_SHORT_COLUMNS | CLINVAR_SV_COLUMNS | {"PHENOTYPES"}
    )


def test_custom_builder_expects_the_combinatorial_columns():
    # default gnomAD_exomes = All + Both + UKB -> the overall AF column, plus the
    # bare `short_name` match column (always emitted by a custom)
    assert _expected(gnomad_exomes=True) == {"gnomAD_exomes", "gnomAD_exomes_AF"}
    # adding the afr ancestry (both) adds its column, via the same builder that
    # writes the config `fields=`
    assert _expected(gnomad_exomes=True, gnomad_exomes_afr=True) == {
        "gnomAD_exomes",
        "gnomAD_exomes_AF",
        "gnomAD_exomes_AF_afr",
    }


def test_one_config_to_many_parse_expects_both():
    # eve config feeds both the eve and popeve parsers
    assert _expected(eve=True) == {
        "EVE_CLASS", "EVE_SCORE",
        "popEVE_SCORE", "popEVE_EVE", "popEVE_mutant",
    }


def test_sub_flagged_plugin_is_excluded():
    # ProtVar has from_option sub-flags (a sub-option can drop a column), so its
    # own columns are excluded entirely — even with pocket off, none is
    # (wrongly) required. HGVSg is there because ProtVar `forces_on` it: the
    # `--hgvsg` line really is emitted, and an allele-scoped flag column is
    # present for every variant, so it is expected like any other.
    assert _expected(protvar=True) == {"HGVSg"}
    assert _expected(protvar=True, protvar_pocket=False) == {"HGVSg"}
    # IntAct (variadic flags) and mutfunc (from_option sub-flags) likewise
    assert _expected(intact=True) == set()
    assert _expected(mutfunc=True) == set()


def test_a_forced_option_contributes_its_columns():
    """A forced option's config line is emitted, so its output must be there.

    If the expectation ignored `forces_on`, the missing-field check would go
    silent for exactly the data the option was forced on to produce — which is
    how ClinVar now rides in on Phenotypes.
    """
    # hgvsg is not selected; ProtVar turns it on, and its column comes with it.
    assert "HGVSg" in _expected(protvar=True)
    assert "HGVSg" not in _expected(protvar=False)


def test_flags_require_only_their_allele_scoped_columns():
    # Flag columns are conditional in general and excluded — EXCEPT allele-scoped
    # ones, which every variant carries: HGVSg (from --hgvsg) and SPDI (--spdi).
    assert _expected(hgvsg=True) == {"HGVSg"}
    assert _expected(spdi=True) == {"SPDI"}
    # protein is a flag with no parse plugin; hgvs (c/p) is transcript-scoped and
    # conditional, so neither requires anything.
    assert _expected(protein=True) == set()
    assert _expected(hgvs=True) == set()
    # combined
    assert _expected(hgvsg=True, spdi=True, protein=True) == {"HGVSg", "SPDI"}


def test_disabled_option_contributes_nothing():
    assert "REVEL" not in _expected(revel=False)


# --- Refs into a nested list, and scales named by the data ------------------ #
#
# ClinVar's conditions table reads three list levels deep: a condition holds the
# classifications its submitters gave, and each of those holds the submissions
# it counted. Only the top level used to be resolvable, so a typo below it
# loaded clean and rendered an empty cell. These probe the real assembled
# document rather than a synthetic one, because the shapes that were unchecked
# are the ones only it has.


def _assembled(name="human_grch38"):
    from app.vep.utils.spec_loader import _assemble_payload

    return _assemble_payload(name)


def _clinvar_conditions_table(doc):
    """The first conditions table of the phenotypes option.

    It iterates `clinvar.records` -- one row per RCV, ClinVar's own
    variant+condition aggregate -- with the CLNDN names joined onto it."""
    for option in doc["display"]["options"]:
        if option["option_id"] != "phenotypes":
            continue
        stack = list(option["blocks"])
        while stack:
            block = stack.pop(0)
            if block.get("kind") == "group":
                stack = list(block.get("blocks", [])) + stack
            elif block.get("from") == "clinvar.records":
                return block
    raise AssertionError("the conditions table is no longer where this expects")


def _classification_column(doc):
    for column in _clinvar_conditions_table(doc)["columns"]:
        if column.get("items"):
            return column
    raise AssertionError("no itemised column in the conditions table")


def test_assembled_spec_still_loads():
    # The control for every mutation below: unmutated, this must pass.
    MergedSpec.model_validate(_assembled())


def test_column_items_unknown_field_raises():
    doc = _assembled()
    _classification_column(doc)["items"]["from"] = "classifcation"  # sic
    with pytest.raises(ValidationError, match="item field 'classifcation'"):
        MergedSpec.model_validate(doc)


def test_column_items_unknown_count_field_raises():
    doc = _assembled()
    _classification_column(doc)["items"]["count_from"] = "bogus"
    with pytest.raises(ValidationError, match="item field 'bogus'"):
        MergedSpec.model_validate(doc)


def test_column_link_from_unknown_field_raises():
    # `link_from` sits on whichever value does the linking — the column itself,
    # or the items it renders one per line. Both are element refs and both are
    # checked, so this probes wherever the shipped spec puts it.
    doc = _assembled()
    changed = False
    for column in _clinvar_conditions_table(doc)["columns"]:
        for holder in (column, column.get("items") or {}):
            if holder.get("link_from"):
                holder["link_from"] = "bogus_url"
                changed = True
    assert changed, "no `link_from` left in the conditions table to probe"
    with pytest.raises(ValidationError, match="item field 'bogus_url'"):
        MergedSpec.model_validate(doc)


def test_expanded_cell_unknown_field_raises():
    # Two levels below the table's own list: a submission's field.
    doc = _assembled()
    _classification_column(doc)["items"]["expand"]["cells"][0]["from"] = "submiter"
    with pytest.raises(ValidationError, match="item field 'submiter'"):
        MergedSpec.model_validate(doc)


def test_expand_over_a_field_that_is_not_a_list_raises():
    doc = _assembled()
    _classification_column(doc)["items"]["expand"]["from"] = "count"
    with pytest.raises(ValidationError, match="not a list|does not produce"):
        MergedSpec.model_validate(doc)


def test_emphasis_unknown_field_raises():
    doc = _assembled()
    _classification_column(doc)["items"]["expand"]["emphasis"]["field"] = "bogus"
    with pytest.raises(ValidationError, match="item field 'bogus'"):
        MergedSpec.model_validate(doc)


def test_row_where_unknown_field_raises():
    # The table's `where` was checked; a stacked row's identical one was not.
    doc = _assembled()
    changed = False
    for option in doc["display"]["options"]:
        blocks = list(option["blocks"])
        while blocks:
            block = blocks.pop()
            blocks += block.get("blocks", []) or []
            for row in block.get("rows", []) or []:
                if row.get("where"):
                    row["where"]["field"] = "bogus"
                    changed = True
    assert changed, "no stacked row with a `where` left to probe"
    with pytest.raises(ValidationError, match="item field 'bogus'"):
        MergedSpec.model_validate(doc)


def test_scale_named_by_a_stack_constant_must_exist():
    # `stars_from` names a field, not a scale, so the display alone cannot check
    # it -- but the field is filled from the parse's `const`, which can be.
    doc = _assembled()
    for plugin in doc["parsing"]["plugins"]:
        if plugin["plugin"] != "clinvar":
            continue
        for target in plugin["targets"]:
            for group in target.get("of", []) or []:
                if "rating_scale" in group.get("const", {}):
                    group["const"]["rating_scale"] = "clinvar_agregate"  # sic
    with pytest.raises(ValidationError, match="clinvar_agregate"):
        MergedSpec.model_validate(doc)


def test_nested_item_format_must_suit_its_type():
    # `humanize` calls string methods; ClinVar's per-classification count is a
    # number the join produced, two lists below the option's own fields.
    doc = _assembled()
    items = _classification_column(doc)["items"]
    items["from"] = "count"
    items["format"] = "humanize"
    with pytest.raises(ValidationError, match="formats .*count.* as 'humanize'"):
        MergedSpec.model_validate(doc)


def test_expanded_cell_format_must_suit_its_type():
    doc = _assembled()
    cell = _classification_column(doc)["items"]["expand"]["cells"][0]
    cell["format"] = "num"
    with pytest.raises(ValidationError, match="formats .*submitter.* as 'num'"):
        MergedSpec.model_validate(doc)


def test_stacked_row_cell_format_must_suit_its_type():
    doc = _assembled()
    found = False
    for option in doc["display"]["options"]:
        blocks = list(option["blocks"])
        while blocks:
            block = blocks.pop()
            blocks += block.get("blocks", []) or []
            for row in block.get("rows", []) or []:
                for cell in (row.get("item") or {}).get("cells", []):
                    if cell.get("from") == "supporting":
                        cell["format"] = "humanize"
                        found = True
    assert found, "no stacked cell over a counted field left to probe"
    with pytest.raises(ValidationError, match="as 'humanize'"):
        MergedSpec.model_validate(doc)
