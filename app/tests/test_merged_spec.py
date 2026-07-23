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
    """expected_csq_columns for a submission with these options set (over the
    bundled spec), using ConfigIniParams so sub-option defaults are realistic."""
    params = ConfigIniParams(genome_id="g", assembly_name="GRCh38", **options)
    return SPEC.expected_csq_columns(params.model_dump())


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
    assert len(spec.config_entries()) == 30
    assert len(spec.parse_plugins()) == 27


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


def test_custom_literal_expects_exact_columns():
    assert _expected(clinvar=True) == {"ClinVar_CLNSIG", "ClinVar_CLNSIGCONF"}


def test_custom_builder_expects_the_combinatorial_columns():
    # default gnomAD_exomes = All + Both + UKB -> the overall AF column
    assert _expected(gnomad_exomes=True) == {"gnomAD_exomes_AF"}
    # adding the afr ancestry (both) adds its column, via the same builder that
    # writes the config `fields=`
    assert _expected(gnomad_exomes=True, gnomad_exomes_afr=True) == {
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
    # ProtVar has from_option sub-flags (a sub-option can drop a column), so it is
    # excluded entirely — even with pocket off, nothing is (wrongly) required.
    assert _expected(protvar=True) == set()
    assert _expected(protvar=True, protvar_pocket=False) == set()
    # IntAct (variadic flags) and mutfunc (from_option sub-flags) likewise
    assert _expected(intact=True) == set()
    assert _expected(mutfunc=True) == set()


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
