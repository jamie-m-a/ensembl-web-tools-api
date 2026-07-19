"""Tests for the config<->parsing consistency check on MergedSpec
(merged_spec_model.py, design section 6.1).

The check runs as a model_validator at load time, so a bad merged document fails
loudly the moment it is loaded rather than producing empty annotations later.
"""

import logging

import pytest
from pydantic import ValidationError

from app.vep.models.merged_spec_model import MergedSpec
from app.vep.utils.spec_loader import load_merged_spec


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
    assert len(spec.config_entries()) == 31
    assert len(spec.parse_plugins()) == 21


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
