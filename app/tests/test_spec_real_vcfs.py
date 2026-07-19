"""Differential test: spec-driven parsing vs the hand-written `_parse_*` bank,
run over the *real* dev-data VCFs rather than synthetic fixtures.

test_spec_interpreter.py proves the interpreter matches the parsers on hand-built
rows. This runs the same comparison over actual pipeline output, which is the
stronger check: it exercises the column layouts, '&'-lists, 'NA' placeholders and
per-ancestry frequency columns that real runs emit and fixtures only approximate.
It is how the interpreter's equivalence was validated while the specs were
written, kept here as a standing gate for the results-time wiring.

The dev-data VCFs are gitignored (local pipeline output), so this skips wherever
they are absent — CI included. It is a local pre-flight, not a CI gate.

Excluded from the comparison, deliberately:
  - pathogenicity: the hand-written `_parse_pathogenicity` bundles many
    annotations into one object and is being dropped. Its data is covered by the
    flat revel/alphamissense/cadd/eve + spliceai/popeve plugins, none of which
    has a standalone hand-written parser to diff against.
  - uniprot, protein_matches, SIFT, PolyPhen: not specced (no plugin spec).
  - utr_annotation's `annotation` field only: the parser keeps the raw string,
    the spec parses it into an order-independent dict via `key_value`, so the two
    representations differ by design (5UTR_annotation key order is not stable per
    record). Every other utr field is compared.
  - a GO term with an empty name: the spec yields name=None, `_parse_go` yields
    name="" (the intended value is None; the parser is deleted at the cutover).
    The empty string is coerced to None for the comparison; any other GO
    difference still fails.
"""

from pathlib import Path

import pytest
import vcfpy

from app.vep.models.parsing_spec_model import ParsingSpec, PluginSpec
from app.vep.utils.csq import get_prediction_index_map
from app.vep.utils.spec_interpreter import apply_plugin_spec
from app.vep.utils.spec_loader import load_merged_spec
from app.vep.utils.vcf_results import (
    _parse_clinvar,
    _parse_dosage_sensitivity,
    _parse_frequencies,
    _parse_go,
    _parse_hgvs,
    _parse_intact,
    _parse_mavedb,
    _parse_mutfunc,
    _parse_open_targets,
    _parse_phenotype_data,
    _parse_popeve,
    _parse_population_frequencies,
    _parse_protvar,
    _parse_riboseq_orfs,
    _parse_spliceai,
    _parse_utr_annotation,
)

DEV_DATA = Path(__file__).resolve().parent.parent.parent / "dev-data"
VCF_FILES = ["output.vcf.gz", "has_utr.vcf.gz", "smaller_test.vcf.gz"]
# The small files run in full; the big one is capped so the test stays quick
# (its records carry ~75 CSQ entries each). Coverage of the plugin surface comes
# from breadth of columns, not depth of records.
RECORD_CAP = {"output.vcf.gz": 400}

SPEC: ParsingSpec = load_merged_spec("human_grch38").parsing


def _dump(model):
    return model.model_dump() if model is not None else None


# Plugins whose interpreter output equals `parser(...).model_dump()` outright.
_SIMPLE_PARSERS = {
    "mutfunc": _parse_mutfunc,
    "mavedb": _parse_mavedb,
    "clinvar": _parse_clinvar,
    "protvar": _parse_protvar,
    "intact": _parse_intact,
    "hgvs": _parse_hgvs,
    "dosage_sensitivity": _parse_dosage_sensitivity,
    "riboseq_orfs": _parse_riboseq_orfs,
    "phenotype_data": _parse_phenotype_data,
    "spliceai": _parse_spliceai,
    "popeve": _parse_popeve,
    "opentargets": _parse_open_targets,
}

# Resolve every PluginSpec once; a per-row SPEC.plugin() lookup would rescan the
# 21-plugin list on every CSQ entry.
_SIMPLE_SPECS: dict[str, PluginSpec] = {}
for _name in list(_SIMPLE_PARSERS) + [
    "gnomad_exomes", "gnomad_genomes", "all_of_us", "go", "utr_annotation"
]:
    _spec = SPEC.plugin(_name)
    assert _spec is not None, f"no spec for {_name}"
    _SIMPLE_SPECS[_name] = _spec


def _oracle_pop_freq(csq, index_map, overall_key, pattern):
    """A gnomAD source: the shared model always carries `max_subpopulation`
    (an All of Us concept) as None, but the gnomAD specs never emit that key."""
    freq = _parse_population_frequencies(csq, index_map, overall_key, pattern)
    if freq is None:
        return None
    data = freq.model_dump()
    data.pop("max_subpopulation", None)
    return data


def _oracle_all_of_us(csq, index_map):
    """All of Us: `max_subpopulation` is attached during composition, so the
    oracle is the composer's `all_of_us`, keeping that key."""
    freq = _parse_frequencies(csq, index_map)
    if freq is None or freq.all_of_us is None:
        return None
    return freq.all_of_us.model_dump()


def _check_row(csq_values: list[str], index_map: dict[str, int]) -> None:
    for name, parser in _SIMPLE_PARSERS.items():
        assert apply_plugin_spec(csq_values, index_map, _SIMPLE_SPECS[name]) == _dump(
            parser(csq_values, index_map)
        ), name

    assert apply_plugin_spec(
        csq_values, index_map, _SIMPLE_SPECS["gnomad_exomes"]
    ) == _oracle_pop_freq(
        csq_values, index_map, "gnomAD_exomes_AF", "gnomAD_exomes_AF_{}"
    ), "gnomad_exomes"
    assert apply_plugin_spec(
        csq_values, index_map, _SIMPLE_SPECS["gnomad_genomes"]
    ) == _oracle_pop_freq(
        csq_values, index_map, "gnomAD_genomes_AF", "gnomAD_genomes_AF_{}"
    ), "gnomad_genomes"
    assert apply_plugin_spec(
        csq_values, index_map, _SIMPLE_SPECS["all_of_us"]
    ) == _oracle_all_of_us(csq_values, index_map), "all_of_us"

    # go: interpreter wraps the term list in `{go_terms: [...]}` (or None when
    # empty); the parser returns the bare list. One accepted divergence: for a GO
    # entry with no term name the spec yields name=None while `_parse_go` yields
    # name="" (~4% of real GO entries). None is the intended value — see
    # test_go_entry_without_a_term_name_is_null_not_empty_string — and `_parse_go`
    # is deleted at the flat cutover, so the "" is coerced to None here rather than
    # tracked as a failure. Any other GO difference still fails.
    interp_go = apply_plugin_spec(csq_values, index_map, _SIMPLE_SPECS["go"])
    interp_terms = interp_go["go_terms"] if interp_go else []
    parser_terms = []
    for term in _parse_go(csq_values, index_map):
        data = term.model_dump()
        data["name"] = data["name"] or None
        parser_terms.append(data)
    assert interp_terms == parser_terms, "go"

    # utr_annotation: every field but `annotation` (raw string vs key_value dict).
    interp_utr = apply_plugin_spec(
        csq_values, index_map, _SIMPLE_SPECS["utr_annotation"]
    )
    parser_utr = _dump(_parse_utr_annotation(csq_values, index_map))
    if interp_utr is None or parser_utr is None:
        assert interp_utr == parser_utr, "utr_annotation presence"
    else:
        fields = [key for key in parser_utr if key != "annotation"]
        assert {key: interp_utr[key] for key in fields} == {
            key: parser_utr[key] for key in fields
        }, "utr_annotation fields"


@pytest.mark.parametrize("vcf_name", VCF_FILES)
def test_interpreter_matches_hand_written_over_real_vcf(vcf_name):
    vcf_path = DEV_DATA / vcf_name
    if not vcf_path.exists():
        pytest.skip(f"dev-data VCF not present: {vcf_path}")

    reader = vcfpy.Reader.from_path(str(vcf_path))
    description = reader.header.get_info_field_info("CSQ").description
    index_map = get_prediction_index_map(description)
    cap = RECORD_CAP.get(vcf_name)

    checked = 0
    for record_number, record in enumerate(reader):
        if cap is not None and record_number >= cap:
            break
        for csq_string in record.INFO.get("CSQ", []):
            _check_row(csq_string.split("|"), index_map)
            checked += 1

    assert checked > 0, f"no CSQ rows checked in {vcf_name}"
