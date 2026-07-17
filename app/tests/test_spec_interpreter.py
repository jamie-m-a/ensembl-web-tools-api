"""Differential tests: spec-driven parsing vs the hand-written `_parse_*` bank.

The hand-written parsers are the oracle. For the same CSQ entry, the interpreter
driven by `parsing_specs/human_grch38.json` must produce exactly what the
corresponding `_parse_*` produces (compared as plain data, via model_dump).

This is what proves the spec vocabulary is sufficient before anything is
rewired, so the fixtures are deliberately shared with test_csq_parsers.
"""

import pytest

from app.tests.test_csq_parsers import EMPTY, INDEX_MAP, row_list
from app.vep.models.parsing_spec_model import ParsingSpec, TargetSpec
from app.vep.utils.csq import get_prediction_index_map
from app.vep.utils.spec_interpreter import apply_plugin_spec
from app.vep.utils.spec_loader import SPEC_DIR, load_spec_file
from app.vep.utils.vcf_results import (
    _parse_clinvar,
    _parse_frequencies,
    _parse_mavedb,
    _parse_mutfunc,
    _parse_population_frequencies,
)

SPEC: ParsingSpec = load_spec_file(SPEC_DIR / "human_grch38.json")


def dump(model):
    """A parser's output as plain data, or None."""
    return model.model_dump() if model is not None else None


def dump_frequencies(model):
    """A gnomAD PopulationFrequencies as plain data.

    `max_subpopulation` is dropped: it is an All of Us concept (the label column
    AoU_gvs_max_subpop), so the gnomAD specs do not produce the key at all,
    while the shared pydantic model always carries it as None. The AoU tests
    compare against _parse_frequencies instead and do not use this helper.
    """
    if model is None:
        return None
    data = model.model_dump()
    data.pop("max_subpopulation", None)
    return data


def index_map_for(*columns: str) -> dict[str, int]:
    return get_prediction_index_map("Format: " + "|".join(columns))


def run(plugin: str, csq_values, index_map=INDEX_MAP):
    spec = SPEC.plugin(plugin)
    assert spec is not None, f"no spec for {plugin}"
    return apply_plugin_spec(csq_values, index_map, spec)


# --- the spec document itself ------------------------------------------------


def test_bundled_spec_validates():
    """The shipped JSON round-trips through the strict model."""
    assert SPEC.spec_version
    assert {p.plugin for p in SPEC.plugins} == {
        "mutfunc",
        "mavedb",
        "clinvar",
        "gnomad_exomes",
        "gnomad_genomes",
        "all_of_us",
    }


def test_unknown_key_is_rejected():
    """extra=forbid: a spec we don't understand fails at load, not at parse."""
    with pytest.raises(Exception):
        ParsingSpec.model_validate(
            {"spec_version": "x", "plugins": [], "surprise": True}
        )


def test_zip_requires_matching_as_entries():
    with pytest.raises(Exception):
        TargetSpec.model_validate(
            {
                "field": "assays",
                "from": ["a", "b"],
                "transform": "zip",
                "as": [{"field": "only_one"}],
            }
        )


# --- mutfunc: four scalars ---------------------------------------------------

MUTFUNC_SCORES = dict(
    mutfunc_motif="0.1", mutfunc_int="0.2", mutfunc_mod="0.3", mutfunc_exp="0.4"
)


def test_mutfunc_matches_hand_written_parser():
    csq = row_list(**MUTFUNC_SCORES)
    assert run("mutfunc", csq) == dump(_parse_mutfunc(csq, INDEX_MAP))


def test_mutfunc_empty_matches():
    assert run("mutfunc", EMPTY) == dump(_parse_mutfunc(EMPTY, INDEX_MAP)) == None


def test_mutfunc_partial_matches():
    """Only some scores present: the rest must come back None, not be dropped."""
    csq = row_list(mutfunc_motif="0.1", mutfunc_exp="0.4")
    assert run("mutfunc", csq) == dump(_parse_mutfunc(csq, INDEX_MAP))


# --- MaveDB: positional zip, the hard case -----------------------------------

MAVEDB_MULTI = dict(
    MaveDB_score="1.5&2.5&NA",
    MaveDB_urn="urn:1&urn:2&urn:3",
    MaveDB_doi="10.1/a&NA&10.1/c",
    MaveDB_nt="c.1A>G&NA",
    MaveDB_pro="p.Lys1Arg&NA",
)


def test_mavedb_multi_assay_matches_hand_written_parser():
    csq = row_list(**MAVEDB_MULTI)
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP))


def test_mavedb_multi_assay_shape_is_as_expected():
    """Pin the actual expected values, so a bug in *both* paths can't pass."""
    result = run("mavedb", row_list(**MAVEDB_MULTI))
    assert result["protein_variant"] == "p.Lys1Arg"
    assert [(a["urn"], a["score"]) for a in result["assays"]] == [
        ("urn:1", 1.5),
        ("urn:2", 2.5),
        ("urn:3", None),  # NA score, but a real urn -> assay kept
    ]


def test_mavedb_empty_matches():
    assert run("mavedb", EMPTY) == dump(_parse_mavedb(EMPTY, INDEX_MAP)) == None


def test_mavedb_uneven_columns_match():
    """Fewer scores than urns: `align: max` must pad, not truncate."""
    csq = row_list(MaveDB_score="1.5", MaveDB_urn="urn:1&urn:2")
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP))


def test_mavedb_protein_variant_only_is_none():
    """pro present but no score/urn -> no assays -> whole annotation is None
    (require_any_output), matching the hand-written parser."""
    csq = row_list(MaveDB_pro="p.Lys1Arg")
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP)) == None


def test_mavedb_all_na_assay_dropped():
    """A position where both score and urn are NA is dropped entirely."""
    csq = row_list(MaveDB_score="1.5&NA", MaveDB_urn="urn:1&NA")
    assert run("mavedb", csq) == dump(_parse_mavedb(csq, INDEX_MAP))


# --- ClinVar: the `when` conditional -----------------------------------------

CONFLICTING = "Conflicting_classifications_of_pathogenicity"


def test_clinvar_conflicting_reads_breakdown():
    csq = row_list(
        ClinVar_CLNSIG=CONFLICTING,
        ClinVar_CLNSIGCONF="Likely_pathogenic_(6)&Benign_(2)",
    )
    assert run("clinvar", csq) == dump(_parse_clinvar(csq, INDEX_MAP))


def test_clinvar_conflicting_breakdown_shape():
    result = run(
        "clinvar",
        row_list(
            ClinVar_CLNSIG=CONFLICTING,
            ClinVar_CLNSIGCONF="Likely_pathogenic_(6)&Benign_(2)",
        ),
    )
    assert result["significance"] == [CONFLICTING]
    assert result["conflicting_breakdown"] == [
        {"significance": "Likely_pathogenic", "count": 6},
        {"significance": "Benign", "count": 2},
    ]


def test_clinvar_non_conflicting_ignores_breakdown():
    """The `when` gate: CLNSIGCONF is present but must not be read, because the
    classification is not conflicting."""
    csq = row_list(
        ClinVar_CLNSIG="Pathogenic",
        ClinVar_CLNSIGCONF="Likely_pathogenic_(6)",
    )
    result = run("clinvar", csq)
    assert result == dump(_parse_clinvar(csq, INDEX_MAP))
    assert result["conflicting_breakdown"] == []


def test_clinvar_when_matches_list_membership_not_substring():
    """A value that merely embeds the conflicting term must not trigger the
    breakdown — the condition is membership of the '&'-split list."""
    csq = row_list(
        ClinVar_CLNSIG="Not_" + CONFLICTING,
        ClinVar_CLNSIGCONF="Benign_(2)",
    )
    result = run("clinvar", csq)
    assert result == dump(_parse_clinvar(csq, INDEX_MAP))
    assert result["conflicting_breakdown"] == []


def test_clinvar_unparseable_breakdown_token_skipped():
    csq = row_list(
        ClinVar_CLNSIG=CONFLICTING,
        ClinVar_CLNSIGCONF="Benign_(2)&garbage_no_count",
    )
    result = run("clinvar", csq)
    assert result == dump(_parse_clinvar(csq, INDEX_MAP))
    assert [b["significance"] for b in result["conflicting_breakdown"]] == ["Benign"]


def test_clinvar_empty_matches():
    assert run("clinvar", EMPTY) == dump(_parse_clinvar(EMPTY, INDEX_MAP)) == None


# --- gnomAD / All of Us: pattern_map -----------------------------------------


def test_gnomad_exomes_pattern_map_matches():
    columns = ["gnomAD_exomes_AF", "gnomAD_exomes_AF_afr", "gnomAD_exomes_AF_nfe_XX"]
    index_map = index_map_for(*columns)
    values = ["0.01", "0.02", "0.03"]

    result = run("gnomad_exomes", values, index_map)
    oracle = _parse_population_frequencies(
        values, index_map, "gnomAD_exomes_AF", "gnomAD_exomes_AF_{}"
    )
    assert result == dump_frequencies(oracle)
    # ancestry columns discovered from the header, not named in the spec
    assert result["populations"] == {"afr": 0.02, "nfe_XX": 0.03}
    assert result["overall"] == 0.01


def test_gnomad_exomes_zero_overall_is_kept():
    """A 0.0 frequency is a real value. require_any_output must not treat it as
    absent (plain truthiness would drop the whole annotation)."""
    columns = ["gnomAD_exomes_AF"]
    index_map = index_map_for(*columns)
    values = ["0.0"]

    result = run("gnomad_exomes", values, index_map)
    oracle = _parse_population_frequencies(
        values, index_map, "gnomAD_exomes_AF", "gnomAD_exomes_AF_{}"
    )
    assert result == dump_frequencies(oracle)
    assert result is not None
    assert result["overall"] == 0.0


def test_gnomad_exomes_absent_matches():
    index_map = index_map_for("Allele")
    assert run("gnomad_exomes", ["A"], index_map) is None


def test_gnomad_exomes_legacy_prefix_ignored():
    """The old gnomADe_ prefix must not match the pattern."""
    columns = ["gnomADe_AF", "gnomADe_afr_AF"]
    index_map = index_map_for(*columns)
    assert run("gnomad_exomes", ["0.1", "0.2"], index_map) is None


def test_gnomad_genomes_pattern_map_matches():
    columns = ["gnomAD_genomes_AF", "gnomAD_genomes_AF_ami", "gnomAD_genomes_AF_grpmax"]
    index_map = index_map_for(*columns)
    values = ["0.10", "0.20", "0.30"]

    result = run("gnomad_genomes", values, index_map)
    oracle = _parse_population_frequencies(
        values, index_map, "gnomAD_genomes_AF", "gnomAD_genomes_AF_{}"
    )
    assert result == dump_frequencies(oracle)
    assert result["populations"] == {"ami": 0.20, "grpmax": 0.30}


def test_all_of_us_pattern_map_with_suffix_matches():
    """AoU's pattern has a suffix (AoU_gvs_{pop}_af), unlike gnomAD's.

    The oracle is _parse_frequencies (the composer), not
    _parse_population_frequencies, because max_subpopulation is attached during
    composition — reproducing that attach is exactly what the spec's
    max_subpopulation target has to do.
    """
    columns = ["AoU_gvs_all_af", "AoU_gvs_afr_af", "AoU_gvs_max_af", "AoU_gvs_max_subpop"]
    index_map = index_map_for(*columns)
    values = ["0.10", "0.20", "0.30", "eur"]

    result = run("all_of_us", values, index_map)
    oracle = _parse_frequencies(values, index_map).all_of_us
    assert result == oracle.model_dump()
    assert result["overall"] == 0.10
    assert result["populations"] == {"afr": 0.20, "max": 0.30}
    assert result["max_subpopulation"] == "eur"
    # the label column is not a frequency and must not appear among populations
    assert "max_subpop" not in result["populations"]


def test_all_of_us_label_without_frequencies_is_none():
    """A max_subpop label with no frequencies is not an annotation — which is
    why max_subpopulation is deliberately absent from require_any_output."""
    index_map = index_map_for("AoU_gvs_all_af", "AoU_gvs_max_subpop")
    values = ["", "eur"]

    assert run("all_of_us", values, index_map) is None
    assert _parse_frequencies(values, index_map) is None
