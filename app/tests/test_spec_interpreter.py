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
    _parse_pathogenicity,
    _parse_dosage_sensitivity,
    _parse_hgvs,
    _parse_intact,
    _parse_phenotype_data,
    _parse_popeve,
    _parse_go,
    _parse_spliceai,
    _parse_open_targets,
    _parse_protvar,
    _parse_protvar_pocket,
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
        "protvar",
        "opentargets",
        "go",
        "spliceai",
        "riboseq_orfs",
        "hgvs",
        "phenotype_data",
        "dosage_sensitivity",
        "intact",
        "popeve",
        "revel",
        "alphamissense",
        "cadd",
        "eve",
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


# --- ProtVar: chunk + positional ---------------------------------------------
#
# The happy path matches the hand-written parser exactly. The edge cases below
# do NOT, deliberately: _parse_protvar_pocket collects only the parts that parse
# as a float and then assigns them in order, so one unparseable item silently
# shifts every later value into the wrong field. `positional` assigns strictly by
# index instead. Those tests document the divergence rather than enshrine it.

PROTVAR_FULL = dict(
    ProtVar_stability="0.42",
    ProtVar_pocket="POCKET1&-5.2&0.3&0.8&0.6&12.5&RES",
    ProtVar_int="PARTNER1&0.9&PARTNER2&0.8",
)


def test_protvar_well_formed_matches_hand_written_parser():
    csq = row_list(**PROTVAR_FULL)
    assert run("protvar", csq) == dump(_parse_protvar(csq, INDEX_MAP))


def test_protvar_shape_is_as_expected():
    result = run("protvar", row_list(**PROTVAR_FULL))
    assert result["structure_stability_score"] == 0.42

    pocket = result["pockets"][0]
    assert pocket["pocket_id"] == "POCKET1"
    assert pocket["energy"] == -5.2
    assert pocket["radius_of_gyration"] == 12.5
    # the trailing residues item is unnamed, so ignored -- but `raw` keeps it
    assert pocket["raw"] == PROTVAR_FULL["ProtVar_pocket"]

    assert [i["partner"] for i in result["interaction_interfaces"]] == [
        "PARTNER1",
        "PARTNER2",
    ]
    assert result["interaction_interfaces"][0]["score"] == 0.9
    assert result["interaction_interfaces"][0]["raw"] == "PARTNER1&0.9"


def test_protvar_odd_interaction_token_count_matches():
    """A trailing partner with no score: still one interface, score null."""
    csq = row_list(ProtVar_int="PARTNER1&0.9&PARTNER3")
    assert run("protvar", csq) == dump(_parse_protvar(csq, INDEX_MAP))


def test_protvar_empty_matches():
    assert run("protvar", EMPTY) == dump(_parse_protvar(EMPTY, INDEX_MAP)) == None


def test_protvar_pocket_missing_middle_value_does_not_shift():
    """An unparseable score empties only its own field, in both paths.

    This used to be a deliberate divergence: `positional` assigned by index
    while the parser compacted, mislabelling score as energy_per_volume and so
    on. The parser has since been fixed, so the two now agree.
    """
    raw = "POCKET1&-5.2&NA&0.8&0.6&12.5&RES"
    spec_pocket = run("protvar", row_list(ProtVar_pocket=raw))["pockets"][0]

    assert spec_pocket["energy"] == -5.2
    assert spec_pocket["energy_per_volume"] is None
    assert spec_pocket["score"] == 0.8
    assert spec_pocket["buriedness"] == 0.6
    assert spec_pocket["radius_of_gyration"] == 12.5

    assert spec_pocket == _parse_protvar_pocket(raw).model_dump()


def test_protvar_interaction_na_partner_is_nulled():
    """DIVERGENCE (spec is more consistent).

    The spec treats 'NA' as absent everywhere. The hand-written parser nulls
    'NA' for MaveDB urns but passes it through verbatim as a ProtVar partner.
    """
    interfaces = run("protvar", row_list(ProtVar_int="NA&0.9"))["interaction_interfaces"]
    assert interfaces[0]["partner"] is None

    parser = _parse_protvar(row_list(ProtVar_int="NA&0.9"), INDEX_MAP)
    assert parser.interaction_interfaces[0].partner == "NA"


# --- OpenTargets: align:min + dedup + sort -----------------------------------

OT_COLS = [
    "OpenTargets_gwasDiseases", "OpenTargets_gwasGeneId",
    "OpenTargets_gwasLocusToGeneScore", "OpenTargets_qtlGeneId",
    "OpenTargets_qtlBiosampleName",
]
OT_INDEX = index_map_for(*OT_COLS)


def ot_row(diseases="", genes="", l2g="", qtl_genes="", qtl_biosamples=""):
    return [diseases, genes, l2g, qtl_genes, qtl_biosamples]


def run_ot(**kwargs):
    return run("opentargets", ot_row(**kwargs), OT_INDEX)


def parse_ot(**kwargs):
    return dump(_parse_open_targets(ot_row(**kwargs), OT_INDEX))


def test_opentargets_sorts_strongest_first_and_matches_parser():
    args = dict(
        diseases="EFO_1&EFO_2&EFO_3",
        genes="ENSG1&ENSG2&ENSG3",
        l2g="0.1&0.9&NA",
    )
    result = run_ot(**args)
    assert result == parse_ot(**args)
    # descending by score; the unscored association goes last
    assert [(a["disease"], a["l2g_score"]) for a in result["gwas_associations"]] == [
        ("EFO_2", 0.9),
        ("EFO_1", 0.1),
        ("EFO_3", None),
    ]


def test_opentargets_dedups_repeated_rows():
    """The plugin emits duplicate rows -- dedup fires on 93% of real records."""
    args = dict(diseases="EFO_1&EFO_1", genes="ENSG1&ENSG1", l2g="0.5&0.5")
    result = run_ot(**args)
    assert result == parse_ot(**args)
    assert len(result["gwas_associations"]) == 1


def test_opentargets_drops_row_without_disease():
    args = dict(diseases="NA&EFO_2", genes="ENSG1&ENSG2", l2g="0.1&0.9")
    result = run_ot(**args)
    assert result == parse_ot(**args)
    assert [a["disease"] for a in result["gwas_associations"]] == ["EFO_2"]


def test_opentargets_misaligned_columns_truncate():
    """align:min. Real data contains ragged columns (3 diseases, 2 genes), so
    the plugin's positional alignment is not guaranteed; zip drops the excess.
    Faithful to the parser -- with ragged input the true pairing is unknowable.
    """
    args = dict(diseases="EFO_1&EFO_2&EFO_3", genes="ENSG1&ENSG2", l2g="0.1&0.9")
    result = run_ot(**args)
    assert result == parse_ot(**args)
    assert len(result["gwas_associations"]) == 2  # EFO_3 dropped


def test_opentargets_qtl_dedups_and_nulls_na_biosample():
    args = dict(qtl_genes="ENSG1&ENSG1&ENSG2", qtl_biosamples="liver&liver&NA")
    result = run_ot(**args)
    assert result == parse_ot(**args)
    assert result["qtl_associations"] == [
        {"gene_id": "ENSG1", "biosample": "liver"},
        {"gene_id": "ENSG2", "biosample": None},
    ]


def test_opentargets_empty_matches():
    assert run_ot() == parse_ot() == None


def test_opentargets_absent_columns_is_none():
    assert run("opentargets", ["A"], index_map_for("Allele")) is None


# --- GO: regex + replace/strip -----------------------------------------------

GO_INDEX = index_map_for("GO")


def test_go_terms_match_hand_written_parser():
    values = ["GO:0001558:regulation_of_cell_growth&GO:0005509:calcium_ion_binding"]
    result = run("go", values, GO_INDEX)
    assert result["go_terms"] == [t.model_dump() for t in _parse_go(values, GO_INDEX)]
    assert result["go_terms"] == [
        {"id": "GO:0001558", "name": "regulation of cell growth"},
        {"id": "GO:0005509", "name": "calcium ion binding"},
    ]


def test_go_entry_without_a_term_name_is_null_not_empty_string():
    """DIVERGENCE (deliberate).

    Real data carries ids with no name at all (38 of 368 distinct GO ids in
    dev-data/output.vcf.gz, e.g. "GO:0050911:"). The parser reports name "";
    the spec reports null, consistently with how it treats every other absent
    value. This is the only difference between the two paths across 182,786
    real GO annotations.
    """
    values = ["GO:0050911:"]
    spec_terms = run("go", values, GO_INDEX)["go_terms"]
    assert spec_terms == [{"id": "GO:0050911", "name": None}]
    assert [t.model_dump() for t in _parse_go(values, GO_INDEX)] == [
        {"id": "GO:0050911", "name": ""}
    ]


def test_go_entry_without_a_name_part_is_skipped():
    """Fewer than three ':'-parts is not a term, in either path."""
    values = ["GO:0001558"]
    assert run("go", values, GO_INDEX) is None
    assert _parse_go(values, GO_INDEX) == []


def test_go_absent_is_none():
    assert run("go", [""], GO_INDEX) is None


# --- SpliceAI: mixed float/int scalars ---------------------------------------

SPLICEAI_COLS = [
    "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL",
]
SPLICEAI_INDEX = index_map_for(*SPLICEAI_COLS)


def test_spliceai_matches_hand_written_parser():
    values = ["BRCA1", "0.01", "0.02", "0.03", "0.04", "-5", "10", "-20", "30"]
    result = run("spliceai", values, SPLICEAI_INDEX)
    assert result == dump(_parse_spliceai(values, SPLICEAI_INDEX))
    assert result["symbol"] == "BRCA1"
    assert result["ds_acceptor_gain"] == 0.01
    assert result["dp_donor_loss"] == 30


def test_spliceai_zero_scores_are_kept():
    """0.00 is the commonest real value (170k+ entries): a real score, not
    absence. require_any_output must not discard it."""
    values = ["BRCA1", "0.00", "", "", "", "", "", "", ""]
    result = run("spliceai", values, SPLICEAI_INDEX)
    assert result == dump(_parse_spliceai(values, SPLICEAI_INDEX))
    assert result is not None
    assert result["ds_acceptor_gain"] == 0.0


def test_spliceai_symbol_alone_is_not_an_annotation():
    values = ["BRCA1", "", "", "", "", "", "", "", ""]
    assert run("spliceai", values, SPLICEAI_INDEX) is None
    assert _parse_spliceai(values, SPLICEAI_INDEX) is None


# --- plain scalar/list plugins -----------------------------------------------
#
# No new vocabulary; all five verified against dev-data/output.vcf.gz with zero
# mismatches (hgvs 210,658 / phenotype_data 382,715 / dosage 393,079 /
# intact 25 / popeve 96,953 CSQ entries).


def test_hgvs_matches_hand_written_parser():
    index_map = index_map_for("HGVSg", "HGVSc", "HGVSp")
    values = ["NC_1:g.100A>G", "ENST1:c.50A>G", "ENSP1:p.Lys1Arg"]
    result = run("hgvs", values, index_map)
    assert result == dump(_parse_hgvs(values, index_map))
    assert result["transcript"] == "ENST1:c.50A>G"


def test_hgvs_partial_matches():
    """This run emits HGVSc/HGVSp but no HGVSg -- the absent one is null."""
    index_map = index_map_for("HGVSc", "HGVSp")
    values = ["ENST1:c.50A>G", ""]
    result = run("hgvs", values, index_map)
    assert result == dump(_parse_hgvs(values, index_map))
    assert result["genomic"] is None


def test_hgvs_empty_is_none():
    index_map = index_map_for("HGVSg", "HGVSc", "HGVSp")
    assert run("hgvs", ["", "", ""], index_map) is None


def test_phenotype_data_matches_hand_written_parser():
    index_map = index_map_for("PHENOTYPES", "CLIN_SIG", "PUBMED")
    values = ["cancer&diabetes", "pathogenic&NA", "123&456"]
    result = run("phenotype_data", values, index_map)
    assert result == dump(_parse_phenotype_data(values, index_map))
    assert result["phenotypes"] == ["cancer", "diabetes"]
    assert result["clinical_significance"] == ["pathogenic"]  # NA dropped


def test_phenotype_data_empty_is_none():
    index_map = index_map_for("PHENOTYPES", "CLIN_SIG", "PUBMED")
    assert run("phenotype_data", ["", "", ""], index_map) is None


def test_dosage_sensitivity_matches_hand_written_parser():
    index_map = index_map_for("pHaplo", "pTriplo")
    values = ["0.98", "0.12"]
    result = run("dosage_sensitivity", values, index_map)
    assert result == dump(_parse_dosage_sensitivity(values, index_map))
    assert result == {"phaplo": 0.98, "ptriplo": 0.12}


def test_dosage_sensitivity_zero_is_kept():
    """0.0 is a real probability, not absence."""
    index_map = index_map_for("pHaplo", "pTriplo")
    values = ["0.0", ""]
    result = run("dosage_sensitivity", values, index_map)
    assert result == dump(_parse_dosage_sensitivity(values, index_map))
    assert result["phaplo"] == 0.0


def test_intact_matches_hand_written_parser():
    """This run emits only three of IntAct's columns; the unselected
    sub-options are absent and come back null."""
    index_map = index_map_for(
        "IntAct_feature_type", "IntAct_interaction_ac", "IntAct_feature_ac"
    )
    values = ["mutation", "EBI-123", "EBI-ac"]
    result = run("intact", values, index_map)
    assert result == dump(_parse_intact(values, index_map))
    assert result["feature_type"] == "mutation"
    assert result["pmid"] is None


def test_popeve_matches_hand_written_parser():
    index_map = index_map_for(
        "popEVE_SCORE", "popEVE_EVE", "popEVE_ESM1v", "popEVE_gene",
        "popEVE_mutant", "popEVE_gap_frequency",
    )
    values = ["-0.5", "-1.2", "-3.4", "BRCA1", "K1R", "0.02"]
    result = run("popeve", values, index_map)
    assert result == dump(_parse_popeve(values, index_map))
    assert result["score"] == -0.5
    assert result["gene"] == "BRCA1"


def test_popeve_empty_is_none():
    index_map = index_map_for("popEVE_SCORE", "popEVE_EVE", "popEVE_mutant")
    assert run("popeve", ["", "", ""], index_map) is None


# --- pathogenicity, dissolved ------------------------------------------------
#
# _parse_pathogenicity groups several unrelated predictors into one object and
# nests spliceai/popeve (themselves standalone plugins) inside it. The grouping
# is not what the flat annotation payload wants, so the spec models each member
# as its own plugin. These tests prove the flat set reproduces the nested object
# field for field -- verified on real data too (revel 11,290 / alphamissense
# 62,384 / cadd 419,210 / eve 14,968 CSQ entries, zero mismatches).

PATH_COLS = ["REVEL", "am_class", "am_pathogenicity", "CADD_PHRED", "CADD_RAW",
             "EVE_CLASS", "EVE_SCORE"]
PATH_INDEX = index_map_for(*PATH_COLS)
PATH_VALUES = ["0.7", "likely_pathogenic", "0.9", "25.1", "3.2", "Pathogenic", "0.85"]


def test_flat_plugins_reproduce_the_nested_pathogenicity_object():
    nested = _parse_pathogenicity(PATH_VALUES, PATH_INDEX)

    assert run("revel", PATH_VALUES, PATH_INDEX)["score"] == nested.revel
    alphamissense = run("alphamissense", PATH_VALUES, PATH_INDEX)
    assert alphamissense["classification"] == nested.alphamissense_class
    assert alphamissense["score"] == nested.alphamissense_score
    cadd = run("cadd", PATH_VALUES, PATH_INDEX)
    assert (cadd["phred"], cadd["raw"]) == (nested.cadd_phred, nested.cadd_raw)
    eve = run("eve", PATH_VALUES, PATH_INDEX)
    assert (eve["classification"], eve["score"]) == (nested.eve_class, nested.eve_score)


def test_flat_pathogenicity_members_are_independent():
    """Only REVEL present: revel is an annotation, the others are absent.

    The nested object cannot express this -- it returns one object carrying a
    revel score and six nulls.
    """
    values = ["0.7", "", "", "", "", "", ""]
    assert run("revel", values, PATH_INDEX) == {"score": 0.7}
    assert run("alphamissense", values, PATH_INDEX) is None
    assert run("cadd", values, PATH_INDEX) is None
    assert run("eve", values, PATH_INDEX) is None


def test_flat_pathogenicity_members_absent_are_none():
    values = ["", "", "", "", "", "", ""]
    for plugin in ("revel", "alphamissense", "cadd", "eve"):
        assert run(plugin, values, PATH_INDEX) is None
    assert _parse_pathogenicity(values, PATH_INDEX) is None


def test_cadd_zero_is_kept():
    values = ["", "", "", "0.0", "0.0", "", ""]
    assert run("cadd", values, PATH_INDEX) == {"phred": 0.0, "raw": 0.0}
