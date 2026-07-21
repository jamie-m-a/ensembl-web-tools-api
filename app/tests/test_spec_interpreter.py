"""Tests for the spec-driven parsing interpreter (`apply_plugin_spec`) against
the parsing half of `specs/human_grch38.json`.

These began life as differential tests: every case compared the interpreter's
output to the corresponding hand-written `_parse_*` function, which was the
oracle that proved the spec vocabulary sufficient before anything was rewired.
The go-flat cutover deleted that bank, so the expected values below are the
frozen outputs from the last run in which both paths agreed — the equivalence
proof is now a set of pinned literals rather than a live comparison.

The fixtures are deliberately shared with test_csq_parsers.
"""

import pytest

from app.tests.test_csq_parsers import EMPTY, INDEX_MAP, row_list
from app.vep.models.parsing_spec_model import ParsingSpec, TargetSpec
from app.vep.utils.csq import get_prediction_index_map
from app.vep.utils.spec_interpreter import apply_plugin_spec
from app.vep.utils.spec_loader import load_merged_spec

SPEC: ParsingSpec = load_merged_spec("human_grch38").parsing


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
        "hgvsg",
        "spdi",
        "loeuf",
        "phenotype_data",
        "dosage_sensitivity",
        "intact",
        "popeve",
        "revel",
        "alphamissense",
        "cadd",
        "eve",
        "utr_annotation",
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


def test_key_value_requires_both_delimiters():
    with pytest.raises(Exception):
        TargetSpec.model_validate(
            {"field": "x", "from": "col", "transform": "key_value", "pair_delimiter": ":"}
        )


# --- mutfunc: four scalars ---------------------------------------------------

MUTFUNC_SCORES = dict(
    mutfunc_motif="0.1", mutfunc_int="0.2", mutfunc_mod="0.3", mutfunc_exp="0.4"
)


def test_mutfunc_all_four_scores():
    assert run("mutfunc", row_list(**MUTFUNC_SCORES)) == {
        "linear_motifs": 0.1,
        "protein_interactions": 0.2,
        "protein_structure": 0.3,
        "protein_structure_experimental": 0.4,
    }


def test_mutfunc_empty_is_none():
    assert run("mutfunc", EMPTY) is None


def test_mutfunc_partial_keeps_absent_scores_as_null():
    """Only some scores present: the rest must come back None, not be dropped."""
    csq = row_list(mutfunc_motif="0.1", mutfunc_exp="0.4")
    assert run("mutfunc", csq) == {
        "linear_motifs": 0.1,
        "protein_interactions": None,
        "protein_structure": None,
        "protein_structure_experimental": 0.4,
    }


# --- MaveDB: positional zip, the hard case -----------------------------------

MAVEDB_MULTI = dict(
    MaveDB_score="1.5&2.5&NA",
    MaveDB_urn="urn:1&urn:2&urn:3",
    MaveDB_doi="10.1/a&NA&10.1/c",
    MaveDB_nt="c.1A>G&NA",
    MaveDB_pro="p.Lys1Arg&NA",
)


def test_mavedb_multi_assay_shape_is_as_expected():
    result = run("mavedb", row_list(**MAVEDB_MULTI))
    assert result["protein_variant"] == "p.Lys1Arg"
    assert [(a["urn"], a["score"]) for a in result["assays"]] == [
        ("urn:1", 1.5),
        ("urn:2", 2.5),
        ("urn:3", None),  # NA score, but a real urn -> assay kept
    ]


def test_mavedb_empty_is_none():
    assert run("mavedb", EMPTY) is None


def test_mavedb_uneven_columns_pad_rather_than_truncate():
    """Fewer scores than urns: `align: max` must pad, not truncate."""
    csq = row_list(MaveDB_score="1.5", MaveDB_urn="urn:1&urn:2")
    assert run("mavedb", csq) == {
        "protein_variant": None,
        "assays": [
            {"score": 1.5, "urn": "urn:1"},
            {"score": None, "urn": "urn:2"},
        ],
    }


def test_mavedb_protein_variant_only_is_none():
    """pro present but no score/urn -> no assays -> whole annotation is None
    (require_any_output)."""
    assert run("mavedb", row_list(MaveDB_pro="p.Lys1Arg")) is None


def test_mavedb_all_na_assay_dropped():
    """A position where both score and urn are NA is dropped entirely."""
    csq = row_list(MaveDB_score="1.5&NA", MaveDB_urn="urn:1&NA")
    assert run("mavedb", csq) == {
        "protein_variant": None,
        "assays": [{"score": 1.5, "urn": "urn:1"}],
    }


# --- ClinVar: the `when` conditional -----------------------------------------

CONFLICTING = "Conflicting_classifications_of_pathogenicity"


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
    assert run("clinvar", csq) == {
        "significance": ["Pathogenic"],
        "conflicting_breakdown": [],
    }


def test_clinvar_when_matches_list_membership_not_substring():
    """A value that merely embeds the conflicting term must not trigger the
    breakdown — the condition is membership of the '&'-split list."""
    csq = row_list(
        ClinVar_CLNSIG="Not_" + CONFLICTING,
        ClinVar_CLNSIGCONF="Benign_(2)",
    )
    assert run("clinvar", csq) == {
        "significance": ["Not_" + CONFLICTING],
        "conflicting_breakdown": [],
    }


def test_clinvar_unparseable_breakdown_token_skipped():
    csq = row_list(
        ClinVar_CLNSIG=CONFLICTING,
        ClinVar_CLNSIGCONF="Benign_(2)&garbage_no_count",
    )
    result = run("clinvar", csq)
    assert [b["significance"] for b in result["conflicting_breakdown"]] == ["Benign"]


def test_clinvar_empty_is_none():
    assert run("clinvar", EMPTY) is None


# --- gnomAD / All of Us: pattern_map -----------------------------------------


def test_gnomad_exomes_pattern_map():
    index_map = index_map_for(
        "gnomAD_exomes_AF", "gnomAD_exomes_AF_afr", "gnomAD_exomes_AF_nfe_XX"
    )
    result = run("gnomad_exomes", ["0.01", "0.02", "0.03"], index_map)
    assert result["overall"] == 0.01
    # ancestry columns discovered from the header, not named in the spec
    assert result["populations"] == {"afr": 0.02, "nfe_XX": 0.03}


def test_gnomad_exomes_zero_overall_is_kept():
    """A 0.0 frequency is a real value. require_any_output must not treat it as
    absent (plain truthiness would drop the whole annotation)."""
    index_map = index_map_for("gnomAD_exomes_AF")
    result = run("gnomad_exomes", ["0.0"], index_map)
    assert result == {"overall": 0.0, "populations": {}}


def test_gnomad_exomes_absent_matches():
    index_map = index_map_for("Allele")
    assert run("gnomad_exomes", ["A"], index_map) is None


def test_gnomad_exomes_legacy_prefix_ignored():
    """The old gnomADe_ prefix must not match the pattern."""
    index_map = index_map_for("gnomADe_AF", "gnomADe_afr_AF")
    assert run("gnomad_exomes", ["0.1", "0.2"], index_map) is None


def test_gnomad_genomes_pattern_map():
    index_map = index_map_for(
        "gnomAD_genomes_AF", "gnomAD_genomes_AF_ami", "gnomAD_genomes_AF_grpmax"
    )
    result = run("gnomad_genomes", ["0.10", "0.20", "0.30"], index_map)
    assert result["overall"] == 0.1
    assert result["populations"] == {"ami": 0.20, "grpmax": 0.30}


def test_all_of_us_pattern_map_with_suffix():
    """AoU's pattern has a suffix (AoU_gvs_{pop}_af), unlike gnomAD's, plus a
    label column (AoU_gvs_max_subpop) naming which subpopulation the max
    frequency came from."""
    index_map = index_map_for(
        "AoU_gvs_all_af", "AoU_gvs_afr_af", "AoU_gvs_max_af", "AoU_gvs_max_subpop"
    )
    result = run("all_of_us", ["0.10", "0.20", "0.30", "eur"], index_map)
    assert result["overall"] == 0.10
    assert result["populations"] == {"afr": 0.20, "max": 0.30}
    assert result["max_subpopulation"] == "eur"
    # the label column is not a frequency and must not appear among populations
    assert "max_subpop" not in result["populations"]


def test_all_of_us_label_without_frequencies_is_none():
    """A max_subpop label with no frequencies is not an annotation — which is
    why max_subpopulation is deliberately absent from require_any_output."""
    index_map = index_map_for("AoU_gvs_all_af", "AoU_gvs_max_subpop")
    assert run("all_of_us", ["", "eur"], index_map) is None


# --- ProtVar: chunk + positional ---------------------------------------------

PROTVAR_FULL = dict(
    ProtVar_stability="0.42",
    ProtVar_pocket="POCKET1&-5.2&0.3&0.8&0.6&12.5&RES",
    ProtVar_int="PARTNER1&0.9&PARTNER2&0.8",
)


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


def test_protvar_odd_interaction_token_count():
    """A trailing partner with no score: still one interface, score null."""
    result = run("protvar", row_list(ProtVar_int="PARTNER1&0.9&PARTNER3"))
    assert result["interaction_interfaces"] == [
        {"partner": "PARTNER1", "score": 0.9, "raw": "PARTNER1&0.9"},
        {"partner": "PARTNER3", "score": None, "raw": "PARTNER3"},
    ]


def test_protvar_empty_is_none():
    assert run("protvar", EMPTY) is None


def test_protvar_pocket_missing_middle_value_does_not_shift():
    """An unparseable score empties only its own field: `positional` assigns
    strictly by index, so a bad item cannot pull the later values forward and
    have them silently reported under the wrong names."""
    raw = "POCKET1&-5.2&NA&0.8&0.6&12.5&RES"
    pocket = run("protvar", row_list(ProtVar_pocket=raw))["pockets"][0]

    assert pocket["energy"] == -5.2
    assert pocket["energy_per_volume"] is None
    assert pocket["score"] == 0.8
    assert pocket["buriedness"] == 0.6
    assert pocket["radius_of_gyration"] == 12.5


def test_protvar_interaction_na_partner_is_nulled():
    """'NA' means absent everywhere in the spec, including as a partner id."""
    interfaces = run("protvar", row_list(ProtVar_int="NA&0.9"))["interaction_interfaces"]
    assert interfaces[0]["partner"] is None


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


def test_opentargets_sorts_strongest_first():
    result = run_ot(
        diseases="EFO_1&EFO_2&EFO_3",
        genes="ENSG1&ENSG2&ENSG3",
        l2g="0.1&0.9&NA",
    )
    # descending by score; the unscored association goes last
    assert [(a["disease"], a["l2g_score"]) for a in result["gwas_associations"]] == [
        ("EFO_2", 0.9),
        ("EFO_1", 0.1),
        ("EFO_3", None),
    ]


def test_opentargets_dedups_repeated_rows():
    """The plugin emits duplicate rows -- dedup fires on 93% of real records."""
    result = run_ot(diseases="EFO_1&EFO_1", genes="ENSG1&ENSG1", l2g="0.5&0.5")
    assert result["gwas_associations"] == [
        {"disease": "EFO_1", "gene_id": "ENSG1", "l2g_score": 0.5}
    ]


def test_opentargets_drops_row_without_disease():
    result = run_ot(diseases="NA&EFO_2", genes="ENSG1&ENSG2", l2g="0.1&0.9")
    assert [a["disease"] for a in result["gwas_associations"]] == ["EFO_2"]


def test_opentargets_misaligned_columns_truncate():
    """align:min. Real data contains ragged columns (3 diseases, 2 genes), so
    the plugin's positional alignment is not guaranteed; zip drops the excess --
    with ragged input the true pairing is unknowable.
    """
    result = run_ot(diseases="EFO_1&EFO_2&EFO_3", genes="ENSG1&ENSG2", l2g="0.1&0.9")
    assert len(result["gwas_associations"]) == 2  # EFO_3 dropped


def test_opentargets_qtl_dedups_and_nulls_na_biosample():
    result = run_ot(qtl_genes="ENSG1&ENSG1&ENSG2", qtl_biosamples="liver&liver&NA")
    assert result["qtl_associations"] == [
        {"gene_id": "ENSG1", "biosample": "liver"},
        {"gene_id": "ENSG2", "biosample": None},
    ]


def test_opentargets_empty_is_none():
    assert run_ot() is None


def test_opentargets_absent_columns_is_none():
    assert run("opentargets", ["A"], index_map_for("Allele")) is None


# --- GO: regex + replace/strip -----------------------------------------------

GO_INDEX = index_map_for("GO")


def test_go_terms_split_id_from_name():
    values = ["GO:0001558:regulation_of_cell_growth&GO:0005509:calcium_ion_binding"]
    assert run("go", values, GO_INDEX)["go_terms"] == [
        {"id": "GO:0001558", "name": "regulation of cell growth"},
        {"id": "GO:0005509", "name": "calcium ion binding"},
    ]


def test_go_entry_without_a_term_name_is_null_not_empty_string():
    """Real data carries ids with no name at all (38 of 368 distinct GO ids in
    dev-data/output.vcf.gz, e.g. "GO:0050911:"). An absent name is null, as
    everywhere else in the spec."""
    assert run("go", ["GO:0050911:"], GO_INDEX)["go_terms"] == [
        {"id": "GO:0050911", "name": None}
    ]


def test_go_entry_without_a_name_part_is_skipped():
    """Fewer than three ':'-parts is not a term."""
    assert run("go", ["GO:0001558"], GO_INDEX) is None


def test_go_absent_is_none():
    assert run("go", [""], GO_INDEX) is None


# --- SpliceAI: mixed float/int scalars ---------------------------------------

SPLICEAI_COLS = [
    "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL",
]
SPLICEAI_INDEX = index_map_for(*SPLICEAI_COLS)


def test_spliceai_all_scores():
    values = ["BRCA1", "0.01", "0.02", "0.03", "0.04", "-5", "10", "-20", "30"]
    assert run("spliceai", values, SPLICEAI_INDEX) == {
        "symbol": "BRCA1",
        "ds_acceptor_gain": 0.01,
        "ds_acceptor_loss": 0.02,
        "ds_donor_gain": 0.03,
        "ds_donor_loss": 0.04,
        "dp_acceptor_gain": -5,
        "dp_acceptor_loss": 10,
        "dp_donor_gain": -20,
        "dp_donor_loss": 30,
    }


def test_spliceai_zero_scores_are_kept():
    """0.00 is the commonest real value (170k+ entries): a real score, not
    absence. require_any_output must not discard it."""
    values = ["BRCA1", "0.00", "", "", "", "", "", "", ""]
    result = run("spliceai", values, SPLICEAI_INDEX)
    assert result is not None
    assert result["ds_acceptor_gain"] == 0.0
    assert result["ds_acceptor_loss"] is None


def test_spliceai_symbol_alone_is_not_an_annotation():
    values = ["BRCA1", "", "", "", "", "", "", "", ""]
    assert run("spliceai", values, SPLICEAI_INDEX) is None


# --- plain scalar/list plugins -----------------------------------------------
#
# No new vocabulary; all five verified against dev-data/output.vcf.gz with zero
# mismatches (hgvs 210,658 / phenotype_data 382,715 / dosage 393,079 /
# intact 25 / popeve 96,953 CSQ entries).


def test_hgvs_three_notations():
    index_map = index_map_for("HGVSg", "HGVSc", "HGVSp")
    values = ["NC_1:g.100A>G", "ENST1:c.50A>G", "ENSP1:p.Lys1Arg"]
    assert run("hgvs", values, index_map) == {
        "genomic": "NC_1:g.100A>G",
        "transcript": "ENST1:c.50A>G",
        "protein": "ENSP1:p.Lys1Arg",
    }


def test_hgvs_partial():
    """This run emits HGVSc/HGVSp but no HGVSg -- the absent one is null."""
    index_map = index_map_for("HGVSc", "HGVSp")
    assert run("hgvs", ["ENST1:c.50A>G", ""], index_map) == {
        "genomic": None,
        "transcript": "ENST1:c.50A>G",
        "protein": None,
    }


def test_hgvs_empty_is_none():
    index_map = index_map_for("HGVSg", "HGVSc", "HGVSp")
    assert run("hgvs", ["", "", ""], index_map) is None


def test_phenotype_data_splits_and_drops_na():
    # The Phenotypes plugin produces a single PHENOTYPES column; CLIN_SIG / PUBMED
    # are co-located-variant fields from unrelated options, not Phenotypes output.
    index_map = index_map_for("PHENOTYPES")
    assert run("phenotype_data", ["cancer&NA&diabetes"], index_map) == {
        "phenotypes": ["cancer", "diabetes"],  # NA dropped
    }


def test_phenotype_data_empty_is_none():
    index_map = index_map_for("PHENOTYPES")
    assert run("phenotype_data", [""], index_map) is None


def test_dosage_sensitivity_probabilities():
    index_map = index_map_for("pHaplo", "pTriplo")
    assert run("dosage_sensitivity", ["0.98", "0.12"], index_map) == {
        "phaplo": 0.98, "ptriplo": 0.12
    }


def test_dosage_sensitivity_zero_is_kept():
    """0.0 is a real probability, not absence."""
    index_map = index_map_for("pHaplo", "pTriplo")
    assert run("dosage_sensitivity", ["0.0", ""], index_map) == {
        "phaplo": 0.0, "ptriplo": None
    }


def test_intact_unselected_sub_options_are_null():
    """This run emits only three of IntAct's columns; the unselected
    sub-options are absent and come back null."""
    index_map = index_map_for(
        "IntAct_feature_type", "IntAct_interaction_ac", "IntAct_feature_ac"
    )
    assert run("intact", ["mutation", "EBI-123", "EBI-ac"], index_map) == {
        "feature_type": "mutation",
        "interaction_ac": "EBI-123",
        "feature_ac": "EBI-ac",
        "feature_short_label": None,
        "feature_annotation": None,
        "ap_ac": None,
        "interaction_participants": None,
        "pmid": None,
    }


def test_popeve_scores():
    index_map = index_map_for(
        "popEVE_SCORE", "popEVE_EVE", "popEVE_ESM1v", "popEVE_gene",
        "popEVE_mutant", "popEVE_gap_frequency",
    )
    values = ["-0.5", "-1.2", "-3.4", "BRCA1", "K1R", "0.02"]
    assert run("popeve", values, index_map) == {
        "score": -0.5,
        "eve": -1.2,
        "esm1v": -3.4,
        "pop_adjusted_eve": None,
        "pop_adjusted_esm1v": None,
        "gene": "BRCA1",
        "protein": None,
        "mutant": "K1R",
        "gap_frequency": 0.02,
    }


def test_popeve_empty_is_none():
    index_map = index_map_for("popEVE_SCORE", "popEVE_EVE", "popEVE_mutant")
    assert run("popeve", ["", "", ""], index_map) is None


# --- pathogenicity, dissolved ------------------------------------------------
#
# The deleted `_parse_pathogenicity` grouped several unrelated predictors into
# one object and nested spliceai/popeve (themselves standalone plugins) inside
# it. The grouping is not what the flat annotation payload wants, so the spec
# models each member as its own plugin. Verified on real data at cutover time
# (revel 11,290 / alphamissense 62,384 / cadd 419,210 / eve 14,968 CSQ entries,
# zero mismatches against the nested object).

PATH_COLS = ["REVEL", "am_class", "am_pathogenicity", "CADD_PHRED", "CADD_RAW",
             "EVE_CLASS", "EVE_SCORE"]
PATH_INDEX = index_map_for(*PATH_COLS)
PATH_VALUES = ["0.7", "likely_pathogenic", "0.9", "25.1", "3.2", "Pathogenic", "0.85"]


def test_flat_plugins_carry_the_former_pathogenicity_members():
    assert run("revel", PATH_VALUES, PATH_INDEX) == {"score": 0.7}
    assert run("alphamissense", PATH_VALUES, PATH_INDEX) == {
        "classification": "likely_pathogenic", "score": 0.9
    }
    assert run("cadd", PATH_VALUES, PATH_INDEX) == {"phred": 25.1, "raw": 3.2}
    assert run("eve", PATH_VALUES, PATH_INDEX) == {
        "classification": "Pathogenic", "score": 0.85
    }


def test_flat_pathogenicity_members_are_independent():
    """Only REVEL present: revel is an annotation, the others are absent.

    The old nested object could not express this -- it returned one object
    carrying a revel score and six nulls.
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


def test_cadd_zero_is_kept():
    values = ["", "", "", "0.0", "0.0", "", ""]
    assert run("cadd", values, PATH_INDEX) == {"phred": 0.0, "raw": 0.0}


# --- the plugins added by the go-flat cutover --------------------------------


def test_loeuf_is_a_transcript_scoped_score():
    spec = SPEC.plugin("loeuf")
    assert spec.scope == "transcript"
    index_map = index_map_for("LOEUF")
    assert run("loeuf", ["0.15"], index_map) == {"score": 0.15}
    assert run("loeuf", [""], index_map) is None


def test_spdi_and_hgvsg_are_allele_scoped():
    """Both are allele-scoped because intergenic variants have no transcript
    rows and must still carry their variant representations."""
    assert SPEC.plugin("spdi").scope == "allele"
    assert SPEC.plugin("hgvsg").scope == "allele"

    spdi_map = index_map_for("SPDI")
    assert run("spdi", ["1:79106:T:C"], spdi_map) == {"spdi": "1:79106:T:C"}
    assert run("spdi", [""], spdi_map) is None

    hgvsg_map = index_map_for("HGVSg")
    assert run("hgvsg", ["1:g.79107T>C"], hgvsg_map) == {"genomic": "1:g.79107T>C"}
    assert run("hgvsg", [""], hgvsg_map) is None


# --- UTRAnnotator ------------------------------------------------------------

UTR_COLS = ["5UTR_consequence", "5UTR_annotation", "Existing_uORFs",
            "Existing_InFrame_oORFs", "Existing_OutOfFrame_oORFs"]
UTR_INDEX = index_map_for(*UTR_COLS)
# A real value from dev-data/has_utr.vcf.gz.
UTR_ANNOTATION = (
    "alt_type=uORF:ref_StartDistanceToCDS=324:ref_type=uORF:KozakStrength=Moderate"
    ":KozakContext=GCGATGC:ref_type_length=15:Evidence=False:alt_type_length=189"
)
UTR_VALUES = ["5_prime_UTR_uORF_frameshift_variant", UTR_ANNOTATION, "5", "0", "0"]


def test_utr_annotation_parses_the_detail_string_into_a_dict():
    """The deleted parser copied `annotation` verbatim as a ':'-delimited
    string; the spec parses it into a dict via `key_value`, which is the point
    of the transform (see the ordering test below)."""
    result = run("utr_annotation", UTR_VALUES, UTR_INDEX)
    assert result["consequence"] == "5_prime_UTR_uORF_frameshift_variant"
    assert result["existing_uorfs"] == "5"
    assert result["existing_inframe_oorfs"] == "0"
    assert result["existing_outofframe_oorfs"] == "0"
    assert result["annotation"] == {
        "alt_type": "uORF",
        "ref_StartDistanceToCDS": "324",
        "ref_type": "uORF",
        "KozakStrength": "Moderate",
        "KozakContext": "GCGATGC",
        "ref_type_length": "15",
        "Evidence": "False",
        "alt_type_length": "189",
    }


def test_utr_annotation_key_value_is_order_independent():
    """The actual bug this fixes: UTRAnnotator emits the same pairs in a
    different order every record (all 9 in has_utr.vcf.gz are one identical
    annotation shuffled 9 ways). The raw string is 9 different values; parsed
    as key_value, all 9 must be the same dict."""
    shuffled = "ref_type=uORF:alt_type=uORF:Evidence=False:KozakStrength=Moderate"
    original = "alt_type=uORF:KozakStrength=Moderate:ref_type=uORF:Evidence=False"
    assert shuffled != original  # different raw strings

    def annotation_of(value):
        row = ["", value, "", "", ""]
        return run("utr_annotation", row, UTR_INDEX)["annotation"]

    assert annotation_of(shuffled) == annotation_of(original)  # same parsed value


def test_utr_annotation_malformed_piece_is_dropped_not_raised():
    """A piece without '=' does not break parsing of the rest of the value."""
    row = ["", "alt_type=uORF:garbage:Evidence=False", "", "", ""]
    result = run("utr_annotation", row, UTR_INDEX)
    assert result["annotation"] == {"alt_type": "uORF", "Evidence": "False"}


def test_utr_annotation_empty_is_none():
    assert run("utr_annotation", ["", "", "", "", ""], UTR_INDEX) is None
