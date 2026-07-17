"""Tests for the CSQ field parsers (app/vep/utils/vcf_results._parse_*).

The existing test_vep.py exercises `_get_alt_allele_details` against an *old*
CSQ header (SIFT/PolyPhen/AF era) with none of the modern plugin columns, so the
new plugin parsers are otherwise unexercised. Here a modern CSQ header (with the
plugin columns) is used to build an index_map; each parser is called with one
populated row and one empty row, and there is an end-to-end allele test for a
transcript row and an intergenic row.
"""

from app.vep.models import vcf_results_model as model
from app.vep.utils.csq import get_prediction_index_map
from app.vep.utils.vcf_results import (
    _get_alt_allele_details,
    _parse_frequencies,
    _parse_protvar,
    _parse_protvar_pocket,
    _parse_intact,
    _parse_mutfunc,
    _parse_mavedb,
    _parse_popeve,
    _parse_dosage_sensitivity,
    _parse_utr_annotation,
    _parse_riboseq_orfs,
    _parse_spliceai,
    _parse_pathogenicity,
    _parse_clinvar,
)

# A modern CSQ header: the columns the new plugin parsers read.
ALL_COLS = [
    # core / allele-level
    "Allele", "AF", "Consequence", "Feature", "Feature_type", "BIOTYPE",
    "CANONICAL", "SYMBOL", "Gene", "STRAND", "ENSP", "Existing_variation",
    "MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL",
    "SPDI", "HGVSg", "HGVSc", "HGVSp", "CADD_PHRED", "CADD_RAW", "LOEUF",
    # ProtVar
    "ProtVar_stability", "ProtVar_pocket", "ProtVar_int",
    # IntAct (base + six sub-option columns)
    "IntAct_feature_type", "IntAct_interaction_ac", "IntAct_feature_ac",
    "IntAct_feature_short_label", "IntAct_feature_annotation", "IntAct_ap_ac",
    "IntAct_interaction_participants", "IntAct_pmid",
    # mutfunc
    "mutfunc_motif", "mutfunc_int", "mutfunc_mod", "mutfunc_exp",
    # MaveDB
    "MaveDB_score", "MaveDB_urn", "MaveDB_doi", "MaveDB_nt", "MaveDB_pro",
    # popEVE
    "popEVE_SCORE", "popEVE_EVE", "popEVE_ESM1v", "popEVE_pop_adjusted_EVE",
    "popEVE_pop_adjusted_ESM1v", "popEVE_gene", "popEVE_protein",
    "popEVE_mutant", "popEVE_gap_frequency",
    # dosage sensitivity
    "pHaplo", "pTriplo",
    # UTRAnnotator
    "5UTR_consequence", "5UTR_annotation", "Existing_uORFs",
    "Existing_InFrame_oORFs", "Existing_OutOfFrame_oORFs",
    # Ribo-seq ORFs
    "RiboseqORFs_id", "RiboseqORFs_consequences", "RiboseqORFs_impact",
    "RiboseqORFs_protein_position", "RiboseqORFs_codons",
    "RiboseqORFs_amino_acids", "RiboseqORFs_publications",
    # SpliceAI
    "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL",
    # pathogenicity
    "REVEL", "am_class", "am_pathogenicity", "EVE_CLASS", "EVE_SCORE",
    # ClinVar
    "ClinVar_CLNSIG", "ClinVar_CLNSIGCONF",
]

HEADER = "Consequence annotations from Ensembl VEP. Format: " + "|".join(ALL_COLS)
INDEX_MAP = get_prediction_index_map(HEADER)


def row_list(**values):
    """A CSQ values list (aligned to ALL_COLS); unset columns are empty."""
    return [str(values.get(col, "")) for col in ALL_COLS]


def row_str(**values):
    """A pipe-joined CSQ string (as it appears in the VCF INFO field)."""
    return "|".join(row_list(**values))


EMPTY = row_list()


# --- ProtVar -----------------------------------------------------------------


def test_parse_protvar_pocket_positional_and_no_residues():
    pocket = _parse_protvar_pocket("POCKET1&-5.2&0.3&0.8&0.6&12.5&RES:12,13")
    assert pocket.pocket_id == "POCKET1"
    assert pocket.energy == -5.2
    assert pocket.energy_per_volume == 0.3
    assert pocket.score == 0.8
    assert pocket.buriedness == 0.6
    assert pocket.radius_of_gyration == 12.5
    # residues are intentionally not captured
    assert not hasattr(pocket, "residues")


def test_parse_protvar_pocket_unparseable_item_does_not_shift_later_values():
    """An unparseable score must empty only its own field.

    This previously collected the parts that parsed as a float and assigned them
    in order, so a gap pulled every later value one field forward — the score
    was reported as energy_per_volume, buriedness as score, and so on.
    """
    pocket = _parse_protvar_pocket("POCKET1&-5.2&NA&0.8&0.6&12.5&RES")
    assert pocket.energy == -5.2
    assert pocket.energy_per_volume is None  # the unparseable one, and only it
    assert pocket.score == 0.8
    assert pocket.buriedness == 0.6
    assert pocket.radius_of_gyration == 12.5


def test_parse_protvar_pocket_truncated_value():
    """Fewer items than fields: the missing tail is empty, not misread."""
    pocket = _parse_protvar_pocket("POCKET1&-5.2&0.3")
    assert pocket.pocket_id == "POCKET1"
    assert pocket.energy == -5.2
    assert pocket.energy_per_volume == 0.3
    assert pocket.score is None
    assert pocket.buriedness is None
    assert pocket.radius_of_gyration is None
    assert pocket.raw == "POCKET1&-5.2&0.3"


def test_parse_protvar_populated():
    result = _parse_protvar(
        row_list(
            ProtVar_stability="0.42",
            ProtVar_pocket="POCKET1&-5.2&0.3&0.8&0.6&12.5&RES",
            ProtVar_int="PARTNER1&0.9&PARTNER2&0.8",
        ),
        INDEX_MAP,
    )
    assert result.structure_stability_score == 0.42
    assert len(result.pockets) == 1
    assert result.pockets[0].pocket_id == "POCKET1"
    assert [i.partner for i in result.interaction_interfaces] == [
        "PARTNER1",
        "PARTNER2",
    ]
    assert result.interaction_interfaces[0].score == 0.9


def test_parse_protvar_empty_is_none():
    assert _parse_protvar(EMPTY, INDEX_MAP) is None


# --- IntAct ------------------------------------------------------------------


def test_parse_intact_base_and_sub_options():
    result = _parse_intact(
        row_list(
            IntAct_feature_type="mutation",
            IntAct_interaction_ac="EBI-123",
            IntAct_feature_ac="EBI-ac",
            IntAct_feature_short_label="short",
            IntAct_feature_annotation="annot",
            IntAct_ap_ac="EBI-ap",
            IntAct_interaction_participants="2",
            IntAct_pmid="123456",
        ),
        INDEX_MAP,
    )
    assert result.feature_type == "mutation"
    assert result.interaction_ac == "EBI-123"
    assert result.feature_ac == "EBI-ac"
    assert result.feature_short_label == "short"
    assert result.feature_annotation == "annot"
    assert result.ap_ac == "EBI-ap"
    assert result.interaction_participants == "2"
    assert result.pmid == "123456"


def test_parse_intact_empty_is_none():
    assert _parse_intact(EMPTY, INDEX_MAP) is None


# --- mutfunc -----------------------------------------------------------------


def test_parse_mutfunc_scores():
    result = _parse_mutfunc(
        row_list(
            mutfunc_motif="0.1",
            mutfunc_int="0.2",
            mutfunc_mod="0.3",
            mutfunc_exp="0.4",
        ),
        INDEX_MAP,
    )
    assert result.linear_motifs == 0.1
    assert result.protein_interactions == 0.2
    assert result.protein_structure == 0.3
    assert result.protein_structure_experimental == 0.4


def test_parse_mutfunc_empty_is_none():
    assert _parse_mutfunc(EMPTY, INDEX_MAP) is None


# --- MaveDB (multi-assay) ----------------------------------------------------


def test_parse_mavedb_multi_assay_pairs_urn_and_score():
    result = _parse_mavedb(
        row_list(
            MaveDB_score="1.5&2.5&NA",
            MaveDB_urn="urn:1&urn:2&urn:3",
            MaveDB_doi="10.1/a&NA&10.1/c",
            MaveDB_nt="c.1A>G&NA",
            MaveDB_pro="p.Lys1Arg&NA",
        ),
        INDEX_MAP,
    )
    # Three assays, each pairing its urn with its (positional) score. The third
    # score is NA -> None, but its urn is present so the assay is kept.
    assert [(a.urn, a.score) for a in result.assays] == [
        ("urn:1", 1.5),
        ("urn:2", 2.5),
        ("urn:3", None),
    ]
    assert result.protein_variant == "p.Lys1Arg"


def test_parse_mavedb_empty_is_none():
    assert _parse_mavedb(EMPTY, INDEX_MAP) is None


# --- popEVE ------------------------------------------------------------------


def test_parse_popeve_populated():
    result = _parse_popeve(
        row_list(
            popEVE_SCORE="0.91",
            popEVE_EVE="0.8",
            popEVE_ESM1v="-0.5",
            popEVE_pop_adjusted_EVE="0.7",
            popEVE_pop_adjusted_ESM1v="-0.4",
            popEVE_gene="TP53",
            popEVE_protein="ENSP1",
            popEVE_mutant="K41R",
            popEVE_gap_frequency="0.02",
        ),
        INDEX_MAP,
    )
    assert result.score == 0.91
    assert result.eve == 0.8
    assert result.esm1v == -0.5
    assert result.pop_adjusted_eve == 0.7
    assert result.pop_adjusted_esm1v == -0.4
    assert result.gene == "TP53"
    assert result.protein == "ENSP1"
    assert result.mutant == "K41R"
    assert result.gap_frequency == 0.02


def test_parse_popeve_empty_is_none():
    assert _parse_popeve(EMPTY, INDEX_MAP) is None


# --- dosage sensitivity ------------------------------------------------------


def test_parse_dosage_sensitivity_populated():
    result = _parse_dosage_sensitivity(
        row_list(pHaplo="0.95", pTriplo="0.12"), INDEX_MAP
    )
    assert result.phaplo == 0.95
    assert result.ptriplo == 0.12


def test_parse_dosage_sensitivity_empty_is_none():
    assert _parse_dosage_sensitivity(EMPTY, INDEX_MAP) is None


# --- UTRAnnotator ------------------------------------------------------------


def test_parse_utr_annotation_populated():
    result = _parse_utr_annotation(
        row_list(
            **{
                "5UTR_consequence": "5_prime_UTR_premature_start_codon_gain_variant",
                "5UTR_annotation": "uORF",
                "Existing_uORFs": "2",
                "Existing_InFrame_oORFs": "1",
                "Existing_OutOfFrame_oORFs": "0",
            }
        ),
        INDEX_MAP,
    )
    assert result.consequence == "5_prime_UTR_premature_start_codon_gain_variant"
    assert result.annotation == "uORF"
    assert result.existing_uorfs == "2"
    assert result.existing_inframe_oorfs == "1"
    assert result.existing_outofframe_oorfs == "0"


def test_parse_utr_annotation_empty_is_none():
    assert _parse_utr_annotation(EMPTY, INDEX_MAP) is None


# --- Ribo-seq ORFs -----------------------------------------------------------


def test_parse_riboseq_orfs_populated():
    result = _parse_riboseq_orfs(
        row_list(
            RiboseqORFs_id="c1orf1",
            RiboseqORFs_consequences="missense_variant&synonymous_variant",
            RiboseqORFs_impact="MODERATE",
        ),
        INDEX_MAP,
    )
    assert result.orf_id == "c1orf1"
    assert result.consequences == ["missense_variant", "synonymous_variant"]
    assert result.impact == "MODERATE"


def test_parse_riboseq_orfs_empty_is_none():
    assert _parse_riboseq_orfs(EMPTY, INDEX_MAP) is None


# --- SpliceAI ----------------------------------------------------------------


def test_parse_spliceai_all_fields():
    result = _parse_spliceai(
        row_list(
            SpliceAI_pred_SYMBOL="TP53",
            SpliceAI_pred_DS_AG="0.01",
            SpliceAI_pred_DS_AL="0.02",
            SpliceAI_pred_DS_DG="0.90",
            SpliceAI_pred_DS_DL="0.03",
            SpliceAI_pred_DP_AG="-5",
            SpliceAI_pred_DP_AL="10",
            SpliceAI_pred_DP_DG="2",
            SpliceAI_pred_DP_DL="-7",
        ),
        INDEX_MAP,
    )
    assert result.symbol == "TP53"
    assert result.ds_acceptor_gain == 0.01
    assert result.ds_acceptor_loss == 0.02
    assert result.ds_donor_gain == 0.90
    assert result.ds_donor_loss == 0.03
    assert result.dp_acceptor_gain == -5
    assert result.dp_acceptor_loss == 10
    assert result.dp_donor_gain == 2
    assert result.dp_donor_loss == -7


def test_parse_spliceai_symbol_only_is_none():
    # a symbol with no delta scores is not a real result
    assert _parse_spliceai(row_list(SpliceAI_pred_SYMBOL="TP53"), INDEX_MAP) is None


# --- pathogenicity (aggregate) -----------------------------------------------


def test_parse_pathogenicity_aggregates_nested_and_flat():
    result = _parse_pathogenicity(
        row_list(
            am_class="likely_pathogenic",
            am_pathogenicity="0.98",
            REVEL="0.75",
            CADD_PHRED="25.3",
            CADD_RAW="4.1",
            EVE_CLASS="Pathogenic",
            EVE_SCORE="0.88",
            SpliceAI_pred_DS_DG="0.90",  # nested spliceai
            popEVE_SCORE="0.91",  # nested popeve
        ),
        INDEX_MAP,
    )
    assert result.alphamissense_class == "likely_pathogenic"
    assert result.alphamissense_score == 0.98
    assert result.revel == 0.75
    assert result.cadd_phred == 25.3
    assert result.cadd_raw == 4.1
    assert result.eve_class == "Pathogenic"
    assert result.eve_score == 0.88
    assert result.spliceai is not None
    assert result.spliceai.ds_donor_gain == 0.90
    assert result.popeve is not None
    assert result.popeve.score == 0.91


def test_parse_pathogenicity_empty_is_none():
    assert _parse_pathogenicity(EMPTY, INDEX_MAP) is None


# --- ClinVar clinical significance -------------------------------------------


def test_parse_clinvar_non_conflicting_ignores_clnsigconf():
    result = _parse_clinvar(
        # CLNSIGCONF present but must be ignored when the class isn't conflicting
        row_list(ClinVar_CLNSIG="Pathogenic", ClinVar_CLNSIGCONF="Benign_(3)"),
        INDEX_MAP,
    )
    assert result.significance == ["Pathogenic"]
    assert result.conflicting_breakdown == []


def test_parse_clinvar_conflicting_extracts_clnsigconf_breakdown():
    result = _parse_clinvar(
        row_list(
            ClinVar_CLNSIG="Conflicting_classifications_of_pathogenicity",
            ClinVar_CLNSIGCONF=(
                "Pathogenic_(10)&Likely_pathogenic_(6)&Uncertain_significance_(2)"
            ),
        ),
        INDEX_MAP,
    )
    assert result.significance == ["Conflicting_classifications_of_pathogenicity"]
    assert [(s.significance, s.count) for s in result.conflicting_breakdown] == [
        ("Pathogenic", 10),
        ("Likely_pathogenic", 6),
        ("Uncertain_significance", 2),
    ]


def test_parse_clinvar_empty_is_none():
    assert _parse_clinvar(EMPTY, INDEX_MAP) is None


# --- allele frequencies (All of Us AoU_ prefix) ------------------------------


def test_parse_frequencies_allofus_uses_aou_prefix():
    # short_name=AoU means the custom columns come back prefixed AoU_gvs_*
    cols = [
        "AoU_gvs_all_af",
        "AoU_gvs_afr_af",
        "AoU_gvs_max_af",
        "AoU_gvs_max_subpop",
    ]
    index_map = get_prediction_index_map("Format: " + "|".join(cols))
    values = ["0.10", "0.20", "0.30", "eur"]

    result = _parse_frequencies(values, index_map)
    assert result is not None
    aou = result.all_of_us
    assert aou.overall == 0.10  # AoU_gvs_all_af
    assert aou.populations["afr"] == 0.20
    assert aou.populations["max"] == 0.30  # AoU_gvs_max_af
    # the label column AoU_gvs_max_subpop is not a frequency and is excluded
    assert "max_subpop" not in aou.populations


def test_parse_frequencies_ignores_legacy_allofus_prefix():
    # the old AllOfUs_ prefix must no longer be picked up
    cols = ["AllOfUs_gvs_all_af", "AllOfUs_gvs_afr_af"]
    index_map = get_prediction_index_map("Format: " + "|".join(cols))
    assert _parse_frequencies(["0.4", "0.5"], index_map) is None


def test_parse_frequencies_gnomad_exomes_uses_custom_prefix():
    # short_name=gnomAD_exomes -> gnomAD_exomes_AF (overall) + AF_<...> variants
    cols = [
        "gnomAD_exomes_AF",
        "gnomAD_exomes_AF_afr",
        "gnomAD_exomes_AF_nfe_XX",
    ]
    index_map = get_prediction_index_map("Format: " + "|".join(cols))
    result = _parse_frequencies(["0.01", "0.02", "0.03"], index_map)
    assert result is not None
    exomes = result.gnomad_exomes
    assert exomes.overall == 0.01
    assert exomes.populations["afr"] == 0.02
    assert exomes.populations["nfe_XX"] == 0.03  # sex-split variant captured


def test_parse_frequencies_gnomad_genomes_uses_custom_prefix():
    cols = [
        "gnomAD_genomes_AF",
        "gnomAD_genomes_AF_ami",
        "gnomAD_genomes_AF_grpmax",
    ]
    index_map = get_prediction_index_map("Format: " + "|".join(cols))
    result = _parse_frequencies(["0.10", "0.20", "0.30"], index_map)
    assert result is not None
    genomes = result.gnomad_genomes
    assert genomes.overall == 0.10
    assert genomes.populations["ami"] == 0.20
    assert genomes.populations["grpmax"] == 0.30


def test_parse_frequencies_ignores_legacy_gnomad_prefixes():
    # the old gnomADe_/gnomADg_ prefixes must no longer match
    cols = ["gnomADe_AF", "gnomADg_AF", "gnomADe_afr_AF"]
    index_map = get_prediction_index_map("Format: " + "|".join(cols))
    assert _parse_frequencies(["0.1", "0.2", "0.3"], index_map) is None


# --- end-to-end allele (modern header) ---------------------------------------


def test_get_alt_allele_details_transcript_row_populates_annotations():
    transcript = row_str(
        Allele="T",
        Consequence="missense_variant",
        Feature="ENST00000269305.9",
        Feature_type="Transcript",
        BIOTYPE="protein_coding",
        CANONICAL="YES",
        SYMBOL="TP53",
        Gene="ENSG00000141510",
        STRAND="1",
        ENSP="ENSP00000269305",
        HGVSc="c.123A>G",
        HGVSp="p.Lys41Arg",
        HGVSg="17:g.7676154A>G",
        SPDI="NC_000017.11:7676153:A:G",
        CADD_PHRED="25.3",
        LOEUF="0.15",
        ProtVar_stability="0.42",
        IntAct_feature_type="mutation",
        mutfunc_motif="0.1",
        MaveDB_score="1.5",
        am_class="likely_pathogenic",
        pHaplo="0.95",
        **{"5UTR_consequence": "uORF_variant"},
        RiboseqORFs_id="c1orf1",
    )

    result = _get_alt_allele_details("C", "T", [transcript], INDEX_MAP)
    assert len(result.predicted_molecular_consequences) == 1
    consequence = result.predicted_molecular_consequences[0]

    # compare by value: vcf_results imports the model as `vep.models...` while
    # this test imports `app.vep.models...`, so the enum *members* differ by
    # identity even though their values match (same reason test_vep.py skips it).
    assert consequence.feature_type.value == "transcript"
    assert consequence.protvar is not None
    assert consequence.intact is not None
    assert consequence.mutfunc is not None
    assert consequence.mavedb is not None
    assert consequence.hgvs is not None
    assert consequence.pathogenicity is not None
    assert consequence.loeuf == 0.15
    assert consequence.dosage_sensitivity is not None
    assert consequence.utr_annotation is not None
    assert consequence.riboseq_orfs is not None

    # allele-level fields captured on the allele
    assert result.spdi == "NC_000017.11:7676153:A:G"
    assert result.hgvsg == "17:g.7676154A>G"
    assert result.cadd_phred == 25.3


def test_get_alt_allele_details_intergenic_surfaces_allele_level_fields():
    intergenic = row_str(
        Allele="A",
        Consequence="intergenic_variant",
        Feature_type="",
        SPDI="NC_000017.11:7676153:T:A",
        HGVSg="17:g.7676154T>A",
        CADD_PHRED="8.2",
    )

    result = _get_alt_allele_details("T", "A", [intergenic], INDEX_MAP)
    assert len(result.predicted_molecular_consequences) == 1
    # no transcript consequence for an intergenic variant
    assert result.predicted_molecular_consequences[0].feature_type is None
    # ...but the allele-level representations / CADD are still surfaced
    assert result.spdi == "NC_000017.11:7676153:T:A"
    assert result.hgvsg == "17:g.7676154T>A"
    assert result.cadd_phred == 8.2
