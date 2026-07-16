""" Module for loading a VCF and parsing it into a VepResultsResponse
object as defined in APISpecification"""

from io import StringIO
import gzip
import json
import re
import struct
import subprocess
import zlib
from pathlib import Path
from typing import Iterator
from pydantic import FilePath
import vcfpy
from vep.models import vcf_results_model as model

TARGET_COLUMNS = [
    "Allele",
    "AF",
    "Consequence",
    "Feature",
    "Feature_type",
    "BIOTYPE",
    "CANONICAL",
    "SYMBOL",
    "Gene",
    "STRAND",
    "IMPACT",
    # MANE (human GRCh38)
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    # Protein & functional
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    "DOMAINS",
    "ProtVar_stability",
    "ProtVar_int",
    "ProtVar_pocket",
    "IntAct_feature_ac",
    "IntAct_feature_type",
    "IntAct_interaction_ac",
    "mutfunc_motif",
    "mutfunc_int",
    "mutfunc_mod",
    "mutfunc_exp",
    "MaveDB_score",
    "MaveDB_urn",
    "MaveDB_doi",
    "MaveDB_nt",
    "MaveDB_pro",
]

META_FILE = "results_meta.json"

# Taken from https://github.com/Ensembl/ensembl-hypsipyle
# main/common/file_model/variant.py#L142
# Needs to be moved into a shared module
def _set_allele_type(alt_one_bp: bool, ref_one_bp: bool, ref_alt_equal_bp: bool) -> tuple[str,str]:
    """Create a allele type for a variant based on Variation
    teams logic using ref and largest alt allele sizes"""
    match [alt_one_bp, ref_one_bp, ref_alt_equal_bp]:
        case [True, True, True]:
            allele_type = "SNV"
            so_term = "SO:0001483"

        case [True, False, False]:
            allele_type = "deletion"
            so_term = "SO:0000159"

        case [False, True, False]:
            allele_type = "insertion"
            so_term = "SO:0000667"

        case [False, False, False]:
            allele_type = "indel"
            so_term = "SO:1000032"

        case [False, False, True]:
            allele_type = "substitution"
            so_term = "SO:1000002"
    return allele_type, so_term

def _get_variant_type(ref: str, alt: str) -> str:
    """Helper function to infer variant type from allele values"""
    if alt=="copy_number_variation":
        return alt
    else:
        return _set_allele_type(len(alt) < 2, len(ref) < 2, len(alt) == len(ref))[0]


def _alt_value(alt) -> str:
    """Return an alt allele's sequence string.

    Simple substitution alts expose `.value`; symbolic and breakend alts
    (e.g. structural variants) do not, so fall back to their serialized VCF
    representation."""
    value = getattr(alt, "value", None)
    if value is not None:
        return value
    serialize = getattr(alt, "serialize", None)
    return serialize() if callable(serialize) else str(alt)


def _get_prediction_index_map(
    csq_header: str, target_columns: list[str] | None = None
) -> dict[str, int]:
    """Creates a dictionary of column indexes from the CSQ info description.

    By default every CSQ column is indexed (so any annotation field can be
    read); pass target_columns to restrict to an allow-list."""
    csq_header = csq_header.split(":")[-1].strip()
    csq_headers = csq_header.split("|")

    return {
        header: index
        for index, header in enumerate(csq_headers)
        if target_columns is None or header in target_columns
    }


def _get_csq_value(
    csq_values: list[str], csq_key: str, default_value: str | None, index_map: dict[str, int]
):
    """Helper method to return CSQ values or a default value
    if either the key or the value is missing"""
    if csq_key in index_map and csq_values[index_map[csq_key]]:
        return csq_values[index_map[csq_key]]
    return default_value


def _has_any_column(index_map: dict[str, int], *columns: str) -> bool:
    """Whether any of `columns` is present in the CSQ header at all. Lets a parser
    skip its work when the plugin that produces its columns wasn't run — the
    header is fixed for the whole file, so a column that is absent here is absent
    for every record (and the parser could only ever return None/empty)."""
    return any(column in index_map for column in columns)


def _to_float(value: str | None) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_uniprot(csq_values, index_map) -> model.UniprotIds | None:
    """Build Uniprot cross-references from the SWISSPROT/TREMBL/UNIPARC/isoform
    CSQ columns; returns None if none are present."""
    if not _has_any_column(
        index_map, "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"
    ):
        return None
    swissprot = _get_csq_value(csq_values, "SWISSPROT", None, index_map)
    trembl = _get_csq_value(csq_values, "TREMBL", None, index_map)
    uniparc = _get_csq_value(csq_values, "UNIPARC", None, index_map)
    isoform = _get_csq_value(csq_values, "UNIPROT_ISOFORM", None, index_map)
    if not any([swissprot, trembl, uniparc, isoform]):
        return None
    return model.UniprotIds(
        swissprot=swissprot, trembl=trembl, uniparc=uniparc, isoform=isoform
    )


def _parse_protein_matches(csq_values, index_map) -> list[model.ProteinMatch]:
    """Parse the DOMAINS CSQ column (e.g. AlphaFold-DB / PDB mappings).
    Multiple matches are '&'-joined; each is 'source:id'."""
    domains = _get_csq_value(csq_values, "DOMAINS", None, index_map)
    if not domains:
        return []
    matches = []
    for item in domains.split("&"):
        if not item:
            continue
        source, sep, identifier = item.partition(":")
        matches.append(
            model.ProteinMatch(
                source=source if sep else "",
                id=identifier if sep else source,
            )
        )
    return matches


def _parse_protvar_pocket(raw: str) -> model.ProtVarPocket:
    """Parse a ProtVar_pocket value, positionally:
    id & energy & energy_per_volume & score & buriedness & radius_of_gyration &
    residues. The leading id (and trailing residue token) are non-numeric; the
    middle parts are the numeric scores captured here. Residues are ignored."""
    parts = raw.split("&")
    pocket_id = parts[0] if parts else ""
    numeric: list[float | None] = []
    for part in parts[1:]:
        number = _to_float(part)
        if number is not None:
            numeric.append(number)
    numeric += [None] * (5 - len(numeric))
    return model.ProtVarPocket(
        pocket_id=pocket_id,
        energy=numeric[0],
        energy_per_volume=numeric[1],
        score=numeric[2],
        buriedness=numeric[3],
        radius_of_gyration=numeric[4],
        raw=raw,
    )


def _parse_protvar(csq_values, index_map) -> model.ProtVarAnnotation | None:
    """Build a ProtVar annotation from the stability/pocket/interaction CSQ
    columns; returns None if none are present."""
    if not _has_any_column(index_map, "ProtVar_stability", "ProtVar_pocket", "ProtVar_int"):
        return None
    stability = _get_csq_value(csq_values, "ProtVar_stability", None, index_map)
    pocket = _get_csq_value(csq_values, "ProtVar_pocket", None, index_map)
    interaction = _get_csq_value(csq_values, "ProtVar_int", None, index_map)
    if not any([stability, pocket, interaction]):
        return None

    pockets = [_parse_protvar_pocket(pocket)] if pocket else []

    interfaces = []
    if interaction:
        tokens = interaction.split("&")
        # partner & score [& partner & score ...]
        for i in range(0, len(tokens), 2):
            partner = tokens[i]
            score = tokens[i + 1] if i + 1 < len(tokens) else None
            interfaces.append(
                model.ProtVarInteractionInterface(
                    partner=partner,
                    score=_to_float(score),
                    raw="&".join(t for t in [partner, score] if t),
                )
            )

    return model.ProtVarAnnotation(
        structure_stability_score=_to_float(stability),
        pockets=pockets,
        interaction_interfaces=interfaces,
    )


def _parse_intact(csq_values, index_map) -> model.IntActAnnotation | None:
    """Build an IntAct annotation; returns None if no IntAct columns present.

    Besides the base feature_type / interaction_ac, the plugin emits one column
    per selected sub-option (IntAct_<flag>), so read those too."""
    if not _has_any_column(
        index_map, "IntAct_feature_type", "IntAct_interaction_ac", "IntAct_feature_ac"
    ):
        return None

    def col(name):
        return _get_csq_value(csq_values, name, None, index_map)

    feature_type = col("IntAct_feature_type")
    interaction_ac = col("IntAct_interaction_ac")
    # Sub-option columns (present when the corresponding flag was selected).
    feature_ac = col("IntAct_feature_ac")
    feature_short_label = col("IntAct_feature_short_label")
    feature_annotation = col("IntAct_feature_annotation")
    ap_ac = col("IntAct_ap_ac")
    interaction_participants = col("IntAct_interaction_participants")
    pmid = col("IntAct_pmid")
    if not any(
        [
            feature_type,
            interaction_ac,
            feature_ac,
            feature_short_label,
            feature_annotation,
            ap_ac,
            interaction_participants,
            pmid,
        ]
    ):
        return None
    return model.IntActAnnotation(
        feature_type=feature_type,
        interaction_ac=interaction_ac,
        feature_ac=feature_ac,
        feature_short_label=feature_short_label,
        feature_annotation=feature_annotation,
        ap_ac=ap_ac,
        interaction_participants=interaction_participants,
        pmid=pmid,
    )


def _parse_mutfunc(csq_values, index_map) -> model.MutfuncAnnotation | None:
    """Build a mutfunc annotation from its per-data-type score columns; returns
    None if none are present."""
    if not _has_any_column(
        index_map, "mutfunc_motif", "mutfunc_int", "mutfunc_mod", "mutfunc_exp"
    ):
        return None
    motif = _get_csq_value(csq_values, "mutfunc_motif", None, index_map)
    interactions = _get_csq_value(csq_values, "mutfunc_int", None, index_map)
    structure = _get_csq_value(csq_values, "mutfunc_mod", None, index_map)
    structure_exp = _get_csq_value(csq_values, "mutfunc_exp", None, index_map)
    if not any([motif, interactions, structure, structure_exp]):
        return None
    return model.MutfuncAnnotation(
        linear_motifs=_to_float(motif),
        protein_interactions=_to_float(interactions),
        protein_structure=_to_float(structure),
        protein_structure_experimental=_to_float(structure_exp),
    )


def _normalise_na(value: str | None) -> str | None:
    """Treat the literal 'NA' (used by some plugins for missing values) as None."""
    return None if value in (None, "NA") else value


def _first_amp(value: str | None) -> str | None:
    """First real (non-empty, non-'NA') item of a '&'-joined CSQ list."""
    for item in _split_amp(value):
        return item
    return None


def _parse_mavedb(csq_values, index_map) -> model.MaveDBAnnotation | None:
    """Build a MaveDB annotation; returns None if no MaveDB columns present.

    The MaveDB plugin reports several assays for one variant as parallel
    '&'-joined lists (score/urn/pro). Zip score and urn positionally into one
    assay each, so every score pairs with its score-set URN."""
    if not _has_any_column(
        index_map, "MaveDB_score", "MaveDB_urn", "MaveDB_doi", "MaveDB_nt", "MaveDB_pro"
    ):
        return None
    score = _get_csq_value(csq_values, "MaveDB_score", None, index_map)
    urn = _get_csq_value(csq_values, "MaveDB_urn", None, index_map)
    pro = _get_csq_value(csq_values, "MaveDB_pro", None, index_map)
    if not any([score, urn, pro]):
        return None
    scores = _raw_amp(score)
    urns = _raw_amp(urn)
    assays: list[model.MaveDBAssay] = []
    for i in range(max(len(scores), len(urns))):
        raw_urn = urns[i] if i < len(urns) else ""
        assay_urn = raw_urn if raw_urn and raw_urn != "NA" else None
        assay_score = _to_float(scores[i]) if i < len(scores) else None
        if assay_urn is None and assay_score is None:
            continue
        assays.append(model.MaveDBAssay(urn=assay_urn, score=assay_score))
    if not assays:
        return None
    return model.MaveDBAnnotation(
        protein_variant=_first_amp(pro),
        assays=assays,
    )


def _parse_popeve(csq_values, index_map) -> model.PopEve | None:
    """Build popEVE scores; returns None if no popEVE columns present."""
    if not _has_any_column(index_map, "popEVE_SCORE", "popEVE_EVE", "popEVE_mutant"):
        return None

    def v(name):
        return _get_csq_value(csq_values, name, None, index_map)

    score = v("popEVE_SCORE")
    if not any([score, v("popEVE_EVE"), v("popEVE_mutant")]):
        return None
    return model.PopEve(
        score=_to_float(score),
        eve=_to_float(v("popEVE_EVE")),
        esm1v=_to_float(v("popEVE_ESM1v")),
        pop_adjusted_eve=_to_float(v("popEVE_pop_adjusted_EVE")),
        pop_adjusted_esm1v=_to_float(v("popEVE_pop_adjusted_ESM1v")),
        gene=v("popEVE_gene"),
        protein=v("popEVE_protein"),
        mutant=v("popEVE_mutant"),
        gap_frequency=_to_float(v("popEVE_gap_frequency")),
    )


def _split_amp(value: str | None) -> list[str]:
    """Split a '&'-delimited CSQ list, dropping empties and 'NA' placeholders."""
    if not value:
        return []
    return [v for v in value.split("&") if v and v != "NA"]


def _parse_dosage_sensitivity(
    csq_values, index_map
) -> model.DosageSensitivity | None:
    """Build gnomAD dosage-sensitivity probabilities (pHaplo/pTriplo);
    returns None if neither column is present."""
    if not _has_any_column(index_map, "pHaplo", "pTriplo"):
        return None
    phaplo = _get_csq_value(csq_values, "pHaplo", None, index_map)
    ptriplo = _get_csq_value(csq_values, "pTriplo", None, index_map)
    if not any([phaplo, ptriplo]):
        return None
    return model.DosageSensitivity(
        phaplo=_to_float(phaplo), ptriplo=_to_float(ptriplo)
    )


def _parse_utr_annotation(
    csq_values, index_map
) -> model.FivePrimeUtrAnnotation | None:
    """Build a UTRAnnotator 5' UTR annotation; returns None if no 5'UTR /
    uORF columns are present."""
    if not _has_any_column(
        index_map, "5UTR_consequence", "5UTR_annotation", "Existing_uORFs",
        "Existing_InFrame_oORFs", "Existing_OutOfFrame_oORFs",
    ):
        return None
    consequence = _get_csq_value(csq_values, "5UTR_consequence", None, index_map)
    annotation = _get_csq_value(csq_values, "5UTR_annotation", None, index_map)
    uorfs = _get_csq_value(csq_values, "Existing_uORFs", None, index_map)
    inframe = _get_csq_value(
        csq_values, "Existing_InFrame_oORFs", None, index_map
    )
    outofframe = _get_csq_value(
        csq_values, "Existing_OutOfFrame_oORFs", None, index_map
    )
    if not any([consequence, annotation, uorfs, inframe, outofframe]):
        return None
    return model.FivePrimeUtrAnnotation(
        consequence=consequence,
        annotation=annotation,
        existing_uorfs=uorfs,
        existing_inframe_oorfs=inframe,
        existing_outofframe_oorfs=outofframe,
    )


def _parse_riboseq_orfs(
    csq_values, index_map
) -> model.RiboseqOrfsAnnotation | None:
    """Build a Ribo-seq ORFs annotation; returns None if no RiboseqORFs
    columns are present."""
    if not _has_any_column(
        index_map, "RiboseqORFs_id", "RiboseqORFs_consequences", "RiboseqORFs_impact",
        "RiboseqORFs_protein_position", "RiboseqORFs_codons",
        "RiboseqORFs_amino_acids", "RiboseqORFs_publications",
    ):
        return None
    orf_id = _get_csq_value(csq_values, "RiboseqORFs_id", None, index_map)
    consequences = _split_amp(
        _get_csq_value(csq_values, "RiboseqORFs_consequences", None, index_map)
    )
    impact = _get_csq_value(csq_values, "RiboseqORFs_impact", None, index_map)
    protein_position = _get_csq_value(
        csq_values, "RiboseqORFs_protein_position", None, index_map
    )
    codons = _get_csq_value(csq_values, "RiboseqORFs_codons", None, index_map)
    amino_acids = _get_csq_value(
        csq_values, "RiboseqORFs_amino_acids", None, index_map
    )
    publications = _split_amp(
        _get_csq_value(csq_values, "RiboseqORFs_publications", None, index_map)
    )
    if not any(
        [orf_id, consequences, impact, protein_position, codons, amino_acids, publications]
    ):
        return None
    return model.RiboseqOrfsAnnotation(
        orf_id=orf_id,
        consequences=consequences,
        impact=impact,
        protein_position=protein_position,
        codons=codons,
        amino_acids=amino_acids,
        publications=publications,
    )


_PREDICTION_RE = re.compile(r"^(?P<prediction>[^(]+)\((?P<score>[-\d.eE]+)\)$")


def _parse_prediction(value: str | None) -> model.PredictionWithScore | None:
    """Parse a 'prediction(score)' CSQ value, e.g. SIFT 'tolerated(0.15)'."""
    if not value:
        return None
    match = _PREDICTION_RE.match(value.strip())
    if match:
        return model.PredictionWithScore(
            prediction=match.group("prediction"),
            score=_to_float(match.group("score")),
        )
    return model.PredictionWithScore(prediction=value, score=None)


def _parse_hgvs(csq_values, index_map) -> model.HgvsNotations | None:
    if not _has_any_column(index_map, "HGVSg", "HGVSc", "HGVSp"):
        return None
    genomic = _get_csq_value(csq_values, "HGVSg", None, index_map)
    transcript = _get_csq_value(csq_values, "HGVSc", None, index_map)
    protein = _get_csq_value(csq_values, "HGVSp", None, index_map)
    if not any([genomic, transcript, protein]):
        return None
    return model.HgvsNotations(
        genomic=genomic, transcript=transcript, protein=protein
    )


def _parse_spliceai(csq_values, index_map) -> model.SpliceAiScores | None:
    if not _has_any_column(index_map, "SpliceAI_pred_DS_AG"):
        return None

    def f(name):
        return _to_float(_get_csq_value(csq_values, name, None, index_map))

    def i(name):
        value = _get_csq_value(csq_values, name, None, index_map)
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    # Compute the delta scores first; only build the model when at least one is
    # present (avoids constructing + model_dump-ing an all-None model per
    # transcript just to discover it is empty).
    values = dict(
        ds_acceptor_gain=f("SpliceAI_pred_DS_AG"),
        ds_acceptor_loss=f("SpliceAI_pred_DS_AL"),
        ds_donor_gain=f("SpliceAI_pred_DS_DG"),
        ds_donor_loss=f("SpliceAI_pred_DS_DL"),
        dp_acceptor_gain=i("SpliceAI_pred_DP_AG"),
        dp_acceptor_loss=i("SpliceAI_pred_DP_AL"),
        dp_donor_gain=i("SpliceAI_pred_DP_DG"),
        dp_donor_loss=i("SpliceAI_pred_DP_DL"),
    )
    if not any(value is not None for value in values.values()):
        return None
    symbol = _get_csq_value(csq_values, "SpliceAI_pred_SYMBOL", None, index_map)
    return model.SpliceAiScores(symbol=symbol, **values)


def _parse_pathogenicity(csq_values, index_map) -> model.PathogenicityPredictions | None:
    if not _has_any_column(
        index_map,
        "SIFT", "PolyPhen", "REVEL", "am_class", "am_pathogenicity",
        "CADD_PHRED", "CADD_RAW", "EVE_CLASS", "EVE_SCORE",
        "SpliceAI_pred_DS_AG", "popEVE_SCORE",
    ):
        return None

    def v(name):
        return _get_csq_value(csq_values, name, None, index_map)

    sift = _parse_prediction(v("SIFT"))
    polyphen = _parse_prediction(v("PolyPhen"))
    spliceai = _parse_spliceai(csq_values, index_map)
    fields = dict(
        sift=sift,
        polyphen=polyphen,
        revel=_to_float(v("REVEL")),
        alphamissense_class=v("am_class"),
        alphamissense_score=_to_float(v("am_pathogenicity")),
        cadd_phred=_to_float(v("CADD_PHRED")),
        cadd_raw=_to_float(v("CADD_RAW")),
        spliceai=spliceai,
        eve_class=v("EVE_CLASS"),
        eve_score=_to_float(v("EVE_SCORE")),
        popeve=_parse_popeve(csq_values, index_map),
    )
    if not any(value is not None for value in fields.values()):
        return None
    return model.PathogenicityPredictions(**fields)


def _parse_population_frequencies(
    csq_values, index_map, overall_key: str, pop_pattern: str
) -> model.PopulationFrequencies | None:
    """Build overall + per-population frequencies. pop_pattern has a single
    '{}' where the population code sits, e.g. 'gnomADe_{}_AF'."""
    overall = _to_float(_get_csq_value(csq_values, overall_key, None, index_map))
    prefix, suffix = pop_pattern.split("{}")
    populations: dict[str, float] = {}
    for column in index_map:
        if (
            column.startswith(prefix)
            and column.endswith(suffix)
            and column != overall_key
        ):
            pop = column[len(prefix): len(column) - len(suffix)]
            value = _to_float(_get_csq_value(csq_values, column, None, index_map))
            if value is not None:
                populations[pop] = value
    if overall is None and not populations:
        return None
    return model.PopulationFrequencies(overall=overall, populations=populations)


def _parse_frequencies(csq_values, index_map) -> model.AlleleFrequencies | None:
    if not _has_any_column(
        index_map, "gnomAD_exomes_AF", "gnomAD_genomes_AF", "AoU_gvs_all_af"
    ):
        return None
    # gnomAD exomes/genomes are added via VEP `custom` tracks with
    # short_name=gnomAD_exomes / gnomAD_genomes, so their columns are prefixed
    # accordingly and the emitted field is `AF` (overall) plus `AF_<...>` variants
    # (ancestry / sex / subset). Everything after the `<track>_AF_` prefix becomes
    # the population key (e.g. afr, nfe_XX, non_ukb_afr, grpmax).
    gnomad_exomes = _parse_population_frequencies(
        csq_values, index_map, "gnomAD_exomes_AF", "gnomAD_exomes_AF_{}"
    )
    gnomad_genomes = _parse_population_frequencies(
        csq_values, index_map, "gnomAD_genomes_AF", "gnomAD_genomes_AF_{}"
    )
    # All of Us is added via a VEP `custom` with short_name=AoU, so its columns
    # are prefixed AoU_ (e.g. AoU_gvs_all_af, AoU_gvs_afr_af). The AoU_gvs_max_af
    # subpopulation frequency is picked up as the "max" population; the companion
    # AoU_gvs_max_subpop column names which subpopulation that max came from
    # (a code, e.g. "eur"), attached so the results view can show it in brackets.
    all_of_us = _parse_population_frequencies(
        csq_values, index_map, "AoU_gvs_all_af", "AoU_gvs_{}_af"
    )
    if all_of_us is not None:
        max_subpop = _get_csq_value(csq_values, "AoU_gvs_max_subpop", None, index_map)
        if max_subpop:
            all_of_us.max_subpopulation = max_subpop
    if not any([gnomad_exomes, gnomad_genomes, all_of_us]):
        return None
    return model.AlleleFrequencies(
        gnomad_exomes=gnomad_exomes,
        gnomad_genomes=gnomad_genomes,
        all_of_us=all_of_us,
    )


def _parse_go(csq_values, index_map) -> list[model.GoTerm]:
    """Gene Ontology terms from the GO plugin's GO column. Each '&'-separated
    entry is `GO:<id>:<term>` (term underscore-delimited), e.g.
    `GO:0001558:regulation_of_cell_growth` -> id 'GO:0001558', name
    'regulation of cell growth'."""
    if not _has_any_column(index_map, "GO"):
        return []
    terms: list[model.GoTerm] = []
    for entry in _split_amp(_get_csq_value(csq_values, "GO", None, index_map)):
        parts = entry.split(":")
        if len(parts) < 3:
            continue
        go_id = ":".join(parts[:2])  # "GO:0001558"
        name = ":".join(parts[2:]).replace("_", " ").strip()
        terms.append(model.GoTerm(id=go_id, name=name))
    return terms


def _parse_phenotype_data(csq_values, index_map) -> model.VariantPhenotypeData | None:
    if not _has_any_column(index_map, "PHENOTYPES", "CLIN_SIG", "PUBMED"):
        return None
    phenotypes = _split_amp(_get_csq_value(csq_values, "PHENOTYPES", None, index_map))
    clin_sig = _split_amp(_get_csq_value(csq_values, "CLIN_SIG", None, index_map))
    pubmed = _split_amp(_get_csq_value(csq_values, "PUBMED", None, index_map))
    if not any([phenotypes, clin_sig, pubmed]):
        return None
    return model.VariantPhenotypeData(
        phenotypes=phenotypes,
        clinical_significance=clin_sig,
        pubmed_ids=pubmed,
    )


_CLNSIG_CONFLICTING = "Conflicting_classifications_of_pathogenicity"
# A CLNSIGCONF token, e.g. "Likely_pathogenic_(6)" -> ("Likely_pathogenic", 6).
_CLNSIGCONF_RE = re.compile(r"^(?P<term>.+)_\((?P<count>\d+)\)$")


def _parse_clinvar(csq_values, index_map) -> model.ClinVarAnnotation | None:
    """ClinVar clinical significance from the ClinVar custom track. CLNSIG is the
    overall classification; CLNSIGCONF (the per-classification submission
    breakdown) is surfaced only when the classification is conflicting, and
    ignored otherwise. Returns None if there is no CLNSIG value."""
    if not _has_any_column(index_map, "ClinVar_CLNSIG"):
        return None
    clnsig = _get_csq_value(csq_values, "ClinVar_CLNSIG", None, index_map)
    if not clnsig:
        return None
    significance = _split_amp(clnsig)

    breakdown: list[model.ClinVarSignificance] = []
    if _CLNSIG_CONFLICTING in significance:
        conf = _get_csq_value(csq_values, "ClinVar_CLNSIGCONF", None, index_map)
        for token in _split_amp(conf):
            match = _CLNSIGCONF_RE.match(token)
            if match:
                breakdown.append(
                    model.ClinVarSignificance(
                        significance=match.group("term"),
                        count=int(match.group("count")),
                    )
                )
    return model.ClinVarAnnotation(
        significance=significance, conflicting_breakdown=breakdown
    )


def _raw_amp(value: str | None) -> list[str]:
    """Split a '&'-delimited CSQ list keeping every position (incl. 'NA'), so
    positionally-aligned subfields can be zipped together."""
    return value.split("&") if value else []


def _parse_open_targets(csq_values, index_map) -> model.OpenTargetsAssociation | None:
    if not _has_any_column(
        index_map, "OpenTargets_gwasDiseases", "OpenTargets_qtlGeneId"
    ):
        return None

    def raw(name):
        return _raw_amp(_get_csq_value(csq_values, name, None, index_map))

    diseases = raw("OpenTargets_gwasDiseases")
    gene_ids = raw("OpenTargets_gwasGeneId")
    l2g = raw("OpenTargets_gwasLocusToGeneScore")
    qtl_genes = raw("OpenTargets_qtlGeneId")
    qtl_biosamples = raw("OpenTargets_qtlBiosampleName")

    # GWAS: disease / gene / L2G are positionally aligned. Keep real (non-NA)
    # rows and de-duplicate (the plugin currently emits duplicate rows).
    gwas_associations = []
    seen_gwas = set()
    for disease, gene_id, score in zip(diseases, gene_ids, l2g):
        if not disease or disease == "NA":
            continue
        key = (disease, gene_id, score)
        if key in seen_gwas:
            continue
        seen_gwas.add(key)
        gwas_associations.append(
            model.OpenTargetsGwasAssociation(
                disease=disease, gene_id=gene_id, l2g_score=_to_float(score)
            )
        )

    # QTL: gene / biosample are positionally aligned. Keep real rows, unique.
    qtl_associations = []
    seen_qtl = set()
    for gene_id, biosample in zip(qtl_genes, qtl_biosamples):
        if not gene_id or gene_id == "NA":
            continue
        key = (gene_id, biosample)
        if key in seen_qtl:
            continue
        seen_qtl.add(key)
        qtl_associations.append(
            model.OpenTargetsQtlAssociation(
                gene_id=gene_id,
                biosample=None if biosample in ("", "NA") else biosample,
            )
        )

    # Strongest associations first (missing scores last).
    gwas_associations.sort(
        key=lambda a: a.l2g_score if a.l2g_score is not None else float("-inf"),
        reverse=True,
    )

    if not gwas_associations and not qtl_associations:
        return None
    return model.OpenTargetsAssociation(
        gwas_associations=gwas_associations,
        qtl_associations=qtl_associations,
    )


def _get_alt_allele_details(
    ref: str, alt: str, csqs: list[str], index_map: dict[str, int]
) -> model.AlternativeVariantAllele:
    """Creates  AlternativeVariantAllele based on
    target alt allele and CSQ entires"""
    frequency = None
    consequences = []
    allele_type = _get_variant_type(ref, alt)
    # Allele-level annotations are identical across all of this allele's CSQ
    # rows, so capture them once (from the first matching row).
    frequencies = None
    colocated_variants: list[str] = []
    phenotype_data = None
    open_targets = None
    spdi = None
    hgvsg = None
    cadd_phred = None
    cadd_raw = None
    clinvar = None
    allele_level_captured = False

    for str_csq in csqs:
        csq_values = str_csq.split("|")

        if csq_values[index_map["Allele"]] != alt:
            continue

        frequency = _get_csq_value(csq_values, "AF", frequency, index_map)

        if not allele_level_captured:
            frequencies = _parse_frequencies(csq_values, index_map)
            colocated_variants = _split_amp(
                _get_csq_value(csq_values, "Existing_variation", None, index_map)
            )
            phenotype_data = _parse_phenotype_data(csq_values, index_map)
            open_targets = _parse_open_targets(csq_values, index_map)
            # Variant-level representations / genome-wide scores. These are the
            # same across an allele's transcripts and are the only annotations
            # available for intergenic variants (which have no transcript rows).
            spdi = _get_csq_value(csq_values, "SPDI", None, index_map)
            hgvsg = _get_csq_value(csq_values, "HGVSg", None, index_map)
            cadd_phred = _to_float(
                _get_csq_value(csq_values, "CADD_PHRED", None, index_map)
            )
            cadd_raw = _to_float(
                _get_csq_value(csq_values, "CADD_RAW", None, index_map)
            )
            clinvar = _parse_clinvar(csq_values, index_map)
            allele_level_captured = True

        cons = _get_csq_value(csq_values, "Consequence", "", index_map)
        if len(cons) == 0:
            cons = []
        else:
            cons = cons.split("&")
        if csq_values[index_map["Feature_type"]] == "Transcript":
            is_canonical = (
                _get_csq_value(csq_values, "CANONICAL", "NO", index_map) == "YES"
            )

            # It looks like for Feature_type = Transcript that we always have a STRAND value
            strand = (
                model.Strand.reverse
                if _get_csq_value(csq_values, "STRAND", "1", index_map) == "-1"
                else model.Strand.forward
            )

            # MANE: depending on the VEP run, either the MANE column carries the
            # label (MANE_Select / MANE_Plus_Clinical) or the MANE_SELECT /
            # MANE_PLUS_CLINICAL columns carry the matched RefSeq id. Handle both.
            mane_label = _get_csq_value(csq_values, "MANE", None, index_map)
            mane_select_refseq = _get_csq_value(
                csq_values, "MANE_SELECT", None, index_map
            )
            mane_plus_clinical = _get_csq_value(
                csq_values, "MANE_PLUS_CLINICAL", None, index_map
            )
            is_mane_select = bool(mane_select_refseq) or mane_label == "MANE_Select"
            is_mane_plus_clinical = (
                bool(mane_plus_clinical) or mane_label == "MANE_Plus_Clinical"
            )

            consequences.append(
                model.PredictedTranscriptConsequence(
                    feature_type=model.FeatureType.transcript,
                    stable_id=_get_csq_value(csq_values, "Feature", "", index_map),
                    gene_stable_id=_get_csq_value(csq_values, "Gene", "", index_map),
                    biotype=_get_csq_value(csq_values, "BIOTYPE", "", index_map),
                    is_canonical=is_canonical,
                    gene_symbol=_get_csq_value(csq_values, "SYMBOL", None, index_map),
                    consequences=cons,
                    strand=strand,
                    # MANE
                    is_mane_select=is_mane_select,
                    is_mane_plus_clinical=is_mane_plus_clinical,
                    mane_select_refseq_id=mane_select_refseq,
                    # Protein & functional annotations
                    ensembl_protein_id=_get_csq_value(
                        csq_values, "ENSP", None, index_map
                    ),
                    uniprot=_parse_uniprot(csq_values, index_map),
                    protein_matches=_parse_protein_matches(csq_values, index_map),
                    protvar=_parse_protvar(csq_values, index_map),
                    intact=_parse_intact(csq_values, index_map),
                    mutfunc=_parse_mutfunc(csq_values, index_map),
                    mavedb=_parse_mavedb(csq_values, index_map),
                    # Variant representations
                    spdi=_get_csq_value(csq_values, "SPDI", None, index_map),
                    # HGVS, pathogenicity, gene constraint
                    hgvs=_parse_hgvs(csq_values, index_map),
                    pathogenicity=_parse_pathogenicity(csq_values, index_map),
                    loeuf=_to_float(
                        _get_csq_value(csq_values, "LOEUF", None, index_map)
                    ),
                    dosage_sensitivity=_parse_dosage_sensitivity(
                        csq_values, index_map
                    ),
                    # Transcript-level predictions
                    utr_annotation=_parse_utr_annotation(csq_values, index_map),
                    riboseq_orfs=_parse_riboseq_orfs(csq_values, index_map),
                    go_terms=_parse_go(csq_values, index_map),
                )
            )
        elif "intergenic_variant" in cons:
            consequences.append(
                model.PredictedIntergenicConsequence(
                    feature_type=None,
                    consequences=["intergenic_variant"],
                )
            )

    return model.AlternativeVariantAllele(
        allele_sequence=("" if alt=="copy_number_variation" else alt),
        allele_type=allele_type,
        representative_population_allele_frequency=frequency,
        spdi=spdi,
        hgvsg=hgvsg,
        cadd_phred=cadd_phred,
        cadd_raw=cadd_raw,
        frequencies=frequencies,
        colocated_variants=colocated_variants,
        phenotype_data=phenotype_data,
        open_targets=open_targets,
        clinvar=clinvar,
        predicted_molecular_consequences=consequences,
    )


# ---------------------------------------------------------------------------
# DEV/LOCAL-ONLY: stale metadata-cache guard
#
# In the full pipeline, every run writes its output VCF (and this metadata
# cache) into its own directory, so the cache can never go stale. This guard
# only matters when a single local output file is parsed repeatedly and then
# regenerated (e.g. re-running the pipeline against the fixed LOCAL_RESULTS_VCF
# path during development): without it the old header/record counts would be
# reused and mis-slice the new file. Safe to remove once outputs always land in
# per-run directories.
# ---------------------------------------------------------------------------
def _is_meta_cache_stale(meta_path: FilePath, vcf_path: FilePath) -> bool:
    """True if the metadata cache exists but predates its VCF."""
    return (
        meta_path.exists()
        and meta_path.stat().st_mtime < vcf_path.stat().st_mtime
    )


def _get_vcf_meta(vcf_path: FilePath) -> model.VcfMetadata:
    """Helper method to manage metainfo for a VCF file"""

    meta_path = vcf_path.with_name(META_FILE)
    if not meta_path.exists() or _is_meta_cache_stale(meta_path, vcf_path):
        variant_count_str = subprocess.check_output(
            f"bcftools stats {vcf_path} | grep 'number of records:'",
            shell=True, text=True
        )
        header_count_str = subprocess.check_output(
            f"bcftools view -h {vcf_path} | wc -l",
            shell=True, text=True
        )
        try:
            vcf_info = model.VcfMetadata(
                variant_count=int(variant_count_str.split(":")[-1]),
                header_count=int(header_count_str)
            )
        except ValueError as e:
            e.args = (
                f"_get_vcf_meta: unexpected bcftools output: variant_count: {variant_count_str} | header_count: {header_count_str}",
                *e.args,
            )
            raise

        with open(meta_path, "w") as meta_file:
            meta_file.write(vcf_info.model_dump_json())
    else:
        with open(meta_path, "r") as meta_file:
            vcf_info = model.VcfMetadata.model_validate_json(meta_file.read())
    return vcf_info


# ---------------------------------------------------------------------------
# BGZF page-index seek path
#
# When the pipeline emits a `<vcf>.pageidx.json` sidecar (see
# pagination-design.md / build_page_index.py), a page can be fetched by seeking
# straight to it instead of scanning from the top with bcftools. The sidecar
# stores, every `stride` records, the packed BGZF virtual offset
# (compressed_block_offset << 16 | within_block_offset) of that record's line.
# ---------------------------------------------------------------------------
PAGE_INDEX_SUFFIX = ".pageidx.json"


def _bc_block_size(extra: bytes) -> int | None:
    """BSIZE from a gzip member's FEXTRA field, or None if no 'BC' subfield."""
    i = 0
    while i + 4 <= len(extra):
        si1, si2, slen = extra[i], extra[i + 1], struct.unpack("<H", extra[i + 2 : i + 4])[0]
        if si1 == 0x42 and si2 == 0x43:  # 'B', 'C'
            return struct.unpack("<H", extra[i + 4 : i + 6])[0]
        i += 4 + slen
    return None


class _BgzfReader:
    """Minimal, dependency-free BGZF reader supporting seek by virtual offset.

    Only the operations the page-index seek needs are implemented: `seek` to a
    packed virtual offset, line-oriented `readline` (stitching lines that cross
    block boundaries), and `tell` (the virtual offset of the next byte).
    """

    def __init__(self, path: str):
        self._file = open(path, "rb")
        self._block_coffset = 0  # compressed offset of the current block
        self._data = b""  # decompressed bytes of the current block
        self._pos = 0  # position within the current block
        self._load_block_at(0)

    def _read_block(self):
        """Read/decompress the block at the current file position -> (coffset,
        data); data is None at EOF, b'' for the empty BGZF EOF-marker block."""
        coffset = self._file.tell()
        header = self._file.read(12)
        if not header:
            return coffset, None
        if len(header) < 12 or header[:2] != b"\x1f\x8b":
            raise ValueError("not a BGZF block")
        xlen = struct.unpack("<H", header[10:12])[0]
        extra = self._file.read(xlen)
        bsize = _bc_block_size(extra)
        if bsize is None:
            raise ValueError("missing BGZF 'BC' subfield")
        cdata = self._file.read((bsize + 1) - 12 - xlen - 8)
        self._file.read(8)  # trailer (CRC32 + ISIZE)
        return coffset, (zlib.decompress(cdata, -15) if cdata else b"")

    def _load_block_at(self, coffset: int) -> None:
        self._file.seek(coffset)
        self._block_coffset, data = self._read_block()
        self._data = data or b""
        self._pos = 0

    def seek(self, voffset: int) -> None:
        self._load_block_at(voffset >> 16)
        self._pos = voffset & 0xFFFF

    def tell(self) -> int:
        return (self._block_coffset << 16) | self._pos

    def _advance_block(self) -> bool:
        """Load the next non-empty block; False at EOF (empty EOF-marker blocks
        are skipped)."""
        while True:
            coffset, data = self._read_block()
            if data is None:
                return False
            self._block_coffset, self._data, self._pos = coffset, data, 0
            if data:
                return True

    def readline(self) -> bytes:
        chunks: list[bytes] = []
        while True:
            if self._pos >= len(self._data):
                if not self._advance_block():
                    break  # EOF
            newline = self._data.find(b"\n", self._pos)
            if newline == -1:
                chunks.append(self._data[self._pos :])
                self._pos = len(self._data)
                continue  # line continues into the next block
            chunks.append(self._data[self._pos : newline + 1])
            self._pos = newline + 1
            break
        # Normalise an end-of-block position to the next block's start, so tell()
        # reports the next line's virtual offset the way the generator recorded
        # it (a line beginning a block has within-block offset 0).
        if self._pos >= len(self._data):
            self._advance_block()
        return b"".join(chunks)

    def close(self) -> None:
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def _load_page_index(vcf_path: FilePath) -> dict | None:
    """The parsed `<vcf>.pageidx.json` sidecar, or None if it doesn't exist."""
    index_path = Path(str(vcf_path) + PAGE_INDEX_SUFFIX)
    if not index_path.exists():
        return None
    return json.loads(index_path.read_text())


def _read_indexed_page(
    vcf_path: FilePath, index: dict, page: int, page_size: int
) -> tuple[str, str]:
    """Return (header_text, page_rows_text) for the requested page by seeking to
    the nearest checkpoint and reading forward. `page` is 1-based; a page past
    the end yields empty rows."""
    total = index["total_records"]
    stride = index["stride"]
    checkpoints = index["checkpoints"]
    header_end = index["header_end_voffset"]
    start = (max(page, 1) - 1) * page_size

    header_lines: list[bytes] = []
    rows: list[bytes] = []
    with _BgzfReader(str(vcf_path)) as reader:
        # Header = every line before the first data record.
        while reader.tell() < header_end:
            line = reader.readline()
            if not line:
                break
            header_lines.append(line)
        # Seek to the checkpoint at/before the page start, skip the remainder.
        if page_size > 0 and start < total:
            checkpoint = start // stride
            reader.seek(checkpoints[checkpoint])
            for _ in range(start - checkpoint * stride):
                reader.readline()
            for _ in range(min(page_size, total - start)):
                line = reader.readline()
                if not line:
                    break
                rows.append(line)

    return b"".join(header_lines).decode(), b"".join(rows).decode()


def get_results_from_path(
    page_size: int, page: int, vcf_path: FilePath
) -> model.VepResultsResponse:
    """Returns a page of VCF data from the given filepath.
    Slices the input VCF file to a smaller one
    and converts it to stream for get_results_from_stream"""

    # Fast path: if the pipeline emitted a page-index sidecar, seek to the page
    # instead of scanning the file / shelling out to bcftools.
    index = _load_page_index(vcf_path)
    if index is not None:
        page = max(page, 1)
        page_size = max(page_size, 0)
        header_text, rows_text = _read_indexed_page(vcf_path, index, page, page_size)
        return get_results_from_stream(
            page_size,
            page,
            index["total_records"],
            StringIO(header_text + rows_text),
            presliced=True,
        )

    # Fallback (no sidecar): scan the file from the top through page*page_size
    # records and shell out to bcftools for the counts. `head` short-circuits so
    # it stops at the offset rather than scanning the whole file, but deep pages
    # get slower and the last page is a full pass. Runs from the pipeline now ship
    # a page-index sidecar (handled above); this remains for older/un-indexed
    # outputs. Longer term, a queryable store (SQLite/Parquet) would also enable
    # sorting/filtering (see pagination-design.md).
    # Fetch a pageful of variant records with headers
    vcf_info = _get_vcf_meta(vcf_path)
    total = vcf_info.variant_count
    page = max(page, 1) # normalize values
    page_size = min(max(page_size, 0), total)
    row_offset = min(page * page_size, total) + vcf_info.header_count
    vcf_headers = subprocess.check_output( # fetch all header rows
        f"bcftools view -h {vcf_path}", shell=True, text=True
    )
    vcf_slice = subprocess.check_output( # fetch subset of variant rows
        f"bcftools view {vcf_path} | head -n{row_offset} | tail -n{page_size}",
        shell=True, text=True
    )
    vcf_stream = StringIO(vcf_headers + vcf_slice)

    return get_results_from_stream(page_size, page, total, vcf_stream)


def get_results_from_stream(
    page_size: int, page: int, total: int, vcf_stream: StringIO,
    presliced: bool = False,
) -> model.VepResultsResponse:
    """Helper method to convert a filestream to VCF records for _get_results_from_vcfpy"""

    # Load vcf
    vcf_records = vcfpy.Reader.from_stream(vcf_stream)
    return _get_results_from_vcfpy(page_size, page, total, vcf_records, presliced)


def _get_results_from_vcfpy(
    page_size: int, page: int, total: int, vcf_records: vcfpy.Reader,
    presliced: bool = False,
) -> model.VepResultsResponse:
    """Generates a page of VCF data in the format described in
    APISpecification.yaml for a given VCFPY reader"""

    # Parse csq header
    csq_header = vcf_records.header.get_info_field_info("CSQ").description
    if not csq_header:
        raise Exception("CSQ header missing")

    prediction_index_map = _get_prediction_index_map(csq_header)
    # Required CSQ column (the rest use fallback values)
    if "Allele" not in prediction_index_map:
        raise Exception("Allele column missing from CSQ header")

    variants = []
    # populate variants page. `presliced` means the stream already contains
    # exactly this page's rows (the index seek path), so the page-bounds guard —
    # which the scan path needs to return empty past the end — is skipped.
    if presliced or page*page_size <= total:
        for record in vcf_records:
            if record is None:
                break
            if record.CHROM.startswith("chr"):
                record.CHROM = record.CHROM[3:]

            # https://github.com/bihealth/vcfpy/blob/697768d032b6b476766fb4c524c91c8d24559330/vcfpy/record.py#L63
            # end does not look like it is implemented.
            # https://github.com/Penghui-Wang/PyVCF/blob/master/vcf/model.py#L190
            # from competing vcf module
            location = model.Location(
                region_name=record.CHROM,
                start=record.POS,
                end=record.POS + len(record.REF),
            )

            if "CSQ" not in record.INFO:
                csq_strings = []
                alt_allele_strings = [_alt_value(alt) for alt in record.ALT]
            else:
                csq_strings = record.INFO["CSQ"]
                alt_allele_strings = list(set([
                    csq_string.split("|")[prediction_index_map["Allele"]]
                    for csq_string in csq_strings
                ]))

            alt_alleles = [
                _get_alt_allele_details(record.REF, alt, csq_strings, prediction_index_map)
                for alt in alt_allele_strings
            ]

            longest_alt = max((_alt_value(a) for a in record.ALT), key=len)

            variants.append(
                model.Variant(
                    name=";".join(record.ID) if len(record.ID) > 0 else ".",
                    location=location,
                    reference_allele=model.ReferenceVariantAllele(
                        allele_sequence=record.REF
                    ),
                    alternative_alleles=alt_alleles,
                    allele_type=_get_variant_type(record.REF, longest_alt),
                )
            )

    return model.VepResultsResponse(
        metadata=model.Metadata(
            pagination=model.PaginationMetadata(
                page=page, per_page=page_size, total=total
            )
        ),
        variants=variants,
    )


# ---------------------------------------------------------------------------
# Flattened "expanded columnar" export (human-readable / spreadsheet-friendly).
#
# One tab-separated row per CSQ entry — every allele x feature, transcript AND
# intergenic, i.e. fully expanded — with a column for the variant location plus
# every CSQ field (all VEP + plugin annotations). Streamed line by line so the
# whole file needn't be held in memory. Reads the VCF directly (bgzip is
# gzip-readable), so it needs neither vcfpy nor bcftools.
# ---------------------------------------------------------------------------
def _parse_csq_format(info_line: str) -> list[str]:
    """Pull the pipe-delimited CSQ field names out of the ##INFO CSQ header."""
    fmt = info_line.split("Format:", 1)[1].split('">')[0].strip()
    return fmt.split("|")


def gzip_text_stream(chunks: Iterator[str], level: int = 6) -> Iterator[bytes]:
    """Gzip-compress a stream of text chunks on the fly, yielding gzip-format
    bytes. Used to serve the flattened TSV download compressed (plain gzip, for
    broad compatibility) without buffering the whole table in memory."""
    # wbits 16 + MAX_WBITS => a standalone gzip container (header + trailer).
    compressor = zlib.compressobj(level, zlib.DEFLATED, 16 + zlib.MAX_WBITS)
    for chunk in chunks:
        compressed = compressor.compress(chunk.encode("utf-8"))
        if compressed:
            yield compressed
    tail = compressor.flush()
    if tail:
        yield tail


def stream_vep_tsv(vcf_path: FilePath) -> Iterator[str]:
    """Yield the VEP output VCF flattened to tab-separated rows (with header)."""
    csq_fields: list[str] | None = None
    header_emitted = False
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = _parse_csq_format(line)
                continue
            if line.startswith("#") or csq_fields is None:
                continue
            if not header_emitted:
                yield (
                    "\t".join(["Uploaded_variation", "Location", "Ref"] + csq_fields)
                    + "\n"
                )
                header_emitted = True
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8:
                continue
            chrom, pos, variant_id, ref = columns[:4]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            location = f"{chrom}:{pos}"
            csq = next(
                (c[4:] for c in columns[7].split(";") if c.startswith("CSQ=")),
                None,
            )
            if not csq:
                continue
            for entry in csq.split(","):
                values = entry.split("|")
                if len(values) < len(csq_fields):
                    values += [""] * (len(csq_fields) - len(values))
                row = [variant_id, location, ref] + values[: len(csq_fields)]
                yield "\t".join(row) + "\n"
