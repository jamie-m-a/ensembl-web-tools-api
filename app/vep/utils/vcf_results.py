""" Module for loading a VCF and parsing it into a VepResultsResponse
object as defined in APISpecification"""

from io import StringIO
import re
import subprocess
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


def _to_float(value: str | None) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_uniprot(csq_values, index_map) -> model.UniprotIds | None:
    """Build Uniprot cross-references from the SWISSPROT/TREMBL/UNIPARC/isoform
    CSQ columns; returns None if none are present."""
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
    residues. The leading id and trailing residue token are non-numeric; the
    middle parts are numeric scores."""
    parts = raw.split("&")
    pocket_id = parts[0] if parts else ""
    numeric: list[float | None] = []
    residues: list[str] = []
    for part in parts[1:]:
        number = _to_float(part)
        if number is not None:
            numeric.append(number)
        elif part:
            # residue token, e.g. "p79p81p82" -> ["79", "81", "82"]
            residues.extend(r for r in part.split("p") if r)
    numeric += [None] * (5 - len(numeric))
    return model.ProtVarPocket(
        pocket_id=pocket_id,
        energy=numeric[0],
        energy_per_volume=numeric[1],
        score=numeric[2],
        buriedness=numeric[3],
        radius_of_gyration=numeric[4],
        residues=residues,
        raw=raw,
    )


def _parse_protvar(csq_values, index_map) -> model.ProtVarAnnotation | None:
    """Build a ProtVar annotation from the stability/pocket/interaction CSQ
    columns; returns None if none are present."""
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
    """Build an IntAct annotation; returns None if no IntAct columns present."""
    feature_ac = _get_csq_value(csq_values, "IntAct_feature_ac", None, index_map)
    feature_type = _get_csq_value(
        csq_values, "IntAct_feature_type", None, index_map
    )
    interaction_ac = _get_csq_value(
        csq_values, "IntAct_interaction_ac", None, index_map
    )
    if not any([feature_ac, feature_type, interaction_ac]):
        return None
    return model.IntActAnnotation(
        feature_ac=feature_ac,
        feature_type=feature_type,
        interaction_ac=interaction_ac,
    )


def _parse_mutfunc(csq_values, index_map) -> model.MutfuncAnnotation | None:
    """Build a mutfunc annotation from its per-data-type score columns; returns
    None if none are present."""
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


def _parse_mavedb(csq_values, index_map) -> model.MaveDBAnnotation | None:
    """Build a MaveDB annotation; returns None if no MaveDB columns present."""
    score = _get_csq_value(csq_values, "MaveDB_score", None, index_map)
    urn = _get_csq_value(csq_values, "MaveDB_urn", None, index_map)
    doi = _get_csq_value(csq_values, "MaveDB_doi", None, index_map)
    nt = _get_csq_value(csq_values, "MaveDB_nt", None, index_map)
    pro = _get_csq_value(csq_values, "MaveDB_pro", None, index_map)
    if not any([score, urn, doi, nt, pro]):
        return None
    return model.MaveDBAnnotation(
        score=_to_float(score),
        urn=_normalise_na(urn),
        doi=_normalise_na(doi),
        nucleotide_variant=_normalise_na(nt),
        protein_variant=_normalise_na(pro),
    )


def _split_amp(value: str | None) -> list[str]:
    """Split a '&'-delimited CSQ list, dropping empties and 'NA' placeholders."""
    if not value:
        return []
    return [v for v in value.split("&") if v and v != "NA"]


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
    genomic = _get_csq_value(csq_values, "HGVSg", None, index_map)
    transcript = _get_csq_value(csq_values, "HGVSc", None, index_map)
    protein = _get_csq_value(csq_values, "HGVSp", None, index_map)
    if not any([genomic, transcript, protein]):
        return None
    return model.HgvsNotations(
        genomic=genomic, transcript=transcript, protein=protein
    )


def _parse_spliceai(csq_values, index_map) -> model.SpliceAiScores | None:
    def f(name):
        return _to_float(_get_csq_value(csq_values, name, None, index_map))

    def i(name):
        value = _get_csq_value(csq_values, name, None, index_map)
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    symbol = _get_csq_value(csq_values, "SpliceAI_pred_SYMBOL", None, index_map)
    scores = model.SpliceAiScores(
        symbol=symbol,
        ds_acceptor_gain=f("SpliceAI_pred_DS_AG"),
        ds_acceptor_loss=f("SpliceAI_pred_DS_AL"),
        ds_donor_gain=f("SpliceAI_pred_DS_DG"),
        ds_donor_loss=f("SpliceAI_pred_DS_DL"),
        dp_acceptor_gain=i("SpliceAI_pred_DP_AG"),
        dp_acceptor_loss=i("SpliceAI_pred_DP_AL"),
        dp_donor_gain=i("SpliceAI_pred_DP_DG"),
        dp_donor_loss=i("SpliceAI_pred_DP_DL"),
    )
    if scores.model_dump(exclude={"symbol"}) == {
        k: None for k in scores.model_dump(exclude={"symbol"})
    }:
        return None
    return scores


def _parse_pathogenicity(csq_values, index_map) -> model.PathogenicityPredictions | None:
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
    gnomad_exomes = _parse_population_frequencies(
        csq_values, index_map, "gnomADe_AF", "gnomADe_{}_AF"
    )
    gnomad_genomes = _parse_population_frequencies(
        csq_values, index_map, "gnomADg_AF", "gnomADg_{}_AF"
    )
    all_of_us = _parse_population_frequencies(
        csq_values, index_map, "AllOfUs_gvs_all_af", "AllOfUs_gvs_{}_af"
    )
    if not any([gnomad_exomes, gnomad_genomes, all_of_us]):
        return None
    return model.AlleleFrequencies(
        gnomad_exomes=gnomad_exomes,
        gnomad_genomes=gnomad_genomes,
        all_of_us=all_of_us,
    )


def _parse_phenotype_data(csq_values, index_map) -> model.VariantPhenotypeData | None:
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


def _parse_open_targets(csq_values, index_map) -> model.OpenTargetsAssociation | None:
    diseases = _split_amp(
        _get_csq_value(csq_values, "OpenTargets_gwasDiseases", None, index_map)
    )
    gene_ids = _split_amp(
        _get_csq_value(csq_values, "OpenTargets_gwasGeneId", None, index_map)
    )
    l2g = [
        score
        for score in (
            _to_float(s)
            for s in _split_amp(
                _get_csq_value(
                    csq_values, "OpenTargets_gwasLocusToGeneScore", None, index_map
                )
            )
        )
        if score is not None
    ]
    qtl_genes = _split_amp(
        _get_csq_value(csq_values, "OpenTargets_qtlGeneId", None, index_map)
    )
    qtl_biosamples = _split_amp(
        _get_csq_value(csq_values, "OpenTargets_qtlBiosampleName", None, index_map)
    )
    if not any([diseases, gene_ids, l2g, qtl_genes, qtl_biosamples]):
        return None
    return model.OpenTargetsAssociation(
        gwas_diseases=diseases,
        gwas_gene_ids=gene_ids,
        gwas_l2g_scores=l2g,
        qtl_gene_ids=qtl_genes,
        qtl_biosamples=qtl_biosamples,
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
        predicted_molecular_consequences=consequences,
    )


def _get_vcf_meta(vcf_path: FilePath) -> model.VcfMetadata:
    """Helper method to manage metainfo for a VCF file"""

    meta_path = vcf_path.with_name(META_FILE)
    if not meta_path.exists():
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


def get_results_from_path(
    page_size: int, page: int, vcf_path: FilePath
) -> model.VepResultsResponse:
    """Returns a page of VCF data from the given filepath.
    Slices the input VCF file to a smaller one
    and converts it to stream for get_results_from_stream"""

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
    page_size: int, page: int, total: int, vcf_stream: StringIO
) -> model.VepResultsResponse:
    """Helper method to convert a filestream to VCF records for _get_results_from_vcfpy"""

    # Load vcf
    vcf_records = vcfpy.Reader.from_stream(vcf_stream)
    return _get_results_from_vcfpy(page_size, page, total, vcf_records)


def _get_results_from_vcfpy(
    page_size: int, page: int, total: int, vcf_records: vcfpy.Reader
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
    # populate variants page
    if page*page_size <= total:
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
