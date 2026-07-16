from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, Field

class VcfMetadata(BaseModel):
    variant_count: int = Field(
        description="Total number of variant records in the VCF file"
    )
    header_count: int = Field(
        description="Number of header rows in the VCF file"
    )

class PaginationMetadata(BaseModel):
    page: int
    per_page: int
    total: int


class PredictedIntergenicConsequence(BaseModel):
    feature_type: Any | None = Field(
        default=None,
        description="The value of this field is always null. The presence of null in this field will serve as a marker that this is a consequence of an intergenic variant.",
    )
    consequences: list[str] = Field(
        default=['intergenic_variant'],
        description="The only expected member of this array will be the string 'intergenic_variant'",
    )


class FeatureType(Enum):
    transcript = "transcript"


class Strand(Enum):
    forward = "forward"
    reverse = "reverse"


# --- Protein & functional annotations (consequence level) -------------------


class UniprotIds(BaseModel):
    """Uniprot cross-references for the protein product (from VEP --uniprot)."""
    swissprot: str | None = None
    trembl: str | None = None
    uniparc: str | None = None
    isoform: str | None = None


class ProteinMatch(BaseModel):
    """A protein structure/domain match (from the VEP DOMAINS field), e.g.
    AlphaFold-DB or PDB mappings."""
    source: str = Field(..., description="e.g. 'AFDB-ENSP_mappings'")
    id: str = Field(..., description="e.g. 'AF-Q9UGM6-F1'")


class ProtVarPocket(BaseModel):
    """A ProtVar predicted ligand-binding pocket.

    NOTE: field order is parsed positionally from the ProtVar_pocket CSQ value
    (id & energy & energy_per_volume & score & buriedness & radius_of_gyration &
    residues). The names are best-effort; `raw` preserves the original value so
    nothing is lost if the ProtVar format differs. Residues are not captured."""
    pocket_id: str
    energy: float | None = None
    energy_per_volume: float | None = None
    score: float | None = None
    buriedness: float | None = None
    radius_of_gyration: float | None = None
    raw: str


class ProtVarInteractionInterface(BaseModel):
    """A ProtVar protein-protein interaction interface (partner & score)."""
    partner: str
    score: float | None = None
    raw: str


class ProtVarAnnotation(BaseModel):
    structure_stability_score: float | None = Field(
        default=None, description="ProtVar_stability"
    )
    pockets: list[ProtVarPocket] = []
    interaction_interfaces: list[ProtVarInteractionInterface] = []


class IntActAnnotation(BaseModel):
    """IntAct molecular interaction annotation (from the IntAct VEP plugin).
    Besides the base feature_type / interaction_ac, each field corresponds to a
    selected IntAct sub-option (emitted only when its flag was set)."""
    feature_type: str | None = None
    interaction_ac: str | None = None
    feature_ac: str | None = None
    feature_short_label: str | None = None
    feature_annotation: str | None = None
    ap_ac: str | None = None  # affected protein AC
    interaction_participants: str | None = None
    pmid: str | None = None


class MutfuncAnnotation(BaseModel):
    """mutfunc predicted effect scores (from the mutfunc VEP plugin). Each is a
    score for one data type, or None when not available for this variant."""
    linear_motifs: float | None = None  # mutfunc_motif
    protein_interactions: float | None = None  # mutfunc_int
    protein_structure: float | None = None  # mutfunc_mod
    protein_structure_experimental: float | None = None  # mutfunc_exp


class MaveDBAssay(BaseModel):
    """One MaveDB multiplexed-assay measurement: a score-set URN and its score."""
    urn: str | None = None
    score: float | None = None


class MaveDBAnnotation(BaseModel):
    """MaveDB multiplexed-assay measurements for the variant. The plugin reports
    several assays as parallel &-joined columns; each (urn, score) pair is one
    assay. `protein_variant` is the (shared) protein change."""
    protein_variant: str | None = None  # MaveDB_pro
    assays: list[MaveDBAssay] = []


# --- HGVS, pathogenicity, conservation (consequence level) ------------------


class HgvsNotations(BaseModel):
    genomic: str | None = None  # HGVSg
    transcript: str | None = None  # HGVSc
    protein: str | None = None  # HGVSp


class PredictionWithScore(BaseModel):
    """A categorical prediction plus its score, e.g. SIFT 'tolerated(0.15)'."""
    prediction: str | None = None
    score: float | None = None


class SpliceAiScores(BaseModel):
    symbol: str | None = None
    ds_acceptor_gain: float | None = None  # DS_AG
    ds_acceptor_loss: float | None = None  # DS_AL
    ds_donor_gain: float | None = None  # DS_DG
    ds_donor_loss: float | None = None  # DS_DL
    dp_acceptor_gain: int | None = None  # DP_AG
    dp_acceptor_loss: int | None = None  # DP_AL
    dp_donor_gain: int | None = None  # DP_DG
    dp_donor_loss: int | None = None  # DP_DL


class PopEve(BaseModel):
    """popEVE scores (population-adjusted EVE/ESM1v), from the EVE plugin's
    popeve_file. One prediction per protein-altering variant."""
    score: float | None = None  # popEVE_SCORE
    eve: float | None = None  # popEVE_EVE
    esm1v: float | None = None  # popEVE_ESM1v
    pop_adjusted_eve: float | None = None  # popEVE_pop_adjusted_EVE
    pop_adjusted_esm1v: float | None = None  # popEVE_pop_adjusted_ESM1v
    gene: str | None = None
    protein: str | None = None
    mutant: str | None = None
    gap_frequency: float | None = None


class PathogenicityPredictions(BaseModel):
    sift: PredictionWithScore | None = None
    polyphen: PredictionWithScore | None = None
    revel: float | None = None
    alphamissense_class: str | None = None
    alphamissense_score: float | None = None
    cadd_phred: float | None = None
    cadd_raw: float | None = None
    spliceai: SpliceAiScores | None = None
    eve_class: str | None = None
    eve_score: float | None = None
    popeve: PopEve | None = None


# --- gene constraint, UTR, Ribo-seq (consequence level) ---------------------


class DosageSensitivity(BaseModel):
    """gnomAD dosage sensitivity probabilities (gene-level, shown per
    transcript). From the DosageSensitivity plugin."""
    phaplo: float | None = None  # pHaplo: haploinsufficiency probability
    ptriplo: float | None = None  # pTriplo: triplosensitivity probability


class FivePrimeUtrAnnotation(BaseModel):
    """UTRAnnotator 5' UTR consequence for variants affecting upstream ORFs."""
    consequence: str | None = None  # 5UTR_consequence
    annotation: str | None = None  # 5UTR_annotation (raw detail string)
    existing_uorfs: str | None = None  # Existing_uORFs
    existing_inframe_oorfs: str | None = None  # Existing_InFrame_oORFs
    existing_outofframe_oorfs: str | None = None  # Existing_OutOfFrame_oORFs


class RiboseqOrfsAnnotation(BaseModel):
    """Ribo-seq ORFs plugin: overlap with translated ORFs identified by
    Ribo-seq (Ensembl/GENCODE phase 2)."""
    orf_id: str | None = None  # RiboseqORFs_id
    consequences: list[str] = []  # RiboseqORFs_consequences (&-split)
    impact: str | None = None  # RiboseqORFs_impact
    protein_position: str | None = None  # RiboseqORFs_protein_position
    codons: str | None = None  # RiboseqORFs_codons
    amino_acids: str | None = None  # RiboseqORFs_amino_acids
    publications: list[str] = []  # RiboseqORFs_publications (&-split)


class GoTerm(BaseModel):
    """A Gene Ontology annotation for a transcript (from the GO plugin's GO
    column): the GO id plus its (human-readable) term name."""
    id: str  # e.g. GO:0001558
    name: str  # e.g. "regulation of cell growth"


# --- frequencies, phenotypes, associations (allele level) -------------------


class PopulationFrequencies(BaseModel):
    """An overall allele frequency plus per-population frequencies."""
    overall: float | None = None
    populations: dict[str, float] = {}
    # For All of Us: the subpopulation code the "max" frequency came from
    # (e.g. "eur"), from the AoU_gvs_max_subpop label column. None otherwise.
    max_subpopulation: str | None = None


class AlleleFrequencies(BaseModel):
    gnomad_exomes: PopulationFrequencies | None = None
    gnomad_genomes: PopulationFrequencies | None = None
    all_of_us: PopulationFrequencies | None = None


class VariantPhenotypeData(BaseModel):
    phenotypes: list[str] = []
    clinical_significance: list[str] = []
    pubmed_ids: list[str] = []


class OpenTargetsGwasAssociation(BaseModel):
    disease: str  # EFO ontology id
    gene_id: str  # Ensembl gene id
    l2g_score: float | None = None  # locus-to-gene confidence


class OpenTargetsQtlAssociation(BaseModel):
    gene_id: str  # Ensembl gene id
    biosample: str | None = None  # affected tissue


class OpenTargetsAssociation(BaseModel):
    # Structured, de-duplicated associations. The parallel CSQ subfields are
    # positionally aligned, so they are zipped together at parse time.
    gwas_associations: list[OpenTargetsGwasAssociation] = []
    qtl_associations: list[OpenTargetsQtlAssociation] = []


class PredictedTranscriptConsequence(BaseModel):
    feature_type: FeatureType
    stable_id: str = Field(..., description="transcript stable id, versioned")
    gene_stable_id: str = Field(..., description="gene stable id, versioned")
    gene_symbol: str | None = None
    biotype: str
    is_canonical: bool
    consequences: list[str]
    strand: Strand
    # MANE (human GRCh38 only): is_mane_select flags the MANE Select transcript;
    # mane_select_refseq_id is its matched RefSeq id (e.g. NM_001242672.3).
    is_mane_select: bool = False
    is_mane_plus_clinical: bool = False
    mane_select_refseq_id: str | None = None
    # Protein & functional annotations (optional; populated when the relevant
    # VEP options/plugins were enabled for the run).
    ensembl_protein_id: str | None = None
    uniprot: UniprotIds | None = None
    protein_matches: list[ProteinMatch] = []
    protvar: ProtVarAnnotation | None = None
    intact: IntActAnnotation | None = None
    mutfunc: MutfuncAnnotation | None = None
    mavedb: MaveDBAnnotation | None = None
    # Variant representations
    spdi: str | None = None  # SPDI (e.g. "1:11021:G:A")
    # HGVS, pathogenicity and gene constraint
    hgvs: HgvsNotations | None = None
    pathogenicity: PathogenicityPredictions | None = None
    loeuf: float | None = None  # gnomAD LOEUF (gene-level, shown per transcript)
    dosage_sensitivity: DosageSensitivity | None = None
    # Transcript-level predictions (populated when the relevant plugins ran).
    utr_annotation: FivePrimeUtrAnnotation | None = None
    riboseq_orfs: RiboseqOrfsAnnotation | None = None
    go_terms: list[GoTerm] = []


class ReferenceVariantAllele(BaseModel):
    allele_sequence: str


class Location(BaseModel):
    region_name: str
    start: int
    end: int


class FilterStat(BaseModel):
    """How many records a single active filter removed (among those that reached
    it in the pipeline). Captured so the filter ordering can be tuned later."""

    field: str
    removed: int


class FilterMetadata(BaseModel):
    """Summary of server-side filtering, present only when filters were applied."""

    unfiltered_total: int = Field(
        description="Total records before any filtering"
    )
    filtered_total: int = Field(
        description="Records remaining after all filters (the paginated total)"
    )
    stats: list[FilterStat] = Field(
        description="Per-filter removed counts, in the order the pipeline ran"
    )


class AfSource(BaseModel):
    """An allele-frequency column available to filter on (i.e. an AF option that
    was selected at input). `population` is empty for the source's overall AF."""

    key: str = Field(description="CSQ column name, e.g. gnomAD_exomes_AF_nfe")
    source: str = Field(
        description="gnomad_exomes | gnomad_genomes | all_of_us"
    )
    population: str = Field(description="Population code, or '' for overall")


class Metadata(BaseModel):
    pagination: PaginationMetadata
    filters: FilterMetadata | None = None
    # AF columns present in this result set (the AF options chosen at input),
    # so the frontend can populate the allele-frequency filter.
    available_af_sources: list[AfSource] = []


class ClinVarSignificance(BaseModel):
    """A ClinVar clinical significance term and how many submissions reported it
    (from the CLNSIGCONF breakdown, e.g. 'Pathogenic_(10)')."""
    significance: str
    count: int


class ClinVarAnnotation(BaseModel):
    """ClinVar clinical significance (CLNSIG). When the overall classification is
    'Conflicting classifications of pathogenicity', `conflicting_breakdown` holds
    the per-classification submission counts (from CLNSIGCONF); it is empty
    otherwise (CLNSIGCONF is ignored for non-conflicting classifications)."""
    significance: list[str]
    conflicting_breakdown: list[ClinVarSignificance] = []


class AlternativeVariantAllele(BaseModel):
    allele_sequence: str
    allele_type: str
    representative_population_allele_frequency: float | None = None
    # Allele-level annotations (same across all transcripts of this allele, and
    # the only annotations available for intergenic variants).
    spdi: str | None = None  # SPDI (e.g. "1:79106:T:C")
    hgvsg: str | None = None  # genomic HGVS (e.g. "1:g.79107T>C")
    cadd_phred: float | None = None  # CADD is genome-wide (per position)
    cadd_raw: float | None = None
    frequencies: AlleleFrequencies | None = None
    colocated_variants: list[str] = []  # Existing_variation
    phenotype_data: VariantPhenotypeData | None = None
    open_targets: OpenTargetsAssociation | None = None
    # ClinVar clinical significance (from the ClinVar custom track).
    clinvar: ClinVarAnnotation | None = None
    predicted_molecular_consequences: list[
        PredictedTranscriptConsequence | PredictedIntergenicConsequence
    ]


class Variant(BaseModel):
    name: str | None = Field(
        default=None,
        description="User's name for the variant; optional"
    )
    allele_type: str
    location: Location
    reference_allele: ReferenceVariantAllele
    alternative_alleles: list[AlternativeVariantAllele]


class VepResultsResponse(BaseModel):
    metadata: Metadata
    variants: list[Variant]
