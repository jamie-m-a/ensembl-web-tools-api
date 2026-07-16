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
    nothing is lost if the ProtVar format differs."""
    pocket_id: str
    energy: float | None = None
    energy_per_volume: float | None = None
    score: float | None = None
    buriedness: float | None = None
    radius_of_gyration: float | None = None
    residues: list[str] = []
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
    """IntAct molecular interaction annotation (from the IntAct VEP plugin)."""
    feature_ac: str | None = None
    feature_type: str | None = None
    interaction_ac: str | None = None


class MutfuncAnnotation(BaseModel):
    """mutfunc predicted effect scores (from the mutfunc VEP plugin). Each is a
    score for one data type, or None when not available for this variant."""
    linear_motifs: float | None = None  # mutfunc_motif
    protein_interactions: float | None = None  # mutfunc_int
    protein_structure: float | None = None  # mutfunc_mod
    protein_structure_experimental: float | None = None  # mutfunc_exp


class MaveDBAnnotation(BaseModel):
    """MaveDB multiplexed assay measurement for the variant."""
    score: float | None = None
    urn: str | None = None
    doi: str | None = None
    nucleotide_variant: str | None = None  # MaveDB_nt
    protein_variant: str | None = None  # MaveDB_pro


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


# --- frequencies, phenotypes, associations (allele level) -------------------


class PopulationFrequencies(BaseModel):
    """An overall allele frequency plus per-population frequencies."""
    overall: float | None = None
    populations: dict[str, float] = {}


class AlleleFrequencies(BaseModel):
    gnomad_exomes: PopulationFrequencies | None = None
    gnomad_genomes: PopulationFrequencies | None = None
    all_of_us: PopulationFrequencies | None = None


class VariantPhenotypeData(BaseModel):
    phenotypes: list[str] = []
    clinical_significance: list[str] = []
    pubmed_ids: list[str] = []


class OpenTargetsAssociation(BaseModel):
    gwas_diseases: list[str] = []
    gwas_gene_ids: list[str] = []
    gwas_l2g_scores: list[float] = []
    qtl_gene_ids: list[str] = []
    qtl_biosamples: list[str] = []


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


class ReferenceVariantAllele(BaseModel):
    allele_sequence: str


class Location(BaseModel):
    region_name: str
    start: int
    end: int


class Metadata(BaseModel):
    pagination: PaginationMetadata


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
