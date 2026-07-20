from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, Field

from vep.models.display_panels_model import DisplayPanel

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


class Annotation(BaseModel):
    """A generic, spec-driven plugin annotation: the payload produced by
    `spec_interpreter.apply_plugin_spec` for one plugin, tagged with its plugin
    id and scope. This is the only source of annotation data on the wire — the
    envelope (variant/allele/consequence) stays typed; the annotations
    themselves are generic (`data`). A plugin absent from the list means "did
    not run / no data"."""

    plugin: str  # spec plugin id, e.g. "mavedb", "gnomad_exomes"
    scope: str  # "allele" | "transcript"
    data: dict[str, Any]


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


class PredictionWithScore(BaseModel):
    """A categorical prediction plus its score, e.g. SIFT 'tolerated(0.15)'."""
    prediction: str | None = None
    score: float | None = None


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
    #
    # TODO: uniprot / protein_matches / sift / polyphen are the unspecced tail
    # of the go-flat cutover — deliberately left typed for now, still to be
    # converted to plugin specs. No sample data carries their columns, so no
    # spec could be validated for them. Everything else moved to `annotations`.
    ensembl_protein_id: str | None = None
    uniprot: UniprotIds | None = None
    protein_matches: list[ProteinMatch] = []
    sift: PredictionWithScore | None = None
    polyphen: PredictionWithScore | None = None
    # Generic spec-driven annotations for this transcript consequence (scope
    # "transcript"). The only source of annotation data beside the tail above.
    annotations: list[Annotation] = []


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
    # The option panels this job was submitted against, pinned at submission.
    # None for jobs that predate the pin — the results view then falls back to
    # the live form-config panels, as it did before.
    display_panels: list[DisplayPanel] | None = None


class AlternativeVariantAllele(BaseModel):
    allele_sequence: str
    allele_type: str
    # TODO: colocated_variants is part of the unspecced tail of the go-flat
    # cutover — deliberately left typed for now, still to be converted to a
    # plugin spec (no sample data carries its column, so no spec could be
    # validated for it).
    colocated_variants: list[str] = []  # Existing_variation
    # Generic spec-driven annotations for this allele (scope "allele"). Same
    # across all of this allele's transcripts, and the only annotations
    # available for intergenic variants (which have no transcript rows).
    annotations: list[Annotation] = []
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
