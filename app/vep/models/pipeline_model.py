import logging, os
from pydantic import (
    BaseModel,
    DirectoryPath,
    FilePath,
    model_serializer,
    Field,
    AliasPath,
    field_serializer,
)
from requests import HTTPError

from core.config import (
    NF_COMPUTE_ENV_ID,
    NF_PIPELINE_URL,
)
from core.logging import InterceptHandler

from vep.models.config_spec_model import ConfigSpec
from vep.utils.config_interpreter import emit_config_lines
from vep.utils.web_metadata import get_vep_support_location

logging.getLogger().handlers = [InterceptHandler()]


class VEPConfigParams(BaseModel):
    vcf: FilePath
    vep_config: FilePath
    outdir: DirectoryPath
    bin_size: int = 3000
    sort: bool = True
    vep_version: str = "115.2"

    @model_serializer
    def vep_config_serialiser(self):
        vcf_str = f'"input": "{self.vcf.as_posix()}"'
        config_str = f'"vep_config": "{self.vep_config.as_posix()}"'
        outdir_str = f'"outdir": "{self.outdir.as_posix()}"'
        bin_str = f'"bin_size": {self.bin_size}'
        sort_str = f'"sort": {"true" if self.sort else "false"}'
        vep_version_str = f'"vep_version": "{self.vep_version}"'
        json_str = (
            "{" + ", ".join([vcf_str, config_str, outdir_str, bin_str, sort_str, vep_version_str]) + "}"
        )
        return json_str


class LaunchParams(BaseModel):
    computeEnvId: str = NF_COMPUTE_ENV_ID
    pipeline: str = NF_PIPELINE_URL
    workDir: DirectoryPath
    revision: str = "release/115"
    pullLatest: bool = True
    configProfiles: list[str] = ["ensembl"]
    paramsText: VEPConfigParams

    @field_serializer("workDir")
    def serialize_workdir(self, workdir: DirectoryPath):
        return workdir.as_posix()


class PipelineParams(BaseModel):
    launch: LaunchParams


# Placeholder root for plugin data files. These paths are not yet resolved
# per-genome; they will be substituted for real locations later.
# TODO (pre-production, required): replace this placeholder with real per-genome
# plugin-data resolution. Every `plugin ...` line below interpolates it, so no
# plugin can run against real data until this is wired up.
PLUGIN_PATH = "/[placeholder_path]"

# The option→config.ini translation — which `plugin …` / `custom …` line each
# selected option emits — now lives in the declarative config spec (the `config`
# section of specs/human_grch38.json), applied by vep.utils.config_interpreter.
# What used to be the hardcoded PLUGIN_CONFIG_LINES / PLUGIN_CONFIG_LINES_BY_ASSEMBLY
# maps and the gnomAD/AoU field builders is gone; create_config_ini_file is now a
# thin runtime over that spec plus the always-on base below.


def base_config_lines(
    *,
    assembly_name: str,
    gff: str,
    fasta: str,
    force_overwrite: int = 1,
    transcript_version: int = 1,
    canonical: int = 1,
) -> list[str]:
    """The always-on VEP config.ini lines — invocation invariants not exposed as
    options. Centralised here rather than scattered through the ini builder
    because, when the option-driven lines move to the declarative config spec,
    these stay in the backend: they are VEP invariants plus the two
    runtime-resolved paths (`gff`/`fasta`), which cannot be static spec data.
    See docs/design/merged-annotation-spec.md §4.5.

    Assembly gating mirrors ConfigIniParams' own prefix checks:
      mane 1                 — human GRCh38 and the mouse reference (GRCm39) only
      assembly               — the human reference assemblies (GRCh38 / GRCh37)
      flag_gencode_primary 1 — human GRCh38 only
    """
    is_human_grch38 = assembly_name.startswith("GRCh38")
    is_human_grch37 = assembly_name.startswith("GRCh37")
    is_mouse_reference = assembly_name.startswith("GRCm39")

    lines = [
        f"force_overwrite {force_overwrite}",
        "numbers 1",
    ]
    # MANE annotations only exist for human GRCh38 and the mouse reference
    # (GRCm39); requesting `mane` for other species has no data.
    if is_human_grch38 or is_mouse_reference:
        lines.append("mane 1")
    # VEP assembly name, always on for the human reference assemblies.
    if is_human_grch38:
        lines.append("assembly GRCh38")
        # GENCODE primary annotation flag — human GRCh38 only.
        lines.append("flag_gencode_primary 1")
    elif is_human_grch37:
        lines.append("assembly GRCh37")
    lines += [
        "symbol 1",
        "biotype 1",
        "gene_version 1",
        f"transcript_version {transcript_version}",
        f"canonical {canonical}",
        # Disable VEP's database connection (cache/plugin-file mode only). A new
        # always-on invariant, not previously emitted anywhere. See design §4.5.
        "database 0",
        f"gff {gff}",
        f"fasta {fasta}",
    ]
    return lines


class ConfigIniParams(BaseModel):
    genome_id: str
    force_overwrite: int = 1
    transcript_version: int = 1
    canonical: int = 1
    # HGVS notations (client-selectable). `hgvs` implies HGVSc + HGVSp; `hgvsg`
    # is the genomic notation and is selected independently.
    hgvs: bool = False  # HGVSc + HGVSp (linked); off by default
    hgvsg: bool = False
    # Annotation plugins (client-selectable); each enabled flag appends its
    # `plugin ...` line via the config spec (config_specs → config_interpreter).
    mavedb: bool = False
    revel: bool = False
    riboseqorfs: bool = False
    alphamissense: bool = False
    opentargets: bool = False
    eve: bool = False
    cadd: bool = False
    spliceai: bool = False
    protvar: bool = False
    # ProtVar sub-features (only used when protvar is on); default all on.
    protvar_stability: bool = True
    protvar_pocket: bool = True
    protvar_int: bool = True
    loeuf: bool = False
    phenotypes: bool = False
    # Gene Ontology terms (GO plugin; human GRCh38 for now).
    go: bool = False
    # SPDI variant notation (flag line, like hgvs/hgvsg).
    spdi: bool = False
    # Protein ID (VEP --protein, adds the Ensembl protein id). Flag line.
    protein: bool = False
    # Distance to TSS (TSSDistance plugin); direction radio defaults to upstream.
    tss_distance: bool = False
    tss_distance_direction: str = "upstream"  # upstream | downstream | both
    # Nearest gene (NearestGene plugin).
    nearest_gene: bool = False
    nearest_gene_both_directions: bool = False
    # Nearest exon junction boundary (NearestExonJB plugin).
    nearest_exon_jb: bool = False
    nearest_exon_jb_max_range: int = 10000  # max search range (bp)
    nearest_exon_jb_intronic: bool = False
    # Geno2MP variant associations (assembly-specific; see *_BY_ASSEMBLY).
    geno2mp: bool = False
    # UTRAnnotator 5' UTR variants (assembly-specific; see *_BY_ASSEMBLY).
    utrannotator: bool = False
    # NMD escape prediction (NMD plugin; no params). Human GRCh37/38.
    nmd: bool = False
    # Dosage sensitivity (DosageSensitivity plugin); `cover` sub-flag.
    dosage_sensitivity: bool = False
    dosage_sensitivity_cover: bool = False
    # IntAct molecular interactions (human GRCh38). Sub-flags default off.
    intact: bool = False
    intact_feature_ac: bool = False
    intact_feature_short_label: bool = False
    intact_feature_annotation: bool = False
    intact_ap_ac: bool = False
    intact_interaction_participants: bool = False
    intact_pmid: bool = False
    # mutfunc (human GRCh38). Sub-flags default off.
    mutfunc: bool = False
    mutfunc_motif: bool = False
    mutfunc_int: bool = False
    mutfunc_mod: bool = False
    mutfunc_exp: bool = False
    # gnomAD exomes v4.1 frequencies (human GRCh38). Rendered as a VEP `custom`
    # line whose `fields` list is built from the selected genetic-ancestry groups
    # x sexes (and whether UK Biobank samples are included). Field grammar:
    # AF[_non_ukb][_<ancestry>][_XX|_XY] (XX=female, XY=male; base = both sexes).
    gnomad_exomes: bool = False
    gnomad_exomes_include_ukb: bool = True  # False -> the _non_ukb subset fields
    # Ancestry toggles ("all" pre-selected so enabling the option yields fields=AF).
    gnomad_exomes_all: bool = True
    gnomad_exomes_afr: bool = False
    gnomad_exomes_amr: bool = False
    gnomad_exomes_asj: bool = False
    gnomad_exomes_eas: bool = False
    gnomad_exomes_fin: bool = False
    gnomad_exomes_mid: bool = False
    gnomad_exomes_nfe: bool = False
    # Per-ancestry sex sub-options (Both defaults on = the base combined-sex field).
    gnomad_exomes_all_both: bool = True
    gnomad_exomes_all_female: bool = False
    gnomad_exomes_all_male: bool = False
    gnomad_exomes_afr_both: bool = True
    gnomad_exomes_afr_female: bool = False
    gnomad_exomes_afr_male: bool = False
    gnomad_exomes_amr_both: bool = True
    gnomad_exomes_amr_female: bool = False
    gnomad_exomes_amr_male: bool = False
    gnomad_exomes_asj_both: bool = True
    gnomad_exomes_asj_female: bool = False
    gnomad_exomes_asj_male: bool = False
    gnomad_exomes_eas_both: bool = True
    gnomad_exomes_eas_female: bool = False
    gnomad_exomes_eas_male: bool = False
    gnomad_exomes_fin_both: bool = True
    gnomad_exomes_fin_female: bool = False
    gnomad_exomes_fin_male: bool = False
    gnomad_exomes_mid_both: bool = True
    gnomad_exomes_mid_female: bool = False
    gnomad_exomes_mid_male: bool = False
    gnomad_exomes_nfe_both: bool = True
    gnomad_exomes_nfe_female: bool = False
    gnomad_exomes_nfe_male: bool = False
    # gnomAD genomes v4.1 frequencies (human GRCh38). Same custom-line grammar as
    # exomes but with no UK Biobank subset (no _non_ukb), extra ancestry groups
    # (Amish, Remaining) and grpmax (max across groups; no XX/XY split).
    gnomad_genomes: bool = False
    # Ancestry toggles ("all" pre-selected so enabling yields fields=AF).
    gnomad_genomes_all: bool = True
    gnomad_genomes_afr: bool = False
    gnomad_genomes_amr: bool = False
    gnomad_genomes_asj: bool = False
    gnomad_genomes_eas: bool = False
    gnomad_genomes_fin: bool = False
    gnomad_genomes_mid: bool = False
    gnomad_genomes_nfe: bool = False
    gnomad_genomes_ami: bool = False
    gnomad_genomes_remaining: bool = False
    gnomad_genomes_grpmax: bool = False  # no sex sub-options (no XX/XY)
    # Per-ancestry sex sub-options (Both on = base combined-sex field). grpmax has
    # none.
    gnomad_genomes_all_both: bool = True
    gnomad_genomes_all_female: bool = False
    gnomad_genomes_all_male: bool = False
    gnomad_genomes_afr_both: bool = True
    gnomad_genomes_afr_female: bool = False
    gnomad_genomes_afr_male: bool = False
    gnomad_genomes_amr_both: bool = True
    gnomad_genomes_amr_female: bool = False
    gnomad_genomes_amr_male: bool = False
    gnomad_genomes_asj_both: bool = True
    gnomad_genomes_asj_female: bool = False
    gnomad_genomes_asj_male: bool = False
    gnomad_genomes_eas_both: bool = True
    gnomad_genomes_eas_female: bool = False
    gnomad_genomes_eas_male: bool = False
    gnomad_genomes_fin_both: bool = True
    gnomad_genomes_fin_female: bool = False
    gnomad_genomes_fin_male: bool = False
    gnomad_genomes_mid_both: bool = True
    gnomad_genomes_mid_female: bool = False
    gnomad_genomes_mid_male: bool = False
    gnomad_genomes_nfe_both: bool = True
    gnomad_genomes_nfe_female: bool = False
    gnomad_genomes_nfe_male: bool = False
    gnomad_genomes_ami_both: bool = True
    gnomad_genomes_ami_female: bool = False
    gnomad_genomes_ami_male: bool = False
    gnomad_genomes_remaining_both: bool = True
    gnomad_genomes_remaining_female: bool = False
    gnomad_genomes_remaining_male: bool = False
    # NIH All of Us frequencies (human GRCh38). A VEP `custom` line whose `fields`
    # list is built from the selected population toggles (no sex split).
    allofus: bool = False
    allofus_all: bool = True  # pre-selected so enabling yields fields=gvs_all_af
    allofus_max: bool = False  # -> gvs_max_af + gvs_max_subpop
    allofus_afr: bool = False
    allofus_amr: bool = False
    allofus_eas: bool = False
    allofus_eur: bool = False
    allofus_mid: bool = False
    allofus_sas: bool = False
    allofus_oth: bool = False
    # ClinVar clinical significance (human GRCh37/GRCh38). A VEP `custom` line
    # surfacing the CLNSIG field; not assembly-specific.
    clinvar: bool = False
    # Assembly name (from the selected species, e.g. "GRCh38"/"GRCh37"); used to
    # pick assembly-specific plugin data files. Defaults to GRCh38.
    assembly_name: str = ""
    # Species taxonomy id of the selected species (e.g. "9606" for human). Sent
    # by the client alongside `assembly_name`, and used only to compute the
    # option panels pinned to this job — the same species/assembly predicates
    # /form_config uses (form_panels.is_human_grch37_or_38 / is_human_grch38).
    # Without it every human-specific panel would be silently dropped from the
    # pin. It emits no config.ini line.
    species_taxonomy_id: str = ""
    gff: str = ""
    fasta: str = ""

    def create_config_ini_file(self, directory, config_spec: ConfigSpec):
        """Write the VEP config.ini for this submission: the always-on base
        (`base_config_lines`) plus the option-driven `plugin …` / `custom …` /
        flag lines the config interpreter emits from `config_spec` — the
        `.config` half of the job's pinned merged spec. A thin runtime over the
        declarative spec; the option→line rules are data, not code here. See
        docs/design/merged-annotation-spec.md."""
        vep_support_location = get_vep_support_location(self.genome_id)
        self.gff = vep_support_location["gff_location"]
        self.fasta = vep_support_location["faa_location"]

        # Assembly of the selected species (e.g. "GRCh38.p14", "GRCh37.p13"),
        # resolved to the value the spec's `by_assembly` params key on (default
        # GRCh38) — the same notion of "which genome" the base config uses.
        assembly_name = self.assembly_name or ""
        assembly = "GRCh37" if assembly_name.startswith("GRCh37") else "GRCh38"

        lines = base_config_lines(
            assembly_name=assembly_name,
            gff=self.gff,
            fasta=self.fasta,
            force_overwrite=self.force_overwrite,
            transcript_version=self.transcript_version,
            canonical=self.canonical,
        )
        lines += emit_config_lines(
            config_spec,
            self.model_dump(),
            assembly=assembly,
            plugin_path=PLUGIN_PATH,
            gff=self.gff,
        )

        config_ini = "\n".join(lines) + "\n"
        try:
            with open(os.path.join(directory, "config.ini"), "w") as ini_file:
                ini_file.write(config_ini)
            return ini_file
        except Exception as e:
            raise RuntimeError(f"Could not create VEP config ini file: {e}")


class PipelineStatus(BaseModel):
    submission_id: str
    status: str = Field(
        validation_alias=AliasPath("status", "workflow", "status"), default="FAILED"
    )

    @field_serializer("status")
    def serialize_status(self, status: str):
        if status == "UNKNOWN":
            status = "FAILED"
            logging.info(
                f"Unknown status was returned for submission {self.submission_id}"
            )
        return status
