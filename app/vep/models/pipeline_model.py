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
PLUGIN_PATH = "/[placeholder_path]"

# Maps an annotation option (also the boolean parameter name sent by the client)
# to the VEP `plugin ...` line it adds to the config when enabled. Order here is
# the order the lines appear in the generated ini.
PLUGIN_CONFIG_LINES: dict[str, str] = {
    "mavedb": f"plugin MaveDB,file={PLUGIN_PATH}/MaveDB_variants.tsv.gz",
    "revel": f"plugin REVEL,file={PLUGIN_PATH}/new_tabbed_revel_grch38.tsv.gz",
    "riboseqorfs": f"plugin RiboseqORFs,file={PLUGIN_PATH}/Ribo-seq_ORFs.phase2.comprehensive.v1like.final.bed.gz",
    "alphamissense": f"plugin AlphaMissense,file={PLUGIN_PATH}/AlphaMissense_hg38.tsv.gz",
    "opentargets": f"plugin OpenTargets,file={PLUGIN_PATH}/open_targets_vep.tsv.gz",
    "eve": f"plugin EVE,file={PLUGIN_PATH}/eve_merged.vcf.gz,popeve_file={PLUGIN_PATH}/grch38_popEVE_ukbb.vcf.gz",
    "cadd": f"plugin CADD,snv={PLUGIN_PATH}/CADD_GRCh38_1.7_InDels.tsv.gz",
    "spliceai": (
        f"plugin SpliceAI,snv={PLUGIN_PATH}/spliceai_scores.masked.snv.hg38.vcf.gz,"
        f"indel={PLUGIN_PATH}/spliceai_scores.masked.indel.hg38.vcf.gz,"
        f"snv_ensembl={PLUGIN_PATH}/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz"
    ),
    "protvar": f"plugin ProtVar,db={PLUGIN_PATH}/ProtVar_data.db,stability=1,pocket=1,int=1",
    "loeuf": f"plugin LOEUF,file={PLUGIN_PATH}/loeuf_dataset_grch38.tsv.gz,match_by=gene",
    "phenotypes": f"plugin Phenotypes,dir={PLUGIN_PATH}/,phenotype_feature=1,exclude_sources=COSMIC\\&HGMD-PUBLIC\\&Cancer_Gene_Census",
}


class ConfigIniParams(BaseModel):
    genome_id: str
    force_overwrite: int = 1
    transcript_version: int = 1
    canonical: int = 1
    # HGVS notations (client-selectable). `hgvs` implies HGVSc + HGVSp; `hgvsg`
    # is the genomic notation and is selected independently.
    hgvs: bool = True  # HGVSc + HGVSp (linked); on by default
    hgvsg: bool = False
    # Annotation plugins (client-selectable); each enabled flag appends its
    # `plugin ...` line (see PLUGIN_CONFIG_LINES).
    mavedb: bool = False
    revel: bool = False
    riboseqorfs: bool = False
    alphamissense: bool = False
    opentargets: bool = False
    eve: bool = False
    cadd: bool = False
    spliceai: bool = False
    protvar: bool = False
    loeuf: bool = False
    phenotypes: bool = False
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
    gff: str = ""
    fasta: str = ""

    def create_config_ini_file(self, directory):
        vep_support_location = get_vep_support_location(self.genome_id)
        self.gff = vep_support_location["gff_location"]
        self.fasta = vep_support_location["faa_location"]

        # Assembly of the selected species (e.g. "GRCh38.p14", "GRCh37.p13",
        # "GRCm39"), used to gate assembly-conditional lines and to pick
        # assembly-specific plugin data files.
        assembly_name = self.assembly_name or ""
        is_human_grch38 = assembly_name.startswith("GRCh38")
        is_human_grch37 = assembly_name.startswith("GRCh37")
        is_mouse_reference = assembly_name.startswith("GRCm39")

        # Always-on defaults (not exposed on the input form).
        lines = [
            f"force_overwrite {self.force_overwrite}",
            "numbers 1",
        ]
        # MANE annotations only exist for human GRCh38 and the mouse reference
        # (GRCm39); requesting `mane` for other species has no data.
        if is_human_grch38 or is_mouse_reference:
            lines.append("mane 1")
        # VEP assembly name, always on for the human reference assemblies.
        if is_human_grch38:
            lines.append("assembly GRCh38")
        elif is_human_grch37:
            lines.append("assembly GRCh37")
        lines += [
            "symbol 1",
            "biotype 1",
            f"transcript_version {self.transcript_version}",
            f"canonical {self.canonical}",
            f"hgvs {1 if self.hgvs else 0}",
            f"hgvsg {1 if self.hgvsg else 0}",
            f"spdi {1 if self.spdi else 0}",
            f"protein {1 if self.protein else 0}",
            f"gff {self.gff}",
            f"fasta {self.fasta}",
        ]
        # Assembly used to pick assembly-specific plugin files (default GRCh38).
        assembly = "GRCh37" if is_human_grch37 else "GRCh38"

        # Append a `plugin ...` line for every enabled plugin option.
        for option, plugin_line in PLUGIN_CONFIG_LINES.items():
            if getattr(self, option):
                lines.append(plugin_line)

        # gff3-based Genes & transcripts plugins, built from their sub-options.
        # gff3 points at the genome's resolved gff (same as the `gff` line above).
        gff3 = f"gff3={self.gff}"
        if self.tss_distance:
            direction = self.tss_distance_direction
            lines.append(
                "plugin TSSDistance,"
                f"upstream={1 if direction == 'upstream' else 0},"
                f"downstream={1 if direction == 'downstream' else 0},"
                f"both={1 if direction == 'both' else 0},"
                f"{gff3}"
            )
        if self.nearest_gene:
            both_directions = 1 if self.nearest_gene_both_directions else 0
            lines.append(
                "plugin NearestGene,"
                f"both_directions={both_directions},"
                f"{gff3}"
            )
        if self.nearest_exon_jb:
            lines.append(
                "plugin NearestExonJB,"
                f"max_range={self.nearest_exon_jb_max_range},"
                f"intronic={1 if self.nearest_exon_jb_intronic else 0},"
                f"{gff3}"
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
