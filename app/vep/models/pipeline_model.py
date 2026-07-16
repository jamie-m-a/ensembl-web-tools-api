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
    "riboseqorfs": f"plugin RiboseqORFs,file={PLUGIN_PATH}/Ribo-seq_ORFs.phase2.comprehensive.v1like.final.bed.gz",
    "opentargets": f"plugin OpenTargets,file={PLUGIN_PATH}/open_targets_vep.tsv.gz",
    "eve": f"plugin EVE,file={PLUGIN_PATH}/eve_merged.vcf.gz,popeve_file={PLUGIN_PATH}/grch38_popEVE_ukbb.vcf.gz",
    "maxentscan": f"plugin MaxEntScan,file={PLUGIN_PATH}/",
    "gnomad_mt": f"plugin gnomADMt,file={PLUGIN_PATH}/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
    # ProtVar's stability/pocket/int sub-flags are filled in from the client's
    # sub-options at build time (see create_config_ini_file).
    "protvar": (
        f"plugin ProtVar,db={PLUGIN_PATH}/ProtVar_data.db,"
        "stability={stability},pocket={pocket},int={int}"
    ),
    # IntAct base line; the selected sub-option flags (or `all=1`) are appended
    # at build time (see create_config_ini_file).
    "intact": (
        f"plugin IntAct,mutation_file={PLUGIN_PATH}/mutations.tsv,"
        f"mapping_file={PLUGIN_PATH}/mutation_gc_map.txt.gz"
    ),
    # mutfunc sub-flags filled in at build time; extended is always on.
    "mutfunc": (
        "plugin mutfunc,motif={motif},int={int},mod={mod},exp={exp},extended=1,"
        f"db={PLUGIN_PATH}/mutfunc_data.db"
    ),
    # DosageSensitivity's `cover` sub-flag is filled in at build time.
    "dosage_sensitivity": (
        f"plugin DosageSensitivity,"
        f"file={PLUGIN_PATH}/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz,"
        "cover={cover}"
    ),
    "phenotypes": f"plugin Phenotypes,dir={PLUGIN_PATH}/,phenotype_feature=1,exclude_sources=COSMIC\\&HGMD-PUBLIC\\&Cancer_Gene_Census",
}

# Assembly-specific plugin lines: option -> {assembly: line}. Resolved against
# the submission's assembly (GRCh37/GRCh38) at build time.
PLUGIN_CONFIG_LINES_BY_ASSEMBLY: dict[str, dict[str, str]] = {
    "alphamissense": {
        "GRCh38": f"plugin AlphaMissense,file={PLUGIN_PATH}/AlphaMissense_hg38.tsv.gz",
        "GRCh37": f"plugin AlphaMissense,file={PLUGIN_PATH}/AlphaMissense_hg19.tsv.gz",
    },
    "cadd": {
        "GRCh38": (
            f"plugin CADD,snv={PLUGIN_PATH}/CADD_GRCh38_1.7_whole_genome_SNVs.tsv.gz,"
            f"indels={PLUGIN_PATH}/CADD_GRCh38_1.7_InDels.tsv.gz"
        ),
        "GRCh37": (
            f"plugin CADD,snv={PLUGIN_PATH}/CADD_GRCh37_1.7_whole_genome_SNVs.tsv.gz,"
            f"indels={PLUGIN_PATH}/CADD_GRCh37_1.7_InDels.tsv.gz"
        ),
    },
    "revel": {
        "GRCh38": f"plugin REVEL,file={PLUGIN_PATH}/new_tabbed_revel_grch38.tsv.gz",
        "GRCh37": f"plugin REVEL,file={PLUGIN_PATH}/new_tabbed_revel_grch37.tsv.gz",
    },
    "loeuf": {
        "GRCh38": f"plugin LOEUF,file={PLUGIN_PATH}/loeuf_dataset_grch38.tsv.gz,match_by=gene",
        "GRCh37": f"plugin LOEUF,file={PLUGIN_PATH}/loeuf_dataset_grch37.tsv.gz,match_by=gene",
    },
    "geno2mp": {
        "GRCh38": f"plugin Geno2MP,file={PLUGIN_PATH}/Geno2MP.variants_GRCh38.vcf.gz",
        "GRCh37": f"plugin Geno2MP,file={PLUGIN_PATH}/Geno2MP.variants_GRCh37.vcf.gz",
    },
    "enformer": {
        "GRCh38": f"plugin Enformer,file={PLUGIN_PATH}/enformer_grch38.vcf.gz",
        "GRCh37": f"plugin Enformer,file={PLUGIN_PATH}/enformer_grch37.vcf.gz",
    },
    "spliceai": {
        "GRCh38": (
            f"plugin SpliceAI,snv={PLUGIN_PATH}/spliceai_scores.masked.snv.hg38.vcf.gz,"
            f"indel={PLUGIN_PATH}/spliceai_scores.masked.indel.hg38.vcf.gz,"
            f"snv_ensembl={PLUGIN_PATH}/spliceai_scores.raw.snv.ensembl_mane.grch38.110.vcf.gz"
        ),
        "GRCh37": (
            f"plugin SpliceAI,snv={PLUGIN_PATH}/spliceai_scores.masked.snv.hg19.vcf.gz,"
            f"indel={PLUGIN_PATH}/spliceai_scores.masked.indel.hg19.vcf.gz"
        ),
    },
    "utrannotator": {
        "GRCh38": f"plugin UTRAnnotator,file={PLUGIN_PATH}/uORF_5UTR_GRCh38_PUBLIC.txt",
        "GRCh37": f"plugin UTRAnnotator,file={PLUGIN_PATH}/uORF_5UTR_GRCh37_PUBLIC.txt",
    },
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
    # ProtVar sub-features (only used when protvar is on); default all on.
    protvar_stability: bool = True
    protvar_pocket: bool = True
    protvar_int: bool = True
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
    # MaxEntScan splicing (MaxEntScan plugin).
    maxentscan: bool = False
    # Geno2MP variant associations (assembly-specific; see *_BY_ASSEMBLY).
    geno2mp: bool = False
    # Enformer non-coding predictions (assembly-specific; see *_BY_ASSEMBLY).
    enformer: bool = False
    # UTRAnnotator 5' UTR variants (assembly-specific; see *_BY_ASSEMBLY).
    utrannotator: bool = False
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
    # gnomAD mitochondrial frequencies (human GRCh38).
    gnomad_mt: bool = False
    # Assembly name (from the selected species, e.g. "GRCh38"/"GRCh37"); used to
    # pick assembly-specific plugin data files. Defaults to GRCh38.
    assembly_name: str = ""
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
            if not getattr(self, option):
                continue
            if option == "protvar":
                plugin_line = plugin_line.format(
                    stability=int(self.protvar_stability),
                    pocket=int(self.protvar_pocket),
                    int=int(self.protvar_int),
                )
            elif option == "dosage_sensitivity":
                plugin_line = plugin_line.format(
                    cover=int(self.dosage_sensitivity_cover)
                )
            elif option == "intact":
                # Append the selected IntAct sub-option flags: none -> base line
                # only; all selected -> `all=1`; otherwise just the chosen ones.
                intact_flags = {
                    "feature_ac": self.intact_feature_ac,
                    "feature_short_label": self.intact_feature_short_label,
                    "feature_annotation": self.intact_feature_annotation,
                    "ap_ac": self.intact_ap_ac,
                    "interaction_participants": self.intact_interaction_participants,
                    "pmid": self.intact_pmid,
                }
                selected = [name for name, on in intact_flags.items() if on]
                if len(selected) == len(intact_flags):
                    plugin_line = f"{plugin_line},all=1"
                elif selected:
                    plugin_line = plugin_line + "," + ",".join(
                        f"{name}=1" for name in selected
                    )
            elif option == "mutfunc":
                plugin_line = plugin_line.format(
                    motif=int(self.mutfunc_motif),
                    int=int(self.mutfunc_int),
                    mod=int(self.mutfunc_mod),
                    exp=int(self.mutfunc_exp),
                )
            lines.append(plugin_line)

        # Assembly-specific plugins (GRCh37/GRCh38 data files).
        for option, by_assembly in PLUGIN_CONFIG_LINES_BY_ASSEMBLY.items():
            if getattr(self, option):
                lines.append(by_assembly.get(assembly, by_assembly["GRCh38"]))

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
