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
    # Phenotypes: species/assembly-specific data file. Only the human GRCh38 file
    # exists for now; other species files follow (then move to a species-keyed
    # map, like PLUGIN_CONFIG_LINES_BY_ASSEMBLY). The form only offers this option
    # for human GRCh38, so the single GRCh38 line is safe here.
    "phenotypes": f"plugin Phenotypes,file={PLUGIN_PATH}/Phenotypes.pm_homo_sapiens_116_GRCh38.gvf.gz",
    # Gene Ontology: species/assembly-specific data file, same GRCh38-only caveat
    # as phenotypes above (see phenotypes-species-todo.md).
    "go": f"plugin GO,file={PLUGIN_PATH}/GO.pm_homo_sapiens_116_GRCh38.gff.gz",
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

# gnomAD exomes v4.1 custom-line field building (human GRCh38). The `fields`
# value is assembled from the selected ancestry x sex combinations.
# Ancestry param suffix -> field code component ("" = the all-ancestries field).
GNOMAD_EXOMES_ANCESTRIES: list[tuple[str, str]] = [
    ("all", ""),
    ("afr", "afr"),
    ("amr", "amr"),
    ("asj", "asj"),
    ("eas", "eas"),
    ("fin", "fin"),
    ("mid", "mid"),
    ("nfe", "nfe"),
]
# gnomAD genomes v4.1: the same sex-split ancestries as exomes plus Amish and
# Remaining ("grpmax" — max across groups — is handled separately: no XX/XY
# split; and genomes has no UK Biobank subset, so no _non_ukb fields).
GNOMAD_GENOMES_ANCESTRIES: list[tuple[str, str]] = GNOMAD_EXOMES_ANCESTRIES + [
    ("ami", "ami"),
    ("remaining", "remaining"),
]

# Sex param suffix -> field code component ("" = both/combined sexes;
# XX = female, XY = male). Shared by gnomAD exomes and genomes.
GNOMAD_SEXES: list[tuple[str, str]] = [
    ("both", ""),
    ("female", "XX"),
    ("male", "XY"),
]

# NIH All of Us (AoU) population toggles (human GRCh38). Param suffix -> the
# custom `fields` code(s) it contributes. "Maximum subpopulation" uniquely emits
# two fields (the max AF and the subpopulation label).
ALLOFUS_POPULATIONS: list[tuple[str, list[str]]] = [
    ("all", ["gvs_all_af"]),
    ("max", ["gvs_max_af", "gvs_max_subpop"]),
    ("afr", ["gvs_afr_af"]),
    ("amr", ["gvs_amr_af"]),
    ("eas", ["gvs_eas_af"]),
    ("eur", ["gvs_eur_af"]),
    ("mid", ["gvs_mid_af"]),
    ("sas", ["gvs_sas_af"]),
    ("oth", ["gvs_oth_af"]),
]


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
    gff: str = ""
    fasta: str = ""

    def _gnomad_exomes_fields(self) -> list[str]:
        """The `fields` values for the gnomAD exomes custom line, one per selected
        ancestry x sex. Grammar: AF[_non_ukb][_<ancestry>][_XX|_XY] — base `AF` is
        all ancestries, both sexes, UKB included; XX=female, XY=male."""
        fields: list[str] = []
        for anc_param, anc_code in GNOMAD_EXOMES_ANCESTRIES:
            if not getattr(self, f"gnomad_exomes_{anc_param}"):
                continue
            for sex_param, sex_code in GNOMAD_SEXES:
                if not getattr(self, f"gnomad_exomes_{anc_param}_{sex_param}"):
                    continue
                parts = ["AF"]
                if not self.gnomad_exomes_include_ukb:
                    parts.append("non_ukb")
                if anc_code:
                    parts.append(anc_code)
                if sex_code:
                    parts.append(sex_code)
                fields.append("_".join(parts))
        return fields

    def _gnomad_genomes_fields(self) -> list[str]:
        """The `fields` values for the gnomAD genomes custom line. Grammar:
        AF[_<ancestry>][_XX|_XY] (no UK Biobank subset). grpmax (max across
        groups) has no sex split, so it contributes a single `AF_grpmax`."""
        fields: list[str] = []
        for anc_param, anc_code in GNOMAD_GENOMES_ANCESTRIES:
            if not getattr(self, f"gnomad_genomes_{anc_param}"):
                continue
            for sex_param, sex_code in GNOMAD_SEXES:
                if not getattr(self, f"gnomad_genomes_{anc_param}_{sex_param}"):
                    continue
                parts = ["AF"]
                if anc_code:
                    parts.append(anc_code)
                if sex_code:
                    parts.append(sex_code)
                fields.append("_".join(parts))
        if self.gnomad_genomes_grpmax:
            fields.append("AF_grpmax")
        return fields

    def _allofus_fields(self) -> list[str]:
        """The `fields` values for the All of Us custom line, in population order.
        "Maximum subpopulation" contributes both gvs_max_af and gvs_max_subpop."""
        fields: list[str] = []
        for pop, codes in ALLOFUS_POPULATIONS:
            if getattr(self, f"allofus_{pop}"):
                fields.extend(codes)
        return fields

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

        # gnomAD exomes v4.1 (human GRCh38): a VEP `custom` line whose `fields`
        # list is built from the selected ancestry x sex combinations. Emitted
        # only when at least one field resolves (nothing selected -> no line).
        if self.gnomad_exomes:
            exome_fields = self._gnomad_exomes_fields()
            if exome_fields:
                lines.append(
                    "custom "
                    f"file={PLUGIN_PATH}/gnomad.exomes.v4.1.sites.chr###CHR###.vcf.bgz,"
                    "short_name=gnomAD_exomes,"
                    f"fields={'%'.join(exome_fields)},"
                    "format=vcf"
                )

        # gnomAD genomes v4.1 (human GRCh38): as gnomAD exomes above.
        if self.gnomad_genomes:
            genome_fields = self._gnomad_genomes_fields()
            if genome_fields:
                lines.append(
                    "custom "
                    f"file={PLUGIN_PATH}/gnomad.genomes.v4.1.sites.chr###CHR###.vcf.bgz,"
                    "short_name=gnomAD_genomes,"
                    f"fields={'%'.join(genome_fields)},"
                    "format=vcf"
                )

        # NIH All of Us (human GRCh38): custom line built from the selected
        # population fields.
        if self.allofus:
            allofus_fields = self._allofus_fields()
            if allofus_fields:
                lines.append(
                    "custom "
                    f"file={PLUGIN_PATH}/AllOfUs_chr###CHR###.vcf.gz,"
                    "short_name=AoU,"
                    f"fields={'%'.join(allofus_fields)},"
                    "format=vcf"
                )

        # ClinVar clinical significance (human GRCh37/GRCh38): custom line
        # surfacing the CLNSIG field, from the assembly-specific ClinVar file.
        if self.clinvar:
            lines.append(
                "custom "
                f"file={PLUGIN_PATH}/clinvar_{assembly}.vcf.gz,"
                "short_name=ClinVar,"
                "fields=CLNSIG%CLNSIGCONF,"
                "format=vcf,"
                "type=exact"
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
