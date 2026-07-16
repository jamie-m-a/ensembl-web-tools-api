"""Tests for VEP config.ini construction (ConfigIniParams.create_config_ini_file).

The builder turns the client's submission booleans into a flat VEP config.ini:
always-on defaults, `key value` flag lines, and `plugin ...` lines (some of
which are assembly-specific or carry sub-option flags). These tests monkeypatch
the metadata lookup so no network call is made, build the ini into a tmp dir,
read it back, and assert on the emitted lines.
"""

import re

import pytest

from app.vep.models.pipeline_model import ConfigIniParams, PLUGIN_PATH

GFF = "/vep_support/test.gff3.gz"
FASTA = "/vep_support/test.fa"


def build_lines(monkeypatch, tmp_path, *, assembly="GRCh38.p14", **kwargs):
    """Build a config.ini for the given assembly/options and return its lines."""
    monkeypatch.setattr(
        "app.vep.models.pipeline_model.get_vep_support_location",
        lambda genome_id: {"gff_location": GFF, "faa_location": FASTA},
    )
    params = ConfigIniParams(
        genome_id="genome-under-test", assembly_name=assembly, **kwargs
    )
    params.create_config_ini_file(str(tmp_path))
    return (tmp_path / "config.ini").read_text().splitlines()


def find_line(lines, needle):
    """The first line containing `needle`, or None."""
    return next((line for line in lines if needle in line), None)


def plugin_lines(lines):
    return [line for line in lines if line.startswith("plugin ")]


# --- 1. always-on defaults ---------------------------------------------------


def test_always_on_defaults(monkeypatch, tmp_path):
    lines = build_lines(monkeypatch, tmp_path)
    for expected in [
        "force_overwrite 1",
        "numbers 1",
        "symbol 1",
        "biotype 1",
        "transcript_version 1",
        "canonical 1",
        f"gff {GFF}",
        f"fasta {FASTA}",
    ]:
        assert expected in lines


# --- 2. flag options ---------------------------------------------------------


def test_flag_options_render_as_one_or_zero(monkeypatch, tmp_path):
    lines = build_lines(
        monkeypatch, tmp_path, hgvs=True, hgvsg=False, spdi=True, protein=False
    )
    assert "hgvs 1" in lines
    assert "hgvsg 0" in lines
    assert "spdi 1" in lines
    assert "protein 0" in lines


# --- 3. mane gating ----------------------------------------------------------


@pytest.mark.parametrize(
    "assembly,expected",
    [
        ("GRCh38.p14", True),  # human reference
        ("GRCm39", True),  # mouse reference
        ("GRCh37.p13", False),  # human GRCh37 has no MANE
        ("T2T-CHM13v2.0", False),  # human T2T
        ("ARS-UCD1.2", False),  # non-human (cow)
    ],
)
def test_mane_gating(monkeypatch, tmp_path, assembly, expected):
    lines = build_lines(monkeypatch, tmp_path, assembly=assembly)
    assert ("mane 1" in lines) is expected


# --- 4. assembly line --------------------------------------------------------


@pytest.mark.parametrize(
    "assembly,expected",
    [
        ("GRCh38.p14", "assembly GRCh38"),
        ("GRCh37.p13", "assembly GRCh37"),
        ("T2T-CHM13v2.0", None),
        ("GRCm39", None),
    ],
)
def test_assembly_line(monkeypatch, tmp_path, assembly, expected):
    lines = build_lines(monkeypatch, tmp_path, assembly=assembly)
    assembly_lines = [line for line in lines if line.startswith("assembly ")]
    if expected is None:
        assert assembly_lines == []
    else:
        assert assembly_lines == [expected]


# --- 5. per-assembly plugin files --------------------------------------------

# option -> (substring only in the GRCh38 file, substring only in the GRCh37 file)
ASSEMBLY_PLUGIN_MARKERS = {
    "alphamissense": ("AlphaMissense_hg38", "AlphaMissense_hg19"),
    "cadd": ("CADD_GRCh38", "CADD_GRCh37"),
    "revel": ("revel_grch38", "revel_grch37"),
    "loeuf": ("loeuf_dataset_grch38", "loeuf_dataset_grch37"),
    "geno2mp": ("Geno2MP.variants_GRCh38", "Geno2MP.variants_GRCh37"),
    "enformer": ("enformer_grch38", "enformer_grch37"),
    "utrannotator": ("GRCh38_PUBLIC", "GRCh37_PUBLIC"),
}


@pytest.mark.parametrize("option,markers", ASSEMBLY_PLUGIN_MARKERS.items())
def test_per_assembly_plugin_files(monkeypatch, tmp_path, option, markers):
    grch38_marker, grch37_marker = markers

    lines38 = build_lines(
        monkeypatch, tmp_path, assembly="GRCh38.p14", **{option: True}
    )
    assert find_line(lines38, grch38_marker) is not None
    assert find_line(lines38, grch37_marker) is None

    lines37 = build_lines(
        monkeypatch, tmp_path, assembly="GRCh37.p13", **{option: True}
    )
    assert find_line(lines37, grch37_marker) is not None
    assert find_line(lines37, grch38_marker) is None


def test_spliceai_grch37_omits_snv_ensembl(monkeypatch, tmp_path):
    line38 = find_line(
        build_lines(monkeypatch, tmp_path, assembly="GRCh38.p14", spliceai=True),
        "plugin SpliceAI",
    )
    line37 = find_line(
        build_lines(monkeypatch, tmp_path, assembly="GRCh37.p13", spliceai=True),
        "plugin SpliceAI",
    )
    assert "snv_ensembl=" in line38
    assert "snv_ensembl=" not in line37


# --- 6. static (assembly-independent) plugins --------------------------------


@pytest.mark.parametrize(
    "option,prefix",
    [
        ("mavedb", "plugin MaveDB,"),
        ("opentargets", "plugin OpenTargets,"),
        ("eve", "plugin EVE,"),
        ("phenotypes", "plugin Phenotypes,"),
        ("gnomad_mt", "plugin gnomADMt,"),
        ("maxentscan", "plugin MaxEntScan,"),
        ("riboseqorfs", "plugin RiboseqORFs,"),
    ],
)
def test_static_plugins_emit_expected_line(monkeypatch, tmp_path, option, prefix):
    lines = build_lines(monkeypatch, tmp_path, **{option: True})
    assert any(line.startswith(prefix) for line in lines)


# --- 7. sub-flag plugins (ProtVar / mutfunc / DosageSensitivity) -------------


def test_protvar_sub_flags(monkeypatch, tmp_path):
    # defaults: all three sub-options on
    default_line = find_line(
        build_lines(monkeypatch, tmp_path, protvar=True), "plugin ProtVar"
    )
    assert "stability=1,pocket=1,int=1" in default_line

    # a partial selection is reflected one-to-one
    partial_line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            protvar=True,
            protvar_stability=False,
            protvar_pocket=True,
            protvar_int=False,
        ),
        "plugin ProtVar",
    )
    assert "stability=0,pocket=1,int=0" in partial_line


def test_mutfunc_sub_flags_and_extended(monkeypatch, tmp_path):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            mutfunc=True,
            mutfunc_motif=True,
            mutfunc_int=False,
            mutfunc_mod=True,
            mutfunc_exp=False,
        ),
        "plugin mutfunc",
    )
    assert "motif=1,int=0,mod=1,exp=0" in line
    assert "extended=1" in line  # always on


@pytest.mark.parametrize("cover,expected", [(False, "cover=0"), (True, "cover=1")])
def test_dosage_sensitivity_cover(monkeypatch, tmp_path, cover, expected):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            dosage_sensitivity=True,
            dosage_sensitivity_cover=cover,
        ),
        "plugin DosageSensitivity",
    )
    assert expected in line


# --- 8. IntAct tri-state -----------------------------------------------------


def test_intact_no_sub_options_is_base_line(monkeypatch, tmp_path):
    line = find_line(
        build_lines(monkeypatch, tmp_path, intact=True), "plugin IntAct"
    )
    assert f"mutation_file={PLUGIN_PATH}/mutations.tsv" in line
    assert "mapping_file=" in line
    assert "all=1" not in line
    assert "=1" not in line.split("mapping_file=")[1]  # no sub-flags appended


def test_intact_all_sub_options_uses_all(monkeypatch, tmp_path):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            intact=True,
            intact_feature_ac=True,
            intact_feature_short_label=True,
            intact_feature_annotation=True,
            intact_ap_ac=True,
            intact_interaction_participants=True,
            intact_pmid=True,
        ),
        "plugin IntAct",
    )
    assert line.endswith(",all=1")


def test_intact_partial_sub_options_lists_selected_in_order(monkeypatch, tmp_path):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            intact=True,
            intact_feature_ac=True,
            intact_pmid=True,
        ),
        "plugin IntAct",
    )
    assert line.endswith(",feature_ac=1,pmid=1")  # fixed order, only the selected
    assert "all=1" not in line


# --- 9. gff3-based Genes & transcripts plugins -------------------------------


@pytest.mark.parametrize(
    "direction,expected",
    [
        ("upstream", "upstream=1,downstream=0,both=0"),
        ("downstream", "upstream=0,downstream=1,both=0"),
        ("both", "upstream=0,downstream=0,both=1"),
    ],
)
def test_tss_distance_direction(monkeypatch, tmp_path, direction, expected):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            tss_distance=True,
            tss_distance_direction=direction,
        ),
        "plugin TSSDistance",
    )
    assert expected in line
    assert line.endswith(f"gff3={GFF}")


def test_nearest_gene_both_directions_and_gff3(monkeypatch, tmp_path):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            nearest_gene=True,
            nearest_gene_both_directions=True,
        ),
        "plugin NearestGene",
    )
    assert "both_directions=1" in line
    assert line.endswith(f"gff3={GFF}")


def test_nearest_exon_jb_range_intronic_and_gff3(monkeypatch, tmp_path):
    line = find_line(
        build_lines(
            monkeypatch,
            tmp_path,
            nearest_exon_jb=True,
            nearest_exon_jb_intronic=True,
        ),
        "plugin NearestExonJB",
    )
    assert "max_range=10000" in line  # default
    assert "intronic=1" in line
    assert line.endswith(f"gff3={GFF}")


# --- 10. disabled options omitted -------------------------------------------


def test_disabled_options_emit_no_plugin_lines(monkeypatch, tmp_path):
    # everything left at its default (off) => no `plugin ...` lines at all
    assert plugin_lines(build_lines(monkeypatch, tmp_path)) == []


# --- 11. gnomAD exomes custom line (ancestry x sex field grammar) ------------


def gnomad_exomes_line(lines):
    return find_line(lines, "short_name=gnomAD_exomes")


def gnomad_exomes_fields(lines):
    line = gnomad_exomes_line(lines)
    if line is None:
        return None
    return re.search(r"fields=([^,]+)", line).group(1)


def test_gnomad_exomes_off_emits_no_line(monkeypatch, tmp_path):
    assert gnomad_exomes_line(build_lines(monkeypatch, tmp_path)) is None


def test_gnomad_exomes_default_is_all_both_ukb_included(monkeypatch, tmp_path):
    # enabling with defaults (All + Both + include UKB) => fields=AF
    line = gnomad_exomes_line(build_lines(monkeypatch, tmp_path, gnomad_exomes=True))
    assert line is not None
    assert (
        f"file={PLUGIN_PATH}/gnomad.exomes.v4.1.sites.chr###CHR###.vcf.bgz" in line
    )
    assert "short_name=gnomAD_exomes" in line
    assert line.endswith("format=vcf")
    assert "fields=AF," in line


def test_gnomad_exomes_all_female_and_male(monkeypatch, tmp_path):
    # spec example 1: All, Female and Male -> AF_XX%AF_XY
    fields = gnomad_exomes_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            gnomad_exomes=True,
            gnomad_exomes_all_both=False,
            gnomad_exomes_all_female=True,
            gnomad_exomes_all_male=True,
        )
    )
    assert fields == "AF_XX%AF_XY"


def test_gnomad_exomes_non_ukb_male(monkeypatch, tmp_path):
    # spec example 2: All, male, excluding UK Biobank -> AF_non_ukb_XY
    fields = gnomad_exomes_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            gnomad_exomes=True,
            gnomad_exomes_include_ukb=False,
            gnomad_exomes_all_both=False,
            gnomad_exomes_all_male=True,
        )
    )
    assert fields == "AF_non_ukb_XY"


def test_gnomad_exomes_multiple_ancestries_ordered(monkeypatch, tmp_path):
    # afr (both) + nfe (female); ancestry order preserved, sex suffix applied
    fields = gnomad_exomes_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            gnomad_exomes=True,
            gnomad_exomes_all=False,
            gnomad_exomes_afr=True,
            gnomad_exomes_nfe=True,
            gnomad_exomes_nfe_both=False,
            gnomad_exomes_nfe_female=True,
        )
    )
    assert fields == "AF_afr%AF_nfe_XX"


def test_gnomad_exomes_enabled_but_nothing_selected_emits_no_line(
    monkeypatch, tmp_path
):
    # on, but "All" (the only default-on ancestry) has no sex selected => no line
    line = gnomad_exomes_line(
        build_lines(
            monkeypatch,
            tmp_path,
            gnomad_exomes=True,
            gnomad_exomes_all_both=False,
        )
    )
    assert line is None


# --- 12. gnomAD genomes custom line (no UKB subset; ami/remaining/grpmax) -----


def gnomad_genomes_line(lines):
    return find_line(lines, "short_name=gnomAD_genomes")


def gnomad_genomes_fields(lines):
    line = gnomad_genomes_line(lines)
    if line is None:
        return None
    return re.search(r"fields=([^,]+)", line).group(1)


def test_gnomad_genomes_off_emits_no_line(monkeypatch, tmp_path):
    assert gnomad_genomes_line(build_lines(monkeypatch, tmp_path)) is None


def test_gnomad_genomes_default_line(monkeypatch, tmp_path):
    line = gnomad_genomes_line(build_lines(monkeypatch, tmp_path, gnomad_genomes=True))
    assert line is not None
    assert (
        f"file={PLUGIN_PATH}/gnomad.genomes.v4.1.sites.chr###CHR###.vcf.bgz" in line
    )
    assert "short_name=gnomAD_genomes" in line
    assert line.endswith("format=vcf")
    assert "fields=AF," in line  # default: All + Both


def test_gnomad_genomes_has_no_ukb_param(monkeypatch, tmp_path):
    # genomes has no UK Biobank subset, so no _non_ukb fields are possible
    assert "gnomad_genomes_include_ukb" not in ConfigIniParams.model_fields


def test_gnomad_genomes_ancestry_and_sex(monkeypatch, tmp_path):
    # Amish (both) + Remaining (male)
    fields = gnomad_genomes_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            gnomad_genomes=True,
            gnomad_genomes_all=False,
            gnomad_genomes_ami=True,
            gnomad_genomes_remaining=True,
            gnomad_genomes_remaining_both=False,
            gnomad_genomes_remaining_male=True,
        )
    )
    assert fields == "AF_ami%AF_remaining_XY"


def test_gnomad_genomes_grpmax_has_no_sex_split(monkeypatch, tmp_path):
    # grpmax on its own -> single AF_grpmax (no XX/XY)
    fields = gnomad_genomes_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            gnomad_genomes=True,
            gnomad_genomes_all=False,
            gnomad_genomes_grpmax=True,
        )
    )
    assert fields == "AF_grpmax"


def test_gnomad_genomes_all_plus_grpmax(monkeypatch, tmp_path):
    # default All (both) plus grpmax -> AF%AF_grpmax
    fields = gnomad_genomes_fields(
        build_lines(monkeypatch, tmp_path, gnomad_genomes=True, gnomad_genomes_grpmax=True)
    )
    assert fields == "AF%AF_grpmax"


# --- 13. NIH All of Us custom line -------------------------------------------


def allofus_line(lines):
    return find_line(lines, "short_name=AoU")


def allofus_fields(lines):
    line = allofus_line(lines)
    if line is None:
        return None
    return re.search(r"fields=([^,]+)", line).group(1)


def test_allofus_off_emits_no_line(monkeypatch, tmp_path):
    assert allofus_line(build_lines(monkeypatch, tmp_path)) is None


def test_allofus_default_line(monkeypatch, tmp_path):
    line = allofus_line(build_lines(monkeypatch, tmp_path, allofus=True))
    assert line is not None
    assert f"file={PLUGIN_PATH}/AllOfUs_chr###CHR###.vcf.gz" in line
    assert "short_name=AoU" in line
    assert line.endswith("format=vcf")
    assert "fields=gvs_all_af," in line  # default: All


def test_allofus_max_emits_two_fields(monkeypatch, tmp_path):
    # "Maximum subpopulation" contributes both gvs_max_af and gvs_max_subpop
    fields = allofus_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            allofus=True,
            allofus_all=False,
            allofus_max=True,
        )
    )
    assert fields == "gvs_max_af%gvs_max_subpop"


def test_allofus_multiple_populations_in_order(monkeypatch, tmp_path):
    # order follows ALLOFUS_POPULATIONS (all, max, afr, ...)
    fields = allofus_fields(
        build_lines(
            monkeypatch,
            tmp_path,
            allofus=True,
            allofus_afr=True,
            allofus_sas=True,
        )
    )
    assert fields == "gvs_all_af%gvs_afr_af%gvs_sas_af"


def test_allofus_enabled_but_nothing_selected_emits_no_line(monkeypatch, tmp_path):
    line = allofus_line(
        build_lines(monkeypatch, tmp_path, allofus=True, allofus_all=False)
    )
    assert line is None
