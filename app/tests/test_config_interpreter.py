"""Differential test: the declarative config interpreter must reproduce exactly
the option-driven lines that the hardcoded create_config_ini_file emits.

Compared as *sets* of lines — a VEP config.ini is order-independent, and the
migration deliberately re-homes the always-on base (which stays in the backend),
so line position is not meaningful. The always-on base lines (identified by their
keyword) are subtracted from the golden output; everything else must match what
the interpreter emits from the same options.

Scope: the representative subset currently in config_specs/human_grch38.json (one
of each emitter kind). The gnomAD/AoU `custom` field builders land in a later
increment and are not exercised here.
"""

import json
from pathlib import Path

import pytest

import vep
from vep.models.config_spec_model import ConfigSpec
from vep.models.pipeline_model import ConfigIniParams, PLUGIN_PATH
from vep.utils.config_interpreter import emit_config_lines

SPEC_PATH = Path(vep.__file__).resolve().parent / "config_specs" / "human_grch38.json"
CONFIG_SPEC = ConfigSpec.model_validate(json.loads(SPEC_PATH.read_text()))

# The always-on base config the backend keeps (not the interpreter's job). Lines
# are matched by their first token.
BASE_KEYWORDS = {
    "force_overwrite", "numbers", "mane", "assembly", "symbol", "biotype",
    "transcript_version", "canonical", "database", "gff", "fasta",
}

FAKE_SUPPORT = {"gff_location": "/data/x.gff", "faa_location": "/data/x.faa"}


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    # create_config_ini_file resolves gff/fasta from a web lookup; pin it.
    monkeypatch.setattr(
        "vep.models.pipeline_model.get_vep_support_location",
        lambda genome_id: FAKE_SUPPORT,
    )


def _resolved_assembly(assembly_name: str) -> str:
    return "GRCh37" if assembly_name.startswith("GRCh37") else "GRCh38"


def _golden_option_lines(params: ConfigIniParams, tmp_path) -> set[str]:
    """The real builder's output with the always-on base removed."""
    params.create_config_ini_file(str(tmp_path))
    lines = [
        ln.strip()
        for ln in (tmp_path / "config.ini").read_text().splitlines()
        if ln.strip()
    ]
    return {ln for ln in lines if ln.split()[0] not in BASE_KEYWORDS}


# Every plugin/custom entry in the spec (the flags are exercised via every case).
PLUGIN_CUSTOM_IDS = [
    e.id for e in CONFIG_SPEC.entries if e.config.emit in ("plugin", "custom")
]

# Each option alone, on both assemblies — hits every line and every by_assembly
# branch (incl. SpliceAI's GRCh38-only snv_ensembl).
SINGLE_CASES = [
    (f"{oid} {asm}", asm, {oid: True})
    for oid in PLUGIN_CUSTOM_IDS
    for asm in ("GRCh38", "GRCh37")
]

# Sub-options / value params, which a plain "option on" leaves at their defaults.
SUBOPTION_CASES = [
    ("all flags on", "GRCh38", {"hgvsg": True, "spdi": True, "protein": True}),
    ("protvar pocket off", "GRCh38", {"protvar": True, "protvar_pocket": False}),
    ("mutfunc some subs", "GRCh38", {"mutfunc": True, "mutfunc_motif": True, "mutfunc_exp": True}),
    ("dosage cover on", "GRCh38", {"dosage_sensitivity": True, "dosage_sensitivity_cover": True}),
    ("intact all -> all=1", "GRCh38", {
        "intact": True, "intact_feature_ac": True, "intact_feature_short_label": True,
        "intact_feature_annotation": True, "intact_ap_ac": True,
        "intact_interaction_participants": True, "intact_pmid": True,
    }),
    ("intact two sub-flags", "GRCh38", {
        "intact": True, "intact_feature_ac": True, "intact_pmid": True,
    }),
    ("tss downstream", "GRCh38", {"tss_distance": True, "tss_distance_direction": "downstream"}),
    ("tss both", "GRCh38", {"tss_distance": True, "tss_distance_direction": "both"}),
    ("nearest_gene both dirs", "GRCh38", {"nearest_gene": True, "nearest_gene_both_directions": True}),
    ("nearest_exon_jb custom", "GRCh38", {
        "nearest_exon_jb": True, "nearest_exon_jb_max_range": 5000,
        "nearest_exon_jb_intronic": True,
    }),
    # gnomAD / All of Us custom `fields=` builders.
    ("gnomad_exomes non_ukb", "GRCh38", {
        "gnomad_exomes": True, "gnomad_exomes_include_ukb": False,
    }),
    ("gnomad_exomes multi anc+sex", "GRCh38", {
        "gnomad_exomes": True, "gnomad_exomes_afr": True,
        "gnomad_exomes_afr_female": True, "gnomad_exomes_nfe": True,
    }),
    ("gnomad_exomes no fields -> omitted", "GRCh38", {
        "gnomad_exomes": True, "gnomad_exomes_all": False,
    }),
    ("gnomad_genomes grpmax", "GRCh38", {
        "gnomad_genomes": True, "gnomad_genomes_grpmax": True,
    }),
    ("gnomad_genomes ami+remaining+sex", "GRCh38", {
        "gnomad_genomes": True, "gnomad_genomes_ami": True,
        "gnomad_genomes_remaining": True, "gnomad_genomes_all_male": True,
    }),
    ("allofus max + several", "GRCh38", {
        "allofus": True, "allofus_max": True, "allofus_afr": True, "allofus_eur": True,
    }),
    ("allofus no fields -> omitted", "GRCh38", {
        "allofus": True, "allofus_all": False,
    }),
]


def _all_on(assembly: str) -> tuple:
    overrides = {oid: True for oid in PLUGIN_CUSTOM_IDS}
    overrides.update({"hgvsg": True, "spdi": True, "protein": True})
    return (f"kitchen sink {assembly}", assembly, overrides)


CASES = (
    [("defaults GRCh38", "GRCh38", {}), ("defaults GRCh37", "GRCh37", {})]
    + SINGLE_CASES
    + SUBOPTION_CASES
    + [_all_on("GRCh38"), _all_on("GRCh37")]
)


@pytest.mark.parametrize("desc, assembly_name, overrides", CASES, ids=[c[0] for c in CASES])
def test_interpreter_matches_builder(desc, assembly_name, overrides, tmp_path):
    params = ConfigIniParams(genome_id="test", assembly_name=assembly_name, **overrides)
    expected = _golden_option_lines(params, tmp_path)

    got = set(
        emit_config_lines(
            CONFIG_SPEC,
            params.model_dump(),
            assembly=_resolved_assembly(assembly_name),
            plugin_path=PLUGIN_PATH,
            gff=params.gff,
        )
    )

    assert got == expected, (
        f"\n{desc}\n  only in builder: {sorted(expected - got)}"
        f"\n  only in interpreter: {sorted(got - expected)}"
    )
