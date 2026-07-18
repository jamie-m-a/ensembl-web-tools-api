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
    "transcript_version", "canonical", "gff", "fasta",
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


# (description, assembly_name, option overrides)
CASES = [
    ("defaults GRCh38", "GRCh38", {}),
    ("defaults GRCh37", "GRCh37", {}),
    ("all flags on", "GRCh38", {"hgvsg": True, "spdi": True, "protein": True}),
    ("mavedb", "GRCh38", {"mavedb": True}),
    ("protvar all sub-flags", "GRCh38", {"protvar": True}),
    ("protvar pocket off", "GRCh38", {"protvar": True, "protvar_pocket": False}),
    ("intact base only", "GRCh38", {"intact": True}),
    ("intact all -> all=1", "GRCh38", {
        "intact": True, "intact_feature_ac": True, "intact_feature_short_label": True,
        "intact_feature_annotation": True, "intact_ap_ac": True,
        "intact_interaction_participants": True, "intact_pmid": True,
    }),
    ("intact two sub-flags", "GRCh38", {
        "intact": True, "intact_feature_ac": True, "intact_pmid": True,
    }),
    ("alphamissense GRCh38", "GRCh38", {"alphamissense": True}),
    ("alphamissense GRCh37", "GRCh37", {"alphamissense": True}),
    ("nearest_exon_jb defaults", "GRCh38", {"nearest_exon_jb": True}),
    ("nearest_exon_jb custom", "GRCh38", {
        "nearest_exon_jb": True, "nearest_exon_jb_max_range": 5000,
        "nearest_exon_jb_intronic": True,
    }),
    ("clinvar GRCh38", "GRCh38", {"clinvar": True}),
    ("clinvar GRCh37", "GRCh37", {"clinvar": True}),
    ("several together", "GRCh38", {
        "mavedb": True, "protvar": True, "protvar_int": False, "clinvar": True,
        "alphamissense": True,
    }),
]


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
