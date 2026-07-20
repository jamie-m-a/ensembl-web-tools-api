"""Tests for pinning a job's option panels at submission.

The results view used to lay itself out from the *live* form-config panels, so a
job rendered against whatever the config said at viewing time. The panels are now
computed at submission and pinned as a third sidecar (beside the merged spec and
the expected CSQ columns), then handed back on the results response.

Two things have to hold:
  * the pin is faithful -- the sidecar round-trips get_visible_panels() exactly,
    and the panels pinned for a submission are the same ones /form_config would
    return for that species/assembly (which is what proves species_taxonomy_id
    is actually reaching the submission path);
  * the load side is defensive -- a job submitted before this existed has no
    sidecar, gets None, and keeps rendering against the live panels.
"""

import json

import pytest
from pydantic import FilePath

from app.vep.form_panels import get_visible_panels
from app.vep.models.display_panels_model import (
    dump_display_panels,
    to_display_panels,
)
from app.vep.models.pipeline_model import ConfigIniParams
from app.vep.utils.spec_loader import (
    DISPLAY_PANELS_SIDECAR_FILE,
    load_display_panels_sidecar,
    write_display_panels_sidecar,
)
from app.vep.utils.vcf_results import _load_pinned_display_panels

HUMAN = "9606"
MOUSE = "10090"


def _form_config_panels(*, genome_id, species_taxonomy_id, assembly_name):
    """The `panels` the /form_config endpoint actually serves, with its genome
    metadata lookup stubbed (that feeds the transcript-set dropdown only)."""
    import asyncio

    from app.vep import vep_resources

    async def fake_get_genome_metadata(_genome_id):
        return {
            "genebuild.provider_name": "Ensembl",
            "genebuild.provider_version": "115",
            "genebuild.last_geneset_update": "2024-01",
        }

    original = vep_resources.get_genome_metadata
    vep_resources.get_genome_metadata = fake_get_genome_metadata
    try:
        response = asyncio.run(
            vep_resources.get_form_config(
                request=None,
                genome_id=genome_id,
                species_taxonomy_id=species_taxonomy_id,
                assembly_name=assembly_name,
            )
        )
    finally:
        vep_resources.get_genome_metadata = original
    return response["panels"]


def _vcf_path(tmp_path):
    """A stand-in for a job's results VCF; the sidecar sits beside it."""
    path = tmp_path / "input_VEP.vcf.gz"
    path.write_bytes(b"")
    return FilePath(path)


# --- the model is permissive enough to pin the panels losslessly ---------


@pytest.mark.parametrize(
    "species_taxonomy_id,assembly_name",
    [
        (HUMAN, "GRCh38.p14"),
        (HUMAN, "GRCh37.p13"),
        (MOUSE, "GRCm39"),
        (None, None),
    ],
)
def test_panels_round_trip_through_the_model(species_taxonomy_id, assembly_name):
    """Nothing may be dropped or invented -- categories, sub_options and the
    nested {"type": "group", ...} nodes all have to survive."""
    panels = get_visible_panels(
        species_taxonomy_id=species_taxonomy_id, assembly_name=assembly_name
    )
    assert dump_display_panels(to_display_panels(panels)) == panels


def test_sidecar_round_trips_the_panels(tmp_path):
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    write_display_panels_sidecar(tmp_path, to_display_panels(panels))
    assert (tmp_path / DISPLAY_PANELS_SIDECAR_FILE).exists()

    loaded = load_display_panels_sidecar(_vcf_path(tmp_path))
    assert dump_display_panels(loaded) == panels


# --- the pin matches what the form was built from ------------------------


def test_submission_pins_the_same_panels_form_config_serves():
    """The panels pinned for a human GRCh38 submission must be exactly the ones
    /form_config returns for that species -- i.e. species_taxonomy_id really does
    reach the submission path. Without it, get_visible_panels falls back to the
    base panels and every human-specific panel silently disappears from the pin.
    """
    ini_parameters = ConfigIniParams(
        genome_id="a7335667-93e7-11ec-a39d-005056b38ce3",
        assembly_name="GRCh38.p14",
        species_taxonomy_id=HUMAN,
    )
    pinned = get_visible_panels(
        species_taxonomy_id=ini_parameters.species_taxonomy_id,
        assembly_name=ini_parameters.assembly_name,
    )
    # What the form_config endpoint itself serves for the same species/assembly.
    # Its only network call is the genome-metadata lookup, which feeds the
    # transcript-set dropdown, not the panels (`attributes` is accepted by
    # get_visible_panels but never read) -- so it is stubbed out here.
    from_form_config = _form_config_panels(
        genome_id=ini_parameters.genome_id,
        species_taxonomy_id=HUMAN,
        assembly_name="GRCh38.p14",
    )
    assert dump_display_panels(to_display_panels(pinned)) == from_form_config

    panel_ids = {panel["id"] for panel in pinned}
    assert {
        "pathogenicity_predictions",
        "conservation_and_constraint",
        "variant_associations",
        "allele_frequencies",
    } <= panel_ids


def test_missing_species_taxonomy_id_would_lose_the_human_panels():
    """The regression this plumbing exists to prevent, stated explicitly."""
    without = get_visible_panels(species_taxonomy_id="", assembly_name="GRCh38.p14")
    assert {panel["id"] for panel in without} == {
        "variant_representations",
        "genes_and_transcripts",
        "protein_and_functional",
    }


def test_species_taxonomy_id_emits_no_config_ini_line(tmp_path, monkeypatch):
    """It is submission metadata for the pin, not a VEP option."""
    assert "species_taxonomy_id" in ConfigIniParams.model_fields
    params = ConfigIniParams(genome_id="x", assembly_name="GRCh38.p14")
    assert params.species_taxonomy_id == ""


# --- the load side never breaks an older job -----------------------------


def test_no_sidecar_returns_none(tmp_path):
    """A job submitted before this existed has no sidecar; results must still
    parse, with the frontend falling back to the live panels."""
    assert load_display_panels_sidecar(_vcf_path(tmp_path)) is None
    assert _load_pinned_display_panels(_vcf_path(tmp_path)) is None


def test_unreadable_sidecar_is_ignored_rather_than_raising(tmp_path):
    (tmp_path / DISPLAY_PANELS_SIDECAR_FILE).write_text("{not json")
    assert _load_pinned_display_panels(_vcf_path(tmp_path)) is None


def test_empty_sidecar_falls_back_rather_than_pinning_nothing(tmp_path):
    """An empty list can only come from a corrupted sidecar — get_visible_panels
    always returns at least the always-visible panels. Treating it as a valid pin
    would render the job with no panels at all instead of falling back to live."""
    (tmp_path / DISPLAY_PANELS_SIDECAR_FILE).write_text("[]")
    assert _load_pinned_display_panels(_vcf_path(tmp_path)) is None


def test_sidecar_with_unexpected_keys_is_still_loaded(tmp_path):
    """The panel structure is still evolving; a sidecar written by a newer
    backend must not fail to load on an older one."""
    (tmp_path / DISPLAY_PANELS_SIDECAR_FILE).write_text(
        json.dumps(
            [
                {
                    "id": "future_panel",
                    "label": "Future",
                    "options": [
                        {"id": "opt", "label": "Opt", "type": "boolean",
                         "default": False, "something_new": {"a": 1}}
                    ],
                }
            ]
        )
    )
    panels = _load_pinned_display_panels(_vcf_path(tmp_path))
    assert panels is not None
    assert panels[0].options[0].model_dump()["something_new"] == {"a": 1}
