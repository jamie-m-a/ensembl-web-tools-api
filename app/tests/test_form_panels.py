"""Tests for form-panel visibility (form_panels.get_visible_panels).

The form_config endpoint returns the panels/options to show for the selected
genome. The set is species/assembly conditional: a common base for every genome,
extra options + panels for human GRCh37/38, and further options for human GRCh38
only. Option (and sub-option) ids double as ConfigIniParams parameter names, so
the form round-trips into the generated config.ini.
"""

from app.vep import form_panels
from app.vep.form_panels import get_visible_panels
from app.vep.models.pipeline_model import ConfigIniParams

HUMAN = "9606"
MOUSE = "10090"

ALWAYS_VISIBLE_PANEL_IDS = {
    "variant_representations",
    "genes_and_transcripts",
    "protein_and_functional",
}
HUMAN_37_38_PANEL_IDS = {
    "pathogenicity_predictions",
    "conservation_and_constraint",
    "variant_associations",
}
GRCH38_ONLY_OPTION_IDS = {
    "eve",
    "intact",
    "mavedb",
    "opentargets",
    "protvar",
    "riboseqorfs",
    "gnomad_mt",
}


def panel_ids(panels):
    return {panel["id"] for panel in panels}


def option_ids(panels, *, include_sub_options=True):
    ids = set()
    for panel in panels:
        for option in panel["options"]:
            ids.add(option["id"])
            if include_sub_options:
                for sub in option.get("sub_options", []):
                    ids.add(sub["id"])
    return ids


def categories(panels, panel_id):
    panel = next(p for p in panels if p["id"] == panel_id)
    return {opt["category"] for opt in panel["options"] if "category" in opt}


# --- 1. human GRCh38 ---------------------------------------------------------


def test_human_grch38_shows_all_panels_and_options():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    ids = panel_ids(panels)

    assert ALWAYS_VISIBLE_PANEL_IDS <= ids
    assert HUMAN_37_38_PANEL_IDS <= ids
    assert "allele_frequencies" in ids  # GRCh38-only panel

    opts = option_ids(panels)
    assert GRCH38_ONLY_OPTION_IDS <= opts
    assert "utrannotator" in opts  # 37/38 option


def test_human_grch38_category_labels():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    assert categories(panels, "pathogenicity_predictions") == {
        "Missense",
        "Splicing",
        "Genome wide",
        "Non-coding",
    }
    assert categories(panels, "protein_and_functional") == {"Protein", "Functional"}


# --- 2. human GRCh37 ---------------------------------------------------------


def test_human_grch37_has_37_38_options_but_not_38_only():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh37.p13"
    )
    ids = panel_ids(panels)

    assert ALWAYS_VISIBLE_PANEL_IDS <= ids
    assert HUMAN_37_38_PANEL_IDS <= ids
    assert "allele_frequencies" not in ids  # GRCh38-only

    opts = option_ids(panels)
    assert "utrannotator" in opts
    assert opts.isdisjoint(GRCH38_ONLY_OPTION_IDS)


# --- 3. non-human / non-GRCh37-38 -------------------------------------------


def test_mouse_shows_only_always_visible_panels():
    panels = get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    assert panel_ids(panels) == ALWAYS_VISIBLE_PANEL_IDS

    genes_opts = option_ids(
        [p for p in panels if p["id"] == "genes_and_transcripts"]
    )
    assert "utrannotator" not in genes_opts
    assert "riboseqorfs" not in genes_opts


def test_human_t2t_is_not_treated_as_grch37_38():
    # human taxonomy but a non-GRCh37/38 assembly gets only the base panels
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="T2T-CHM13v2.0"
    )
    assert panel_ids(panels) == ALWAYS_VISIBLE_PANEL_IDS


def test_no_species_info_defaults_to_base_panels():
    assert panel_ids(get_visible_panels()) == ALWAYS_VISIBLE_PANEL_IDS


# --- 4. deep-copy isolation (guards the earlier shared-reference bug) --------


def test_calls_return_equal_but_independent_structures():
    a = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    b = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")

    assert a == b
    assert a is not b

    # mutating one result must not affect the other
    a[0]["options"].append({"id": "injected"})
    assert a != b


def test_module_constants_are_not_mutated_between_calls():
    # a GRCh38 call mutates its *copy* of the always-visible panels (adds
    # categories, extra options); the module constants must be untouched.
    get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")

    assert len(form_panels._ALWAYS_VISIBLE_PANELS) == 3
    protein_panel = next(
        p
        for p in form_panels._ALWAYS_VISIBLE_PANELS
        if p["id"] == "protein_and_functional"
    )
    protein_option = protein_panel["options"][0]
    assert protein_option["id"] == "protein"
    # GRCh38 adds category="Protein" to the copy, not the constant
    assert "category" not in protein_option
    # ...and does not leave ProtVar/MaveDB/etc. on the constant
    assert [o["id"] for o in protein_panel["options"]] == ["protein"]


# --- 5. id contract: option ids are ConfigIniParams parameters --------------


def test_option_ids_are_valid_config_ini_parameters():
    # The GRCh38 set is the superset of every option/sub-option.
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    config_fields = set(ConfigIniParams.model_fields)

    # locked_children (hgvs_c/hgvs_p) are display-only, not parameters, so they
    # are excluded here; every real option/sub-option id must be a parameter.
    for option_id in option_ids(panels, include_sub_options=True):
        assert option_id in config_fields, f"{option_id} is not a ConfigIniParams field"
