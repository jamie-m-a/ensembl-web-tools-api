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
}


def panel_ids(panels):
    return {panel["id"] for panel in panels}


def option_ids(panels, *, include_sub_options=True):
    ids = set()

    def add_option(option):
        ids.add(option["id"])
        if not include_sub_options:
            return
        for sub in option.get("sub_options", []):
            # A 'group' sub-option has no id of its own; recurse into its nested
            # options (e.g. gnomAD exomes' ancestry toggles + their sex options).
            if sub.get("type") == "group":
                for nested in sub["options"]:
                    add_option(nested)
            else:
                ids.add(sub["id"])

    for panel in panels:
        for option in panel["options"]:
            add_option(option)
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
    }
    assert categories(panels, "protein_and_functional") == {"Protein", "Functional"}


def test_maxentscan_and_enformer_are_not_offered():
    # Removed entirely: enabled but never parsed/displayed, so dropped from the
    # human GRCh37/38 pathogenicity panel.
    for assembly in ("GRCh38.p14", "GRCh37.p13"):
        opts = option_ids(
            get_visible_panels(species_taxonomy_id=HUMAN, assembly_name=assembly)
        )
        assert "maxentscan" not in opts
        assert "enformer" not in opts


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


def test_gnomad_exomes_structure_grch38():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    exomes = next(o for o in af["options"] if o["id"] == "gnomad_exomes")

    sub_ids = [s.get("id") for s in exomes["sub_options"]]
    assert "gnomad_exomes_include_ukb" in sub_ids

    group = next(s for s in exomes["sub_options"] if s.get("type") == "group")
    assert group["label"] == "Genetic ancestry group"
    assert [o["id"] for o in group["options"]] == [
        f"gnomad_exomes_{a}"
        for a in ["all", "afr", "amr", "asj", "eas", "fin", "mid", "nfe"]
    ]

    all_ancestry = group["options"][0]
    assert all_ancestry["default"] is True  # "All" pre-selected
    assert [s["id"] for s in all_ancestry["sub_options"]] == [
        "gnomad_exomes_all_both",
        "gnomad_exomes_all_female",
        "gnomad_exomes_all_male",
    ]
    both, female, male = all_ancestry["sub_options"]
    assert both["default"] is True  # combined sexes on by default
    assert female["default"] is False and male["default"] is False


def test_gnomad_exomes_absent_below_grch38():
    for assembly in ("GRCh37.p13", "T2T-CHM13v2.0"):
        panels = get_visible_panels(
            species_taxonomy_id=HUMAN, assembly_name=assembly
        )
        assert "gnomad_exomes" not in option_ids(panels)


def test_gnomad_genomes_structure_grch38():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    genomes = next(o for o in af["options"] if o["id"] == "gnomad_genomes")

    # no UK Biobank toggle for genomes; only the ancestry group
    assert all(s.get("type") == "group" for s in genomes["sub_options"])
    group = genomes["sub_options"][0]
    assert [o["id"] for o in group["options"]] == [
        f"gnomad_genomes_{a}"
        for a in [
            "all", "afr", "amr", "asj", "eas", "fin", "mid", "nfe",
            "ami", "remaining", "grpmax",
        ]
    ]

    # grpmax is a plain toggle: no sex sub-options
    grpmax = next(o for o in group["options"] if o["id"] == "gnomad_genomes_grpmax")
    assert "sub_options" not in grpmax

    # the other ancestries carry Both/Female/Male
    ami = next(o for o in group["options"] if o["id"] == "gnomad_genomes_ami")
    assert [s["id"] for s in ami["sub_options"]] == [
        "gnomad_genomes_ami_both",
        "gnomad_genomes_ami_female",
        "gnomad_genomes_ami_male",
    ]


def test_gnomad_genomes_absent_below_grch38():
    for assembly in ("GRCh37.p13", "T2T-CHM13v2.0"):
        panels = get_visible_panels(
            species_taxonomy_id=HUMAN, assembly_name=assembly
        )
        assert "gnomad_genomes" not in option_ids(panels)


def test_allofus_structure_grch38():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    allofus = next(o for o in af["options"] if o["id"] == "allofus")

    group = allofus["sub_options"][0]
    assert group["type"] == "group"
    assert "label" not in group  # no heading
    assert [o["id"] for o in group["options"]] == [
        f"allofus_{p}"
        for p in ["all", "max", "afr", "amr", "eas", "eur", "mid", "sas", "oth"]
    ]
    # population toggles are plain booleans (no sex sub-options)
    assert all("sub_options" not in o for o in group["options"])
    # "All" pre-selected
    all_pop = next(o for o in group["options"] if o["id"] == "allofus_all")
    assert all_pop["default"] is True


def test_allofus_absent_below_grch38():
    for assembly in ("GRCh37.p13", "T2T-CHM13v2.0"):
        panels = get_visible_panels(
            species_taxonomy_id=HUMAN, assembly_name=assembly
        )
        assert "allofus" not in option_ids(panels)


def test_clinvar_in_variant_associations_for_grch37_and_grch38():
    for assembly in ("GRCh37.p13", "GRCh38.p14"):
        panels = get_visible_panels(
            species_taxonomy_id=HUMAN, assembly_name=assembly
        )
        va = next(p for p in panels if p["id"] == "variant_associations")
        assert "clinvar" in [o["id"] for o in va["options"]]


def test_clinvar_absent_for_non_human():
    panels = get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    assert "clinvar" not in option_ids(panels)


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


# ---------------------------------------------------------------------------
# AF population-code -> form-label decoders (af_population_label /
# af_max_subpopulation_label). These decode the population codes the results
# parser emits (see results_filters.af_source_descriptor) back to their form
# labels, reusing the option tuples above as the single source of truth. Ported
# from the frontend's former frequencyPopulationLabels util, which now reads the
# decoded label off the response rather than keeping its own copy of the tables.
# ---------------------------------------------------------------------------

label = form_panels.af_population_label


def test_af_label_gnomad_bare_ancestry():
    assert label("gnomad_exomes", "afr") == "African & African-American"
    assert label("gnomad_exomes", "nfe") == "Non-Finnish European"


def test_af_label_gnomad_sex_suffix():
    assert label("gnomad_exomes", "nfe_XX") == "Non-Finnish European · Female"
    assert label("gnomad_genomes", "afr_XY") == "African & African-American · Male"


def test_af_label_gnomad_bare_sex_is_all_that_sex():
    assert label("gnomad_exomes", "XX") == "All · Female"
    assert label("gnomad_exomes", "XY") == "All · Male"


def test_af_label_gnomad_non_ukb_subset():
    assert label("gnomad_exomes", "non_ukb") == "All · excl. UK Biobank"
    assert label("gnomad_exomes", "non_ukb_afr") == (
        "African & African-American · excl. UK Biobank"
    )
    assert label("gnomad_exomes", "non_ukb_nfe_XX") == (
        "Non-Finnish European · Female · excl. UK Biobank"
    )


def test_af_label_gnomad_grpmax_and_genomes_only():
    assert label("gnomad_genomes", "grpmax") == "Maximum across all groups"
    assert label("gnomad_genomes", "ami") == "Amish"
    assert label("gnomad_genomes", "remaining") == "Remaining"


def test_af_label_allofus_flat_codes():
    assert label("all_of_us", "afr") == "African"
    assert label("all_of_us", "amr") == "Latino/Ad Mixed American"
    assert label("all_of_us", "eur") == "European"
    assert label("all_of_us", "max") == "Maximum subpopulation"


def test_af_label_overall_is_all_for_every_source():
    assert label("gnomad_exomes", "") == "All"
    assert label("gnomad_genomes", "") == "All"
    assert label("all_of_us", "") == "All"


def test_af_label_unrecognised_falls_back_to_code():
    assert label("gnomad_exomes", "zzz") == "zzz"
    assert label("all_of_us", "zzz") == "zzz"


def test_af_max_subpopulation_label_single_and_joined():
    assert form_panels.af_max_subpopulation_label("eur") == "European"
    assert form_panels.af_max_subpopulation_label("eur&afr") == "European / African"
