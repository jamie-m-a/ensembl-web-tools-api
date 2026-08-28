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
from app.vep.submission_options import submittable_options

HUMAN = "9606"
MOUSE = "10090"

ALWAYS_VISIBLE_PANEL_IDS = {
    "variant_representations",
    "genes_and_transcripts",
    "protein_and_functional",
}
HUMAN_37_38_PANEL_IDS = {
    "variant_impact_predictions",
    "phenotype_and_disease_associations",
}
GRCH38_ONLY_OPTION_IDS = {
    "eve",
    "gerp",
    "mavedb",
    "opentargets",
    "protvar",
    "riboseqorfs",
}

# Reuse-tier options available for human GRCh37 as well as GRCh38 (data exists
# for both). go/phenotypes are really multi-species, gated on human for now.
HUMAN_37_38_EXTRA_OPTION_IDS = {"go", "phenotypes", "intact", "clinvar_sv"}


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
    assert categories(panels, "variant_impact_predictions") == {
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
    # allele_frequencies is present for GRCh37 too now (gnomAD v2 exomes/genomes),
    # but only those two sources — not the GRCh38-only allofus/cnv.
    assert "allele_frequencies" in ids

    opts = option_ids(panels)
    assert "utrannotator" in opts
    # the reuse-tier extras (go / phenotypes / intact / clinvar_sv) are offered
    assert HUMAN_37_38_EXTRA_OPTION_IDS <= opts
    # gnomAD v2 exomes/genomes/SV, but not the GRCh38-only AF sources (allofus/cnv)
    assert {"gnomad_exomes", "gnomad_genomes", "gnomad_sv"} <= opts
    assert opts.isdisjoint({"allofus", "gnomad_cnv"})
    # but none of the GRCh38-only ones
    assert opts.isdisjoint(GRCH38_ONLY_OPTION_IDS)


# --- 3. non-human / non-GRCh37-38 -------------------------------------------


def test_mouse_gets_the_base_panels_plus_its_own_data_options():
    """Mouse carries GO and Phenotypes data files, so it is offered those two on
    top of the always-visible panels — and none of the human-only options."""
    panels = get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    assert panel_ids(panels) == ALWAYS_VISIBLE_PANEL_IDS | {"phenotype_and_disease_associations"}

    genes_opts = option_ids(
        [p for p in panels if p["id"] == "genes_and_transcripts"]
    )
    assert "go" in genes_opts
    assert "utrannotator" not in genes_opts
    assert "riboseqorfs" not in genes_opts

    associations = option_ids([p for p in panels if p["id"] == "phenotype_and_disease_associations"])
    assert associations == {"phenotypes"}  # not geno2mp / clinvar / opentargets


def test_a_species_with_no_data_files_shows_only_the_always_visible_panels():
    panels = get_visible_panels(species_taxonomy_id="1", assembly_name="Wibble_v1")
    assert panel_ids(panels) == ALWAYS_VISIBLE_PANEL_IDS
    assert "go" not in option_ids(panels)
    assert "phenotypes" not in option_ids(panels)


def test_a_go_only_species_is_not_offered_phenotypes():
    # platypus has a GO file but no phenotypes file
    panels = get_visible_panels(species_taxonomy_id="9258", assembly_name="mOrnAna1.p.v1")
    ids = option_ids(panels)
    assert "go" in ids and "phenotypes" not in ids
    assert "phenotype_and_disease_associations" not in panel_ids(panels)


# Every dataset a species row can name; the option id is the dataset name.
DATASET_OPTION_IDS = {"go", "phenotypes", "cadd"}


def test_form_options_match_the_spec_a_submission_would_get():
    """The form and the spec loader read the same table, so what a species is
    offered and what its submission actually configures cannot drift."""
    from app.vep.utils.spec_loader import _species_annotations, resolve_merged_spec

    for row in _species_annotations()["species"]:
        panels = get_visible_panels(
            species_taxonomy_id=row["species_taxonomy_id"],
            assembly_name=row["assembly"],
        )
        offered = {i for i in option_ids(panels) if i in DATASET_OPTION_IDS}
        configured = {
            e.id
            for e in resolve_merged_spec(row["assembly"]).config.entries
            if e.id in DATASET_OPTION_IDS
        }
        assert offered == configured == set(row["datasets"]), row["assembly"]


def test_updownstream_distance_available_for_all_species():
    # A base (all-species) Genes & transcripts option with a bounded numeric field.
    for taxon, assembly in [(HUMAN, "GRCh38.p14"), (MOUSE, "GRCm39"), (HUMAN, "GRCh37.p13")]:
        panels = get_visible_panels(species_taxonomy_id=taxon, assembly_name=assembly)
        genes = next(p for p in panels if p["id"] == "genes_and_transcripts")
        option = next(o for o in genes["options"] if o["id"] == "updownstream_distance")
        field = next(
            s for s in option["sub_options"] if s["id"] == "updownstream_distance_bp"
        )
        assert field["type"] == "number"
        assert (field["default"], field["min"], field["max"]) == (5000, 0, 1000000)


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
    """The panel registry must not accumulate options.

    It used to be possible to mutate the shared literals — a GRCh38 call added
    ProtVar and friends to its *copy*, and a missing `deepcopy` would have left
    them on the constant for the next species. The registry now carries no
    options at all and each call builds its own lists, so the bug is gone by
    construction; this is what says so.
    """
    get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")

    assert all("options" not in panel for panel in form_panels._PANELS)
    assert [panel["id"] for panel in form_panels._PANELS][:2] == [
        "variant_representations",
        "variant_impact_predictions",
    ]


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


def test_gnomad_exomes_absent_for_unspecced_assembly():
    # v4 on GRCh38, v2 on GRCh37; nothing on an assembly without a spec.
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="T2T-CHM13v2.0")
    assert "gnomad_exomes" not in option_ids(panels)


def test_gnomad_v2_exomes_structure_grch37():
    # GRCh37 gets the v2 shape: a Subset group + a Genetic-ancestry group with a
    # plain popmax row, sex splits, and NFE/EAS sub-populations. No UK-Biobank.
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh37.p13")
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    exomes = next(o for o in af["options"] if o["id"] == "gnomad_exomes")

    sub_ids = [s.get("id") for s in exomes["sub_options"]]
    assert "gnomad_exomes_include_ukb" not in sub_ids  # v2 has no UK-Biobank toggle

    subset_grp = next(s for s in exomes["sub_options"] if s.get("label") == "Subset")
    assert [o["id"] for o in subset_grp["options"]] == [
        f"gnomad_exomes_subset_{s}"
        for s in ["full", "controls", "non_neuro", "non_topmed", "non_cancer"]
    ]
    assert subset_grp["options"][0]["default"] is True  # Full dataset pre-selected

    anc_grp = next(
        s for s in exomes["sub_options"] if s.get("label") == "Genetic ancestry group"
    )
    anc_ids = [o["id"] for o in anc_grp["options"]]
    assert anc_ids == [
        f"gnomad_exomes_{a}"
        for a in ["all", "popmax", "afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
    ]

    all_anc = anc_grp["options"][0]
    assert all_anc["default"] is True  # "All" pre-selected -> fields=AF baseline
    popmax = anc_grp["options"][1]
    assert "sub_options" not in popmax  # popmax is a plain toggle

    nfe = next(o for o in anc_grp["options"] if o["id"] == "gnomad_exomes_nfe")
    subpop_grp = next(s for s in nfe["sub_options"] if s.get("type") == "group")
    assert [o["id"] for o in subpop_grp["options"]] == [
        f"gnomad_exomes_nfe_{sp}"
        for sp in ["seu", "bgr", "onf", "swe", "nwe", "est"]
    ]


def test_gnomad_genomes_structure_grch38():
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    genomes = next(o for o in af["options"] if o["id"] == "gnomad_genomes")

    # no UK Biobank toggle for genomes; only the ancestry group
    assert all(s.get("type") == "group" for s in genomes["sub_options"])
    group = genomes["sub_options"][0]
    # grpmax sits directly under "All": the two summary figures together, ahead
    # of the individual ancestries.
    assert [o["id"] for o in group["options"]] == [
        f"gnomad_genomes_{a}"
        for a in [
            "all", "grpmax", "afr", "amr", "asj", "eas", "fin", "mid", "nfe",
            "ami", "remaining",
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


def test_gnomad_genomes_absent_for_unspecced_assembly():
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="T2T-CHM13v2.0")
    assert "gnomad_genomes" not in option_ids(panels)


def test_gnomad_v2_genomes_structure_grch37():
    # Genomes v2 drops SAS, the EAS sub-populations, non_cancer, and two of the
    # NFE sub-populations relative to exomes.
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh37.p13")
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    genomes = next(o for o in af["options"] if o["id"] == "gnomad_genomes")

    subset_grp = next(s for s in genomes["sub_options"] if s.get("label") == "Subset")
    assert [o["id"] for o in subset_grp["options"]] == [
        f"gnomad_genomes_subset_{s}"
        for s in ["full", "controls", "non_neuro", "non_topmed"]  # no non_cancer
    ]

    anc_grp = next(
        s for s in genomes["sub_options"] if s.get("label") == "Genetic ancestry group"
    )
    anc_ids = [o["id"] for o in anc_grp["options"]]
    assert "gnomad_genomes_sas" not in anc_ids  # no South Asian
    assert anc_ids == [
        f"gnomad_genomes_{a}"
        for a in ["all", "popmax", "afr", "amr", "asj", "eas", "fin", "nfe", "oth"]
    ]

    eas = next(o for o in anc_grp["options"] if o["id"] == "gnomad_genomes_eas")
    assert all(s.get("type") != "group" for s in eas["sub_options"])  # no EAS sub-pops
    nfe = next(o for o in anc_grp["options"] if o["id"] == "gnomad_genomes_nfe")
    subpop_grp = next(s for s in nfe["sub_options"] if s.get("type") == "group")
    assert [o["id"] for o in subpop_grp["options"]] == [
        f"gnomad_genomes_nfe_{sp}" for sp in ["seu", "onf", "nwe", "est"]
    ]


def test_gnomad_sv_v2_structure_grch37():
    # GRCh37 SV is v2: 5 continental populations (afr/amr/eas/eur/oth) — broader
    # and fewer than v4.1's 11 — plain toggles + the overlap-cutoff select, no sex.
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh37.p13")
    af = next(p for p in panels if p["id"] == "allele_frequencies")
    sv = next(o for o in af["options"] if o["id"] == "gnomad_sv")
    assert sv["label"] == "gnomAD SV v2.1"

    assert any(
        s.get("id") == "gnomad_sv_overlap_cutoff" and s["type"] == "select"
        for s in sv["sub_options"]
    )
    group = next(s for s in sv["sub_options"] if s.get("type") == "group")
    assert [o["id"] for o in group["options"]] == [
        "gnomad_sv_af",
        *(f"gnomad_sv_af_{p}" for p in ["afr", "amr", "eas", "eur", "oth"]),
    ]
    assert all("sub_options" not in o for o in group["options"])  # no sex splits
    # overall AF pre-selected, populations opt-in; the v4-only pops aren't offered
    af_opts = {o["id"]: o for o in group["options"]}
    assert af_opts["gnomad_sv_af"]["default"] is True
    assert af_opts["gnomad_sv_af_eur"]["default"] is False
    ids = option_ids(panels)
    assert "gnomad_sv_af_eur" in ids and "gnomad_sv_af_oth" in ids
    assert ids.isdisjoint({"gnomad_sv_af_ami", "gnomad_sv_af_nfe", "gnomad_sv_af_sas"})


def test_allele_frequency_sources_are_categorised():
    # The AF sources carry category headings so the form can lay them out side by
    # side under "Short variants" / "Structural variants" instead of stacking
    # every source in one very tall column.
    expected = {
        "GRCh38.p14": {
            "gnomad_exomes": "Short variants",
            "gnomad_genomes": "Short variants",
            "allofus": "Short variants",
            "gnomad_sv": "Structural variants",
            "gnomad_cnv": "Structural variants",
        },
        "GRCh37.p13": {
            "gnomad_exomes": "Short variants",
            "gnomad_genomes": "Short variants",
            "gnomad_sv": "Structural variants",
        },
    }
    for assembly, categories in expected.items():
        panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name=assembly)
        af = next(p for p in panels if p["id"] == "allele_frequencies")
        assert {o["id"]: o["category"] for o in af["options"]} == categories
        # ...and they stay in category order, so each heading owns one run of
        # options (groupByCategory preserves first-seen order).
        seen = [o["category"] for o in af["options"]]
        assert seen == sorted(seen, key=["Short variants", "Structural variants"].index)


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


def _clinvar_option(assembly):
    panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name=assembly)
    va = next(p for p in panels if p["id"] == "phenotype_and_disease_associations")
    return next(o for o in va["options"] if o["id"] == "clinvar_sv")


def test_clinvar_structural_is_its_own_option_not_a_master_and_a_child():
    """It used to sit under a "Clinical Significance (ClinVar)" master. Once the
    germline data moved to Phenotypes, that master gated nothing else and
    ticking it ran nothing, so the one real control stands on its own."""
    for assembly in ("GRCh37.p13", "GRCh38.p14"):
        assert _clinvar_option(assembly)["label"] == "ClinVar structural variants"
        panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name=assembly)
        assert "clinvar" not in option_ids(panels)


def test_clinvar_germline_has_no_form_control():
    """The germline data is served under Phenotypes, which turns its custom on
    behind the scenes — so there is nothing to tick for it here. Offering a
    control that Phenotypes overrides would only misreport what ran."""
    for assembly in ("GRCh38.p14", "GRCh37.p13"):
        panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name=assembly)
        assert "clinvar_short" not in option_ids(panels)


def test_clinvar_structural_available_for_both_human_assemblies():
    for assembly in ("GRCh38.p14", "GRCh37.p13"):
        assert "sub_options" not in _clinvar_option(assembly)


def test_clinvar_absent_for_non_human():
    panels = get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    assert "clinvar" not in option_ids(panels)


def test_regulatory_panel_is_grch38_only():
    g38 = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    assert "regulatory" in panel_ids(g38)
    assert "gencode_promoters" in option_ids(g38)
    # Not for human GRCh37 (no spec) nor other species.
    assert "regulatory" not in panel_ids(
        get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh37")
    )
    assert "regulatory" not in panel_ids(
        get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    )


def test_gnomad_sv_option_is_grch38_allele_frequency():
    g38 = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    af = next(p for p in g38 if p["id"] == "allele_frequencies")
    sv = next((o for o in af["options"] if o["id"] == "gnomad_sv"), None)
    assert sv is not None
    ids = option_ids(g38)
    assert {
        "gnomad_sv",
        "gnomad_sv_af",
        "gnomad_sv_af_afr",
        "gnomad_sv_af_sas",
        "gnomad_sv_overlap_cutoff",
    } <= ids
    # overall AF is pre-selected; the populations are opt-in
    af_opts = {
        opt["id"]: opt
        for grp in sv["sub_options"]
        if grp.get("type") == "group"
        for opt in grp["options"]
    }
    assert af_opts["gnomad_sv_af"]["default"] is True
    assert af_opts["gnomad_sv_af_afr"]["default"] is False
    # GRCh38-only (mouse has no allele_frequencies panel)
    assert "gnomad_sv" not in option_ids(
        get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    )


def test_gnomad_cnv_option_is_grch38_allele_frequency():
    g38 = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    ids = option_ids(g38)
    # SF prefix, "remaining" spelled out, no ami
    assert {
        "gnomad_cnv",
        "gnomad_cnv_sf",
        "gnomad_cnv_sf_remaining",
        "gnomad_cnv_overlap_cutoff",
    } <= ids
    assert "gnomad_cnv_sf_ami" not in ids
    assert "gnomad_cnv" not in option_ids(
        get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    )


def test_option_ids_round_trip_into_a_submission():
    """Every id the form offers is one a submission can set, and comes back out.

    This used to assert each id was a `ConfigIniParams` field — one per option,
    199 of them. The options are now a map validated against the spec for the
    submission's assembly, so the contract is the round trip rather than the
    field list, and it can be checked per assembly instead of against one flat
    set that accepted everything for everyone.
    """
    for assembly in ("GRCh38.p14", "GRCh37.p13"):
        panels = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name=assembly)
        # locked_children (hgvs_c/hgvs_p) are display-only, not parameters.
        offered = set(option_ids(panels, include_sub_options=True))

        known = submittable_options(
            species_taxonomy_id=HUMAN, assembly_name=assembly
        )
        submitted = {
            option_id: (
                True if known[option_id].type is bool else known[option_id].default
            )
            for option_id in offered
        }
        params = ConfigIniParams(
            genome_id="g",
            assembly_name=assembly,
            species_taxonomy_id=HUMAN,
            options=submitted,
        )
        assert offered <= set(params.options), assembly
        assert all(params.options[option_id] == submitted[option_id] for option_id in offered)


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
    assert label("gnomad_exomes", "nfe_XX") == "Non-Finnish European · XX"
    assert label("gnomad_genomes", "afr_XY") == "African & African-American · XY"


def test_af_label_gnomad_bare_sex_is_all_that_sex():
    assert label("gnomad_exomes", "XX") == "All · XX"
    assert label("gnomad_exomes", "XY") == "All · XY"


def test_af_label_gnomad_non_ukb_subset():
    assert label("gnomad_exomes", "non_ukb") == "All · excl. UK Biobank"
    assert label("gnomad_exomes", "non_ukb_afr") == (
        "African & African-American · excl. UK Biobank"
    )
    assert label("gnomad_exomes", "non_ukb_nfe_XX") == (
        "Non-Finnish European · XX · excl. UK Biobank"
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


def test_an_unsent_option_comes_back_at_the_form_s_declared_default():
    """The form default and the submission default are now one statement.

    They used to be two. The form declared `default`, `ConfigIniParams` declared
    a field default beside it, and nothing kept them in step — the All of Us
    "Maximum subpopulation" default was set in the form and had no effect until
    someone noticed. The 199 option fields are gone, so an unsent option is now
    filled from the same `default` the form renders
    (`ConfigIniParams._resolve_options` -> `submission_options.option_values`).

    So the guard moves with the mechanism: read what the form declares, submit
    nothing, and require the option map to come back saying exactly that. The
    complement of `test_option_ids_round_trip_into_a_submission`, which sends
    every option — this one sends none.
    """
    for species, assembly in (
        (HUMAN, "GRCh38.p14"),
        (HUMAN, "GRCh37.p13"),
        (MOUSE, "GRCm39"),
    ):
        declared = {}
        undeclared = []

        def walk(option):
            # A 'group' is a heading around nested controls, with no id or value
            # of its own — and the only node whose `options` are controls. A
            # select's `options` are its choices ({label, value}), not controls.
            if option.get("type") == "group":
                for nested in option["options"]:
                    walk(nested)
                return
            if "default" in option:
                declared[option["id"]] = option["default"]
            else:
                undeclared.append(option["id"])
            for sub in option.get("sub_options", []):
                walk(sub)

        for panel in get_visible_panels(
            species_taxonomy_id=species, assembly_name=assembly
        ):
            for option in panel["options"]:
                walk(option)

        # What this test replaces read its expectations off a field list that
        # later emptied out, so it passed over nothing for a while. If the walk
        # finds no controls, everything below is vacuous too.
        assert declared, f"no options walked out of the form for {assembly}"
        assert not undeclared, (
            f"{assembly}: options the form declares no default for, so a "
            "submission that omits them carries None: " + ", ".join(sorted(undeclared))
        )

        submitted = ConfigIniParams(
            genome_id="g",
            assembly_name=assembly,
            species_taxonomy_id=species,
            options={},
        ).options

        # Compare the type alongside the value: `True == 1` in Python, so a
        # boolean widened to an int would otherwise slip through.
        def typed(value):
            return type(value).__name__, value

        mismatched = [
            f"{option_id}: form={declared[option_id]!r}, submitted="
            + (
                repr(submitted[option_id])
                if option_id in submitted
                else "<not in the map>"
            )
            for option_id in sorted(declared)
            if option_id not in submitted
            or typed(submitted[option_id]) != typed(declared[option_id])
        ]
        assert not mismatched, (
            f"{assembly}: a submission that sends nothing does not come back at "
            "the defaults the form declares:\n  " + "\n  ".join(mismatched)
        )
        assert set(submitted) == set(declared), assembly


def test_panels_come_back_in_the_agreed_order():
    """One order, stated once, for both surfaces: the input form and the results
    annotation panel each render the list this returns."""
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    assert [panel["label"] for panel in panels] == [
        "Variant representations",
        "Variant impact predictions",
        "Allele frequencies",
        "Genes & transcripts",
        "Protein & functional",
        "Regulatory",
        "Phenotype & disease associations",
    ]


def test_a_species_with_fewer_panels_keeps_the_same_relative_order():
    """Ordering is applied to whatever a species actually has, so a genome
    missing several panels is not silently reshuffled."""
    panels = get_visible_panels(assembly_name="GRCg6a")  # chicken: no AF panel
    assert [panel["label"] for panel in panels] == [
        "Variant representations",
        "Variant impact predictions",
        "Genes & transcripts",
        "Protein & functional",
        "Phenotype & disease associations",
    ]


def test_opentargets_follows_phenotypes_where_both_exist():
    """Phenotypes leads because every genome that has associations has it;
    OpenTargets is GRCh38-only and slots in after."""
    grch38 = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    associations = next(
        panel
        for panel in grch38
        if panel["id"] == "phenotype_and_disease_associations"
    )
    ids = [option["id"] for option in associations["options"]]
    assert ids.index("phenotypes") < ids.index("opentargets")

    # GRCh37 has Phenotypes and no OpenTargets, and is unaffected
    grch37 = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh37.p13")
    associations_37 = next(
        panel
        for panel in grch37
        if panel["id"] == "phenotype_and_disease_associations"
    )
    ids_37 = [option["id"] for option in associations_37["options"]]
    assert "phenotypes" in ids_37 and "opentargets" not in ids_37


def test_genes_and_transcripts_options_are_grouped_into_three_categories():
    """Every option in the panel carries a category, so it renders as three
    sub-headings rather than an unlabelled run followed by one heading.

    Constraint used to be a panel of its own ("Conservation & constraint"); it is
    now one of the three, grouped the way Variant impact predictions groups
    Missense / Splicing / Genome wide. It lost the "Conservation" half of its
    name when GERP moved to Genome wide: what is left measures a gene's
    tolerance to variation, not a position's conservation.
    """
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    assert "conservation_and_constraint" not in {panel["id"] for panel in panels}

    genes = next(panel for panel in panels if panel["id"] == "genes_and_transcripts")

    # Grouped by category in first-seen order, which is what the frontend's
    # `groupByCategory` does for both the form and the results panel.
    grouped: dict[str, list[str]] = {}
    for option in genes["options"]:
        grouped.setdefault(option.get("category"), []).append(option["id"])

    assert None not in grouped, "every option should be categorised"
    assert list(grouped) == ["Locations", "Additional molecular consequence predictions", "Constraint"]
    assert grouped["Locations"] == [
        "tss_distance",
        "nearest_gene",
        "nearest_exon_jb",
        "updownstream_distance",
    ]
    assert grouped["Additional molecular consequence predictions"] == [
        "utrannotator",
        "nmd",
        "riboseqorfs",
        "go",
    ]
    # pLI is GRCh38-only, so it arrives with the GRCh38 additions and lands
    # after the two the 37/38 tier contributes — same group, appended.
    assert grouped["Constraint"] == ["loeuf", "dosage_sensitivity", "pli"]


def test_gerp_sits_with_cadd_under_genome_wide():
    """GERP scores a position's conservation across species, which is a
    genome-wide measure like CADD — not a gene's tolerance to variation, which
    is what the Constraint group it used to sit in now means."""
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    impact = next(
        panel for panel in panels if panel["id"] == "variant_impact_predictions"
    )
    grouped: dict[str, list[str]] = {}
    for option in impact["options"]:
        grouped.setdefault(option.get("category"), []).append(option["id"])

    assert grouped["Genome wide"] == ["avi", "cadd", "gerp"]
    genes = next(panel for panel in panels if panel["id"] == "genes_and_transcripts")
    assert "gerp" not in {option["id"] for option in genes["options"]}


def test_a_species_with_fewer_options_keeps_the_same_category_names():
    """The categories are declared across three tiers (base / human 37-38 /
    GRCh38-only), so a species that gets only some of the options must still see
    those it does get under the same headings — and no empty ones."""
    panels = get_visible_panels(species_taxonomy_id=MOUSE, assembly_name="GRCm39")
    genes = next(panel for panel in panels if panel["id"] == "genes_and_transcripts")
    grouped: dict[str, list[str]] = {}
    for option in genes["options"]:
        grouped.setdefault(option.get("category"), []).append(option["id"])

    assert None not in grouped
    # Mouse has GO but none of the human-only annotations or conservation data.
    assert list(grouped) == ["Locations", "Additional molecular consequence predictions"]
    assert grouped["Additional molecular consequence predictions"] == ["go"]


def test_a_non_human_species_gets_neither_the_panel_nor_the_options():
    """Conservation data is human-only, and moving it inside another panel must
    not leak it to species that have none."""
    panels = get_visible_panels(assembly_name="GRCg6a")  # chicken
    genes = next(panel for panel in panels if panel["id"] == "genes_and_transcripts")
    option_ids = {option["id"] for option in genes["options"]}
    assert "loeuf" not in option_ids and "dosage_sensitivity" not in option_ids
    assert "conservation_and_constraint" not in {panel["id"] for panel in panels}


# --- options declared by their own config entry ----------------------------
#
# Two entries carry a `form` block and are placed from it; form_panels still
# writes out the other 33. See docs/form-panels-to-json.md.

AF_OPTION_IDS = {
    "gnomad_exomes",
    "gnomad_genomes",
    "allofus",
    "gnomad_sv",
    "gnomad_cnv",
}

GOLDEN_CASES = {
    "human_grch38": {"species_taxonomy_id": HUMAN, "assembly_name": "GRCh38.p14"},
    "human_grch37": {"species_taxonomy_id": HUMAN, "assembly_name": "GRCh37.p13"},
    "mouse": {"species_taxonomy_id": MOUSE, "assembly_name": "GRCm39"},
    "unlisted": {"species_taxonomy_id": "7955", "assembly_name": "GRCz11"},
}


def _golden() -> dict:
    import json
    from pathlib import Path

    return json.loads(
        (Path(__file__).parent / "form_panels.golden.json").read_text()
    )


def test_the_spec_really_is_declaring_those_options():
    """Guards the golden comparison from passing for the wrong reason.

    If the `form` blocks stopped being read the options would simply vanish, and
    the golden test would catch it — but it would report "an option is missing"
    rather than "the spec is not being read", which is a slower thing to
    diagnose. This says which.
    """
    declared = {option_id for option_id, _, _ in _declared_ids("GRCh38.p14")}
    # Every option the form shows except the allele frequencies, which are still
    # generated from their ancestry tables (see docs/form-panels-to-json.md).
    # Counted from the golden rather than written down, so it cannot go stale.
    shown = {
        option["id"]
        for panel in _golden()["human_grch38"]
        for option in panel["options"]
    }
    # Every one of them now — the allele frequencies included, their sub-option
    # trees grown from the same `fields=` tables that write their config lines.
    assert declared == shown
    assert {"hgvs", "pli", "mutfunc", "protein"} | AF_OPTION_IDS <= declared

    # `pli` and the rest of the GRCh38-only set are declared only by
    # human_grch38.json, so GRCh37 does not see them — the same selection that
    # drops their parse plugins and display options.
    grch37 = {i for i, _, _ in _declared_ids("GRCh37.p13")}
    assert "pli" not in grch37 and "eve" not in grch37
    assert {"hgvs", "nearest_exon_jb", "loeuf"} <= grch37


def _declared_ids(assembly: str):
    """(id, panel, order) for each entry that declares a control."""
    from app.vep.utils.spec_loader import resolve_merged_spec

    return [
        (entry.id, entry.form.panel, entry.form.order)
        for entry in resolve_merged_spec(assembly).config.entries
        if entry.form is not None
    ]


def test_panels_are_unchanged_by_placing_options_from_the_spec():
    """The whole point: an option declared by its config entry lands exactly
    where the hand-written list used to put it.

    Against a golden file captured before the two options moved, across the four
    paths through `get_visible_panels` — human GRCh38, human GRCh37, a species
    with its own annotation data, and one with none.
    """
    golden = _golden()
    for name, kwargs in GOLDEN_CASES.items():
        assert get_visible_panels(**kwargs) == golden[name], name


def test_a_declared_option_sits_between_the_coded_ones():
    """`nearest_exon_jb` is third of four in Locations — not first or last, so
    its `order` has to place it *between* options this module still writes.

    Appending would have been enough for `pli` alone, which is last; this is the
    case that proves the ordering actually orders.
    """
    genes = next(
        panel
        for panel in get_visible_panels(
            species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
        )
        if panel["id"] == "genes_and_transcripts"
    )
    ids = [option["id"] for option in genes["options"]]
    assert ids.index("nearest_gene") < ids.index("nearest_exon_jb")
    assert ids.index("nearest_exon_jb") < ids.index("updownstream_distance")
    # ...and pli is last, after the GRCh38 additions.
    assert ids[-1] == "pli"


def test_an_option_naming_an_unshown_panel_is_an_error(monkeypatch):
    """A control whose panel this genome does not show would otherwise vanish
    without a word — the silent-drop failure this spec keeps having to guard
    against."""
    import pytest

    monkeypatch.setattr(
        form_panels,
        "_spec_form_options",
        lambda _assembly: [("no_such_panel", 10, {"id": "x"})],
    )
    with pytest.raises(ValueError, match="no_such_panel"):
        get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")


def test_an_entry_with_no_control_declares_no_form_block():
    """`clinvar_short` is forced on by Phenotypes and `hgvsg` is hidden pending
    chromosome synonyms — neither has a control, and the absence of a `form`
    block is how that is said. A stray one would put an unusable toggle on the
    form."""
    from app.vep.utils.spec_loader import resolve_merged_spec

    entries = {e.id: e for e in resolve_merged_spec("GRCh38.p14").config.entries}
    assert entries["clinvar_short"].form is None
    assert entries["hgvsg"].form is None


# --- the allele-frequency tables now have one source ------------------------


def test_af_labels_come_from_the_same_tables_that_write_the_config_line():
    """The ancestry/population labels are read off each entry's `fields=`
    builder, not a second copy of the list.

    That copy was real and had already drifted — the spec carried `grpmax` and
    the form's table did not — which is what this whole move is for.
    """
    from app.vep.utils.spec_loader import load_merged_spec

    entry = next(
        e
        for e in load_merged_spec("human_grch38").config.entries
        if e.id == "gnomad_genomes"
    )
    ancestries = {a.code: a.label for a in entry.config.fields.ancestries if a.code}
    assert ancestries["afr"] == "African & African-American"
    assert ancestries["grpmax"] == "Maximum across all groups"

    # The decode the results metadata uses reads exactly those.
    assert form_panels.af_population_label("gnomad_genomes", "afr") == (
        "African & African-American"
    )
    assert form_panels.af_population_label("gnomad_genomes", "grpmax") == (
        "Maximum across all groups"
    )
    # ...and the form draws the same label for the same code.
    panels = get_visible_panels(
        species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14"
    )
    labels = {}

    def walk(option):
        if "id" in option:
            labels[option["id"]] = option.get("label")
        for child in option.get("sub_options", []) + option.get("options", []):
            walk(child)

    for panel in panels:
        for option in panel["options"]:
            walk(option)
    assert labels["gnomad_genomes_afr"] == ancestries["afr"]
    assert labels["gnomad_genomes_grpmax"] == ancestries["grpmax"]


def test_the_two_gnomad_sv_sources_keep_their_own_labels():
    """GRCh38's SV populations and GRCh37's share option ids but not labels, and
    the parser reports different codes for them (`afr` against `AFR`).

    Harvesting them into one table silently gave GRCh38 v2.1's wording, which is
    why the population code is stated on the entry rather than derived from the
    option id.
    """
    assert form_panels.af_population_label("gnomad_sv", "afr") == (
        "African & African-American"  # v4.1, GRCh38
    )
    assert form_panels.af_population_label("gnomad_sv", "AFR") == "African"  # v2.1
