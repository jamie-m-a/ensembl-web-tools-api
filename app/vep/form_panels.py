"""Which option panels/options are visible on the VEP input form for a genome.

Returned by the form_config endpoint (called on species selection). For now the
set is the same for every species ("always visible"); species-conditional rules
will be layered on later by inspecting the genome metadata attributes.

Option (and sub-option) `id`s match the ConfigIniParams parameter names, so the
form's selections round-trip back into the generated config.ini. Options may
carry a `category` label which the form uses to group them within a panel.
"""

import copy

# Always-visible panels/options.
_ALWAYS_VISIBLE_PANELS: list[dict] = [
    {
        "id": "variant_representations",
        "label": "Variant representations",
        "options": [
            # HGVS renders as a single control with linked HGVSc/HGVSp (the
            # `hgvs` param) and a separate HGVSg (the `hgvsg` param), both off by
            # default. The frontend builds the linked UI; the panel just carries
            # the `hgvs` option (default off).
            {"id": "hgvs", "label": "HGVS", "type": "boolean", "default": False},
            {"id": "spdi", "label": "SPDI", "type": "boolean", "default": False},
        ],
    },
    {
        "id": "genes_and_transcripts",
        "label": "Genes & transcripts",
        "options": [
            {
                "id": "tss_distance",
                "label": "Distance to TSS",
                "type": "boolean",
                "default": False,
            },
            {
                "id": "nearest_gene",
                "label": "Nearest gene",
                "type": "boolean",
                "default": False,
                "sub_options": [
                    {
                        "id": "nearest_gene_both_directions",
                        "label": "Both directions",
                        "type": "boolean",
                        "default": False,
                    },
                ],
            },
            {
                "id": "nearest_exon_jb",
                "label": "Nearest exon junction boundary",
                "type": "boolean",
                "default": False,
                "sub_options": [
                    {
                        "id": "nearest_exon_jb_max_range",
                        "label": "Max search range (bp)",
                        "type": "number",
                        "default": 10000,
                    },
                    {
                        "id": "nearest_exon_jb_intronic",
                        "label": "Intronic",
                        "type": "boolean",
                        "default": False,
                    },
                ],
            },
            {
                # Up/downstream distance for consequence calling (VEP `distance`).
                # A toggle revealing a numeric field (bp) that overrides VEP's
                # default of 5000; no output, nothing parsed.
                "id": "updownstream_distance",
                "label": "Up/downstream distance",
                "type": "boolean",
                "default": False,
                "sub_options": [
                    {
                        "id": "updownstream_distance_bp",
                        "label": "Distance (bp)",
                        "type": "number",
                        "default": 5000,
                        "min": 0,
                        "max": 1000000,
                    },
                ],
            },
        ],
    },
    {
        "id": "protein_and_functional",
        "label": "Protein & functional",
        "options": [
            {"id": "protein", "label": "Protein ID", "type": "boolean", "default": False},
        ],
    },
]


# Options added to the existing Genes & transcripts panel for human GRCh37/38.
_HUMAN_37_38_GENES_OPTIONS: list[dict] = [
    {"id": "utrannotator", "label": "UTRAnnotator", "type": "boolean", "default": False},
    {"id": "nmd", "label": "NMD", "type": "boolean", "default": False},
]

# Extra panels shown only for human GRCh37/38. Variant-impact-prediction options
# carry a `category` label used to group them within the panel.
_HUMAN_37_38_PANELS: list[dict] = [
    {
        "id": "pathogenicity_predictions",
        "label": "Variant Impact Predictions",
        "options": [
            {"id": "alphamissense", "label": "AlphaMissense", "type": "boolean", "default": False, "category": "Missense"},
            {"id": "revel", "label": "REVEL", "type": "boolean", "default": False, "category": "Missense"},
            {"id": "clinpred", "label": "ClinPred", "type": "boolean", "default": False, "category": "Missense"},
            {"id": "spliceai", "label": "SpliceAI", "type": "boolean", "default": False, "category": "Splicing"},
            {"id": "cadd", "label": "CADD", "type": "boolean", "default": False, "category": "Genome wide"},
        ],
    },
    {
        "id": "conservation_and_constraint",
        "label": "Conservation & constraint",
        "options": [
            {"id": "loeuf", "label": "LOEUF", "type": "boolean", "default": False},
            {
                "id": "dosage_sensitivity",
                "label": "Dosage sensitivity",
                "type": "boolean",
                "default": False,
                "sub_options": [
                    {"id": "dosage_sensitivity_cover", "label": "Coverage", "type": "boolean", "default": False},
                ],
            },
        ],
    },
    {
        "id": "variant_associations",
        "label": "Variant associations",
        "options": [
            {"id": "geno2mp", "label": "Geno2MP", "type": "boolean", "default": False},
            # ClinVar master: two independent sub-option toggles. "Short variants"
            # (the original ClinVar custom, human 37/38) is here; "Structural
            # variants" (the ClinVar_SV custom) is GRCh38-only and appended in
            # _add_human_grch38_options. Both default off, so enabling the master
            # alone runs nothing until one is picked; each gates its own custom.
            {
                "id": "clinvar",
                "label": "Clinical Significance (ClinVar)",
                "type": "boolean",
                "default": False,
                "sub_options": [
                    {"id": "clinvar_short", "label": "Short variants", "type": "boolean", "default": False},
                ],
            },
        ],
    },
]


def is_human_grch37_or_38(
    species_taxonomy_id: str | None, assembly_name: str | None
) -> bool:
    """True for human GRCh37 / GRCh38."""
    return species_taxonomy_id == "9606" and (assembly_name or "").startswith(
        ("GRCh37", "GRCh38")
    )


def is_human_grch38(
    species_taxonomy_id: str | None, assembly_name: str | None
) -> bool:
    """True for human GRCh38."""
    return species_taxonomy_id == "9606" and (assembly_name or "").startswith(
        "GRCh38"
    )


# Human GRCh38-only sub-options.
_PROTVAR_SUBOPTIONS = [
    {"id": "protvar_stability", "label": "Protein Structure Stability", "type": "boolean", "default": True},
    {"id": "protvar_pocket", "label": "Protein Pockets", "type": "boolean", "default": True},
    {"id": "protvar_int", "label": "Protein-Protein Interaction Interface", "type": "boolean", "default": True},
]
_INTACT_SUBOPTIONS = [
    {"id": "intact_feature_ac", "label": "Feature AC", "type": "boolean", "default": False},
    {"id": "intact_feature_short_label", "label": "Feature short label", "type": "boolean", "default": False},
    {"id": "intact_feature_annotation", "label": "Feature annotation", "type": "boolean", "default": False},
    {"id": "intact_ap_ac", "label": "Affected protein AC", "type": "boolean", "default": False},
    {"id": "intact_interaction_participants", "label": "Interaction participants", "type": "boolean", "default": False},
    {"id": "intact_pmid", "label": "PubMed ID", "type": "boolean", "default": False},
]
_MUTFUNC_SUBOPTIONS = [
    {"id": "mutfunc_motif", "label": "Linear motifs", "type": "boolean", "default": False},
    {"id": "mutfunc_int", "label": "Protein interactions", "type": "boolean", "default": False},
    {"id": "mutfunc_mod", "label": "Protein structure", "type": "boolean", "default": False},
    {"id": "mutfunc_exp", "label": "Protein structure (exp.)", "type": "boolean", "default": False},
]


# gnomAD exomes v4.1 (human GRCh38): a master toggle revealing an "Include UK
# Biobank samples" switch and a "Genetic ancestry group" of ancestry toggles,
# each with Both / Female / Male sub-options. Option/sub-option ids match the
# ConfigIniParams parameter names, so selections round-trip into the ini.
_GNOMAD_EXOMES_ANCESTRIES = [
    ("all", "All"),
    ("afr", "African & African-American"),
    ("amr", "Admixed American"),
    ("asj", "Ashkenazi Jewish"),
    ("eas", "East Asian"),
    ("fin", "Finnish"),
    ("mid", "Middle Eastern"),
    ("nfe", "Non-Finnish European"),
]


def _gnomad_sex_suboptions(option_id: str) -> list[dict]:
    """Both / Female / Male toggles for one ancestry option (Both on = combined
    sexes). `option_id` is the ancestry option's id, e.g. `gnomad_exomes_afr`."""
    return [
        {"id": f"{option_id}_both", "label": "Combined", "type": "boolean", "default": True},
        {"id": f"{option_id}_female", "label": "XX", "type": "boolean", "default": False},
        {"id": f"{option_id}_male", "label": "XY", "type": "boolean", "default": False},
    ]


def _gnomad_exomes_option() -> dict:
    """The gnomAD Exomes v4.1.1 option (freshly built so callers can mutate it)."""
    ancestry_options = [
        {
            "id": f"gnomad_exomes_{anc}",
            "label": label,
            "type": "boolean",
            "default": anc == "all",  # "All" pre-selected -> fields=AF baseline
            "sub_options": _gnomad_sex_suboptions(f"gnomad_exomes_{anc}"),
        }
        for anc, label in _GNOMAD_EXOMES_ANCESTRIES
    ]
    return {
        "id": "gnomad_exomes",
        "label": "gnomAD Exomes v4.1.1",
        "type": "boolean",
        "default": False,
        "sub_options": [
            {"id": "gnomad_exomes_include_ukb", "label": "Include UK BioBank samples", "type": "boolean", "default": True},
            {"type": "group", "label": "Genetic ancestry group", "options": ancestry_options},
        ],
    }


# gnomAD genomes v4.1 (human GRCh38): as exomes but no UK Biobank toggle, plus
# Amish / Remaining, and "Maximum across all groups" (grpmax) which has no sex
# split (a plain toggle).
_GNOMAD_GENOMES_ANCESTRIES = [
    ("all", "All"),
    ("afr", "African & African-American"),
    ("amr", "Admixed American"),
    ("asj", "Ashkenazi Jewish"),
    ("eas", "East Asian"),
    ("fin", "Finnish"),
    ("mid", "Middle Eastern"),
    ("nfe", "Non-Finnish European"),
    ("ami", "Amish"),
    ("remaining", "Remaining"),
]


def _gnomad_genomes_option() -> dict:
    """The gnomAD Genomes v4.1.1 option (freshly built so callers can mutate it)."""
    ancestry_options = [
        {
            "id": f"gnomad_genomes_{anc}",
            "label": label,
            "type": "boolean",
            "default": anc == "all",  # "All" pre-selected -> fields=AF baseline
            "sub_options": _gnomad_sex_suboptions(f"gnomad_genomes_{anc}"),
        }
        for anc, label in _GNOMAD_GENOMES_ANCESTRIES
    ]
    # grpmax (max across groups) has no XX/XY split -> a plain toggle.
    ancestry_options.append({
        "id": "gnomad_genomes_grpmax",
        "label": "Maximum across all groups",
        "type": "boolean",
        "default": False,
    })
    return {
        "id": "gnomad_genomes",
        "label": "gnomAD Genomes v4.1.1",
        "type": "boolean",
        "default": False,
        "sub_options": [
            {"type": "group", "label": "Genetic ancestry group", "options": ancestry_options},
        ],
    }


# NIH All of Us (human GRCh38): a flat list of population toggles (no sex split).
# "Maximum subpopulation" contributes two fields (gvs_max_af + gvs_max_subpop);
# that is handled by the ini builder, not the form.
_ALLOFUS_POPULATIONS = [
    ("all", "All"),
    ("max", "Maximum subpopulation"),
    ("afr", "African"),
    ("amr", "Latino/Ad Mixed American"),
    ("eas", "East Asian"),
    ("eur", "European"),
    ("mid", "Middle Eastern"),
    ("sas", "South Asian"),
    ("oth", "Other"),
]


def _allofus_option() -> dict:
    """The NIH All of Us option (freshly built so callers can mutate it)."""
    population_options = [
        {
            "id": f"allofus_{pop}",
            "label": label,
            "type": "boolean",
            "default": pop == "all",  # "All" pre-selected -> fields=gvs_all_af
        }
        for pop, label in _ALLOFUS_POPULATIONS
    ]
    return {
        "id": "allofus",
        "label": "NIH All of Us",
        "type": "boolean",
        "default": False,
        # A label-less group keeps the population list full-width (reusing the
        # nested-group renderer) without adding a heading.
        "sub_options": [
            {"type": "group", "options": population_options},
        ],
    }


# gnomAD SV v4.1 (human GRCh38): a flat list of AF toggles (no sex split). The
# SV id (`gnomAD_SV`) and `gnomAD_SV_SVTYPE` are always returned; these gate the
# per-population AF columns. Population code -> label; "" is the overall AF.
_GNOMAD_SV_POPULATIONS = [
    ("", "All"),
    ("afr", "African & African-American"),
    ("ami", "Amish"),
    ("amr", "Admixed American"),
    ("asj", "Ashkenazi Jewish"),
    ("eas", "East Asian"),
    ("fin", "Finnish"),
    ("mid", "Middle Eastern"),
    ("nfe", "Non-Finnish European"),
    ("rmi", "Remaining"),
    ("sas", "South Asian"),
]


def _gnomad_sv_af_option_id(code: str) -> str:
    """Form option id for a gnomAD SV AF population (`""` = overall)."""
    return "gnomad_sv_af" if code == "" else f"gnomad_sv_af_{code}"


def _gnomad_sv_option() -> dict:
    """The gnomAD SV v4.1 option: an overlap-cutoff select plus per-population AF
    toggles (overall AF pre-selected). SVTYPE + the SV id ride along always."""
    population_options = [
        {
            "id": _gnomad_sv_af_option_id(code),
            "label": label,
            "type": "boolean",
            "default": code == "",  # overall AF pre-selected
        }
        for code, label in _GNOMAD_SV_POPULATIONS
    ]
    return {
        "id": "gnomad_sv",
        "label": "gnomAD SV v4.1",
        "type": "boolean",
        "default": False,
        "sub_options": [
            {
                "id": "gnomad_sv_overlap_cutoff",
                "label": "Overlap cutoff",
                "type": "select",
                "default": "100",
                "options": [
                    {"label": "80%", "value": "80"},
                    {"label": "90%", "value": "90"},
                    {"label": "100%", "value": "100"},
                ],
            },
            {"type": "group", "options": population_options},
        ],
    }


# gnomAD CNV v4.1 (human GRCh38): like gnomAD SV, but *sample* frequencies (SF)
# and a slightly different population set (no Amish; "remaining" spelled out).
_GNOMAD_CNV_POPULATIONS = [
    ("", "All"),
    ("afr", "African & African-American"),
    ("amr", "Admixed American"),
    ("asj", "Ashkenazi Jewish"),
    ("eas", "East Asian"),
    ("fin", "Finnish"),
    ("mid", "Middle Eastern"),
    ("nfe", "Non-Finnish European"),
    ("sas", "South Asian"),
    ("remaining", "Remaining"),
]


def _gnomad_cnv_sf_option_id(code: str) -> str:
    """Form option id for a gnomAD CNV SF population (`""` = overall)."""
    return "gnomad_cnv_sf" if code == "" else f"gnomad_cnv_sf_{code}"


def _gnomad_cnv_option() -> dict:
    """The gnomAD CNV v4.1 option: an overlap-cutoff select plus per-population
    SF toggles (overall SF pre-selected). SVTYPE + the CNV id ride along always."""
    population_options = [
        {
            "id": _gnomad_cnv_sf_option_id(code),
            "label": label,
            "type": "boolean",
            "default": code == "",  # overall SF pre-selected
        }
        for code, label in _GNOMAD_CNV_POPULATIONS
    ]
    return {
        "id": "gnomad_cnv",
        "label": "gnomAD CNV v4.1",
        "type": "boolean",
        "default": False,
        "sub_options": [
            {
                "id": "gnomad_cnv_overlap_cutoff",
                "label": "Overlap cutoff",
                "type": "select",
                "default": "100",
                "options": [
                    {"label": "80%", "value": "80"},
                    {"label": "90%", "value": "90"},
                    {"label": "100%", "value": "100"},
                ],
            },
            {"type": "group", "options": population_options},
        ],
    }


# --------------------------------------------------------------------------- #
# AF population-code -> form-label decoders
#
# The results parser emits the same population codes these option ids are built
# from (see results_filters.af_source_descriptor), so a served AF column can be
# labelled by decoding its code back to the form label. That decode lives here,
# beside the option tuples it reuses, so the label vocabulary is defined exactly
# once — the frontend reads the decoded label off the response rather than
# keeping its own copy of these tables.
# --------------------------------------------------------------------------- #

# gnomAD ancestry codes -> labels. Genomes is a superset of exomes, so its list
# covers both sources; grpmax (genomes-only, added separately as a plain toggle)
# is folded in. The form's "all" ancestry is the overall AF, which the parser
# reports as an empty population code, so it is left out of the map.
_GNOMAD_ANCESTRY_LABELS = {
    code: label for code, label in _GNOMAD_GENOMES_ANCESTRIES if code != "all"
} | {"grpmax": "Maximum across all groups"}

# All of Us population codes -> labels ("all" is the overall AF -> empty code).
_ALLOFUS_POPULATION_LABELS = {
    code: label for code, label in _ALLOFUS_POPULATIONS if code != "all"
}

# gnomAD SV population codes -> labels ("" is the overall AF -> "All").
_GNOMAD_SV_POPULATION_LABELS = {
    code: label for code, label in _GNOMAD_SV_POPULATIONS if code != ""
}

# gnomAD CNV population codes -> labels ("" is the overall SF -> "All").
_GNOMAD_CNV_POPULATION_LABELS = {
    code: label for code, label in _GNOMAD_CNV_POPULATIONS if code != ""
}

# gnomAD sex-split suffixes — shown as the chromosomal notation (matching the
# form's XX/XY sub-option labels).
_SEX_LABELS = {"XX": "XX", "XY": "XY"}


def _gnomad_population_label(code: str) -> str:
    """Decode a compound gnomAD population code, e.g. `nfe_XX` -> 'Non-Finnish
    European · Female', `non_ukb_afr` -> 'African & African-American · excl. UK
    Biobank'. A code is an optional `non_ukb` subset prefix, then an ancestry (or
    a bare sex code = all ancestries), then an optional `_XX`/`_XY` sex suffix."""
    rest = code

    exclude_ukb = False
    if rest == "non_ukb" or rest.startswith("non_ukb_"):
        exclude_ukb = True
        rest = rest[len("non_ukb"):].lstrip("_")

    sex = None
    if rest.endswith(("_XX", "_XY")):
        sex = _SEX_LABELS[rest[-2:]]
        rest = rest[:-3]
    elif rest in ("XX", "XY"):
        sex = _SEX_LABELS[rest]
        rest = ""

    ancestry = "All" if rest == "" else _GNOMAD_ANCESTRY_LABELS.get(rest, rest)
    parts = [ancestry]
    if sex:
        parts.append(sex)
    if exclude_ukb:
        parts.append("excl. UK Biobank")
    return " · ".join(parts)


def af_population_label(source: str, code: str) -> str:
    """The form label for an AF population `code` within an AF `source`
    (`gnomad_exomes` / `gnomad_genomes` / `all_of_us`). The empty code is the
    source's overall AF ("All"); an unrecognised code falls back to itself.

    Reused by the results metadata (each `AfSource.label`) and the All of Us
    `max_subpopulation` decode, so the label vocabulary stays defined once."""
    if code == "":
        return "All"
    if source == "all_of_us":
        return _ALLOFUS_POPULATION_LABELS.get(code, code)
    if source == "gnomad_sv":
        return _GNOMAD_SV_POPULATION_LABELS.get(code, code)
    if source == "gnomad_cnv":
        return _GNOMAD_CNV_POPULATION_LABELS.get(code, code)
    return _gnomad_population_label(code)


def af_max_subpopulation_label(raw: str) -> str:
    """Decode All of Us `max_subpopulation` — the subpopulation(s) the max AF came
    from, given as `&`-joined population codes — to ` / `-joined labels."""
    return " / ".join(
        af_population_label("all_of_us", part) for part in raw.split("&")
    )


def _add_human_grch38_options(panels: list[dict]) -> None:
    """Layer the human GRCh38-only options onto the (already human 37/38) panels.

    Mutates `panels` in place. Assumes the GRCh37/38 panels are already present.
    """
    by_id = {panel["id"]: panel for panel in panels}

    # Genes & transcripts: RiboSeqORFs + Gene Ontology.
    by_id["genes_and_transcripts"]["options"].extend([
        {"id": "riboseqorfs", "label": "RiboSeqORFs", "type": "boolean", "default": False},
        # GO plugin (human GRCh38 for now; other species to follow).
        {"id": "go", "label": "Gene Ontology", "type": "boolean", "default": False},
    ])

    # Protein & functional: Protein (protein + ProtVar) / Functional (MaveDB,
    # IntAct, mutfunc).
    protein_panel = by_id["protein_and_functional"]
    for option in protein_panel["options"]:
        if option["id"] == "protein":
            option["category"] = "Protein"
    protein_panel["options"].extend([
        {"id": "protvar", "label": "ProtVar", "type": "boolean", "default": False,
         "category": "Protein", "sub_options": copy.deepcopy(_PROTVAR_SUBOPTIONS)},
        {"id": "mavedb", "label": "MaveDB", "type": "boolean", "default": False, "category": "Functional"},
        {"id": "intact", "label": "IntAct", "type": "boolean", "default": False,
         "category": "Functional", "sub_options": copy.deepcopy(_INTACT_SUBOPTIONS)},
        {"id": "mutfunc", "label": "mutfunc", "type": "boolean", "default": False,
         "category": "Functional", "sub_options": copy.deepcopy(_MUTFUNC_SUBOPTIONS)},
    ])

    # Variant Impact Predictions panel: EVE (Missense).
    if "pathogenicity_predictions" in by_id:
        by_id["pathogenicity_predictions"]["options"].append(
            {"id": "eve", "label": "EVE", "type": "boolean", "default": False, "category": "Missense"}
        )

    # Variant associations: OpenTargets + Phenotypes.
    if "variant_associations" in by_id:
        by_id["variant_associations"]["options"].extend([
            {"id": "opentargets", "label": "OpenTargets", "type": "boolean", "default": False},
            # Phenotypes plugin (human GRCh38 for now; other species to follow).
            {"id": "phenotypes", "label": "Phenotypes", "type": "boolean", "default": False},
        ])
        # ClinVar's "Structural variants" sub-option (the ClinVar_SV custom) —
        # GRCh38-only, so it joins the master's sub-options only here.
        clinvar = next(
            (o for o in by_id["variant_associations"]["options"] if o["id"] == "clinvar"),
            None,
        )
        if clinvar is not None:
            clinvar.setdefault("sub_options", []).append(
                {"id": "clinvar_sv", "label": "Structural variants", "type": "boolean", "default": False}
            )

    # Regulatory: GENCODE promoter windows (a gff-overlap custom annotation).
    panels.append({
        "id": "regulatory",
        "label": "Regulatory",
        "options": [
            {"id": "gencode_promoters", "label": "GENCODE promoter", "type": "boolean", "default": False},
        ],
    })

    # Allele frequencies: gnomAD exomes/genomes v4.1, NIH All of Us.
    panels.append({
        "id": "allele_frequencies",
        "label": "Allele frequencies",
        "options": [
            _gnomad_exomes_option(),
            _gnomad_genomes_option(),
            _allofus_option(),
            _gnomad_sv_option(),
            _gnomad_cnv_option(),
        ],
    })


def get_visible_panels(
    attributes: dict | None = None,
    *,
    species_taxonomy_id: str | None = None,
    assembly_name: str | None = None,
) -> list[dict]:
    """Return the panels/options to show for a genome.

    `attributes` is the genome metadata (genebuild.* etc.). `species_taxonomy_id`
    and `assembly_name` are passed by the client (from the selected species) so
    visibility can depend on species/assembly — e.g. human GRCh37/38.
    """
    panels = copy.deepcopy(_ALWAYS_VISIBLE_PANELS)

    if is_human_grch37_or_38(species_taxonomy_id, assembly_name):
        # UTRAnnotator extends the existing Genes & transcripts panel.
        for panel in panels:
            if panel["id"] == "genes_and_transcripts":
                panel["options"].extend(copy.deepcopy(_HUMAN_37_38_GENES_OPTIONS))
        # Variant-impact / conservation / associations panels are human-only.
        panels.extend(copy.deepcopy(_HUMAN_37_38_PANELS))

    if is_human_grch38(species_taxonomy_id, assembly_name):
        _add_human_grch38_options(panels)

    return panels
