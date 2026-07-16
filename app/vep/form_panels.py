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
            # `hgvs` param, on by default) and a separate HGVSg (the `hgvsg`
            # param, off by default). The frontend builds the linked UI; the
            # panel just carries the `hgvs` option (default on).
            {"id": "hgvs", "label": "HGVS", "type": "boolean", "default": True},
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
                "sub_options": [
                    {
                        "id": "tss_distance_direction",
                        "type": "select",
                        "default": "upstream",
                        "options": [
                            {"label": "Upstream", "value": "upstream"},
                            {"label": "Downstream", "value": "downstream"},
                            {"label": "Both", "value": "both"},
                        ],
                    }
                ],
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
]

# Extra panels shown only for human GRCh37/38. Pathogenicity options carry a
# `category` label used to group them within the panel.
_HUMAN_37_38_PANELS: list[dict] = [
    {
        "id": "pathogenicity_predictions",
        "label": "Pathogenicity predictions",
        "options": [
            {"id": "alphamissense", "label": "AlphaMissense", "type": "boolean", "default": False, "category": "Missense"},
            {"id": "revel", "label": "Revel", "type": "boolean", "default": False, "category": "Missense"},
            {"id": "spliceai", "label": "SpliceAI", "type": "boolean", "default": False, "category": "Splicing"},
            {"id": "maxentscan", "label": "MaxEntScan", "type": "boolean", "default": False, "category": "Splicing"},
            {"id": "cadd", "label": "CADD", "type": "boolean", "default": False, "category": "Genome wide"},
            {"id": "enformer", "label": "Enformer", "type": "boolean", "default": False, "category": "Non-coding"},
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


def _add_human_grch38_options(panels: list[dict]) -> None:
    """Layer the human GRCh38-only options onto the (already human 37/38) panels.

    Mutates `panels` in place. Assumes the GRCh37/38 panels are already present.
    """
    by_id = {panel["id"]: panel for panel in panels}

    # Genes & transcripts: RiboSeqORFs.
    by_id["genes_and_transcripts"]["options"].append(
        {"id": "riboseqorfs", "label": "RiboSeqORFs", "type": "boolean", "default": False}
    )

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

    # Pathogenicity predictions: EVE (Missense).
    if "pathogenicity_predictions" in by_id:
        by_id["pathogenicity_predictions"]["options"].append(
            {"id": "eve", "label": "EVE", "type": "boolean", "default": False, "category": "Missense"}
        )

    # Variant associations: OpenTargets.
    if "variant_associations" in by_id:
        by_id["variant_associations"]["options"].append(
            {"id": "opentargets", "label": "OpenTargets", "type": "boolean", "default": False}
        )

    # Allele frequencies: gnomAD mitochondrial.
    panels.append({
        "id": "allele_frequencies",
        "label": "Allele frequencies",
        "options": [
            {"id": "gnomad_mt", "label": "gnomAD mitochondrial", "type": "boolean", "default": False},
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
        # Pathogenicity / conservation / associations panels are human-only.
        panels.extend(copy.deepcopy(_HUMAN_37_38_PANELS))

    if is_human_grch38(species_taxonomy_id, assembly_name):
        _add_human_grch38_options(panels)

    return panels
