"""Which option panels/options are visible on the VEP input form for a genome.

Returned by the form_config endpoint (called on species selection). For now the
set is the same for every species ("always visible"); species-conditional rules
will be layered on later by inspecting the genome metadata attributes.

Option (and sub-option) `id`s match the ConfigIniParams parameter names, so the
form's selections round-trip back into the generated config.ini.
"""

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


def is_human_grch37_or_38(
    species_taxonomy_id: str | None, assembly_name: str | None
) -> bool:
    """True for human GRCh37 / GRCh38."""
    return species_taxonomy_id == "9606" and (assembly_name or "").startswith(
        ("GRCh37", "GRCh38")
    )


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
    panels = list(_ALWAYS_VISIBLE_PANELS)

    # Species-conditional panels/options are layered on here, e.g.:
    # if is_human_grch37_or_38(species_taxonomy_id, assembly_name):
    #     panels += _HUMAN_37_38_PANELS
    return panels
