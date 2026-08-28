"""The help behind an option's (?) button, served from the spec's `help` section.

An option given no help simply shows no tooltip, and nothing downstream errors —
so these tests assert that the text arrived rather than trusting that it did.
"""

import json

import pytest

from app.vep.form_panels import get_visible_panels
from app.vep.models.merged_spec_model import MergedSpec
from app.vep.models.option_help_model import HelpSpec, OptionHelp, OptionHelpLink
from app.vep.utils.spec_loader import load_merged_spec

HUMAN = "9606"


def _options(assembly: str = "GRCh38.p14", taxon: str = HUMAN) -> dict[str, dict]:
    panels = get_visible_panels(species_taxonomy_id=taxon, assembly_name=assembly)
    return {o["id"]: o for panel in panels for o in panel["options"]}


def test_the_library_carries_help():
    spec = load_merged_spec("human_grch38")
    assert spec.help is not None
    assert len(spec.help.options) == 35


def test_a_served_option_carries_its_help():
    hgvs = _options()["hgvs"]
    assert hgvs["help"]["description"].startswith("HGVS")
    assert hgvs["help"]["links"] == [{"href": "https://hgvs-nomenclature.org/stable/"}]


def test_an_option_with_no_help_says_nothing():
    """A missing key, not an empty one: the frontend falls back on `undefined`,
    and an empty description would render an empty tooltip instead."""
    assert "help" not in _options()["updownstream_distance"]


def test_help_follows_the_option_to_every_genome_that_offers_it():
    """CADD is declared by both human documents; the sentence is written once."""
    assert _options("GRCh38.p14")["cadd"]["help"] == _options("GRCh37.p13")["cadd"]["help"]


def test_a_version_placeholder_is_served_uninterpolated():
    """`{version}` is resolved from the option's rendered label, which is the
    frontend's business — the backend must not guess at it."""
    assert "{version}" in _options()["gnomad_exomes"]["help"]["description"]


def test_a_version_scoped_link_keeps_the_frontend_s_spelling():
    """`majorVersion`, not `major_version`: the payload is handed straight to the
    frontend, whose OptionHelpLink names the field that way. Snake_case would
    drop the filter silently — every link would show, and one assembly's help
    would cite the other's release announcement."""
    links = _options()["gnomad_sv"]["help"]["links"]
    assert [link.get("majorVersion") for link in links] == ["4", "2"]
    assert not any("major_version" in link for link in links)


def test_a_link_without_a_label_omits_the_key():
    """The frontend supplies its own default label. A null would override it."""
    for option in _options().values():
        for link in option.get("help", {}).get("links", []):
            assert "label" not in link or link["label"]


def test_a_spec_pinned_before_help_existed_still_loads():
    """`help` is optional so a spec without it still loads: `_load_pinned_spec`
    swallows a validation error and returns None, which would render those
    results with no annotations at all."""
    doc = json.loads(
        load_merged_spec("human_grch38").model_dump_json(by_alias=True, exclude_none=True)
    )
    doc.pop("help")
    assert MergedSpec.model_validate(doc).help is None


def test_an_unknown_field_in_a_help_entry_is_an_error():
    """The models are `extra="forbid"` throughout; help is authored by hand and
    a typo should fail loudly rather than silently serve less."""
    with pytest.raises(Exception, match="[Ee]xtra"):
        HelpSpec.model_validate(
            {"options": [{"option_id": "x", "description": "y", "linkz": []}]}
        )


def test_the_payload_drops_the_key_it_is_keyed_by():
    help_ = OptionHelp(
        option_id="x",
        description="y",
        links=[OptionHelpLink(href="https://example.org", major_version="4")],
    )
    assert help_.as_payload() == {
        "description": "y",
        "links": [{"href": "https://example.org", "majorVersion": "4"}],
    }


# --- wording ------------------------------------------------------------------
#
# Not spell-checks: each is a distinction someone had to think about, and the
# kind of thing a well-meaning tidy-up would flatten.


def test_every_allele_frequency_option_notes_the_source_naming():
    """Population codes are the source's own, not harmonised across sources, and
    the help says so on all five — the short-variant sources and both
    structural-variant ones."""
    opts = _options()
    for option_id in (
        "gnomad_exomes",
        "gnomad_genomes",
        "allofus",
        "gnomad_sv",
        "gnomad_cnv",
    ):
        assert "Populations are named as at source." in opts[option_id]["help"][
            "description"
        ], option_id


def test_gnomad_cnv_is_described_as_sample_frequencies():
    """gnomAD reports CNVs as the fraction of samples carrying the call, not as
    an allele count. The wording is deliberate, not a slip."""
    description = _options()["gnomad_cnv"]["help"]["description"]
    assert "Sample frequencies for copy number variants" in description
    assert "Allele frequencies" not in description


def test_gnomad_sv_is_described_as_allele_frequencies():
    """Its sibling genuinely is an allele frequency, which is what makes the
    CNV wording above worth pinning."""
    assert "Allele frequencies for structural variants" in (
        _options()["gnomad_sv"]["help"]["description"]
    )
