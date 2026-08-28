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


# --- version resolution ------------------------------------------------------
#
# One authored entry serves an option offered at a different version per
# assembly. Both halves resolve against the option's own label, and what is
# served is finished: no placeholder, and only the links that belong with it.


def test_the_placeholder_is_resolved_from_the_label():
    description = _options()["gnomad_exomes"]["help"]["description"]
    assert "{version}" not in description
    assert "(gnomAD) v4.1.1." in description


def test_the_other_assembly_gets_its_own_version_from_the_same_entry():
    description = _options("GRCh37.p13")["gnomad_exomes"]["help"]["description"]
    assert "(gnomAD) v2.1.1." in description
    assert "4.1.1" not in description


def test_no_option_serves_a_placeholder():
    for assembly in ("GRCh38.p14", "GRCh37.p13"):
        for option_id, option in _options(assembly).items():
            description = option.get("help", {}).get("description", "")
            assert "{version}" not in description, f"{assembly} {option_id}"


def test_a_label_with_no_version_leaves_no_gap():
    """The placeholder takes the space before it, so the sentence still reads."""
    help_ = OptionHelp(option_id="x", description="Frequencies from gnomAD{version}.")
    assert help_.as_payload("gnomAD Genomes")["description"] == (
        "Frequencies from gnomAD."
    )


def test_each_assembly_is_served_only_its_own_release():
    """gnomAD SV is v4.1 on GRCh38 and v2.1 on GRCh37, and the v4 release
    announcement does not describe the v2 callset."""
    grch38 = _options()["gnomad_sv"]["help"]["links"]
    grch37 = _options("GRCh37.p13")["gnomad_sv"]["help"]["links"]
    assert len(grch38) == 1 and "v4-structural-variants" in grch38[0]["href"]
    assert grch37 == [{"href": "https://europepmc.org/article/MED/32461652"}]


def test_a_point_release_keeps_its_link():
    """Matching on the major alone: v4.2.1 still describes the v4 callset."""
    help_ = OptionHelp(
        option_id="x",
        description="d",
        links=[OptionHelpLink(href="https://example.org/v4", major_version="4")],
    )
    assert help_.as_payload("gnomAD SV v4.2.1")["links"] == [
        {"href": "https://example.org/v4"}
    ]


def test_a_version_specific_link_is_dropped_when_the_label_has_no_version():
    """Citing the wrong release is worse than citing none."""
    help_ = OptionHelp(
        option_id="x",
        description="d",
        links=[OptionHelpLink(href="https://example.org/v4", major_version="4")],
    )
    assert "links" not in help_.as_payload("gnomAD SV")


def test_a_link_without_a_version_is_always_served():
    help_ = OptionHelp(
        option_id="x",
        description="d",
        links=[OptionHelpLink(href="https://gnomad.broadinstitute.org/")],
    )
    for label in ("gnomAD Exomes v2.1.1", "gnomAD Exomes"):
        assert help_.as_payload(label)["links"] == [
            {"href": "https://gnomad.broadinstitute.org/"}
        ]


def test_major_version_is_spent_not_sent():
    """It chooses the links; the frontend has no use for it afterwards."""
    for option in _options().values():
        for link in option.get("help", {}).get("links", []):
            assert "majorVersion" not in link
            assert "major_version" not in link


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
    assert help_.as_payload("Some option v4.1") == {
        "description": "y",
        "links": [{"href": "https://example.org"}],
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


# --- form-only help -----------------------------------------------------------


def test_form_only_help_is_served_to_the_form():
    """It is ordinary help there — the flag only decides where it is shown."""
    assert _options()["avi"]["help"]["description"].startswith("AlphaGenome")


def test_the_flag_itself_is_never_served():
    """It says where the help goes, not what it says."""
    assert "form_only" not in json.dumps(_options())


def test_the_results_panels_drop_it():
    """An option whose rows carry their own help would otherwise show two (?)
    buttons on the row the option's help anchors to."""
    from app.vep.models.display_panels_model import (
        dump_display_panels,
        to_display_panels,
    )
    from app.vep.utils.vcf_results import _drop_form_only_help
    from app.vep.form_panels import get_visible_panels

    merged = load_merged_spec("human_grch38")
    raw = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    panels = _drop_form_only_help(to_display_panels(raw), merged)
    served = {o["id"]: o for p in dump_display_panels(panels) for o in p["options"]}
    assert "help" not in served["avi"]
    # and only that option
    assert "help" in served["cadd"]


def test_an_unflagged_option_is_untouched_by_the_pass():
    from app.vep.models.display_panels_model import (
        dump_display_panels,
        to_display_panels,
    )
    from app.vep.utils.vcf_results import _drop_form_only_help
    from app.vep.form_panels import get_visible_panels

    merged = load_merged_spec("human_grch38")
    raw = get_visible_panels(species_taxonomy_id=HUMAN, assembly_name="GRCh38.p14")
    before = {o["id"]: dict(o) for p in raw for o in p["options"] if o["id"] != "avi"}
    panels = _drop_form_only_help(to_display_panels(raw), merged)
    after = {o["id"]: o for p in dump_display_panels(panels) for o in p["options"]}
    for option_id, option in before.items():
        assert after[option_id] == option, option_id
