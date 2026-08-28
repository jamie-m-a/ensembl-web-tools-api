"""URLs a plugin builds from the variant rather than from a CSQ column.

A resource that keys its pages on the variant has no column for the parse to
read. `variant_links` names a template per field, resolved against the variant
the annotation belongs to, so what reaches the frontend is an ordinary parsed
field that a display row links to with `link_from`.
"""

from app.vep.utils.vcf_results import _resolve_variant_links

TOKENS = {
    "chromosome": "17",
    "position": "42872062",
    "reference": "G",
    "alternative": "A",
}


def test_a_template_is_filled_from_the_variant():
    assert _resolve_variant_links(
        {"link": "https://x.org/variant/{chromosome}-{position}-{reference}-{alternative}"},
        TOKENS,
    ) == {"link": "https://x.org/variant/17-42872062-G-A"}


def test_a_template_naming_an_unknown_token_yields_nothing():
    """Rather than a URL with a literal `{gene}` in it."""
    assert _resolve_variant_links({"link": "https://x.org/{gene}"}, TOKENS) == {}


def test_no_variant_means_no_link():
    """The row then renders as plain text, as it does for a plugin that emitted
    no URL — never a link to a variant that cannot be named."""
    assert _resolve_variant_links({"link": "https://x.org/{chromosome}"}, None) == {}


def test_one_bad_template_does_not_lose_the_others():
    resolved = _resolve_variant_links(
        {"good": "https://x.org/{chromosome}", "bad": "https://x.org/{nope}"},
        TOKENS,
    )
    assert resolved == {"good": "https://x.org/17"}
