"""The CSQ access vocabulary: reading values out of a VEP CSQ entry.

A VEP output VCF packs every annotation into one pipe-delimited CSQ string per
allele x feature, whose column layout is declared once in the
``##INFO=<ID=CSQ ...>`` header. These helpers turn that layout into an index map
and read values through it. They know about the CSQ *format* and nothing about
what any particular annotation means.

Shared by both parsing paths — the hand-written parsers in vcf_results and the
spec-driven interpreter — which is why these are public: they are the vocabulary
those modules are written in, not internals of either.
"""

import re


def get_prediction_index_map(csq_header: str) -> dict[str, int]:
    """Creates a dictionary of column indexes from the CSQ info description.

    Every CSQ column is indexed, so any annotation field can be read."""
    csq_header = csq_header.split(":")[-1].strip()
    csq_headers = csq_header.split("|")

    return {header: index for index, header in enumerate(csq_headers)}


def csq_index_map_from_header(header_lines: list[str]) -> dict[str, int]:
    """CSQ column -> index, parsed from the raw ##INFO=<ID=CSQ ...> header line.
    Used by the filter scan, which reads raw text rather than via vcfpy."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            match = re.search(r'Description="([^"]*)"', line)
            if match:
                return get_prediction_index_map(match.group(1))
    return {}


def get_csq_value(
    csq_values: list[str], csq_key: str, default_value: str | None, index_map: dict[str, int]
):
    """Helper method to return CSQ values or a default value
    if either the key or the value is missing"""
    if csq_key in index_map and csq_values[index_map[csq_key]]:
        return csq_values[index_map[csq_key]]
    return default_value


def has_any_column(index_map: dict[str, int], *columns: str) -> bool:
    """Whether any of `columns` is present in the CSQ header at all. Lets a parser
    skip its work when the plugin that produces its columns wasn't run — the
    header is fixed for the whole file, so a column that is absent here is absent
    for every record (and the parser could only ever return None/empty)."""
    return any(column in index_map for column in columns)


def to_float(value: str | None) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def split_amp(value: str | None) -> list[str]:
    """Split a '&'-delimited CSQ list, dropping empties and 'NA' placeholders."""
    if not value:
        return []
    return [v for v in value.split("&") if v and v != "NA"]


def raw_amp(value: str | None) -> list[str]:
    """Split a '&'-delimited CSQ list keeping every position (incl. 'NA'), so
    positionally-aligned subfields can be zipped together."""
    return value.split("&") if value else []


def first_amp(value: str | None) -> str | None:
    """First real (non-empty, non-'NA') item of a '&'-joined CSQ list."""
    for item in split_amp(value):
        return item
    return None
