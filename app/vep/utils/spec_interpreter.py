"""Applies a parsing spec to a CSQ entry, in place of a hand-written `_parse_*`.

This is the generic half of the planned annotation-API work: the spec says what
to read and how to shape it, this module does it. Output is a plain dict (the
generic annotation payload), not a per-plugin pydantic model.

Currently additive — the hand-written parsers in vcf_results are still the ones
wired into the response. This runs alongside them so the two can be compared
over the same CSQ fixtures (see tests/test_spec_interpreter.py).

NB it imports the CSQ primitives from vcf_results. Those are the vocabulary the
spec-driven path shares with the hand-written one; when this path takes over,
they want extracting into their own `csq` module.
"""

from vep.models.parsing_spec_model import PluginSpec, TargetSpec
from vep.utils.vcf_results import (
    _first_amp,
    _get_csq_value,
    _has_any_column,
    _raw_amp,
    _split_amp,
    _to_float,
)

# Some plugins write a literal 'NA' for "no value here".
_NULLISH = ("", "NA")


def _coerce(raw: str | None, value_type: str):
    """A raw CSQ value as `value_type`, or None if absent/'NA'/unparseable."""
    if raw is None or raw in _NULLISH:
        return None
    if value_type == "float":
        return _to_float(raw)
    if value_type == "int":
        try:
            return int(raw)
        except (TypeError, ValueError):
            return None
    return raw


def _column(csq_values: list[str], name: str, index_map: dict[str, int]) -> str | None:
    return _get_csq_value(csq_values, name, None, index_map)


def _apply_zip(csq_values, index_map, target: TargetSpec) -> list[dict]:
    """N positionally-aligned '&'-lists -> a list of objects.

    Uses the position-preserving split: an 'NA' still occupies a slot, which is
    what keeps the columns aligned with each other.
    """
    columns = [_raw_amp(_column(csq_values, name, index_map)) for name in target.source]
    lengths = [len(column) for column in columns]
    length = (max(lengths) if target.align == "max" else min(lengths)) if lengths else 0

    rows: list[dict] = []
    for i in range(length):
        row = {
            field_spec.field: _coerce(
                column[i] if i < len(column) else None, field_spec.type
            )
            for column, field_spec in zip(columns, target.as_fields)
        }
        if target.drop_when == "all_null" and all(v is None for v in row.values()):
            continue
        rows.append(row)
    return rows


def _apply_target(csq_values, index_map, target: TargetSpec):
    if target.transform == "zip":
        return _apply_zip(csq_values, index_map, target)

    raw = _column(csq_values, target.source, index_map)
    if target.transform == "scalar":
        return _coerce(raw, target.type)
    if target.transform == "list":
        return _split_amp(raw)
    if target.transform == "first":
        return _coerce(_first_amp(raw), target.type)
    raise ValueError(f"unknown transform: {target.transform}")


def apply_plugin_spec(
    csq_values: list[str], index_map: dict[str, int], spec: PluginSpec
) -> dict | None:
    """One plugin's annotation for this CSQ entry, or None if there is nothing.

    None means "no annotation", matching the hand-written parsers: either the
    plugin's columns are absent from the header (it never ran), or they are
    present but this record has no values in them.
    """
    if not _has_any_column(index_map, *spec.csq_fields):
        return None

    # Raw presence, deliberately: a literal 'NA' counts as present here, which
    # is what the hand-written parsers do.
    if spec.require_any_input and not any(
        _column(csq_values, column, index_map) for column in spec.require_any_input
    ):
        return None

    output = {
        target.field: _apply_target(csq_values, index_map, target)
        for target in spec.targets
    }

    if spec.require_any_output and not any(
        output.get(field) for field in spec.require_any_output
    ):
        return None

    return output
