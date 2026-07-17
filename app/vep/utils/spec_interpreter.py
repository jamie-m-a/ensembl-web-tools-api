"""Applies a parsing spec to a CSQ entry, in place of a hand-written `_parse_*`.

This is the generic half of the planned annotation-API work: the spec says what
to read and how to shape it, this module does it. Output is a plain dict (the
generic annotation payload), not a per-plugin pydantic model.

Currently additive — the hand-written parsers in vcf_results are still the ones
wired into the response. This runs alongside them so the two can be compared
over the same CSQ fixtures (see tests/test_spec_interpreter.py).
"""

import re

from vep.models.parsing_spec_model import PluginSpec, TargetSpec, WhenSpec
from vep.utils.csq import (
    first_amp,
    get_csq_value,
    has_any_column,
    raw_amp,
    split_amp,
    to_float,
)

# Some plugins write a literal 'NA' for "no value here".
_NULLISH = ("", "NA")

_PLACEHOLDER_RE = re.compile(r"\{[^}]*\}")


def _is_present(value) -> bool:
    """Whether a built output actually carries something.

    Deliberately not plain truthiness: an allele frequency of 0.0 is a real
    value, and `not 0.0` would throw it away.
    """
    if value is None:
        return False
    if isinstance(value, (str, list, dict, tuple)):
        return len(value) > 0
    return True


def _coerce(raw: str | None, value_type: str, field_spec=None):
    """A raw CSQ value as `value_type`, or None if absent/'NA'/unparseable."""
    if raw is None or raw in _NULLISH:
        return None
    if value_type == "float":
        return to_float(raw)
    if value_type == "int":
        try:
            return int(raw)
        except (TypeError, ValueError):
            return None
    if field_spec is not None:
        for find, replacement in (field_spec.replace or {}).items():
            raw = raw.replace(find, replacement)
        if field_spec.strip:
            raw = raw.strip()
    return raw


def _column(csq_values: list[str], name: str, index_map: dict[str, int]) -> str | None:
    return get_csq_value(csq_values, name, None, index_map)


def _should_drop(row: dict, drop_when) -> bool:
    if drop_when is None:
        return False
    if drop_when.all_null:
        return all(value is None for value in row.values())
    return row.get(drop_when.null) is None


def _apply_post(rows: list[dict], post) -> list[dict]:
    """Whole-list operations, in the order the spec lists them."""
    for operation in post or []:
        if operation.op == "dedup":
            seen = set()
            unique = []
            for row in rows:
                key = tuple(row.values())
                if key in seen:
                    continue
                seen.add(key)
                unique.append(row)
            rows = unique
        elif operation.op == "sort":
            nulls_last = operation.nulls == "last"
            # A null key sorts to whichever end `nulls` asks for, whatever `desc`
            # does to the rest.
            sentinel = float("-inf") if nulls_last == operation.desc else float("inf")
            rows = sorted(
                rows,
                key=lambda row: (
                    row[operation.by] if row[operation.by] is not None else sentinel
                ),
                reverse=operation.desc,
            )
    return rows


def _apply_zip(csq_values, index_map, target: TargetSpec) -> list[dict]:
    """N positionally-aligned '&'-lists -> a list of objects.

    Uses the position-preserving split: an 'NA' still occupies a slot, which is
    what keeps the columns aligned with each other.
    """
    columns = [raw_amp(_column(csq_values, name, index_map)) for name in target.source]
    lengths = [len(column) for column in columns]
    length = (max(lengths) if target.align == "max" else min(lengths)) if lengths else 0

    rows: list[dict] = []
    for i in range(length):
        row = {
            field_spec.field: _coerce(
                column[i] if i < len(column) else None, field_spec.type, field_spec
            )
            for column, field_spec in zip(columns, target.as_fields)
        }
        if _should_drop(row, target.drop_when):
            continue
        rows.append(row)
    return _apply_post(rows, target.post)


def _apply_regex(csq_values, index_map, target: TargetSpec):
    """Named regex groups -> object(s). Non-matching items are skipped."""
    raw = _column(csq_values, target.source, index_map)
    compiled = re.compile(target.pattern)
    items = split_amp(raw) if target.each else ([raw] if raw else [])

    rows = []
    for item in items:
        match = compiled.match(item)
        if not match:
            continue
        rows.append(
            {
                field_spec.field: _coerce(
                    match.group(field_spec.field), field_spec.type, field_spec
                )
                for field_spec in target.as_fields
            }
        )
    if target.each:
        return rows
    return rows[0] if rows else None


def _apply_pattern_map(csq_values, index_map, target: TargetSpec) -> dict:
    """Columns matching `from_pattern` -> {wildcard: value}.

    The columns are discovered from the CSQ header, so whichever ancestries a
    run actually emitted come through without being named in the spec.
    """
    placeholder = _PLACEHOLDER_RE.search(target.from_pattern)
    prefix = target.from_pattern[: placeholder.start()]
    suffix = target.from_pattern[placeholder.end() :]
    excluded = set(target.exclude or [])

    values: dict = {}
    for column in index_map:
        if column in excluded:
            continue
        if not (column.startswith(prefix) and column.endswith(suffix)):
            continue
        key = column[len(prefix) : len(column) - len(suffix)]
        value = _coerce(_column(csq_values, column, index_map), target.type)
        if value is not None:
            values[key] = value
    return values


def _build_object(tokens: list[str], field_specs, source_text: str) -> dict:
    """Assign `tokens` to `field_specs` strictly by index.

    A `raw`-typed field takes `source_text` and consumes no slot. Slots past the
    end of `tokens` are null: a missing item leaves *its own* field empty and
    does not shift its neighbours along.
    """
    built: dict = {}
    position = 0
    for field_spec in field_specs:
        if field_spec.type == "raw":
            built[field_spec.field] = source_text
            continue
        token = tokens[position] if position < len(tokens) else None
        built[field_spec.field] = _coerce(token, field_spec.type, field_spec)
        position += 1
    return built


def _apply_chunk(csq_values, index_map, target: TargetSpec) -> list[dict]:
    """Fixed-size groups of '&'-items -> a list of objects."""
    raw = _column(csq_values, target.source, index_map)
    tokens = raw.split("&") if raw else []

    rows = []
    for start in range(0, len(tokens), target.size):
        group = tokens[start : start + target.size]
        row = _build_object(group, target.as_fields, "&".join(t for t in group if t))
        if _should_drop(row, target.drop_when):
            continue
        rows.append(row)
    return _apply_post(rows, target.post)


def _apply_positional(csq_values, index_map, target: TargetSpec):
    """'&'-items -> one object, assigned by index."""
    raw = _column(csq_values, target.source, index_map)
    if not raw:
        return [] if target.wrap == "list" else None
    built = _build_object(raw.split("&"), target.as_fields, raw)
    return [built] if target.wrap == "list" else built


def _apply_key_value(csq_values, index_map, target: TargetSpec) -> dict:
    """A ':'-delimited 'k=v' string -> {k: v}.

    Order-independent by construction, unlike a plain scalar copy of the same
    string. A piece without `kv_delimiter` is dropped rather than raising: a
    malformed piece should not break parsing of an otherwise-good value.
    """
    raw = _column(csq_values, target.source, index_map)
    if not raw:
        return {}
    values: dict = {}
    for piece in raw.split(target.pair_delimiter):
        if target.kv_delimiter not in piece:
            continue
        key, value = piece.split(target.kv_delimiter, 1)
        values[key] = value
    return values


def _when_holds(csq_values, index_map, when: WhenSpec | None) -> bool:
    if when is None:
        return True
    return when.includes in split_amp(_column(csq_values, when.field, index_map))


def _empty_value(target: TargetSpec):
    """What a target yields when its `when` condition does not hold."""
    if target.transform in ("list", "zip", "chunk"):
        return []
    if target.transform == "regex":
        return [] if target.each else None
    if target.transform in ("pattern_map", "key_value"):
        return {}
    if target.transform == "positional":
        return [] if target.wrap == "list" else None
    return None


def _apply_target(csq_values, index_map, target: TargetSpec):
    if not _when_holds(csq_values, index_map, target.when):
        return _empty_value(target)

    if target.transform == "zip":
        return _apply_zip(csq_values, index_map, target)
    if target.transform == "regex":
        return _apply_regex(csq_values, index_map, target)
    if target.transform == "pattern_map":
        return _apply_pattern_map(csq_values, index_map, target)
    if target.transform == "chunk":
        return _apply_chunk(csq_values, index_map, target)
    if target.transform == "positional":
        return _apply_positional(csq_values, index_map, target)
    if target.transform == "key_value":
        return _apply_key_value(csq_values, index_map, target)

    raw = _column(csq_values, target.source, index_map)
    if target.transform == "scalar":
        return _coerce(raw, target.type)
    if target.transform == "list":
        return split_amp(raw)
    if target.transform == "first":
        return _coerce(first_amp(raw), target.type)
    raise ValueError(f"unknown transform: {target.transform}")


def apply_plugin_spec(
    csq_values: list[str], index_map: dict[str, int], spec: PluginSpec
) -> dict | None:
    """One plugin's annotation for this CSQ entry, or None if there is nothing.

    None means "no annotation", matching the hand-written parsers: either the
    plugin's columns are absent from the header (it never ran), or they are
    present but this record has no values in them.
    """
    if not has_any_column(index_map, *spec.csq_fields):
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
        _is_present(output.get(field)) for field in spec.require_any_output
    ):
        return None

    return output
