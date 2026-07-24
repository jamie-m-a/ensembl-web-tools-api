"""Server-side filtering of VEP result records.

Filters are applied during a sequential scan of the results VCF. The page-index
fast path (see vcf_results.get_results_from_path) can't be used once records are
filtered, because a filtered record's position no longer maps to a fixed page
offset — so a filtered request scans and applies predicates.

Filters operate at the CSQ-entry (per transcript/consequence) level: a record's
matching entries are those satisfying every condition, and the record is kept if
at least one entry matches. Crucially, the *non-matching entries are pruned* from
the kept record, so a filtered variant only carries the transcripts that actually
match (e.g. filtering by a consequence hides the variant's other-consequence
transcripts). Allele-level annotations are identical across an allele's CSQ rows,
so they survive as long as one entry for that allele remains.

Filters run as an ordered *pipeline*: for each record the entry set is narrowed
by each filter in turn, short-circuiting the moment nothing survives. We tally
how many records each filter removed (among those that reached it), so the
ordering can later be tuned to run the cheapest / highest-yield filters first.
We don't rank them yet — this just captures the numbers to inform that later.

This module is deliberately self-contained: it owns the request model, the
predicate compilation and the pipeline, so the whole feature can be removed in
one piece.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Iterable, Iterator

from pydantic import BaseModel

from vep.form_panels import af_population_label
from vep.utils.spec_interpreter import pattern_affixes

if TYPE_CHECKING:
    from vep.models.parsing_spec_model import ParsingSpec

# Field identifiers understood by the query builder. Must match the frontend's
# filter `field` values (see the results filters UI).
CONSEQUENCE_FIELD = "consequence"
TRANSCRIPT_FIELD = "transcript"
GENE_SYMBOL_FIELD = "gene_symbol"
GENE_ID_FIELD = "gene_id"
TRANSCRIPT_GROUP_FIELD = "transcript_group"
ALLELE_FREQUENCY_FIELD = "allele_frequency"

# Operators understood by the query builder.
OPERATOR_IN = "in"  # "is any of"
# Numeric comparisons (allele frequency): <=, ==, >=.
OPERATOR_LE = "le"
OPERATOR_EQ = "eq"
OPERATOR_GE = "ge"

# Allele-frequency match modes.
AF_MATCH_ANY = "any"  # keep if any tested AF meets the comparison
AF_MATCH_ALL = "all"  # keep only if every tested AF meets the comparison


def _strip_version(feature_id: str) -> str:
    """The stable-id portion of a versioned feature id (drop the '.version')."""
    return feature_id.split(".", 1)[0]


class ResultsFilter(BaseModel):
    """One condition from the results query builder. Conditions are AND-combined
    (a CSQ entry must satisfy every condition to be kept)."""

    field: str
    operator: str
    values: list[str] = []
    # Allele-frequency filter only: the comparison threshold (0-1) and whether to
    # match any/all of the tested AF columns.
    threshold: float | None = None
    match: str | None = None


class FilterError(ValueError):
    """Raised for a malformed filters payload or an unsupported field/operator."""


@dataclass
class CompiledFilter:
    """A ready-to-run predicate for one condition. `keep_entry` takes a single
    CSQ entry (already split on '|') and returns whether that entry matches.

    `line_prefilter`, when set, is a cheap *necessary condition* on the raw,
    unsplit data line: if it returns False the record cannot possibly match, so
    it is dropped without splitting the (large, many-column) CSQ. It must never
    reject a record that would match — false positives are fine (they fall
    through to the exact `keep_entry` check), false negatives are not. Only
    literal-substring membership filters can supply one (see `_membership_prefilter`)."""

    field: str
    keep_entry: Callable[[list[str]], bool]
    line_prefilter: Callable[[str], bool] | None = None


@dataclass
class FilterStat:
    """How many records a filter removed, among those that reached it."""

    field: str
    removed: int


@dataclass
class FilterOutcome:
    """The result of one streaming filter pass over a record source.

    `page` holds only the survivors in the requested `[start, start+count)` slice
    (so memory is bounded by the page, not the file); `matched_total` and
    `scanned_total` are the full counts needed for pagination and metadata."""

    page: list[str]
    matched_total: int
    scanned_total: int
    stats: list[FilterStat]


def parse_filters(raw: str | None) -> list[ResultsFilter]:
    """Parse the `filters` query param (a JSON array of conditions) into models.

    Returns [] when the param is absent or empty. Raises FilterError on malformed
    JSON or shape, so the route can turn it into a 4xx rather than a 500."""
    if not raw:
        return []
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise FilterError(f"filters is not valid JSON: {exc}") from exc
    if not isinstance(data, list):
        raise FilterError("filters must be a JSON array")
    try:
        return [ResultsFilter.model_validate(item) for item in data]
    except Exception as exc:  # pydantic ValidationError et al.
        raise FilterError(f"invalid filter condition: {exc}") from exc


def _split_line(line: str) -> tuple[list[str], bool]:
    """Split a VCF data line into columns, remembering whether it had a newline."""
    has_newline = line.endswith("\n")
    return (line[:-1] if has_newline else line).split("\t"), has_newline


def _find_csq(columns: list[str]) -> tuple[int, int, list[list[str]]] | None:
    """Locate the CSQ annotation in a record's columns.

    Returns (info_part_index, info_field_count_marker, entries) where entries is a
    list of CSQ entries each split into its '|' subfields; or None if the record
    has no CSQ (or too few columns to have an INFO field)."""
    if len(columns) < 8:
        return None
    info_parts = columns[7].split(";")
    for part_index, part in enumerate(info_parts):
        if part.startswith("CSQ="):
            entries = [entry.split("|") for entry in part[4:].split(",")]
            return part_index, len(info_parts), entries
    return None


def extract_csq_entries(line: str) -> list[list[str]]:
    """The CSQ entries of a raw data line, each split into '|' subfields; empty
    when the line has no CSQ. (Convenience for tests / callers that only read.)"""
    columns, _ = _split_line(line)
    found = _find_csq(columns)
    return found[2] if found else []


def _rebuild_line(
    columns: list[str],
    csq_part_index: int,
    kept_entries: list[list[str]],
    has_newline: bool,
) -> str:
    """Rebuild a data line with its CSQ narrowed to `kept_entries`."""
    info_parts = columns[7].split(";")
    info_parts[csq_part_index] = "CSQ=" + ",".join(
        "|".join(entry) for entry in kept_entries
    )
    columns = list(columns)
    columns[7] = ";".join(info_parts)
    return "\t".join(columns) + ("\n" if has_newline else "")


def _require_operator(f: "ResultsFilter", *allowed: str) -> None:
    if f.operator not in allowed:
        raise FilterError(
            f"unsupported operator {f.operator!r} for field {f.field!r}"
        )


def _membership_prefilter(tokens: set[str]) -> Callable[[str], bool] | None:
    """A raw-line necessary-condition test for a literal-membership filter: the
    line must contain at least one selected token as a substring, since a matching
    CSQ entry carries that token verbatim in the line. A cheap C-level `in` run
    before the record is split — it rejects the (usually vast) majority that can't
    match without touching the many-column CSQ. None for an empty token set."""
    if not tokens:
        return None
    values = tuple(tokens)
    return lambda line: any(token in line for token in values)


def _compile_consequence(f: ResultsFilter, index_map: dict[str, int]) -> CompiledFilter | None:
    """A CSQ entry matches if its Consequence carries one of the selected terms.
    VEP '&'-joins co-occurring terms within an entry."""
    _require_operator(f, OPERATOR_IN)
    if not f.values:
        return None
    consequence_index = index_map.get("Consequence")
    if consequence_index is None:
        raise FilterError("Consequence column missing from CSQ header")
    selected = set(f.values)

    def keep_entry(entry: list[str]) -> bool:
        if consequence_index >= len(entry):
            return False
        return bool(selected.intersection(entry[consequence_index].split("&")))

    return CompiledFilter(
        field=CONSEQUENCE_FIELD,
        keep_entry=keep_entry,
        line_prefilter=_membership_prefilter(selected),
    )


def _compile_transcript(f: ResultsFilter, index_map: dict[str, int]) -> CompiledFilter | None:
    """A CSQ entry matches if its Feature (transcript) stable id is one of the
    selected ids. Match is version-insensitive: the '.version' suffix is ignored
    on both sides, so 'ENST0000012345' and 'ENST0000012345.7' are equivalent."""
    _require_operator(f, OPERATOR_IN)
    selected = {_strip_version(value) for value in f.values if value}
    if not selected:
        return None
    feature_index = index_map.get("Feature")
    if feature_index is None:
        raise FilterError("Feature column missing from CSQ header")

    def keep_entry(entry: list[str]) -> bool:
        if feature_index >= len(entry):
            return False
        feature = entry[feature_index]
        if not feature:
            return False
        return _strip_version(feature) in selected

    # The stripped id is a substring of the versioned id as it appears in the line
    # (e.g. "ENST00000012345" within "ENST00000012345.7"), so it is a valid
    # necessary-condition token.
    return CompiledFilter(
        field=TRANSCRIPT_FIELD,
        keep_entry=keep_entry,
        line_prefilter=_membership_prefilter(selected),
    )


def _compile_gene_symbol(f: ResultsFilter, index_map: dict[str, int]) -> CompiledFilter | None:
    """A CSQ entry matches if its SYMBOL (gene name) is one of the selected names.
    Match is case-insensitive (human symbols are conventionally upper-case)."""
    _require_operator(f, OPERATOR_IN)
    selected = {value.upper() for value in f.values if value}
    if not selected:
        return None
    symbol_index = index_map.get("SYMBOL")
    if symbol_index is None:
        raise FilterError("SYMBOL column missing from CSQ header")

    def keep_entry(entry: list[str]) -> bool:
        if symbol_index >= len(entry):
            return False
        symbol = entry[symbol_index]
        return bool(symbol) and symbol.upper() in selected

    return CompiledFilter(field=GENE_SYMBOL_FIELD, keep_entry=keep_entry)


def _compile_gene_id(f: ResultsFilter, index_map: dict[str, int]) -> CompiledFilter | None:
    """A CSQ entry matches if its Gene (Ensembl gene stable id) is one of the
    selected ids. Version-insensitive, like the transcript filter."""
    _require_operator(f, OPERATOR_IN)
    selected = {_strip_version(value) for value in f.values if value}
    if not selected:
        return None
    gene_index = index_map.get("Gene")
    if gene_index is None:
        raise FilterError("Gene column missing from CSQ header")

    def keep_entry(entry: list[str]) -> bool:
        if gene_index >= len(entry):
            return False
        gene = entry[gene_index]
        return bool(gene) and _strip_version(gene) in selected

    return CompiledFilter(
        field=GENE_ID_FIELD,
        keep_entry=keep_entry,
        line_prefilter=_membership_prefilter(selected),
    )


def _entry_value(entry: list[str], name: str, index_map: dict[str, int]) -> str | None:
    """A CSQ subfield value for an entry, or None if the column is absent/empty."""
    i = index_map.get(name)
    if i is None or i >= len(entry) or not entry[i]:
        return None
    return entry[i]


# Transcript-group id -> a per-entry test. Mirrors the canonical / MANE detection
# in vcf_results._get_alt_allele_details, so the filter agrees with what the
# results view labels each transcript.
_TRANSCRIPT_GROUP_TESTS: dict[str, Callable[[list[str], dict[str, int]], bool]] = {
    "canonical": lambda entry, index_map: _entry_value(entry, "CANONICAL", index_map)
    == "YES",
    "mane_select": lambda entry, index_map: bool(
        _entry_value(entry, "MANE_SELECT", index_map)
    )
    or _entry_value(entry, "MANE", index_map) == "MANE_Select",
    "mane_plus_clinical": lambda entry, index_map: bool(
        _entry_value(entry, "MANE_PLUS_CLINICAL", index_map)
    )
    or _entry_value(entry, "MANE", index_map) == "MANE_Plus_Clinical",
    "gencode_primary": lambda entry, index_map: _entry_value(
        entry, "GENCODE_PRIMARY", index_map
    )
    == "1",
}


def _compile_transcript_group(f: ResultsFilter, index_map: dict[str, int]) -> CompiledFilter | None:
    """A CSQ entry matches if it belongs to any of the selected transcript groups
    (canonical / MANE Select / MANE Plus Clinical / GENCODE primary). Which groups
    are offered is a frontend, species-dependent concern; unknown group ids are
    rejected here."""
    _require_operator(f, OPERATOR_IN)
    if not f.values:
        return None
    tests = []
    for group in f.values:
        test = _TRANSCRIPT_GROUP_TESTS.get(group)
        if test is None:
            raise FilterError(f"unsupported transcript group: {group!r}")
        tests.append(test)

    def keep_entry(entry: list[str]) -> bool:
        return any(test(entry, index_map) for test in tests)

    return CompiledFilter(field=TRANSCRIPT_GROUP_FIELD, keep_entry=keep_entry)


def af_columns(index_map: dict[str, int], spec: ParsingSpec | None = None) -> list[str]:
    """The AF-bearing CSQ columns present, in header order. These are exactly the
    allele-frequency options selected at input (VEP only emits chosen ones).

    With a parsing `spec`, a column is AF iff some `frequencies.*` plugin claims
    it (see `af_source_descriptor`) — so the gnomAD v2 grammar (a subset prefix
    *before* `AF`), or any future source, is discovered without this list knowing
    the shape. Without one (the filter path, which has no spec) it falls back to a
    prefix match: any `gnomAD_(exomes|genomes)_…AF…` column, the SV/CNV `AF`/`SF`
    columns, and All of Us `…_af` (the `AoU_gvs_max_subpop` label column excluded).
    """
    if spec is not None:
        columns = [name for name in index_map if af_source_descriptor(name, spec)]
    else:
        columns = [
            name
            for name in index_map
            if (
                name.startswith(("gnomAD_exomes_", "gnomAD_genomes_"))
                and "AF" in name.split("_")
            )
            or name.startswith(("gnomAD_SV_AF", "gnomAD_CNV_SF"))
            or (name.startswith("AoU_gvs_") and name.endswith("_af"))
        ]
    return sorted(columns, key=lambda name: index_map[name])


def _af_source_specs(spec: ParsingSpec) -> list[tuple]:
    """Per allele-frequency plugin (output under `frequencies.`):
    `(source, overall_column, prefix, suffix, exclude)`. `overall_column` is the
    `field=="overall"` scalar's column; `prefix`/`suffix` bracket the populations
    `pattern_map`'s placeholder, so a matched column's key is exactly what sits
    between them — the same key the parse stores the value under."""
    specs: list[tuple] = []
    for plugin in spec.plugins:
        if not plugin.output.startswith("frequencies."):
            continue
        source = plugin.output.split(".")[-1]
        overall = next(
            (
                target.source
                for target in plugin.targets
                if target.transform == "scalar"
                and target.field == "overall"
                and isinstance(target.source, str)
            ),
            None,
        )
        pattern = next(
            (
                target
                for target in plugin.targets
                if target.transform == "pattern_map" and target.from_pattern
            ),
            None,
        )
        if pattern is None:
            prefix = suffix = None
            exclude: set[str] = set()
        else:
            prefix, suffix = pattern_affixes(pattern.from_pattern)
            exclude = set(pattern.exclude or [])
        specs.append((source, overall, prefix, suffix, exclude))
    return specs


def _af_descriptor(column: str, source: str, population: str) -> dict:
    return {
        "key": column,
        "source": source,
        "population": population,
        "label": af_population_label(source, population),
    }


def af_source_descriptor(column: str, spec: ParsingSpec | None = None) -> dict | None:
    """Split an AF column into {key, source, population, label} for the results
    metadata (population "" = the source's overall AF, labelled "All"). The
    `population` code is exactly the key the parse stores the value under, so the
    frontend can join a source's populations against each allele's parsed
    frequencies; the `label` is the human population name (decoded once here, from
    form_panels). None for a non-AF column.

    With a `spec`, source + population are derived from the `frequencies.*`
    plugins — handling gnomAD v2's subset-before-`AF` grammar (population
    `controls_AF_afr`, not `afr`) and any future source. Without one (older jobs)
    it falls back to the fixed gnomAD v4 / SV / CNV / All of Us layout."""
    if spec is not None:
        for source, overall, prefix, suffix, exclude in _af_source_specs(spec):
            if overall is not None and column == overall:
                return _af_descriptor(column, source, "")
            if (
                prefix is not None
                and column not in exclude
                and column.startswith(prefix)
                and column.endswith(suffix)
                and len(column) > len(prefix) + len(suffix)
            ):
                population = column[len(prefix): len(column) - len(suffix)]
                return _af_descriptor(column, source, population)
        return None

    # Legacy (spec-less) path: the fixed gnomAD v4 / SV / CNV / All of Us layout.
    if column.startswith("gnomAD_exomes_AF"):
        source = "gnomad_exomes"
        population = column[len("gnomAD_exomes_AF"):].lstrip("_")
    elif column.startswith("gnomAD_genomes_AF"):
        source = "gnomad_genomes"
        population = column[len("gnomAD_genomes_AF"):].lstrip("_")
    elif column.startswith("gnomAD_SV_AF"):
        source = "gnomad_sv"
        population = column[len("gnomAD_SV_AF"):].lstrip("_")
    elif column.startswith("gnomAD_CNV_SF"):
        source = "gnomad_cnv"
        population = column[len("gnomAD_CNV_SF"):].lstrip("_")
    elif column.startswith("AoU_gvs_") and column.endswith("_af"):
        source = "all_of_us"
        raw = column[len("AoU_gvs_"):-len("_af")]
        population = "" if raw == "all" else raw
    else:
        return None
    return _af_descriptor(column, source, population)


def _compile_allele_frequency(f: ResultsFilter, index_map: dict[str, int]) -> CompiledFilter | None:
    """Keep an allele whose AF meets the comparison. AF is allele-level (identical
    across an allele's CSQ rows), so this reads like an entry test but effectively
    keeps/drops whole alleles. Tests either the specified AF columns (`values`) or,
    when none are given, all AF columns present. `match` = any/all across them.

    NO-DATA: missing/empty AF values are ignored (not tested). If an allele has no
    AF data at all for the tested columns it is currently dropped. The no-data
    semantics need revisiting (see results-filtering-notes.md)."""
    _require_operator(f, OPERATOR_LE, OPERATOR_EQ, OPERATOR_GE)
    if f.threshold is None:
        return None  # nothing to compare against -> no-op
    threshold = f.threshold

    requested = [c for c in f.values if c]
    columns = requested or af_columns(index_map)
    indices = [index_map[c] for c in columns if c in index_map]
    if not indices:
        return None  # no AF columns to test (e.g. AF not run) -> no-op

    if f.operator == OPERATOR_LE:
        compare = lambda value: value <= threshold
    elif f.operator == OPERATOR_GE:
        compare = lambda value: value >= threshold
    else:
        compare = lambda value: value == threshold
    match_all = f.match == AF_MATCH_ALL

    def keep_entry(entry: list[str]) -> bool:
        values: list[float] = []
        for i in indices:
            if i < len(entry) and entry[i] not in ("", ".", None):
                try:
                    values.append(float(entry[i]))
                except ValueError:
                    pass  # non-numeric -> treat as no-data
        if not values:  # all no-data -> excluded for now (revisit)
            return False
        return all(map(compare, values)) if match_all else any(map(compare, values))

    return CompiledFilter(field=ALLELE_FREQUENCY_FIELD, keep_entry=keep_entry)


# Field id -> builder that compiles a ResultsFilter into a CompiledFilter (or None
# for a no-op). Adding a filter is a matter of adding an entry here plus its
# builder.
_BUILDERS: dict[
    str, Callable[[ResultsFilter, dict[str, int]], CompiledFilter | None]
] = {
    CONSEQUENCE_FIELD: _compile_consequence,
    TRANSCRIPT_FIELD: _compile_transcript,
    GENE_SYMBOL_FIELD: _compile_gene_symbol,
    GENE_ID_FIELD: _compile_gene_id,
    TRANSCRIPT_GROUP_FIELD: _compile_transcript_group,
    ALLELE_FREQUENCY_FIELD: _compile_allele_frequency,
}


def compile_filters(
    filters: list[ResultsFilter], index_map: dict[str, int]
) -> list[CompiledFilter]:
    """Validate and compile request filters against the file's CSQ layout. A
    builder returns None for a no-op condition (e.g. empty values), which is
    skipped."""
    compiled: list[CompiledFilter] = []
    for f in filters:
        builder = _BUILDERS.get(f.field)
        if builder is None:
            raise FilterError(f"unsupported filter field: {f.field!r}")
        cf = builder(f, index_map)
        if cf is not None:
            compiled.append(cf)
    return compiled


# What a surviving record carries forward so its line can be rebuilt lazily (only
# the page slice is ever rebuilt): the split columns, the CSQ part index, the
# surviving entries, and whether the line had a trailing newline.
_Survivor = tuple[list[str], int, list[list[str]], bool]


def _evaluate_record(
    line: str, compiled: list[CompiledFilter], removed: list[int]
) -> _Survivor | None:
    """Run the ordered pipeline over one record. Returns what's needed to rebuild
    the kept line (its CSQ narrowed to surviving entries), or None if a filter
    dropped it — crediting the removal to that filter in `removed`.

    Rebuilding is left to the caller so only the requested page slice pays for it.
    Each filter's cheap `line_prefilter` runs first, on the raw unsplit line, so a
    record that can't match is rejected before the expensive CSQ split."""
    # Necessary-condition rejection on the raw line, before any splitting.
    for i, cf in enumerate(compiled):
        if cf.line_prefilter is not None and not cf.line_prefilter(line):
            removed[i] += 1
            return None

    columns, has_newline = _split_line(line)
    found = _find_csq(columns)
    if found is None:
        # No CSQ to match against — a consequence-style filter can't keep it.
        if compiled:
            removed[0] += 1
        return None
    csq_part_index, _, entries = found

    for i, cf in enumerate(compiled):
        survivors = [entry for entry in entries if cf.keep_entry(entry)]
        if not survivors:
            removed[i] += 1
            return None
        entries = survivors

    return columns, csq_part_index, entries, has_newline


def filter_records(
    data_lines: Iterable[str],
    compiled: list[CompiledFilter],
    *,
    start: int = 0,
    count: int | None = None,
) -> FilterOutcome:
    """Stream the ordered filter pipeline over raw VCF data lines in a single pass.

    For each record the CSQ entry set is narrowed by each filter in turn; a filter
    that leaves no surviving entry drops the record (and is credited with the
    removal). Kept records are rebuilt with only their surviving entries.

    Only the survivors in the half-open `[start, start+count)` slice are retained
    (rebuilt) — so a page can be served from a huge result set without holding
    every match in memory. `count=None` keeps every survivor (small callers /
    tests). The full `matched_total` / `scanned_total` counts are always tallied,
    since pagination needs the total regardless of which page is asked for.

    `data_lines` may be a lazy iterator (e.g. a gzip line stream), so the whole
    file need never be materialised."""
    removed = [0] * len(compiled)
    scanned = 0
    matched = 0
    page: list[str] = []
    stop = None if count is None else start + count
    for line in data_lines:
        scanned += 1
        survivor = _evaluate_record(line, compiled, removed)
        if survivor is None:
            continue
        # Rebuild only the survivors that fall in the requested window; the rest
        # are counted but never reassembled.
        if matched >= start and (stop is None or matched < stop):
            page.append(_rebuild_line(*survivor))
        matched += 1

    stats = [
        FilterStat(field=cf.field, removed=removed[i])
        for i, cf in enumerate(compiled)
    ]
    return FilterOutcome(
        page=page, matched_total=matched, scanned_total=scanned, stats=stats
    )


def apply_filter_pipeline(
    data_lines: Iterable[str], compiled: list[CompiledFilter]
) -> tuple[list[str], list[FilterStat]]:
    """Run the ordered filter pipeline over raw VCF data lines, returning every
    kept (rebuilt) line in order plus per-filter removal counts.

    A thin wrapper over `filter_records` that keeps all survivors — retained for
    callers and tests that want the full kept list. The paginated results path
    uses `filter_records` directly with a `start`/`count` window so it never holds
    the whole match set in memory."""
    outcome = filter_records(data_lines, compiled)
    return outcome.page, outcome.stats


def stream_filtered_lines(
    data_lines: Iterable[str], compiled: list[CompiledFilter]
) -> Iterator[str]:
    """Lazily yield every kept (rebuilt) VCF data line in order — its CSQ narrowed
    to the entries surviving all filters — dropping records with no survivor.

    The streaming counterpart to `apply_filter_pipeline`: it accumulates nothing
    and tallies no stats, yielding each survivor as it is found, so a filtered
    download stays bounded in memory no matter how many records match."""
    discard = [0] * len(compiled)
    for line in data_lines:
        survivor = _evaluate_record(line, compiled, discard)
        if survivor is not None:
            yield _rebuild_line(*survivor)
