""" Module for loading a VCF and parsing it into a VepResultsResponse
object as defined in APISpecification"""

from collections import deque, OrderedDict
from dataclasses import dataclass
from io import StringIO
from typing import Iterable, Iterator
import gzip
import itertools
import json
import logging
import re
import subprocess
import threading
from pathlib import Path
from pydantic import FilePath
from vep.models import vcf_results_model as model
from vep.form_panels import af_max_subpopulation_label
from vep.utils import results_filters
from vep.utils.bgzf import _BgzfReader
from vep.utils.csq import (
    csq_index_map_from_header,
    get_csq_value,
    get_prediction_index_map,
    has_any_column,
    split_amp,
    to_float,
)
from vep.utils.vcf_meta import _get_vcf_meta
from vep.utils.vcf_reader import read_records
from vep.utils.spec_loader import (
    load_display_panels_sidecar,
    load_expected_columns_sidecar,
    load_spec_sidecar,
    resolve_merged_spec,
)
from vep.models.display_panels_model import DisplayPanel
from vep.models.display_spec_model import DisplayPayload
from vep.models.merged_spec_model import MergedSpec
from vep.utils.spec_interpreter import (
    apply_plugin_spec,
    compile_parsing_spec,
    pattern_affixes,
)
from vep.models.parsing_spec_model import ParsingSpec

# Taken from https://github.com/Ensembl/ensembl-hypsipyle
# main/common/file_model/variant.py#L142
# Needs to be moved into a shared module
def _set_allele_type(alt_one_bp: bool, ref_one_bp: bool, ref_alt_equal_bp: bool) -> tuple[str,str]:
    """Create a allele type for a variant based on Variation
    teams logic using ref and largest alt allele sizes"""
    match [alt_one_bp, ref_one_bp, ref_alt_equal_bp]:
        case [True, True, True]:
            allele_type = "SNV"
            so_term = "SO:0001483"

        case [True, False, False]:
            allele_type = "deletion"
            so_term = "SO:0000159"

        case [False, True, False]:
            allele_type = "insertion"
            so_term = "SO:0000667"

        case [False, False, False]:
            allele_type = "indel"
            so_term = "SO:1000032"

        case [False, False, True]:
            allele_type = "substitution"
            so_term = "SO:1000002"
    return allele_type, so_term

def _get_variant_type(ref: str, alt: str) -> str:
    """Helper function to infer variant type from allele values"""
    if alt=="copy_number_variation":
        return alt
    else:
        return _set_allele_type(len(alt) < 2, len(ref) < 2, len(alt) == len(ref))[0]


def _alt_value(alt) -> str:
    """Return an alt allele's sequence string.

    Simple substitution alts expose `.value`; symbolic and breakend alts
    (e.g. structural variants) do not, so fall back to their serialized VCF
    representation."""
    value = getattr(alt, "value", None)
    if value is not None:
        return value
    serialize = getattr(alt, "serialize", None)
    return serialize() if callable(serialize) else str(alt)


# VCF SVTYPE (or a symbolic allele's leading `<ID>` token) -> the word shown in
# the results "variant" column. Subtype-proof: <DEL:ME:ALU> is SVTYPE=DEL. See
# the VCF 4.2 spec, symbolic + breakend alternate alleles (§1.4.5, §5.3-5.4).
_SV_TYPE_WORDS = {
    "DEL": "deletion",
    "INS": "insertion",
    "DUP": "duplication",
    "INV": "inversion",
    "CNV": "copy_number_variation",
    "BND": "breakend",
}


def _first_info(value):
    """An INFO value as a scalar — vcfpy hands some back as single-element lists."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _alt_serialized(alt) -> str:
    """The alt allele's VCF text: `<DEL>` / `<DEL:ME:ALU>` for symbolic alleles,
    the `t[p[` notation for breakends, or the sequence for a substitution."""
    serialize = getattr(alt, "serialize", None)
    if callable(serialize):
        return serialize()
    return str(getattr(alt, "value", alt))


_BREAKEND_MATE_RE = re.compile(r"[\[\]]([^\[\]]+)[\[\]]")


def _breakend_junction(record, serialized: str) -> str:
    """A breakend's two loci as `<chrom:pos> ↔ <mate chrom:pos>` — this record's
    position plus the mate parsed from the ALT (e.g. `G[17:198982[` -> the mate is
    `17:198982`). The mate's base lives in the paired record (via MATEID), not this
    one, so only positions are shown."""
    start = f"{record.CHROM}:{record.POS}"
    mate = _BREAKEND_MATE_RE.search(serialized)
    return f"{start} ↔ {mate.group(1)}" if mate else start


def _sv_span(record) -> int | None:
    """A structural variant's length in bases: `abs(SVLEN)` when present and
    non-zero (deletion/insertion/duplication), else `END - POS` (inversion/CNV,
    whose SVLEN is 0 or absent). None when neither is available."""
    svlen = _first_info(record.INFO.get("SVLEN"))
    if svlen not in (None, "", "."):
        try:
            n = abs(int(str(svlen).split(",")[0]))
            if n:
                return n
        except (TypeError, ValueError):
            pass
    end = _first_info(record.INFO.get("END"))
    if end not in (None, "", "."):
        try:
            return abs(int(end) - int(record.POS))
        except (TypeError, ValueError):
            pass
    return None


def _structural_info(record) -> dict | None:
    """Structural-variant display details for a record, or None for a simple
    substitution/indel.

    VEP writes a type word (e.g. `deletion`) into the CSQ `Allele` column, so the
    length heuristic in `_get_variant_type` mis-reads it as an insertion and the
    UI renders the word as a sequence. Instead classify SVs from the VCF record:
    the type from `SVTYPE` (or the symbolic `<ID>`), the displayed allele from the
    symbolic ALT, and a `detail` line (span in bases for sized SVs). Breakends are
    shown as `<BND>` (like the other symbolic alleles) with their two loci as the
    detail, and the type `breakend`.
    """
    if not record.ALT:
        return None
    serialized = _alt_serialized(record.ALT[0])
    is_symbolic = serialized.startswith("<")
    is_breakend = ("[" in serialized) or ("]" in serialized)
    if not (is_symbolic or is_breakend):
        return None
    svtype = _first_info(record.INFO.get("SVTYPE"))
    if is_breakend:
        key = svtype or "BND"
        allele = "<BND>"
        detail = _breakend_junction(record, serialized)
    else:
        # Prefer the declared SVTYPE; fall back to the `<ID>` leading token.
        key = svtype or serialized[1:-1].split(":")[0]
        allele = serialized
        span = _sv_span(record)
        detail = f"{span} bp" if span else None
    return {
        "type_word": _SV_TYPE_WORDS.get(str(key).upper(), str(key).lower()),
        "allele": allele,
        "detail": detail,
    }






def _parse_uniprot(csq_values, index_map) -> model.UniprotIds | None:
    """Build Uniprot cross-references from the SWISSPROT/TREMBL/UNIPARC/isoform
    CSQ columns; returns None if none are present."""
    if not has_any_column(
        index_map, "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"
    ):
        return None
    swissprot = get_csq_value(csq_values, "SWISSPROT", None, index_map)
    trembl = get_csq_value(csq_values, "TREMBL", None, index_map)
    uniparc = get_csq_value(csq_values, "UNIPARC", None, index_map)
    isoform = get_csq_value(csq_values, "UNIPROT_ISOFORM", None, index_map)
    if not any([swissprot, trembl, uniparc, isoform]):
        return None
    return model.UniprotIds(
        swissprot=swissprot, trembl=trembl, uniparc=uniparc, isoform=isoform
    )


def _parse_protein_matches(csq_values, index_map) -> list[model.ProteinMatch]:
    """Parse the DOMAINS CSQ column (e.g. AlphaFold-DB / PDB mappings).
    Multiple matches are '&'-joined; each is 'source:id'."""
    domains = get_csq_value(csq_values, "DOMAINS", None, index_map)
    if not domains:
        return []
    matches = []
    for item in domains.split("&"):
        if not item:
            continue
        source, sep, identifier = item.partition(":")
        matches.append(
            model.ProteinMatch(
                source=source if sep else "",
                id=identifier if sep else source,
            )
        )
    return matches


_PREDICTION_RE = re.compile(r"^(?P<prediction>[^(]+)\((?P<score>[-\d.eE]+)\)$")


def _parse_prediction(value: str | None) -> model.PredictionWithScore | None:
    """Parse a 'prediction(score)' CSQ value, e.g. SIFT 'tolerated(0.15)'."""
    if not value:
        return None
    match = _PREDICTION_RE.match(value.strip())
    if match:
        return model.PredictionWithScore(
            prediction=match.group("prediction"),
            score=to_float(match.group("score")),
        )
    return model.PredictionWithScore(prediction=value, score=None)


def _resolve_variant_links(
    templates: dict[str, str], tokens: dict[str, str] | None
) -> dict[str, str]:
    """Each template filled from the variant, dropping any it cannot complete.

    A missing token means the URL would point at a variant we cannot name, so
    the field is omitted and the row renders as plain text — the same outcome as
    a plugin that emitted no URL.
    """
    if not tokens:
        return {}
    resolved = {}
    for field, template in templates.items():
        try:
            resolved[field] = template.format(**tokens)
        except KeyError:
            continue
    return resolved


def _spec_annotations(
    csq_values: list[str],
    index_map: dict[str, int],
    spec: ParsingSpec | None,
    scope: str,
    cache: dict | None = None,
    plans: dict | None = None,
    variant_tokens: dict[str, str] | None = None,
) -> list[model.Annotation]:
    """The generic annotations for one CSQ entry at the given scope, driving each
    matching spec plugin through `apply_plugin_spec`. Additive to the typed
    fields; when there is no pinned spec this is empty and nothing changes.

    `plans` is every plugin resolved against this file's CSQ header
    (`compile_parsing_spec`). Optional: without it each plugin resolves itself
    per row, which is correct but is the cost the plans exist to remove."""
    if spec is None:
        return []
    annotations: list[model.Annotation] = []
    for plugin in spec.plugins:
        if plugin.scope != scope:
            continue
        plan = plans.get(plugin.plugin) if plans else None
        # A plugin whose columns the header does not carry never ran, and that
        # is a property of the file — skip it here rather than per row.
        if plan is not None and not plan.runnable:
            continue
        data = apply_plugin_spec(csq_values, index_map, plugin, cache, plan)
        if data is not None and plugin.variant_links:
            data = {
                **data,
                **_resolve_variant_links(plugin.variant_links, variant_tokens),
            }
        if data is not None:
            annotations.append(
                model.Annotation(plugin=plugin.plugin, scope=scope, data=data)
            )
    return annotations


def _pool_annotations(variant: model.Variant) -> None:
    """Move this variant's annotations into one pool, referenced by index.

    VEP repeats a plugin's value on every CSQ row it applies to, so the same
    payload was serialised once per transcript consequence: on a 50-variant
    page ClinVar alone was 421 copies of 14 distinct values, 72% of the whole
    response. The payload therefore grew with annotations *times* transcripts,
    which is a multiplier rather than a constant — this removes it.

    Identity is the serialised payload, not the plugin: `hgvs` is genuinely
    different per transcript (744 distinct of 864), so deduplicating by plugin
    would lose data. Equal payloads collapse; different ones stay.

    Mutates in place, after the variant is built, so nothing upstream has to
    know about pooling.
    """
    pool: list[model.Annotation] = []
    seen: dict[str, int] = {}

    def refs(annotations: list[model.Annotation]) -> list[int]:
        out = []
        for annotation in annotations:
            key = annotation.model_dump_json()
            index = seen.get(key)
            if index is None:
                index = len(pool)
                seen[key] = index
                pool.append(annotation)
            out.append(index)
        return out

    for allele in variant.alternative_alleles:
        allele.annotation_refs = refs(allele.annotations)
        for consequence in allele.predicted_molecular_consequences:
            # An intergenic consequence is a different model and carries none.
            if hasattr(consequence, "annotations"):
                consequence.annotation_refs = refs(consequence.annotations)
    variant.annotation_pool = pool


def _get_alt_allele_details(
    ref: str,
    alt: str,
    csqs: list[str],
    index_map: dict[str, int],
    spec: ParsingSpec | None = None,
    sv: dict | None = None,
    plans: dict | None = None,
    chromosome: str | None = None,
    position: str | None = None,
) -> model.AlternativeVariantAllele:
    """Creates  AlternativeVariantAllele based on
    target alt allele and CSQ entires.

    `sv` (from `_structural_info`) overrides the type + displayed allele for
    structural variants: `alt` stays VEP's CSQ `Allele` word for matching the CSQ
    rows below, but the allele is shown as its symbolic/breakend form."""
    consequences = []
    # Resolving the header is meant to happen once per file; a caller that did
    # not do it still gets it once per allele rather than once per plugin per
    # CSQ row, which is what `apply_plugin_spec`'s own fallback would cost.
    if plans is None and spec is not None:
        plans = compile_parsing_spec(index_map, spec)
    allele_type = sv["type_word"] if sv else _get_variant_type(ref, alt)
    # Allele-level annotations are identical across all of this allele's CSQ
    # rows, so capture them once (from the first matching row). They are also
    # the only annotations available for intergenic variants (no transcript
    # rows).
    colocated_variants: list[str] = []
    allele_annotations: list[model.Annotation] = []
    allele_level_captured = False
    # A plugin reads only its own CSQ columns and VEP repeats those on every row
    # of the variant, so a transcript-scoped plugin produces the same annotation
    # for all of them. This lets it be parsed once and reused; the per-row gate
    # (`applies_to`) still runs for each, which is what makes ClinVar attach to
    # one gene and not its neighbour. Per allele, since that is the widest scope
    # over which the columns are guaranteed identical.
    parse_cache: dict = {}
    # This allele, for any plugin declaring `variant_links`. Built per allele
    # rather than per record: a multi-allelic variant would otherwise point all
    # of its alleles at the first one's page. `alt` is VEP's CSQ `Allele` word,
    # which is the VCF's own alternative for the substitutions these links are
    # for; an indel's trimmed form would need the record's ALT instead.
    variant_tokens = (
        {
            "chromosome": chromosome,
            "position": position,
            "reference": ref,
            "alternative": alt,
        }
        if chromosome is not None and position is not None
        else None
    )

    for str_csq in csqs:
        csq_values = str_csq.split("|")

        if csq_values[index_map["Allele"]] != alt:
            continue

        if not allele_level_captured:
            colocated_variants = split_amp(
                get_csq_value(csq_values, "Existing_variation", None, index_map)
            )
            allele_annotations = _spec_annotations(
                csq_values, index_map, spec, "allele", parse_cache, plans,
                variant_tokens,
            )
            allele_level_captured = True

        cons = get_csq_value(csq_values, "Consequence", "", index_map)
        if len(cons) == 0:
            cons = []
        else:
            cons = cons.split("&")
        if csq_values[index_map["Feature_type"]] == "Transcript":
            is_canonical = (
                get_csq_value(csq_values, "CANONICAL", "NO", index_map) == "YES"
            )

            # It looks like for Feature_type = Transcript that we always have a STRAND value
            strand = (
                model.Strand.reverse
                if get_csq_value(csq_values, "STRAND", "1", index_map) == "-1"
                else model.Strand.forward
            )

            # MANE: depending on the VEP run, either the MANE column carries the
            # label (MANE_Select / MANE_Plus_Clinical) or the MANE_SELECT /
            # MANE_PLUS_CLINICAL columns carry the matched RefSeq id. Handle both.
            mane_label = get_csq_value(csq_values, "MANE", None, index_map)
            mane_select_refseq = get_csq_value(
                csq_values, "MANE_SELECT", None, index_map
            )
            mane_plus_clinical = get_csq_value(
                csq_values, "MANE_PLUS_CLINICAL", None, index_map
            )
            is_mane_select = bool(mane_select_refseq) or mane_label == "MANE_Select"
            is_mane_plus_clinical = (
                bool(mane_plus_clinical) or mane_label == "MANE_Plus_Clinical"
            )

            # GENCODE primary: the GENCODE_PRIMARY column (from flag_gencode_primary,
            # human GRCh38 only) carries "1" for the primary transcript, else empty.
            is_gencode_primary = (
                get_csq_value(csq_values, "GENCODE_PRIMARY", None, index_map) == "1"
            )

            consequences.append(
                model.PredictedTranscriptConsequence(
                    feature_type=model.FeatureType.transcript,
                    stable_id=get_csq_value(csq_values, "Feature", "", index_map),
                    gene_stable_id=get_csq_value(csq_values, "Gene", "", index_map),
                    biotype=get_csq_value(csq_values, "BIOTYPE", "", index_map),
                    is_canonical=is_canonical,
                    gene_symbol=get_csq_value(csq_values, "SYMBOL", None, index_map),
                    consequences=cons,
                    strand=strand,
                    # MANE
                    is_mane_select=is_mane_select,
                    is_mane_plus_clinical=is_mane_plus_clinical,
                    mane_select_refseq_id=mane_select_refseq,
                    # GENCODE primary
                    is_gencode_primary=is_gencode_primary,
                    # Protein & functional annotations (ENSP is now the `protein`
                    # parse plugin, in `annotations`).
                    uniprot=_parse_uniprot(csq_values, index_map),
                    protein_matches=_parse_protein_matches(csq_values, index_map),
                    sift=_parse_prediction(
                        get_csq_value(csq_values, "SIFT", None, index_map)
                    ),
                    polyphen=_parse_prediction(
                        get_csq_value(csq_values, "PolyPhen", None, index_map)
                    ),
                    # Generic spec-driven annotations: everything else.
                    annotations=_spec_annotations(
                        csq_values, index_map, spec, "transcript", parse_cache, plans,
                        variant_tokens,
                    ),
                )
            )
        elif "intergenic_variant" in cons:
            consequences.append(
                model.PredictedIntergenicConsequence(
                    feature_type=None,
                    consequences=["intergenic_variant"],
                )
            )

    if sv:
        allele_sequence = sv["allele"]
    elif alt == "copy_number_variation":
        allele_sequence = ""
    else:
        allele_sequence = alt
    return model.AlternativeVariantAllele(
        allele_sequence=allele_sequence,
        allele_type=allele_type,
        structural_variant_detail=sv["detail"] if sv else None,
        colocated_variants=colocated_variants,
        annotations=allele_annotations,
        predicted_molecular_consequences=consequences,
    )


# ---------------------------------------------------------------------------
# BGZF page-index seek path
#
# When the pipeline emits a `<vcf>.pageidx.json` sidecar (see
# pagination-design.md / build_page_index.py), a page can be fetched by seeking
# straight to it (via the _BgzfReader in bgzf.py) instead of scanning from the
# top with bcftools. The sidecar stores, every `stride` records, the packed BGZF
# virtual offset (compressed_block_offset << 16 | within_block_offset) of that
# record's line.
# ---------------------------------------------------------------------------
PAGE_INDEX_SUFFIX = ".pageidx.json"


def _load_page_index(vcf_path: FilePath) -> dict | None:
    """The parsed `<vcf>.pageidx.json` sidecar, or None if it doesn't exist."""
    index_path = Path(str(vcf_path) + PAGE_INDEX_SUFFIX)
    if not index_path.exists():
        return None
    return json.loads(index_path.read_text())


def _read_indexed_page(
    vcf_path: FilePath, index: dict, page: int, page_size: int
) -> tuple[str, str]:
    """Return (header_text, page_rows_text) for the requested page by seeking to
    the nearest checkpoint and reading forward. `page` is 1-based; a page past
    the end yields empty rows."""
    total = index["total_records"]
    stride = index["stride"]
    checkpoints = index["checkpoints"]
    header_end = index["header_end_voffset"]
    start = (max(page, 1) - 1) * page_size

    header_lines: list[bytes] = []
    rows: list[bytes] = []
    with _BgzfReader(str(vcf_path)) as reader:
        # Header = every line before the first data record.
        while reader.tell() < header_end:
            line = reader.readline()
            if not line:
                break
            header_lines.append(line)
        # Seek to the checkpoint at/before the page start, skip the remainder.
        if page_size > 0 and start < total:
            checkpoint = start // stride
            reader.seek(checkpoints[checkpoint])
            for _ in range(start - checkpoint * stride):
                reader.readline()
            for _ in range(min(page_size, total - start)):
                line = reader.readline()
                if not line:
                    break
                rows.append(line)

    return b"".join(header_lines).decode(), b"".join(rows).decode()




# ---------------------------------------------------------------------------
# Filtered-scan cache
#
# Pagination needs the match total, so a filtered page has to scan the whole
# file — and did so again for every page: on a 999k-record output an
# allele-frequency filter cost ~22s per page click, page 200 exactly as much as
# page 1. The scan's expensive half is evaluating the filters, not reading the
# file, so the ordinals of the matching records are remembered and later pages
# of the same filter set skip straight to rebuilding what they need.
#
# Keyed on the file's identity *and* mtime, so a regenerated output (the dev
# harness rewrites one fixed path) can never be served against a stale match
# set. Bounded to a handful of entries: a filter set over a million records is
# ~1M ints worst case, a few MB, and the common case is far smaller.
# ---------------------------------------------------------------------------
_SCAN_CACHE_MAX_ENTRIES = 4
_scan_cache: "OrderedDict[tuple, _ScanResult]" = OrderedDict()
_scan_cache_lock = threading.Lock()


@dataclass
class _ScanResult:
    """What a full filtered pass learned, minus the page itself."""

    matches: list[int]
    scanned_total: int
    stats: list[results_filters.FilterStat]


def _file_identity(vcf_path: FilePath) -> tuple | None:
    """This file as it is right now, or None if it cannot be stat'd.

    None is the "do not cache" answer, deliberately: without a stat there is no
    way to notice the file changing under us, and the dev harness rewrites one
    fixed path — so serving a stale answer is the failure to avoid, not a miss.
    """
    try:
        stat = Path(vcf_path).stat()
    except OSError:
        return None
    return (str(vcf_path), stat.st_mtime_ns, stat.st_size)


def _scan_cache_key(
    vcf_path: FilePath, filters: list[results_filters.ResultsFilter]
) -> tuple | None:
    """Identity of (this file as it is now, this exact filter set), or None if the
    file cannot be stat'd — in which case nothing is cached rather than risking a
    stale answer."""
    identity = _file_identity(vcf_path)
    if identity is None:
        return None
    # Every field a filter's behaviour depends on has to be in here. `match` and
    # `include_missing` both change which records pass without changing the
    # field or the threshold, so leaving either out serves one filter's results
    # for the other — which is exactly what happened to `include_missing`.
    condition = tuple(
        (f.field, f.operator, tuple(f.values), f.threshold, f.match, f.include_missing)
        for f in filters
    )
    return identity + (condition,)


def _scan_cache_get(key: tuple | None) -> _ScanResult | None:
    if key is None:
        return None
    with _scan_cache_lock:
        result = _scan_cache.get(key)
        if result is not None:
            _scan_cache.move_to_end(key)
        return result


def _scan_cache_put(key: tuple | None, result: _ScanResult) -> None:
    if key is None:
        return
    with _scan_cache_lock:
        _scan_cache[key] = result
        _scan_cache.move_to_end(key)
        while len(_scan_cache) > _SCAN_CACHE_MAX_ENTRIES:
            _scan_cache.popitem(last=False)


def clear_scan_cache() -> None:
    """Drop every cached scan (tests; and anything that rewrites outputs)."""
    with _scan_cache_lock:
        _scan_cache.clear()


def _get_filtered_results(
    page_size: int,
    page: int,
    vcf_path: FilePath,
    filters: list[results_filters.ResultsFilter],
    spec: ParsingSpec | None = None,
) -> model.VepResultsResponse:
    """Stream the results VCF applying the filter pipeline, retaining only the
    requested page of survivors. The page-index fast path can't be used once
    records are filtered (positions shift), so this is a full sequential pass —
    but the file is read as a lazy line stream and only the page slice is held, so
    memory is bounded by `page_size` rather than the (multi-GB) file. Attaches
    per-filter removed counts to the response metadata and logs them.

    Note: pagination needs the total match count, so every page still scans the
    whole file. A filtered-index cache keyed by the filter set would remove the
    rescan for later pages (see pagination-design.md); memory is no longer the
    constraint it was."""
    page = max(page, 1)
    page_size = max(page_size, 0)
    start = (page - 1) * page_size
    cache_key = _scan_cache_key(vcf_path, filters)

    header_lines: list[str] = []
    with gzip.open(vcf_path, "rt") as handle:
        # The header is every '#' line, and all of them precede the data records;
        # stop at (and keep) the first data line, then stream the rest lazily so
        # the whole file is never materialised.
        first_data_line: str | None = None
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                first_data_line = line
                break

        index_map = csq_index_map_from_header(header_lines)
        compiled = results_filters.compile_filters(filters, index_map, spec)

        data_lines = (
            handle
            if first_data_line is None
            else itertools.chain((first_data_line,), handle)
        )
        cached = _scan_cache_get(cache_key)
        if cached is not None:
            # The answer is already known; only this page's records need
            # rebuilding, and no filter is evaluated for the rest.
            page_lines = results_filters.replay_matches(
                data_lines, compiled, cached.matches, start=start, count=page_size
            )
            outcome = results_filters.FilterOutcome(
                page=page_lines,
                matched_total=len(cached.matches),
                scanned_total=cached.scanned_total,
                stats=cached.stats,
            )
        else:
            matches: list[int] = []
            outcome = results_filters.filter_records(
                data_lines,
                compiled,
                start=start,
                count=page_size,
                record_matches=matches,
            )
            _scan_cache_put(
                cache_key,
                _ScanResult(
                    matches=matches,
                    scanned_total=outcome.scanned_total,
                    stats=outcome.stats,
                ),
            )

    stream = StringIO("".join(header_lines) + "".join(outcome.page))
    response = get_results_from_stream(
        page_size, page, outcome.matched_total, stream, presliced=True, spec=spec
    )
    response.metadata.filters = model.FilterMetadata(
        unfiltered_total=outcome.scanned_total,
        filtered_total=outcome.matched_total,
        stats=[
            model.FilterStat(field=stat.field, removed=stat.removed)
            for stat in outcome.stats
        ],
    )
    logging.info(
        "VEP results filtered: %d -> %d records (%s)",
        outcome.scanned_total,
        outcome.matched_total,
        ", ".join(f"{stat.field} removed {stat.removed}" for stat in outcome.stats)
        or "no active filters",
    )
    return response


def stream_filtered_vcf_text(
    vcf_path: FilePath,
    filters: list[results_filters.ResultsFilter],
) -> Iterator[str]:
    """A lazy text-line stream of the results VCF reduced to just the CSQ entries
    (and records) that pass `filters`: every header line unchanged, then each kept
    record rebuilt with its CSQ narrowed to the surviving entries.

    The download counterpart of `_get_filtered_results` — same compile-once,
    stream-the-file pipeline, but yielding the whole matched set as VCF text
    rather than parsing a single page into the response model. Feeds both the
    filtered VCF download (gzip this directly) and the filtered TSV download
    (flatten it first).

    Filters are compiled eagerly against the file's CSQ header, so an invalid
    filter raises `results_filters.FilterError` before any streaming begins (the
    download endpoint maps that to a 400). Memory stays bounded: the file is read
    as a lazy line stream and survivors are yielded one at a time, never
    collected."""
    header_lines: list[str] = []
    with gzip.open(vcf_path, "rt") as handle:
        # Every '#' line is header and all precede the data records; stop at the
        # first data line (the header is tiny — a cheap eager read for compile).
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break
    index_map = csq_index_map_from_header(header_lines)
    spec = _load_pinned_spec(vcf_path)
    compiled = results_filters.compile_filters(filters, index_map, spec)

    def generate() -> Iterator[str]:
        yield from header_lines
        with gzip.open(vcf_path, "rt") as handle:
            data_lines = (line for line in handle if not line.startswith("#"))
            yield from results_filters.stream_filtered_lines(data_lines, compiled)

    return generate()


# ---------------------------------------------------------------------------
# The pinned spec sidecar, cached per file.
#
# It is a job's own frozen copy of the spec, so for a given output it can only
# change if the file itself is rewritten — yet reading it means parsing and
# validating a large JSON document, and a single request did that *twice*
# (`_load_pinned_spec` goes through this one) before paging did it all again for
# every page. Same key discipline as the scan cache: the file's identity now, so
# a regenerated output is never served against a stale spec.
#
# A None result is cached too. "This output has no sidecar" is an answer worth
# keeping — otherwise the pre-pin jobs are the ones that re-read on every page.
# ---------------------------------------------------------------------------
_SPEC_CACHE_MAX_ENTRIES = 4
_spec_cache: "OrderedDict[tuple, MergedSpec | None]" = OrderedDict()
_spec_cache_lock = threading.Lock()


def clear_spec_cache() -> None:
    """Drop every cached pinned spec (tests; and anything that rewrites outputs)."""
    with _spec_cache_lock:
        _spec_cache.clear()


def _load_pinned_merged_spec(vcf_path: FilePath) -> MergedSpec | None:
    """The whole merged spec document pinned to this job, loaded defensively.

    Never raises: an output with no sidecar (pre-dating the pin) or an
    unreadable one returns None.
    """
    key = _file_identity(vcf_path)
    if key is not None:
        with _spec_cache_lock:
            if key in _spec_cache:
                _spec_cache.move_to_end(key)
                return _spec_cache[key]

    try:
        merged = load_spec_sidecar(vcf_path)
    except Exception as exc:
        logging.warning(
            "Ignoring unreadable spec sidecar for %s: %s", vcf_path, exc
        )
        # Not cached: an unreadable sidecar is a fault that may be repaired
        # without the output changing, and retrying it costs one read.
        return None
    if merged is None:
        logging.debug(
            "No spec sidecar for %s; no annotations will be emitted", vcf_path
        )

    if key is not None:
        with _spec_cache_lock:
            _spec_cache[key] = merged
            _spec_cache.move_to_end(key)
            while len(_spec_cache) > _SPEC_CACHE_MAX_ENTRIES:
                _spec_cache.popitem(last=False)
    return merged


def _load_pinned_spec(vcf_path: FilePath) -> ParsingSpec | None:
    """The parsing spec pinned to this job at submission, loaded defensively.

    Since the go-flat cutover this spec is the sole source of annotation data:
    every plugin payload on the response comes from driving it through
    spec_interpreter.apply_plugin_spec.

    Never raises: an output with no sidecar (pre-dating the pin) or an
    unreadable one still parses, just with no annotations, so both fall back to
    None.

    The pinned sidecar is now the whole merged document; the parsing half is what
    the results path needs, so that is what this returns.
    """
    merged = _load_pinned_merged_spec(vcf_path)
    if merged is None:
        return None
    spec = merged.parsing
    logging.info("Loaded pinned parsing spec %s for %s", spec.spec_version, vcf_path)
    return spec


def _read_csq_columns(vcf_path: FilePath) -> set[str] | None:
    """The CSQ column names declared in the output VCF header — the fixed layout
    for the whole file, so a set is enough to check presence. Reads only the
    header (stops at the first data line). None if there is no CSQ header line or
    the file can't be read."""
    header_lines: list[str] = []
    try:
        with gzip.open(vcf_path, "rt") as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                header_lines.append(line)
    except OSError:
        return None
    index_map = csq_index_map_from_header(header_lines)
    return set(index_map) or None


def _load_expected_columns(vcf_path: FilePath) -> set[str] | None:
    """The CSQ columns this job's submitted options require, pinned at submission
    (`expected_columns.json`). Defensive: a job predating the pin (no sidecar) or
    an unreadable one returns None. Shared by the missing-column check and the
    AF-source gating, so the sidecar is read once."""
    try:
        return load_expected_columns_sidecar(vcf_path)
    except Exception as exc:
        logging.warning(
            "Ignoring unreadable expected-columns sidecar for %s: %s", vcf_path, exc
        )
        return None


def _check_expected_columns(vcf_path: FilePath, expected: set[str] | None) -> None:
    """Warn if any CSQ column the submitted options require is missing from the
    output header (the runtime missing-expected-field check, design §6.2). A
    missing expected column is a real contract breach — a plugin the user enabled
    produced no column — while extra columns are always tolerated.

    Missing columns currently log warnings and never fail results; a missing pin
    (output predating this) or an unreadable header is a no-op. A retry or
    failure policy requires an explicit workflow contract.
    """
    if not expected:
        return
    actual = _read_csq_columns(vcf_path)
    if actual is None:
        logging.warning(
            "No CSQ header to check expected columns against for %s", vcf_path
        )
        return
    missing = expected - actual
    if missing:
        logging.warning(
            "VEP output %s is missing %d expected CSQ column(s): %s",
            vcf_path, len(missing), ", ".join(sorted(missing)),
        )


def _load_pinned_display_panels(vcf_path: FilePath) -> list[DisplayPanel] | None:
    """The option panels pinned to this job at submission, loaded defensively.

    Never raises: a job submitted before this pin existed (no sidecar), or an
    unreadable one, returns None — the results view then falls back to the live
    form-config panels, exactly as it did before pinning.
    """
    try:
        panels = load_display_panels_sidecar(vcf_path)
    except Exception as exc:
        logging.warning(
            "Ignoring unreadable display-panels sidecar for %s: %s", vcf_path, exc
        )
        return None
    if not panels:
        # None (no sidecar) and [] are both "nothing usable pinned". An empty
        # list can only come from a corrupted sidecar — get_visible_panels always
        # returns at least the always-visible panels — and treating it as a valid
        # pin would render a job with no panels at all rather than falling back.
        logging.debug(
            "No display-panels sidecar for %s; results will use the live panels",
            vcf_path,
        )
        return None
    return panels


def _resolve_display_payload(spec: MergedSpec | None) -> DisplayPayload | None:
    """The display layout to render this job's annotations with.

    Normally the job's *pinned* spec owns it, like everything else about the
    job. But every job submitted before the display section existed has a pinned
    spec with no `display` key, and would otherwise render its twelve
    spec-driven options blank. For those — and only those — fall back to the
    current genome's display spec, resolved from the pinned spec's own genome.

    This deliberately reintroduces a little staleness, narrowly: only for
    pre-display-section jobs, and only for labels/formats. *Which* options ran,
    and how their columns were parsed, still come from the pin.
    """
    if spec is None:
        return None
    payload = spec.display_payload()
    if payload is not None:
        return payload
    assembly = (spec.genome or {}).get("assembly", "")
    try:
        current = resolve_merged_spec(assembly)
    except Exception as exc:
        logging.warning(
            "No current spec to supply a display section for a pinned spec "
            "without one (assembly %r): %s", assembly, exc
        )
        return None
    logging.debug(
        "Pinned spec has no display section; using the current %r display spec",
        assembly,
    )
    current_payload = current.display_payload()
    if current_payload is None:
        return None
    # The scopes must still describe the *pinned* parsers, since those are what
    # produced this job's annotations; only the layout comes from the current
    # spec. Copied rather than rebuilt field by field: listing the fields here
    # meant a new one silently did not reach these jobs, and the rating scales
    # had already been lost that way.
    return current_payload.model_copy(update={"plugin_scopes": spec.plugin_scopes()})


def _drop_form_only_help(
    panels: list[DisplayPanel] | None, merged: MergedSpec | None
) -> list[DisplayPanel] | None:
    """The panels with form-only help removed from their options.

    The pin stays lossless — it records what was submitted — so this happens on
    the way out rather than on the way in. An option marked `form_only` keeps
    its (?) on the input form and loses it here, where its rows already carry
    help of their own.
    """
    if panels is None or merged is None or merged.help is None:
        return panels
    form_only = merged.help.form_only_options()
    if not form_only:
        return panels
    for panel in panels:
        for option in panel.options:
            if option.id in form_only:
                option.help = None
    return panels


def _with_display_panels(
    response: model.VepResultsResponse,
    panels: list[DisplayPanel] | None,
    display: DisplayPayload | None = None,
    *,
    spec: ParsingSpec | None = None,
    expected_columns: set[str] | None = None,
) -> model.VepResultsResponse:
    """Attach the pinned panels and display layout to a response built by the
    parsing path (which knows about neither). None leaves the field absent.

    Also gate the allele-frequency data to the AF columns the submission actually
    *selected* (the pinned expected columns), not merely whatever the output VCF
    carries — two facets of the same full-cache leak: `available_af_sources`
    (the filter's availability) and each AF annotation's populations
    (`_gate_af_populations`). Without a pin (older jobs) both are left as the VCF
    reported them.
    """
    response.metadata.display_panels = panels
    response.metadata.display = display
    # AF is allele-scoped, so its annotations hang off the alt alleles.
    alleles = [
        allele
        for variant in response.variants
        for allele in variant.alternative_alleles
    ]
    if expected_columns is not None:
        response.metadata.available_af_sources = [
            source
            for source in response.metadata.available_af_sources
            if source.key in expected_columns
        ]
        # Same full-cache leak as the AF sources: the VCF may carry impact-score
        # columns this submission never selected. The gate is the score's
        # sentinel column, not the column the filter tests — a plugin's
        # `csq_fields` (and so `expected_columns`) is deliberately
        # under-declared, e.g. SpliceAI declares only its AG column while
        # emitting all four (see results_filters.ScoreSpec).
        response.metadata.available_scores = [
            field
            for field in response.metadata.available_scores
            if results_filters.SCORE_SPECS[field].gate in expected_columns
        ]
        _gate_af_columns(alleles, spec, expected_columns)
    # Decode each All of Us annotation's max-subpopulation code(s) to a label,
    # after any gating (a gated job nulls an unselected max, so a leaked one is
    # not relabelled). Serve-time metadata only — not part of the spec digest.
    _label_af_max_subpopulation(alleles)
    return response


def _gate_af_columns(
    alleles: Iterable[model.AlternativeVariantAllele],
    spec: ParsingSpec | None,
    expected_columns: set[str],
) -> None:
    """Drop allele-frequency values the submission didn't select.

    An AF plugin reads a `pattern_map` of per-population columns plus scalar
    columns (the overall AF, and All of Us's max-subpopulation), matching every
    column *present in the VCF* — a full-cache dev output carries them all, so
    without this the served annotation shows them even when only a few were
    selected. (Production emits only the selected columns, so it never leaks
    there.) An AF plugin is one with a `pattern_map` target; for those, keep a
    population only when its `f"{prefix}{key}{suffix}"` column is in the pinned
    expected set, and null a scalar field whose source column isn't. Every non-AF
    plugin, and the frequency envelope itself, is untouched.
    """
    if spec is None:
        return
    # plugin id -> (population gates, scalar gates), for the AF plugins only (a
    # `pattern_map` target marks one). A population key maps to its column as
    # `f"{prefix}{key}{suffix}"`; a scalar field maps to its `from` column.
    gates: dict[str, tuple[list, list]] = {}
    for plugin in spec.plugins:
        population_gates = [
            (target.field, *pattern_affixes(target.from_pattern))
            for target in plugin.targets
            if target.transform == "pattern_map" and target.from_pattern
        ]
        if not population_gates:
            continue
        scalar_gates = [
            (target.field, target.source)
            for target in plugin.targets
            if target.transform == "scalar" and target.source
        ]
        gates[plugin.plugin] = (population_gates, scalar_gates)
    if not gates:
        return
    for allele in alleles:
        for annotation in allele.annotations:
            plugin_gates = gates.get(annotation.plugin)
            if plugin_gates is None:
                continue
            population_gates, scalar_gates = plugin_gates
            for field, prefix, suffix in population_gates:
                populations = annotation.data.get(field)
                if isinstance(populations, dict):
                    annotation.data[field] = {
                        key: value
                        for key, value in populations.items()
                        if f"{prefix}{key}{suffix}" in expected_columns
                    }
            for field, source in scalar_gates:
                if source not in expected_columns:
                    annotation.data[field] = None


def _label_af_max_subpopulation(
    alleles: Iterable[model.AlternativeVariantAllele],
) -> None:
    """Attach a decoded `max_subpopulation_label` to each AF annotation carrying a
    `max_subpopulation` code — the All of Us subpopulation(s) the max AF came
    from (only that plugin emits the field). The frontend renders the label;
    decoding here keeps the population vocabulary defined only in
    `form_panels.py`. A null/absent code (unselected, or gated away) is skipped."""
    for allele in alleles:
        for annotation in allele.annotations:
            raw = annotation.data.get("max_subpopulation")
            if raw:
                annotation.data["max_subpopulation_label"] = (
                    af_max_subpopulation_label(raw)
                )


def _read_row_slice(vcf_path: FilePath, row_offset: int, page_size: int) -> str:
    """The `page_size` lines ending at `row_offset` of `bcftools view` output.

    Replaces a `bcftools view … | head -nX | tail -nY` shell pipeline. That
    interpolated the VCF path into a shell string, and the path derives from the
    uploaded file's name — so a filename like `a$(…).vcf` ran whatever it liked.

    The pipeline's one virtue is kept: `head` stopped bcftools early rather than
    reading the whole file, so this stops reading at `row_offset` and closes the
    pipe, which SIGPIPEs bcftools exactly as before.
    """
    if page_size <= 0:
        return ""
    keep: deque[str] = deque(maxlen=page_size)
    process = subprocess.Popen(
        ["bcftools", "view", str(vcf_path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        for line_number, line in enumerate(process.stdout, start=1):
            keep.append(line)
            if line_number >= row_offset:
                break
    finally:
        process.stdout.close()
        process.wait()
    return "".join(keep)


def get_results_from_path(
    page_size: int,
    page: int,
    vcf_path: FilePath,
    filters: list[results_filters.ResultsFilter] | None = None,
) -> model.VepResultsResponse:
    """Returns a page of VCF data from the given filepath.
    Slices the input VCF file to a smaller one
    and converts it to stream for get_results_from_stream"""

    # Load the spec pinned to this job at submission. It drives the generic
    # `annotations` on every allele and transcript consequence (threaded down to
    # _get_alt_allele_details). A missing or unreadable pin -> None -> no
    # annotations, never failing results.
    merged = _load_pinned_merged_spec(vcf_path)
    spec = _load_pinned_spec(vcf_path)
    # The CSQ columns this job's submitted options require (pinned at submission).
    # Drives the missing-column check below and gates which AF sources the filter
    # UI is offered — so an AF filter is only available when AF was actually
    # selected, not merely present in the output VCF.
    expected_columns = _load_expected_columns(vcf_path)
    # Runtime missing-expected-field check: warn if the pipeline output is missing
    # a CSQ column the submitted options required. Non-fatal (dev warns only).
    _check_expected_columns(vcf_path, expected_columns)
    # The option panels this job was submitted against (None for older jobs).
    display_panels = _drop_form_only_help(
        _load_pinned_display_panels(vcf_path), merged
    )
    # How those options lay out, from the pin (or the current spec for jobs
    # pinned before the display section existed).
    display = _resolve_display_payload(merged)

    # Filtered requests can't use the page index (filtering shifts record
    # positions), so they take a dedicated scan-and-filter path.
    if filters:
        return _with_display_panels(
            _get_filtered_results(page_size, page, vcf_path, filters, spec),
            display_panels,
            display,
            spec=spec,
            expected_columns=expected_columns,
        )

    # Fast path: if the pipeline emitted a page-index sidecar, seek to the page
    # instead of scanning the file / shelling out to bcftools.
    index = _load_page_index(vcf_path)
    if index is not None:
        page = max(page, 1)
        page_size = max(page_size, 0)
        header_text, rows_text = _read_indexed_page(vcf_path, index, page, page_size)
        return _with_display_panels(
            get_results_from_stream(
                page_size,
                page,
                index["total_records"],
                StringIO(header_text + rows_text),
                presliced=True,
                spec=spec,
            ),
            display_panels,
            display,
            spec=spec,
            expected_columns=expected_columns,
        )

    # Fallback (no sidecar): scan the file from the top through page*page_size
    # records and shell out to bcftools for the counts. `head` short-circuits so
    # it stops at the offset rather than scanning the whole file, but deep pages
    # get slower and the last page is a full pass. Runs from the pipeline now ship
    # a page-index sidecar (handled above); this remains for older/un-indexed
    # outputs. Longer term, a queryable store (SQLite/Parquet) would also enable
    # sorting/filtering (see pagination-design.md).
    # Fetch a pageful of variant records with headers
    vcf_info = _get_vcf_meta(vcf_path)
    total = vcf_info.variant_count
    page = max(page, 1) # normalize values
    page_size = min(max(page_size, 0), total)
    # This page's own bounds. `_read_row_slice` returns the lines *ending* at
    # `row_offset`, so the count it is given is what decides where the page
    # starts — and a last page shorter than `page_size` must ask for only the
    # records that remain. Asking for a full one returned a full page ending at
    # the final record, i.e. records the previous page had already shown.
    first_record = (page - 1) * page_size
    last_record = min(page * page_size, total)
    row_offset = last_record + vcf_info.header_count
    vcf_headers = subprocess.check_output(  # fetch all header rows
        ["bcftools", "view", "-h", str(vcf_path)], text=True
    )
    vcf_slice = _read_row_slice(vcf_path, row_offset, last_record - first_record)
    vcf_stream = StringIO(vcf_headers + vcf_slice)

    return _with_display_panels(
        get_results_from_stream(page_size, page, total, vcf_stream, spec=spec),
        display_panels,
        display,
        spec=spec,
        expected_columns=expected_columns,
    )


def get_results_from_stream(
    page_size: int, page: int, total: int, vcf_stream: StringIO,
    presliced: bool = False, spec: ParsingSpec | None = None,
) -> model.VepResultsResponse:
    """Helper method to split a filestream into header and records.

    The header is kept as raw lines: the only thing read from it is the CSQ
    column list, which `csq_index_map_from_header` takes straight out of the
    `##INFO=<ID=CSQ` line."""
    header_lines: list[str] = []
    data_lines: list[str] = []
    for line in vcf_stream:
        (header_lines if line.startswith("#") else data_lines).append(line)
    return _get_results_from_records(
        page_size, page, total, header_lines, data_lines, presliced, spec
    )


def _get_results_from_records(
    page_size: int, page: int, total: int,
    header_lines: list[str], data_lines: list[str],
    presliced: bool = False, spec: ParsingSpec | None = None,
) -> model.VepResultsResponse:
    """Generates a page of VCF data in the format described in
    APISpecification.yaml for a given VCFPY reader"""

    # Parse csq header
    prediction_index_map = csq_index_map_from_header(header_lines)
    if not prediction_index_map:
        raise Exception("CSQ header missing")
    # Required CSQ column (the rest use fallback values)
    if "Allele" not in prediction_index_map:
        raise Exception("Allele column missing from CSQ header")

    # Resolve every plugin against this file's CSQ header once. The header is
    # fixed for the file, so which columns a plugin reads, whether it ran at
    # all, and which columns a pattern_map matches are all answerable here
    # instead of on every CSQ row. See PluginPlan.
    plans = compile_parsing_spec(prediction_index_map, spec) if spec else None

    variants = []
    # populate variants page. `presliced` means the stream already contains
    # exactly this page's rows (the index seek path), so the page-bounds guard —
    # which the scan path needs to return empty past the end — is skipped.
    #
    # The guard asks where the page *starts*, not where it ends. Asking
    # `page * page_size <= total` treated a page that merely runs past the last
    # record as being past the end, so the final partial page came back empty:
    # 21 records at 5 a page served four pages and lost the 21st. It is the
    # first record that decides whether a page exists at all.
    if presliced or (page - 1) * page_size < total:
        for record in read_records(data_lines):
            if record.CHROM.startswith("chr"):
                record.CHROM = record.CHROM[3:]

            # https://github.com/bihealth/vcfpy/blob/697768d032b6b476766fb4c524c91c8d24559330/vcfpy/record.py#L63
            # end does not look like it is implemented.
            # https://github.com/Penghui-Wang/PyVCF/blob/master/vcf/model.py#L190
            # from competing vcf module
            location = model.Location(
                region_name=record.CHROM,
                start=record.POS,
                end=record.POS + len(record.REF),
            )

            if "CSQ" not in record.INFO:
                csq_strings = []
                alt_allele_strings = [_alt_value(alt) for alt in record.ALT]
            else:
                csq_strings = record.INFO["CSQ"]
                # `dict.fromkeys`, not `set`: both deduplicate, but a set's
                # iteration order depends on string hashes, which Python
                # randomises per process. The alleles came back in a different
                # order in every worker and after every restart, and this order
                # is the order they are displayed in — so a variant's alleles
                # visibly rearranged themselves for no reason, and anything
                # diffing or caching a response saw changes that were not there.
                # First appearance in the CSQ rows is the natural order, and it
                # matches the no-CSQ branch above, which follows record.ALT.
                alt_allele_strings = list(dict.fromkeys(
                    csq_string.split("|")[prediction_index_map["Allele"]]
                    for csq_string in csq_strings
                ))

            sv = _structural_info(record)

            alt_alleles = [
                _get_alt_allele_details(
                    record.REF, alt, csq_strings, prediction_index_map, spec, sv,
                    plans,
                    str(record.CHROM),
                    str(record.POS),
                )
                for alt in alt_allele_strings
            ]

            alt_values = [_alt_value(a) for a in record.ALT]
            longest_alt = max(alt_values, key=len) if alt_values else ""

            variant = model.Variant(
                name=";".join(record.ID) if len(record.ID) > 0 else ".",
                location=location,
                reference_allele=model.ReferenceVariantAllele(
                    allele_sequence=record.REF
                ),
                alternative_alleles=alt_alleles,
                allele_type=(
                    sv["type_word"]
                    if sv
                    else _get_variant_type(record.REF, longest_alt)
                ),
            )
            # One pool per variant, once it holds every allele and consequence.
            _pool_annotations(variant)
            variants.append(variant)

    available_af_sources = [
        model.AfSource(**descriptor)
        for descriptor in (
            results_filters.af_source_descriptor(column, spec)
            for column in results_filters.af_columns(prediction_index_map, spec)
        )
        if descriptor
    ]

    # Which impact scores this output carries: a field is offered when any of
    # the columns its predicate reads is in the CSQ header. That is stage one of
    # two — selection is applied later, with the AF gate, since only the pinned
    # expected columns know what the submission actually chose.
    available_scores = [
        field
        for field, spec in results_filters.SCORE_SPECS.items()
        if any(column in prediction_index_map for column in spec.columns)
    ]

    return model.VepResultsResponse(
        metadata=model.Metadata(
            pagination=model.PaginationMetadata(
                page=page, per_page=page_size, total=total
            ),
            available_af_sources=available_af_sources,
            available_scores=available_scores,
        ),
        variants=variants,
    )
