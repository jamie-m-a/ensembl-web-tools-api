""" Module for loading a VCF and parsing it into a VepResultsResponse
object as defined in APISpecification"""

from io import StringIO
import gzip
import json
import logging
import re
import subprocess
from pathlib import Path
from pydantic import FilePath
import vcfpy
from vep.models import vcf_results_model as model
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
from vep.utils.spec_loader import (
    load_display_panels_sidecar,
    load_expected_columns_sidecar,
    load_spec_sidecar,
)
from vep.models.display_panels_model import DisplayPanel
from vep.utils.spec_interpreter import apply_plugin_spec
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


def _spec_annotations(
    csq_values: list[str],
    index_map: dict[str, int],
    spec: ParsingSpec | None,
    scope: str,
) -> list[model.Annotation]:
    """The generic annotations for one CSQ entry at the given scope, driving each
    matching spec plugin through `apply_plugin_spec`. Additive to the typed
    fields; when there is no pinned spec this is empty and nothing changes."""
    if spec is None:
        return []
    annotations: list[model.Annotation] = []
    for plugin in spec.plugins:
        if plugin.scope != scope:
            continue
        data = apply_plugin_spec(csq_values, index_map, plugin)
        if data is not None:
            annotations.append(
                model.Annotation(plugin=plugin.plugin, scope=scope, data=data)
            )
    return annotations


def _get_alt_allele_details(
    ref: str,
    alt: str,
    csqs: list[str],
    index_map: dict[str, int],
    spec: ParsingSpec | None = None,
) -> model.AlternativeVariantAllele:
    """Creates  AlternativeVariantAllele based on
    target alt allele and CSQ entires"""
    consequences = []
    allele_type = _get_variant_type(ref, alt)
    # Allele-level annotations are identical across all of this allele's CSQ
    # rows, so capture them once (from the first matching row). They are also
    # the only annotations available for intergenic variants (no transcript
    # rows).
    colocated_variants: list[str] = []
    allele_annotations: list[model.Annotation] = []
    allele_level_captured = False

    for str_csq in csqs:
        csq_values = str_csq.split("|")

        if csq_values[index_map["Allele"]] != alt:
            continue

        if not allele_level_captured:
            colocated_variants = split_amp(
                get_csq_value(csq_values, "Existing_variation", None, index_map)
            )
            allele_annotations = _spec_annotations(
                csq_values, index_map, spec, "allele"
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
                    # Protein & functional annotations
                    ensembl_protein_id=get_csq_value(
                        csq_values, "ENSP", None, index_map
                    ),
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
                        csq_values, index_map, spec, "transcript"
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

    return model.AlternativeVariantAllele(
        allele_sequence=("" if alt=="copy_number_variation" else alt),
        allele_type=allele_type,
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



def _get_filtered_results(
    page_size: int,
    page: int,
    vcf_path: FilePath,
    filters: list[results_filters.ResultsFilter],
    spec: ParsingSpec | None = None,
) -> model.VepResultsResponse:
    """Scan the whole results VCF applying the filter pipeline, then paginate the
    filtered records. The page-index fast path can't be used once records are
    filtered (positions shift), so this is a full sequential pass. Attaches
    per-filter removed counts to the response metadata and logs them.

    Note: this loads the kept records into memory and rescans per request. Fine
    for current result sizes; a filtered-index cache keyed by the filter set
    would remove the rescan later (see pagination-design.md)."""
    header_lines: list[str] = []
    data_lines: list[str] = []
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            (header_lines if line.startswith("#") else data_lines).append(line)

    index_map = csq_index_map_from_header(header_lines)
    compiled = results_filters.compile_filters(filters, index_map)
    kept, stats = results_filters.apply_filter_pipeline(data_lines, compiled)

    filtered_total = len(kept)
    page = max(page, 1)
    page_size = max(page_size, 0)
    start = (page - 1) * page_size
    page_rows = kept[start : start + page_size] if page_size > 0 else []

    stream = StringIO("".join(header_lines) + "".join(page_rows))
    response = get_results_from_stream(
        page_size, page, filtered_total, stream, presliced=True, spec=spec
    )
    response.metadata.filters = model.FilterMetadata(
        unfiltered_total=len(data_lines),
        filtered_total=filtered_total,
        stats=[
            model.FilterStat(field=stat.field, removed=stat.removed)
            for stat in stats
        ],
    )
    logging.info(
        "VEP results filtered: %d -> %d records (%s)",
        len(data_lines),
        filtered_total,
        ", ".join(f"{stat.field} removed {stat.removed}" for stat in stats)
        or "no active filters",
    )
    return response


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
    try:
        merged = load_spec_sidecar(vcf_path)
    except Exception as exc:
        logging.warning(
            "Ignoring unreadable spec sidecar for %s: %s", vcf_path, exc
        )
        return None
    if merged is None:
        logging.debug(
            "No spec sidecar for %s; no annotations will be emitted", vcf_path
        )
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


def _check_expected_columns(vcf_path: FilePath) -> None:
    """Warn if any CSQ column the submitted options require is missing from the
    output header (the runtime missing-expected-field check, design §6.2). A
    missing expected column is a real contract breach — a plugin the user enabled
    produced no column — while extra columns are always tolerated.

    Dev only warns and never fails results; a missing pin (output predating this)
    or an unreadable header is a no-op. Production would rerun the pipeline to
    regenerate the headers, capped at 3 retries (decision 15); that path needs
    the real pipeline and is not wired into this dev loop.
    """
    try:
        expected = load_expected_columns_sidecar(vcf_path)
    except Exception as exc:
        logging.warning(
            "Ignoring unreadable expected-columns sidecar for %s: %s", vcf_path, exc
        )
        return
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


def _with_display_panels(
    response: model.VepResultsResponse, panels: list[DisplayPanel] | None
) -> model.VepResultsResponse:
    """Attach the pinned panels to a response built by the parsing path (which
    knows nothing about them). None leaves the field absent."""
    response.metadata.display_panels = panels
    return response


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
    spec = _load_pinned_spec(vcf_path)
    # Runtime missing-expected-field check: warn if the pipeline output is missing
    # a CSQ column the submitted options required. Non-fatal (dev warns only).
    _check_expected_columns(vcf_path)
    # The option panels this job was submitted against (None for older jobs).
    display_panels = _load_pinned_display_panels(vcf_path)

    # Filtered requests can't use the page index (filtering shifts record
    # positions), so they take a dedicated scan-and-filter path.
    if filters:
        return _with_display_panels(
            _get_filtered_results(page_size, page, vcf_path, filters, spec),
            display_panels,
        )

    # Fast path: if the pipeline emitted a page-index sidecar, seek to the page
    # instead of scanning the file / shelling out to bcftools.
    index = _load_page_index(vcf_path)
    if index is not None:
        page = max(page, 1)
        page_size = max(page_size, 0)
        header_text, rows_text = _read_indexed_page(vcf_path, index, page, page_size)
        return _with_display_panels(get_results_from_stream(
            page_size,
            page,
            index["total_records"],
            StringIO(header_text + rows_text),
            presliced=True,
            spec=spec,
        ), display_panels)

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
    row_offset = min(page * page_size, total) + vcf_info.header_count
    vcf_headers = subprocess.check_output( # fetch all header rows
        f"bcftools view -h {vcf_path}", shell=True, text=True
    )
    vcf_slice = subprocess.check_output( # fetch subset of variant rows
        f"bcftools view {vcf_path} | head -n{row_offset} | tail -n{page_size}",
        shell=True, text=True
    )
    vcf_stream = StringIO(vcf_headers + vcf_slice)

    return _with_display_panels(
        get_results_from_stream(page_size, page, total, vcf_stream, spec=spec),
        display_panels,
    )


def get_results_from_stream(
    page_size: int, page: int, total: int, vcf_stream: StringIO,
    presliced: bool = False, spec: ParsingSpec | None = None,
) -> model.VepResultsResponse:
    """Helper method to convert a filestream to VCF records for _get_results_from_vcfpy"""

    # Load vcf
    vcf_records = vcfpy.Reader.from_stream(vcf_stream)
    return _get_results_from_vcfpy(page_size, page, total, vcf_records, presliced, spec)


def _get_results_from_vcfpy(
    page_size: int, page: int, total: int, vcf_records: vcfpy.Reader,
    presliced: bool = False, spec: ParsingSpec | None = None,
) -> model.VepResultsResponse:
    """Generates a page of VCF data in the format described in
    APISpecification.yaml for a given VCFPY reader"""

    # Parse csq header
    csq_header = vcf_records.header.get_info_field_info("CSQ").description
    if not csq_header:
        raise Exception("CSQ header missing")

    prediction_index_map = get_prediction_index_map(csq_header)
    # Required CSQ column (the rest use fallback values)
    if "Allele" not in prediction_index_map:
        raise Exception("Allele column missing from CSQ header")

    variants = []
    # populate variants page. `presliced` means the stream already contains
    # exactly this page's rows (the index seek path), so the page-bounds guard —
    # which the scan path needs to return empty past the end — is skipped.
    if presliced or page*page_size <= total:
        for record in vcf_records:
            if record is None:
                break
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
                alt_allele_strings = list(set([
                    csq_string.split("|")[prediction_index_map["Allele"]]
                    for csq_string in csq_strings
                ]))

            alt_alleles = [
                _get_alt_allele_details(
                    record.REF, alt, csq_strings, prediction_index_map, spec
                )
                for alt in alt_allele_strings
            ]

            longest_alt = max((_alt_value(a) for a in record.ALT), key=len)

            variants.append(
                model.Variant(
                    name=";".join(record.ID) if len(record.ID) > 0 else ".",
                    location=location,
                    reference_allele=model.ReferenceVariantAllele(
                        allele_sequence=record.REF
                    ),
                    alternative_alleles=alt_alleles,
                    allele_type=_get_variant_type(record.REF, longest_alt),
                )
            )

    available_af_sources = [
        model.AfSource(**descriptor)
        for descriptor in (
            results_filters.af_source_descriptor(column)
            for column in results_filters.af_columns(prediction_index_map)
        )
        if descriptor
    ]

    return model.VepResultsResponse(
        metadata=model.Metadata(
            pagination=model.PaginationMetadata(
                page=page, per_page=page_size, total=total
            ),
            available_af_sources=available_af_sources,
        ),
        variants=variants,
    )
