"""Tests for the filtered-download streaming path:

- results_filters.stream_filtered_lines — lazy filtered VCF data-line stream;
- tsv_export.flatten_vcf_lines — flatten any VCF line iterator to TSV rows;
- vcf_results.stream_filtered_vcf_text — header + filtered records from a file,
  the shared source for the filtered VCF and TSV downloads.

Filtered downloads reuse the same compile-once/stream-the-file pipeline as the
paginated results view, but stream the whole matched set rather than one page.
"""

import gzip

import pytest
from pydantic import FilePath

from app.vep.utils import results_filters as rf
from app.vep.utils.tsv_export import flatten_vcf_lines
from app.vep.utils.vcf_results import stream_filtered_vcf_text

# CSQ layout used across these tests: Consequence + the transcript-group columns.
CSQ_COLS = ["Allele", "Consequence", "Feature", "CANONICAL", "GENCODE_PRIMARY"]
CSQ_DESC = "Consequence annotations from Ensembl VEP. Format: " + "|".join(CSQ_COLS)
INDEX_MAP = {name: i for i, name in enumerate(CSQ_COLS)}


def _entry(allele, consequence, feature, canonical="", gencode=""):
    return "|".join([allele, consequence, feature, canonical, gencode])


def _record(pos, entries, chrom="chr1"):
    info = "CSQ=" + ",".join(entries)
    return f"{chrom}\t{pos}\t.\tC\tT\t.\t.\t{info}\n"


def _canonical_filter():
    return rf.ResultsFilter(
        field=rf.TRANSCRIPT_GROUP_FIELD, operator=rf.OPERATOR_IN, values=["canonical"]
    )


# --- stream_filtered_lines --------------------------------------------------


def test_stream_filtered_lines_narrows_entries_and_drops_records():
    compiled = rf.compile_filters([_canonical_filter()], INDEX_MAP)
    lines = [
        # keeps only the canonical entry (ENST_A), drops the non-canonical one
        _record(
            100,
            [
                _entry("T", "missense_variant", "ENST_A", canonical="YES"),
                _entry("T", "missense_variant", "ENST_B"),
            ],
        ),
        # no canonical entry -> whole record dropped
        _record(200, [_entry("T", "intron_variant", "ENST_C")]),
    ]

    kept = list(rf.stream_filtered_lines(iter(lines), compiled))

    assert len(kept) == 1
    entries = rf.extract_csq_entries(kept[0])
    assert [e[INDEX_MAP["Feature"]] for e in entries] == ["ENST_A"]


def test_stream_filtered_lines_no_filters_keeps_everything():
    lines = [_record(100, [_entry("T", "missense_variant", "ENST_A")])]
    kept = list(rf.stream_filtered_lines(iter(lines), []))
    assert kept == lines


# --- flatten_vcf_lines ------------------------------------------------------


def _vcf_lines(data_lines):
    return [
        "##fileformat=VCFv4.2\n",
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="{CSQ_DESC}">\n',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        *data_lines,
    ]


def test_flatten_vcf_lines_emits_header_and_one_row_per_entry():
    lines = _vcf_lines(
        [
            _record(
                100,
                [
                    _entry("T", "missense_variant", "ENST_A", canonical="YES"),
                    _entry("T", "synonymous_variant", "ENST_B"),
                ],
            )
        ]
    )

    rows = [r.rstrip("\n").split("\t") for r in flatten_vcf_lines(lines)]

    assert rows[0] == ["Uploaded_variation", "Location", "Ref"] + CSQ_COLS
    # two CSQ entries -> two data rows, location strips the "chr" prefix
    assert [r[1] for r in rows[1:]] == ["1:100", "1:100"]
    assert [r[5] for r in rows[1:]] == ["ENST_A", "ENST_B"]  # Feature column


def test_flatten_vcf_lines_over_a_prefiltered_line_emits_fewer_rows():
    # A line whose CSQ has already been narrowed (as stream_filtered_lines would
    # leave it) flattens to just those surviving entries.
    lines = _vcf_lines(
        [_record(100, [_entry("T", "missense_variant", "ENST_A", canonical="YES")])]
    )
    rows = list(flatten_vcf_lines(lines))
    assert len(rows) == 2  # header + one surviving entry


# --- stream_filtered_vcf_text (end to end, from a gzipped file) --------------


def _write_vcf(path, data_lines):
    with gzip.open(path, "wt") as handle:
        handle.writelines(_vcf_lines(data_lines))


def test_stream_filtered_vcf_text_preserves_header_and_filters_records(tmp_path):
    vcf = tmp_path / "results.vcf.gz"
    _write_vcf(
        vcf,
        [
            _record(
                100,
                [
                    _entry("T", "missense_variant", "ENST_A", canonical="YES"),
                    _entry("T", "missense_variant", "ENST_B"),
                ],
            ),
            _record(200, [_entry("T", "intron_variant", "ENST_C")]),  # dropped
        ],
    )

    out = list(stream_filtered_vcf_text(FilePath(vcf), [_canonical_filter()]))

    # every header line is preserved verbatim
    header = [line for line in out if line.startswith("#")]
    assert any(line.startswith("##fileformat") for line in header)
    assert any(line.startswith("##INFO=<ID=CSQ") for line in header)
    assert any(line.startswith("#CHROM") for line in header)

    data = [line for line in out if not line.startswith("#")]
    assert len(data) == 1  # the intron-only record was dropped
    entries = rf.extract_csq_entries(data[0])
    assert [e[INDEX_MAP["Feature"]] for e in entries] == ["ENST_A"]


def test_stream_filtered_vcf_text_feeds_the_tsv_flattener(tmp_path):
    # The filtered TSV download is flatten_vcf_lines(stream_filtered_vcf_text(...)):
    # only surviving entries reach the table.
    vcf = tmp_path / "results.vcf.gz"
    _write_vcf(
        vcf,
        [
            _record(
                100,
                [
                    _entry("T", "missense_variant", "ENST_A", canonical="YES"),
                    _entry("T", "missense_variant", "ENST_B"),
                ],
            )
        ],
    )

    rows = list(flatten_vcf_lines(stream_filtered_vcf_text(FilePath(vcf), [_canonical_filter()])))
    features = [r.rstrip("\n").split("\t")[5] for r in rows[1:]]
    assert features == ["ENST_A"]


def test_stream_filtered_vcf_text_raises_eagerly_on_invalid_filter(tmp_path):
    vcf = tmp_path / "results.vcf.gz"
    _write_vcf(vcf, [_record(100, [_entry("T", "missense_variant", "ENST_A")])])

    bad = rf.ResultsFilter(
        field=rf.TRANSCRIPT_GROUP_FIELD, operator=rf.OPERATOR_IN, values=["nonsense"]
    )
    # Compiled eagerly, so the error surfaces on the call, before streaming. Match
    # by message rather than class: vcf_results imports results_filters as
    # `vep.utils...` while this test imports `app.vep.utils...`, so the FilterError
    # classes differ by identity even though they are the same source (same reason
    # test_vep.py compares by value).
    with pytest.raises(Exception, match="unsupported transcript group"):
        stream_filtered_vcf_text(FilePath(vcf), [bad])
