"""Tests for the BGZF page-index seek path (app/vep/utils/vcf_results.py).

When a `<vcf>.pageidx.json` sidecar is present, get_results_from_path seeks to
the page via packed BGZF virtual offsets instead of scanning with bcftools.

These tests build their own BGZF fixtures in pure Python (a tiny block-splitting
writer that also reports the ground-truth virtual offset of every line), so they
run with no external tools/deps and deliberately span BGZF block boundaries.
"""

import gzip
import json
import struct
import subprocess
import zlib

import pytest
from pydantic import FilePath

from app.vep.utils.bgzf import _BgzfReader
from app.vep.utils.vcf_results import (
    _load_page_index,
    _read_indexed_page,
    get_results_from_path,
)

# --- pure-Python BGZF writer (fixtures) --------------------------------------

# The standard 28-byte BGZF end-of-file marker (an empty block).
BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def _bgzf_block(payload: bytes) -> bytes:
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = compressor.compress(payload) + compressor.flush()
    bsize = 12 + 6 + len(cdata) + 8 - 1  # total block size - 1
    header = (
        b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<H", 2)
        + struct.pack("<H", bsize)
    )
    trailer = struct.pack("<I", zlib.crc32(payload) & 0xFFFFFFFF) + struct.pack(
        "<I", len(payload) & 0xFFFFFFFF
    )
    return header + cdata + trailer


def write_bgzf(path, text: str, block_bytes: int = 64) -> list[int]:
    """Write `text` as BGZF, splitting the *uncompressed* stream into blocks of
    `block_bytes` (small, so lines cross block boundaries). Returns the packed
    virtual offset of the start of each line (ground truth, computed from block
    layout — independent of the reader under test)."""
    data = text.encode()
    blocks = [data[i : i + block_bytes] for i in range(0, len(data), block_bytes)] or [b""]

    out = bytearray()
    block_coffsets = []
    for block in blocks:
        block_coffsets.append(len(out))
        out += _bgzf_block(block)
    out += BGZF_EOF
    path.write_bytes(bytes(out))

    def voffset(global_uncompressed_offset: int) -> int:
        block = global_uncompressed_offset // block_bytes
        within = global_uncompressed_offset % block_bytes
        return (block_coffsets[block] << 16) | within

    line_voffsets, cursor = [], 0
    for line in text.splitlines(keepends=True):
        line_voffsets.append(voffset(cursor))
        cursor += len(line.encode())
    return line_voffsets


def write_indexed_vcf(tmp_path, text: str, *, stride: int, block_bytes: int = 64):
    """Write text as BGZF + its `.pageidx.json` sidecar; return the vcf path."""
    vcf_path = tmp_path / "results.vcf.gz"
    line_voffsets = write_bgzf(vcf_path, text, block_bytes=block_bytes)

    lines = text.splitlines(keepends=True)
    data_voffsets = [
        vo for line, vo in zip(lines, line_voffsets) if not line.startswith("#")
    ]
    index = {
        "version": 1,
        "vcf": vcf_path.name,
        "total_records": len(data_voffsets),
        "header_end_voffset": data_voffsets[0] if data_voffsets else 0,
        "stride": stride,
        "checkpoints": [data_voffsets[k] for k in range(0, len(data_voffsets), stride)],
    }
    (tmp_path / "results.vcf.gz.pageidx.json").write_text(json.dumps(index))
    return vcf_path


# --- fixtures ----------------------------------------------------------------

CSQ_DESC = (
    "Consequence annotations from Ensembl VEP. Format: "
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE"
)


def make_vcf(n_records: int) -> str:
    header = (
        "##fileformat=VCFv4.2\n"
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="{CSQ_DESC}">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    rows = []
    for i in range(1, n_records + 1):
        csq = f"T|missense_variant|MODERATE|GENE{i}|ENSG{i}|Transcript|ENST{i}|protein_coding"
        rows.append(f"chr1\t{100 + i}\tid_{i:02d}\tC\tT\t.\t.\tCSQ={csq}\n")
    return header + "".join(rows)


@pytest.fixture
def indexed_vcf(tmp_path):
    # 12 records, checkpoints every 4 records, tiny blocks -> multi-block file.
    return write_indexed_vcf(tmp_path, make_vcf(12), stride=4, block_bytes=64)


# --- _read_indexed_page ------------------------------------------------------


def test_read_indexed_page_matches_sequential(tmp_path):
    text = make_vcf(12)
    vcf_path = write_indexed_vcf(tmp_path, text, stride=4)
    index = _load_page_index(FilePath(vcf_path))

    header_seq, records_seq = [], []
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            (header_seq if line.startswith("#") else records_seq).append(line)

    for page, page_size in [(1, 5), (2, 5), (3, 5), (4, 5), (1, 1), (12, 1), (2, 4)]:
        header, rows = _read_indexed_page(FilePath(vcf_path), index, page, page_size)
        start = (page - 1) * page_size
        assert header.splitlines(keepends=True) == header_seq
        assert rows.splitlines(keepends=True) == records_seq[start : start + page_size]


def test_reader_seeks_across_block_boundaries(tmp_path):
    text = make_vcf(12)
    vcf_path = write_indexed_vcf(tmp_path, text, stride=4, block_bytes=48)
    index = _load_page_index(FilePath(vcf_path))
    # tiny blocks -> the file really is multi-block (more than one checkpoint slot)
    assert len(index["checkpoints"]) == 3

    records_seq = [
        l for l in gzip.open(vcf_path, "rt") if not l.startswith("#")
    ]
    with _BgzfReader(str(vcf_path)) as reader:
        # seek straight to checkpoint 2 (record index 8) and read it back
        reader.seek(index["checkpoints"][2])
        assert reader.readline().decode() == records_seq[8]


# --- get_results_from_path (end to end, via the index) -----------------------


def test_get_results_page_via_index(indexed_vcf):
    result = get_results_from_path(5, 2, FilePath(indexed_vcf))
    assert [v.name for v in result.variants] == [f"id_{i:02d}" for i in range(6, 11)]
    assert result.metadata.pagination.page == 2
    assert result.metadata.pagination.per_page == 5
    assert result.metadata.pagination.total == 12
    # "chr" prefix stripped, location parsed
    assert result.variants[0].location.region_name == "1"


def test_get_results_last_partial_page_via_index(indexed_vcf):
    result = get_results_from_path(5, 3, FilePath(indexed_vcf))
    assert [v.name for v in result.variants] == ["id_11", "id_12"]
    assert result.metadata.pagination.total == 12


def test_get_results_beyond_end_is_empty(indexed_vcf):
    result = get_results_from_path(5, 10, FilePath(indexed_vcf))
    assert result.variants == []
    assert result.metadata.pagination.total == 12


def test_index_path_does_not_shell_out(monkeypatch, indexed_vcf):
    # with a sidecar present, no bcftools subprocess should be invoked
    def boom(*args, **kwargs):
        raise AssertionError("subprocess called despite page index present")

    monkeypatch.setattr(subprocess, "check_output", boom)
    result = get_results_from_path(5, 1, FilePath(indexed_vcf))
    assert [v.name for v in result.variants] == [f"id_{i:02d}" for i in range(1, 6)]


def test_no_sidecar_returns_none(tmp_path):
    plain = tmp_path / "plain.vcf.gz"
    write_bgzf(plain, make_vcf(3))
    assert _load_page_index(FilePath(plain)) is None
