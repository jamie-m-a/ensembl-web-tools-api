"""Minimal, dependency-free BGZF reader supporting seek by virtual offset.

Extracted from vcf_results. Used by the page-index seek path (see
vcf_results._read_indexed_page): when the pipeline emits a `<vcf>.pageidx.json`
sidecar, a page can be fetched by seeking straight to a record's packed BGZF
virtual offset (compressed_block_offset << 16 | within_block_offset) instead of
scanning from the top of the file.
"""

import struct
import zlib


def _bc_block_size(extra: bytes) -> int | None:
    """BSIZE from a gzip member's FEXTRA field, or None if no 'BC' subfield."""
    i = 0
    while i + 4 <= len(extra):
        si1, si2, slen = extra[i], extra[i + 1], struct.unpack("<H", extra[i + 2 : i + 4])[0]
        if si1 == 0x42 and si2 == 0x43:  # 'B', 'C'
            return struct.unpack("<H", extra[i + 4 : i + 6])[0]
        i += 4 + slen
    return None


class _BgzfReader:
    """Minimal, dependency-free BGZF reader supporting seek by virtual offset.

    Only the operations the page-index seek needs are implemented: `seek` to a
    packed virtual offset, line-oriented `readline` (stitching lines that cross
    block boundaries), and `tell` (the virtual offset of the next byte).
    """

    def __init__(self, path: str):
        self._file = open(path, "rb")
        self._block_coffset = 0  # compressed offset of the current block
        self._data = b""  # decompressed bytes of the current block
        self._pos = 0  # position within the current block
        self._load_block_at(0)

    def _read_block(self):
        """Read/decompress the block at the current file position -> (coffset,
        data); data is None at EOF, b'' for the empty BGZF EOF-marker block."""
        coffset = self._file.tell()
        header = self._file.read(12)
        if not header:
            return coffset, None
        if len(header) < 12 or header[:2] != b"\x1f\x8b":
            raise ValueError("not a BGZF block")
        xlen = struct.unpack("<H", header[10:12])[0]
        extra = self._file.read(xlen)
        bsize = _bc_block_size(extra)
        if bsize is None:
            raise ValueError("missing BGZF 'BC' subfield")
        cdata = self._file.read((bsize + 1) - 12 - xlen - 8)
        self._file.read(8)  # trailer (CRC32 + ISIZE)
        return coffset, (zlib.decompress(cdata, -15) if cdata else b"")

    def _load_block_at(self, coffset: int) -> None:
        self._file.seek(coffset)
        self._block_coffset, data = self._read_block()
        self._data = data or b""
        self._pos = 0

    def seek(self, voffset: int) -> None:
        self._load_block_at(voffset >> 16)
        self._pos = voffset & 0xFFFF

    def tell(self) -> int:
        return (self._block_coffset << 16) | self._pos

    def _advance_block(self) -> bool:
        """Load the next non-empty block; False at EOF (empty EOF-marker blocks
        are skipped)."""
        while True:
            coffset, data = self._read_block()
            if data is None:
                return False
            self._block_coffset, self._data, self._pos = coffset, data, 0
            if data:
                return True

    def readline(self) -> bytes:
        chunks: list[bytes] = []
        while True:
            if self._pos >= len(self._data):
                if not self._advance_block():
                    break  # EOF
            newline = self._data.find(b"\n", self._pos)
            if newline == -1:
                chunks.append(self._data[self._pos :])
                self._pos = len(self._data)
                continue  # line continues into the next block
            chunks.append(self._data[self._pos : newline + 1])
            self._pos = newline + 1
            break
        # Normalise an end-of-block position to the next block's start, so tell()
        # reports the next line's virtual offset the way the generator recorded
        # it (a line beginning a block has within-block offset 0).
        if self._pos >= len(self._data):
            self._advance_block()
        return b"".join(chunks)

    def close(self) -> None:
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
