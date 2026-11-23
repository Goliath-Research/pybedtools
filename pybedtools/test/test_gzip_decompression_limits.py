"""Security contract for gzip interval files.

These tests fail on purpose while decompression is unbounded. They pass once
gzip interval input is rejected with ``pybedtools.helpers.pybedtoolsError``
when any of the ceilings below is exceeded, and still accept gzip interval
files that stay within every ceiling.

Ceilings (inclusive; one byte past is an error):

* ``MAX_GZIP_COMPRESSED_BYTES`` — size of the file on disk
* ``MAX_GZIP_DECOMPRESSED_BYTES`` — total decompressed bytes of the member
* ``MAX_GZIP_LINE_BYTES`` — one decompressed line, including its newline

The error must be raised before a caller-supplied output file is created and
before a new temp file is registered. The message must mention gzip,
decompression, a limit, or a size.

Operations that read through a different opener are covered separately:
``BedTool`` iteration (``count``, ``len``, ``str``, equality, indexing,
``file_type``), ``head``, ``saveas``, ``cat``, and the C++ ``IntervalFile``
reader (``intervals``, direct iteration, ``all_hits``).

Not given their own fixtures, because they call the same ``BedTool`` iterator
once that iterator enforces the ceilings: ``cut``, ``each``, ``filter``,
``random_subset``, ``to_dataframe``, and the other Python-side feature loops.
``isBAM`` reads at most four decompressed bytes and is not a full expansion.
Subprocess calls to an external bedtools binary are outside this library's
reader.
"""

import gzip
import os

import pytest

import pybedtools
from pybedtools.helpers import pybedtoolsError

# Inclusive ceilings. A gzip interval file at exactly these sizes is valid.
MAX_GZIP_COMPRESSED_BYTES = 128 * 1024
MAX_GZIP_DECOMPRESSED_BYTES = 256 * 1024
MAX_GZIP_LINE_BYTES = 64 * 1024

_LIMIT_MATCH = r"(?i)(gzip|decompress|limit|exceed|maximum|too (long|large)|size)"

_BED_PREFIX = b"chr1\t0\t1\t"


def teardown_module():
    pybedtools.cleanup()


def _bed_line(total_len):
    """One valid BED line of exactly ``total_len`` bytes, including newline."""
    overhead = len(_BED_PREFIX) + 1
    if total_len < overhead:
        raise ValueError("line shorter than a minimal BED record: %s" % total_len)
    return _BED_PREFIX + (b"n" * (total_len - overhead)) + b"\n"


def _repeated_lines(total_len, line_len):
    """Valid BED payload of exactly ``total_len`` bytes, each line ``line_len``."""
    if total_len % line_len != 0:
        raise ValueError("total_len must be a multiple of line_len")
    line = _bed_line(line_len)
    return line * (total_len // line_len)


def _write_gzip(directory, name, payload):
    path = directory / name
    path.write_bytes(gzip.compress(payload))
    return path


def _write_over_compressed(directory):
    """A valid BED gzip whose on-disk size is one byte past the compressed cap.

    The first record is a short line, so a reader that only capped line length
    or that only counted bytes it had already decompressed would still return
    that record.
    """
    chunks = []
    blob = b""
    while len(blob) <= MAX_GZIP_COMPRESSED_BYTES:
        token = os.urandom(32).hex().encode("ascii")
        chunks.append(_BED_PREFIX + token + b"\n")
        blob = gzip.compress(b"".join(chunks), compresslevel=1)
    path = directory / "over-compressed.bed.gz"
    path.write_bytes(blob)
    assert path.stat().st_size > MAX_GZIP_COMPRESSED_BYTES
    with gzip.open(path, "rt") as handle:
        first = handle.readline()
    assert len(first.encode("utf-8")) <= MAX_GZIP_LINE_BYTES
    return path


def _must_reject(callback):
    with pytest.raises(pybedtoolsError, match=_LIMIT_MATCH):
        callback()


def _consume(op, path):
    """Run one public operation. Full-stream ops read every record."""
    if op == "intervalfile":
        return next(iter(pybedtools.IntervalFile(str(path))))
    if op == "all_hits":
        query = pybedtools.Interval("chr1", 0, 1)
        return pybedtools.IntervalFile(str(path)).all_hits(query)

    tool = pybedtools.BedTool(str(path))
    if op == "count":
        return tool.count()
    if op == "len":
        return len(tool)
    if op == "str":
        return str(tool)
    if op == "iter_first":
        return next(iter(tool))
    if op == "iter_all":
        return list(tool)
    if op == "getitem":
        return tool[0]
    if op == "eq":
        other = pybedtools.BedTool("chr1\t0\t1\n", from_string=True)
        return tool == other
    if op == "file_type":
        return tool.file_type
    if op == "head":
        return tool.head(n=1)
    if op == "intervals":
        return next(iter(tool.intervals))
    raise AssertionError("unknown op %s" % op)


# Readers that materialize the whole member, so a decompressed-size cap applies.
_FULL_READS = ["count", "len", "str", "eq", "iter_all", "all_hits"]

# Readers that must inspect the first record, so a line cap applies.
_FIRST_RECORD = [
    "count",
    "len",
    "str",
    "eq",
    "iter_first",
    "iter_all",
    "getitem",
    "file_type",
    "head",
    "intervals",
    "intervalfile",
    "all_hits",
]

# File size is known before any record is parsed.
_COMPRESSED_SIZE = _FIRST_RECORD


@pytest.mark.parametrize("op", _FIRST_RECORD)
def test_line_over_ceiling_is_rejected(tmp_path, op):
    path = _write_gzip(tmp_path, "long.bed.gz", _bed_line(MAX_GZIP_LINE_BYTES + 1))
    assert path.stat().st_size <= MAX_GZIP_COMPRESSED_BYTES
    _must_reject(lambda: _consume(op, path))


@pytest.mark.parametrize("op", _FULL_READS)
def test_decompressed_over_ceiling_is_rejected(tmp_path, op):
    line_len = 1024
    payload = _repeated_lines(MAX_GZIP_DECOMPRESSED_BYTES, line_len) + _bed_line(line_len)
    assert len(payload) > MAX_GZIP_DECOMPRESSED_BYTES
    path = _write_gzip(tmp_path, "expanded.bed.gz", payload)
    assert path.stat().st_size <= MAX_GZIP_COMPRESSED_BYTES
    _must_reject(lambda: _consume(op, path))


@pytest.mark.parametrize("op", _COMPRESSED_SIZE)
def test_compressed_over_ceiling_is_rejected(tmp_path, op):
    path = _write_over_compressed(tmp_path)
    _must_reject(lambda: _consume(op, path))


def test_saveas_rejects_long_line_without_creating_output(tmp_path):
    src = _write_gzip(tmp_path, "long.bed.gz", _bed_line(MAX_GZIP_LINE_BYTES + 1))
    dest = tmp_path / "out.bed"
    _must_reject(lambda: pybedtools.BedTool(str(src)).saveas(fn=str(dest)))
    assert not dest.exists()


def test_saveas_rejects_expanded_member_without_creating_output(tmp_path):
    line_len = 1024
    payload = _repeated_lines(MAX_GZIP_DECOMPRESSED_BYTES, line_len) + _bed_line(line_len)
    src = _write_gzip(tmp_path, "expanded.bed.gz", payload)
    dest = tmp_path / "out.bed"
    _must_reject(lambda: pybedtools.BedTool(str(src)).saveas(fn=str(dest)))
    assert not dest.exists()


def test_saveas_rejects_over_compressed_file_without_creating_output(tmp_path):
    src = _write_over_compressed(tmp_path)
    dest = tmp_path / "out.bed"
    _must_reject(lambda: pybedtools.BedTool(str(src)).saveas(fn=str(dest)))
    assert not dest.exists()


def test_head_rejects_long_line_without_writing_it(tmp_path, capsys):
    path = _write_gzip(tmp_path, "long.bed.gz", _bed_line(MAX_GZIP_LINE_BYTES + 1))
    _must_reject(lambda: pybedtools.BedTool(str(path)).head(n=1))
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_cat_rejects_long_line_without_creating_a_temp_file(tmp_path):
    src = _write_gzip(tmp_path, "long.bed.gz", _bed_line(MAX_GZIP_LINE_BYTES + 1))
    other = tmp_path / "other.bed"
    other.write_bytes(b"chr1\t0\t1\n")
    before = set(pybedtools.BedTool.TEMPFILES)
    _must_reject(
        lambda: pybedtools.BedTool(str(src)).cat(str(other), postmerge=False)
    )
    created = [path for path in pybedtools.BedTool.TEMPFILES if path not in before]
    assert created == []


def test_line_at_ceiling_is_accepted(tmp_path):
    payload = _bed_line(MAX_GZIP_LINE_BYTES)
    path = _write_gzip(tmp_path, "exact-line.bed.gz", payload)
    assert path.stat().st_size <= MAX_GZIP_COMPRESSED_BYTES
    tool = pybedtools.BedTool(str(path))
    assert tool.count() == 1
    assert len(tool) == 1
    assert tool.file_type == "bed"
    assert len(next(iter(tool)).name) == MAX_GZIP_LINE_BYTES - len(_BED_PREFIX) - 1
    interval = next(iter(pybedtools.IntervalFile(str(path))))
    assert interval.chrom == "chr1"
    assert interval.start == 0
    assert interval.end == 1


def test_decompressed_at_ceiling_is_accepted(tmp_path):
    line_len = 1024
    payload = _repeated_lines(MAX_GZIP_DECOMPRESSED_BYTES, line_len)
    path = _write_gzip(tmp_path, "exact-expanded.bed.gz", payload)
    assert path.stat().st_size <= MAX_GZIP_COMPRESSED_BYTES
    tool = pybedtools.BedTool(str(path))
    assert tool.count() == MAX_GZIP_DECOMPRESSED_BYTES // line_len
    hits = pybedtools.IntervalFile(str(path)).all_hits(pybedtools.Interval("chr1", 0, 1))
    assert len(hits) == MAX_GZIP_DECOMPRESSED_BYTES // line_len


def test_small_gzip_round_trip_still_works(tmp_path):
    payload = b"chr1\t0\t1\tfeature\nchr1\t2\t3\tother\n"
    path = _write_gzip(tmp_path, "small.bed.gz", payload)
    tool = pybedtools.BedTool(str(path))
    assert tool.count() == 2
    assert "feature" in str(tool)
    dest = tmp_path / "copy.bed"
    copied = tool.saveas(fn=str(dest))
    assert dest.exists()
    assert copied.count() == 2
