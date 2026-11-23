"""Oversized interval lines must not kill the process.

``Tokenize`` (pybedtools/include/lineFileUtilities.h) copies every input line
into a stack array sized from that line's length. ``BedFile::GetNextBed``
reads the line onto the heap and then calls ``Tokenize``. These tests drive
that copy through every public helper that reaches it.

Each case runs in a child whose stack soft limit is ``STACK_BYTES``. The
hostile file is one line of that same length and contains no newline. While
the stack copy is present the child dies with SIGSEGV and these tests fail.
They pass when the parser no longer allocates that array on the stack: either
the line is parsed, or the parser rejects it with a Python exception and
returns control to the caller.

Not covered here, and why:

* ``Tokenize`` for ``vector<int>`` has the same stack copy, but nothing in
  this tree calls it.
* ``BedTool.introns`` reaches ``IntervalFile.search`` only after an external
  ``bedtools sort``. That in-process copy is the ``search`` case below.
* ``Window`` only pulls from the iterable it is given. The copy is
  ``IntervalFile`` iteration, covered below.
* ``same_strand=True`` uses a different overlap function, but only after
  ``loadBedFileIntoMap``, which is the crashing call in the unstranded cases.
"""

import gzip
import os
import subprocess
import sys
import textwrap

import pytest

STACK_BYTES = 2 * 1024 * 1024

CHILD = textwrap.dedent(
    r"""
    import sys
    import resource
    from pybedtools import BedTool, Interval, IntervalFile

    hostile, normal, op = sys.argv[1], sys.argv[2], sys.argv[3]
    soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
    print("stack-soft %s" % soft, flush=True)

    query = Interval("chr1", 0, 10, strand="+")

    def exercise(path):
        if op == "intervalfile-iter":
            for _feature in IntervalFile(path):
                pass
        elif op == "intervalfile-file-type":
            IntervalFile(path).file_type
        elif op == "intervalfile-load":
            IntervalFile(path).loadIntoMap()
        elif op == "intervalfile-all-hits":
            IntervalFile(path).all_hits(query)
        elif op == "intervalfile-search":
            IntervalFile(path).search(query)
        elif op == "intervalfile-any-hits":
            IntervalFile(path).any_hits(query)
        elif op == "intervalfile-count-hits":
            IntervalFile(path).count_hits(query)
        elif op == "bedtool-all-hits":
            BedTool(path).all_hits(query)
        elif op == "bedtool-any-hits":
            BedTool(path).any_hits(query)
        elif op == "bedtool-count-hits":
            BedTool(path).count_hits(query)
        elif op == "bedtool-as-intervalfile":
            for _feature in BedTool(path).as_intervalfile():
                pass
        elif op == "bedtool-intervals":
            for _feature in BedTool(path).intervals:
                pass
        else:
            raise SystemExit("unknown op %s" % op)

    exercise(normal)
    print("normal-ok", flush=True)
    try:
        exercise(hostile)
    except Exception as exc:
        print("rejected %s" % type(exc).__name__, flush=True)
    print("returned", flush=True)
    """
)


def _hostile_line(kind, size):
    if kind == "bed":
        prefix = b"chr1\t1\t2\t"
    elif kind == "gff":
        prefix = b"chr1\t.\texon\t1\t2\t.\t+\t.\t"
    elif kind == "vcf":
        prefix = b"chr1\t1\t.\tA\tC\t.\tPASS\t"
    else:
        raise AssertionError(kind)
    if size <= len(prefix):
        raise AssertionError("stack limit %s is shorter than a %s record" % (size, kind))
    line = prefix + b"N" * (size - len(prefix))
    if b"\n" in line or b"\r" in line:
        raise AssertionError("hostile line contains a newline")
    return line


def _write_hostile(directory, kind, gzipped):
    line = _hostile_line(kind, STACK_BYTES)
    if gzipped:
        path = os.path.join(directory, "hostile.%s.gz" % kind)
        with gzip.open(path, "wb") as handle:
            handle.write(line)
    else:
        suffix = {"bed": "bed", "gff": "gff", "vcf": "vcf"}[kind]
        path = os.path.join(directory, "hostile.%s" % suffix)
        with open(path, "wb") as handle:
            handle.write(line)
        assert os.path.getsize(path) == STACK_BYTES
    return path


def _run(tmp_path, op, kind, gzipped):
    normal = tmp_path / "normal.bed"
    normal.write_bytes(b"chr1\t1\t2\tname\t0\t+\n")
    hostile = _write_hostile(str(tmp_path), kind, gzipped)
    proc = subprocess.run(
        [
            "prlimit",
            "--stack=%s:" % STACK_BYTES,
            "--",
            sys.executable,
            "-c",
            CHILD,
            hostile,
            str(normal),
            op,
        ],
        check=False,
        capture_output=True,
        timeout=60,
    )
    stdout = proc.stdout.decode("utf-8", "replace")
    stderr = proc.stderr.decode("utf-8", "replace")
    assert "stack-soft %s" % STACK_BYTES in stdout, (stdout, stderr, proc.returncode)
    assert proc.returncode == 0, (
        "%s on a %s-byte %s line died or failed (returncode=%s). "
        "A negative code is a fatal signal; 11 is SIGSEGV from the stack copy.\n"
        "stdout:\n%s\nstderr:\n%s"
        % (op, STACK_BYTES, kind, proc.returncode, stdout, stderr)
    )
    assert "normal-ok" in stdout, stdout
    assert "returned" in stdout, stdout


BED_OPS = [
    "intervalfile-iter",
    "intervalfile-file-type",
    "intervalfile-load",
    "intervalfile-all-hits",
    "intervalfile-search",
    "intervalfile-any-hits",
    "intervalfile-count-hits",
    "bedtool-all-hits",
    "bedtool-any-hits",
    "bedtool-count-hits",
    "bedtool-as-intervalfile",
    "bedtool-intervals",
]


@pytest.mark.parametrize("op", BED_OPS)
def test_bed_line_of_stack_length_does_not_crash(tmp_path, op):
    _run(tmp_path, op, "bed", gzipped=False)


@pytest.mark.parametrize("kind", ["gff", "vcf"])
def test_gff_and_vcf_line_of_stack_length_does_not_crash(tmp_path, kind):
    _run(tmp_path, "intervalfile-iter", kind, gzipped=False)


@pytest.mark.parametrize(
    "op",
    [
        "intervalfile-iter",
        "intervalfile-all-hits",
        "intervalfile-any-hits",
        "intervalfile-count-hits",
        "bedtool-all-hits",
        "bedtool-any-hits",
        "bedtool-count-hits",
    ],
)
def test_gzip_line_of_stack_length_does_not_crash(tmp_path, op):
    _run(tmp_path, op, "bed", gzipped=True)
