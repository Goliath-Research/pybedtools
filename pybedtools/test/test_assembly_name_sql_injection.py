"""Assembly names must not be interpolated into the chromosome-size SQL.

These tests describe the safe behavior. They fail while an assembly name can
change the statement executed by the mysql client.
"""

import subprocess

import pytest

import pybedtools
from pybedtools import filenames, helpers
from pybedtools.contrib import bigbed, bigwig
from pybedtools.contrib.intersection_matrix import IntersectionMatrix
from pybedtools.parallel import parallel_apply
from pybedtools.scripts import intersection_matrix

# Not a credential. Marker row stored in a table the genome account can read.
LEAKED_CHROM = "LEAKED"
LEAKED_SIZE = 4242
LOCAL_HOST = "127.0.0.1"
MISSING_FETCH = "/no/such/fetchChromSizes"
GENERAL_LOG = "/tmp/mariadb-general.log"

# Closes the intended "ASSEMBLY.chromInfo" reference, unions another table,
# and comments out the trailing ".chromInfo" suffix.
INJECTED_ASSEMBLY = (
    "zzsafe.chromInfo UNION SELECT chrom, size FROM otherdata.notes #"
)


def _log_text():
    return subprocess.check_output(
        ["sudo", "cat", GENERAL_LOG], text=True, errors="replace"
    )


def _assert_rejected(call):
    """The call must refuse the name before mysql runs and before any file write."""
    log_before = _log_text()
    temp_before = len(filenames.TEMPFILES)
    try:
        result = call()
    except ValueError as exc:
        message = str(exc).lower()
        assert "invalid assembly name" in message
        assert "internet" not in message
        new_log = _log_text()[len(log_before) :]
        assert "otherdata" not in new_log
        assert "UNION" not in new_log.upper()
        for path in filenames.TEMPFILES[temp_before:]:
            with open(path) as handle:
                assert LEAKED_CHROM not in handle.read()
        return
    raise AssertionError(
        "assembly name reached the chromosome-size query; got %r" % (result,)
    )


@pytest.fixture(scope="module", autouse=True)
def local_genome_database():
    """Database the product's mysql client can read as user genome, with no password."""
    subprocess.check_call(
        [
            "sudo",
            "mysql",
            "--socket=/run/mysqld/mysqld.sock",
            "-e",
            """
            SET GLOBAL general_log_file='/tmp/mariadb-general.log';
            SET GLOBAL log_output='FILE';
            SET GLOBAL general_log='ON';
            CREATE DATABASE IF NOT EXISTS zzsafe;
            CREATE TABLE IF NOT EXISTS zzsafe.chromInfo (
                chrom VARCHAR(255) NOT NULL,
                size BIGINT NOT NULL
            );
            DELETE FROM zzsafe.chromInfo;
            INSERT INTO zzsafe.chromInfo VALUES ('chr1', 1000), ('chr2', 2000);
            CREATE DATABASE IF NOT EXISTS otherdata;
            CREATE TABLE IF NOT EXISTS otherdata.notes (
                chrom VARCHAR(255) NOT NULL,
                size BIGINT NOT NULL
            );
            DELETE FROM otherdata.notes;
            INSERT INTO otherdata.notes VALUES ('%s', %d);
            CREATE USER IF NOT EXISTS 'genome'@'127.0.0.1';
            GRANT SELECT ON zzsafe.* TO 'genome'@'127.0.0.1';
            GRANT SELECT ON otherdata.* TO 'genome'@'127.0.0.1';
            FLUSH PRIVILEGES;
            """
            % (LEAKED_CHROM, LEAKED_SIZE),
        ]
    )
    original = helpers.get_chromsizes_from_ucsc.__defaults__
    # chromsizes() does not take a host. Point its default at this local server
    # so the real function body runs without contacting the public UCSC host.
    helpers.get_chromsizes_from_ucsc.__defaults__ = original[:-1] + (LOCAL_HOST,)
    yield
    helpers.get_chromsizes_from_ucsc.__defaults__ = original


def test_registered_safe_name_still_returns_database_lengths():
    assert "zzsafe" not in dir(pybedtools.genome_registry)
    sizes = helpers.get_chromsizes_from_ucsc(
        "zzsafe",
        host_url=LOCAL_HOST,
        fetchchromsizes=MISSING_FETCH,
        timeout=5,
    )
    assert sizes == {"chr1": (0, 1000), "chr2": (0, 2000)}
    assert helpers.chromsizes("zzsafe") == sizes


def test_get_chromsizes_from_ucsc_rejects_injected_assembly(tmp_path):
    saveas = tmp_path / "sizes.txt"

    def call():
        return helpers.get_chromsizes_from_ucsc(
            INJECTED_ASSEMBLY,
            saveas=str(saveas),
            host_url=LOCAL_HOST,
            fetchchromsizes=MISSING_FETCH,
            timeout=5,
        )

    _assert_rejected(call)
    assert not saveas.exists()


def test_chromsizes_helper_rejects_unregistered_injected_assembly():
    assert not hasattr(pybedtools.genome_registry, INJECTED_ASSEMBLY)

    def call():
        return helpers.chromsizes(INJECTED_ASSEMBLY)

    _assert_rejected(call)


def test_chromsizes_to_file_rejects_injected_assembly(tmp_path):
    saveas = tmp_path / "from-string.genome"

    def call():
        return helpers.chromsizes_to_file(INJECTED_ASSEMBLY, fn=str(saveas))

    _assert_rejected(call)
    assert not saveas.exists()


def _bed():
    return pybedtools.BedTool("chr1\t10\t20\n", from_string=True)


def test_set_chromsizes_rejects_injected_assembly():
    bed = _bed()

    def call():
        return bed.set_chromsizes(INJECTED_ASSEMBLY)

    _assert_rejected(call)
    assert LEAKED_CHROM not in getattr(bed, "chromsizes", {})


def test_truncate_to_chrom_rejects_injected_assembly():
    _assert_rejected(lambda: _bed().truncate_to_chrom(INJECTED_ASSEMBLY))


@pytest.mark.parametrize(
    "operation",
    [
        "complement",
        "shuffle",
        "slop",
        "flank",
        "shift",
        "sort",
        "genome_coverage",
        "random",
        "window_maker",
        "to_bam",
        "bedpe_to_bam",
    ],
)
def test_genome_operation_rejects_injected_assembly(operation):
    kwargs = {"genome": INJECTED_ASSEMBLY}
    if operation == "random":
        target = pybedtools.BedTool()
    else:
        target = _bed()

    def call():
        return getattr(target, operation)(**kwargs)

    _assert_rejected(call)


def test_multi_intersect_rejects_injected_assembly_when_genome_is_required():
    _assert_rejected(
        lambda: pybedtools.BedTool().multi_intersect(
            i=[_bed().fn], empty=True, genome=INJECTED_ASSEMBLY
        )
    )


def test_fisher_rejects_injected_assembly():
    left = _bed()
    right = _bed()
    _assert_rejected(lambda: left.fisher(right, genome=INJECTED_ASSEMBLY))


def test_relative_distance_rejects_injected_assembly():
    _assert_rejected(
        lambda: list(_bed().relative_distance(_bed(), genome=INJECTED_ASSEMBLY))
    )


def test_parallel_apply_rejects_injected_assembly():
    def call():
        return list(
            parallel_apply(
                _bed(),
                method="intersect",
                genome=INJECTED_ASSEMBLY,
                method_args=(_bed(),),
                method_kwargs={"u": True},
                iterations=1,
                processes=1,
                debug=True,
                reduce_func=lambda result: 0,
            )
        )

    _assert_rejected(call)


def test_intersection_matrix_script_rejects_injected_assembly(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        [
            "intersection_matrix",
            "--enrichment",
            "--genome",
            INJECTED_ASSEMBLY,
            "--iterations",
            "1",
            pybedtools.example_filename("a.bed"),
        ],
    )
    _assert_rejected(intersection_matrix.main)


def test_contrib_intersection_matrix_rejects_injected_assembly():
    example = pybedtools.example_filename("a.bed")
    matrix = IntersectionMatrix([example], INJECTED_ASSEMBLY, iterations=1)
    _assert_rejected(lambda: matrix.run_and_insert(example, example))


def test_bigbed_rejects_injected_assembly(tmp_path):
    _assert_rejected(
        lambda: bigbed.bigbed(_bed(), INJECTED_ASSEMBLY, str(tmp_path / "out.bb"))
    )


@pytest.mark.parametrize(
    "fn",
    [bigwig.bedgraph_to_bigwig, bigwig.wig_to_bigwig, bigwig.bam_to_bigwig],
)
def test_bigwig_helpers_reject_injected_assembly(fn, tmp_path):
    _assert_rejected(lambda: fn(_bed().fn, INJECTED_ASSEMBLY, str(tmp_path / "out.bw")))
