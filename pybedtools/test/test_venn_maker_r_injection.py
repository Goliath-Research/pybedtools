"""
Regression tests for R code injection in pybedtools.contrib.venn_maker.

Interval coordinates, set names, and the figure filename are data. The helper
places them in double-quoted R literals. A quote or backslash in those values
must not change the syntax of the generated script, and running that script
must not execute text that arrived in those values.

`additional_args` is documented as verbatim R and is not part of this contract.
The two-, three-, and four-set overlap helpers all feed the same literal
builder; each count is covered because a chromosome has to survive that
cleanup before it is copied into the script.
"""
import os
import subprocess
import tempfile

import pytest

import pybedtools
from pybedtools.contrib.venn_maker import venn_maker


_COLLECT_R = r"""
collect <- function(x) {
  if (is.character(x)) return(x)
  if (is.symbol(x) || is.null(x) || is.atomic(x)) return(character())
  out <- character()
  if (is.call(x) || is.expression(x) || is.pairlist(x)) {
    nms <- names(x)
    if (!is.null(nms)) out <- c(out, nms[nzchar(nms)])
    for (el in as.list(x)) out <- c(out, collect(el))
    return(out)
  }
  character()
}
args <- commandArgs(trailingOnly = TRUE)
exprs <- parse(args[[1]])
vals <- collect(exprs)
con <- file(args[[2]], "wb")
writeLines(as.character(length(vals)), con)
for (s in vals) {
  b <- charToRaw(enc2utf8(s))
  writeLines(paste(as.integer(b), collapse = ","), con)
}
close(con)
"""


def _write_bed(directory, name, chrom, start=10, stop=20):
    path = os.path.join(directory, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("%s\t%d\t%d\n" % (chrom, start, stop))
    return path


def _r_string_constants(script_path):
    """Parse `script_path` in R and return its character constants.

    The generated script is not evaluated. A syntax error is reported to the
    caller instead of being treated as an empty program.
    """
    collector = script_path + ".collect.R"
    output = script_path + ".strings"
    with open(collector, "w", encoding="utf-8") as handle:
        handle.write(_COLLECT_R)
    completed = subprocess.run(
        ["Rscript", collector, script_path, output],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise AssertionError("generated script is not valid R:\n" + detail)
    with open(output, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    count = int(lines[0])
    values = []
    for line in lines[1:]:
        if line == "":
            values.append("")
            continue
        values.append(bytes(int(part) for part in line.split(",")).decode("utf-8"))
    if len(values) != count:
        raise AssertionError("R string collector returned a truncated list")
    return values


def _run_venn(directory, beds, names, figure_filename):
    script_path = os.path.join(directory, "venn.R")
    try:
        venn_maker(
            beds,
            names=names,
            figure_filename=figure_filename,
            script_filename=script_path,
            run=True,
        )
    except ValueError as exc:
        if "install R" not in str(exc):
            raise
    return script_path


def _assert_literals_are_data(script_path, marker_path, expected_strings):
    problems = []
    if os.path.exists(marker_path):
        problems.append(
            "running the generated script wrote %s from text that was supplied as data"
            % marker_path
        )
    try:
        constants = _r_string_constants(script_path)
    except AssertionError as exc:
        problems.append(str(exc))
        constants = []
    for expected in expected_strings:
        if expected not in constants:
            problems.append(
                "expected this data to be one R string constant:\n%r" % (expected,)
            )
    assert not problems, "\n".join(problems)


def _other_chrom(index):
    return "chr%d" % (index + 2)


@pytest.mark.parametrize("set_count", [2, 3, 4])
def test_chromosome_text_stays_inside_an_r_string(set_count):
    with tempfile.TemporaryDirectory() as directory:
        marker = os.path.join(directory, "chromosome-marker")
        chrom = 'chr1", writeLines("ran", "%s"), "' % marker
        beds = [_write_bed(directory, "set-0.bed", chrom, start=10, stop=20)]
        for index in range(1, set_count):
            beds.append(
                _write_bed(
                    directory,
                    "set-%d.bed" % index,
                    _other_chrom(index),
                    start=100 + index,
                    stop=110 + index,
                )
            )
        names = ["set%d" % index for index in range(set_count)]
        figure = os.path.join(directory, "figure.tiff")
        script = _run_venn(directory, beds, names, figure)
        expected = ['%s|10|20\n' % chrom]
        expected.extend(
            "%s|%d|%d\n" % (_other_chrom(index), 100 + index, 110 + index)
            for index in range(1, set_count)
        )
        _assert_literals_are_data(script, marker, expected)
        pybedtools.cleanup()


def test_set_name_stays_inside_an_r_string():
    with tempfile.TemporaryDirectory() as directory:
        marker = os.path.join(directory, "name-marker")
        malicious = 'alpha", writeLines("ran", "%s"), "beta' % marker
        beds = [
            _write_bed(directory, "a.bed", "chr1", start=10, stop=20),
            _write_bed(directory, "b.bed", "chr2", start=30, stop=40),
        ]
        figure = os.path.join(directory, "figure.tiff")
        script = _run_venn(directory, beds, ["safe", malicious], figure)
        _assert_literals_are_data(
            script,
            marker,
            ["safe", malicious, "chr1|10|20\n", "chr2|30|40\n", figure],
        )
        pybedtools.cleanup()


def test_figure_filename_stays_inside_an_r_string():
    with tempfile.TemporaryDirectory() as directory:
        marker = os.path.join(directory, "filename-marker")
        figure = os.path.join(
            directory, 'figure.tiff"); writeLines("ran", "%s"); "' % marker
        )
        beds = [
            _write_bed(directory, "a.bed", "chr1", start=10, stop=20),
            _write_bed(directory, "b.bed", "chr2", start=30, stop=40),
        ]
        script = _run_venn(directory, beds, ["alpha", "beta"], figure)
        _assert_literals_are_data(
            script,
            marker,
            [figure, "alpha", "beta", "chr1|10|20\n", "chr2|30|40\n"],
        )
        pybedtools.cleanup()


@pytest.mark.parametrize(
    "chrom_prefix",
    [
        "chr1\\",
        "chr1\\\\",
    ],
)
def test_backslash_in_chromosome_text_stays_inside_an_r_string(chrom_prefix):
    with tempfile.TemporaryDirectory() as directory:
        marker = os.path.join(directory, "backslash-marker")
        chrom = '%s", writeLines("ran", "%s"), "' % (chrom_prefix, marker)
        beds = [
            _write_bed(directory, "a.bed", chrom, start=10, stop=20),
            _write_bed(directory, "b.bed", "chr2", start=30, stop=40),
        ]
        figure = os.path.join(directory, "figure.tiff")
        script = _run_venn(directory, beds, ["alpha", "beta"], figure)
        _assert_literals_are_data(
            script,
            marker,
            ["%s|10|20\n" % chrom, "chr2|30|40\n", figure, "alpha", "beta"],
        )
        pybedtools.cleanup()
