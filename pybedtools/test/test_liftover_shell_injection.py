"""Security tests for BedTool.liftover.

The helper builds one command string and passes it to the system shell.
These tests call that helper with the inputs a caller can supply and require
that none of those inputs are executed as shell syntax.

They are expected to fail while the shell string is still used. A fix that
runs liftOver from an argument list, without the shell, should pass even
when liftOver itself is not installed.
"""

import os
import tempfile
from pathlib import Path

import pybedtools
from pybedtools import BedTool
from pybedtools import settings

INTERVAL = "chr1\t0\t10\n"


def _watch(path: Path, expected: str):
    return path, expected


def _assert_shell_did_not_run(marker: Path, watched):
    problems = []
    if marker.exists():
        problems.append(
            "created %s owned by uid %s (process uid %s)"
            % (marker, marker.stat().st_uid, os.getuid())
        )
    for path, expected in watched:
        if not path.exists():
            problems.append("removed %s" % path)
        else:
            actual = path.read_text()
            if actual != expected:
                problems.append("changed %s from %r to %r" % (path, expected, actual))
    assert not problems, "liftover passed input to the system shell: " + "; ".join(
        problems
    )


def _call_liftover(bedtool, **kwargs):
    # The helper may raise after the shell has already run, for example when
    # it tries to open an output file the shell removed or when liftOver is
    # missing. The assertions are about side effects, not the exception.
    try:
        bedtool.liftover(**kwargs)
    except Exception:
        return


def _write_interval(path: Path) -> Path:
    path.write_text(INTERVAL)
    return path


def _snippet(source_name: str, marker_name: str) -> str:
    # Relative names, so this snippet can also be placed inside a filename.
    # `cp` reads the interval data and `printf` changes it. `#` is added by
    # the caller when later shell words must not be parsed as arguments.
    return "cp %s %s; printf Z >> %s" % (source_name, marker_name, source_name)


def test_chainfile_semicolon_reads_and_changes_interval_file(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    chain = "missing.chain; %s; #" % _snippet(str(source), str(marker))
    _call_liftover(BedTool(str(source)), chainfile=chain)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_chainfile_newline_is_not_a_command_separator(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    chain = "missing.chain\n%s\n" % _snippet(str(source), str(marker))
    _call_liftover(BedTool(str(source)), chainfile=chain)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_chainfile_command_substitution_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    chain = "$(%s)" % _snippet(str(source), str(marker))
    _call_liftover(BedTool(str(source)), chainfile=chain)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_chainfile_backticks_are_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    chain = "`%s`" % _snippet(str(source), str(marker))
    _call_liftover(BedTool(str(source)), chainfile=chain)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_existing_chainfile_name_is_not_executed(tmp_path):
    """An existence check on the chain path must not leave the path shell-parsed."""
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    chain = tmp_path / ("chain; %s; #" % _snippet("intervals.bed", "marker"))
    chain.write_text("not a chain file\n")
    previous = os.getcwd()
    os.chdir(tmp_path)
    try:
        _call_liftover(BedTool(str(source)), chainfile=str(chain))
    finally:
        os.chdir(previous)
    _assert_shell_did_not_run(
        marker,
        [_watch(source, INTERVAL), _watch(chain, "not a chain file\n")],
    )


def test_liftover_args_semicolon_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    extra = "; %s; #" % _snippet(str(source), str(marker))
    _call_liftover(BedTool(str(source)), chainfile="missing.chain", liftover_args=extra)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_liftover_args_command_substitution_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    extra = "$(%s)" % _snippet(str(source), str(marker))
    _call_liftover(BedTool(str(source)), chainfile="missing.chain", liftover_args=extra)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_unmapped_path_semicolon_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    unmapped = "%s; %s; #" % (tmp_path / "unmapped.txt", _snippet(str(source), str(marker)))
    _call_liftover(
        BedTool(str(source)),
        chainfile="missing.chain",
        unmapped=unmapped,
    )
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_interval_filename_semicolon_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    nasty = tmp_path / ("a; %s; #.bed" % _snippet("intervals.bed", "marker"))
    nasty.write_text(INTERVAL)
    previous = os.getcwd()
    os.chdir(tmp_path)
    try:
        _call_liftover(BedTool(str(nasty)), chainfile="missing.chain")
    finally:
        os.chdir(previous)
    _assert_shell_did_not_run(
        marker,
        [_watch(source, INTERVAL), _watch(nasty, INTERVAL)],
    )


def test_generated_output_suffix_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    previous_cwd = os.getcwd()
    previous_suffix = settings.tempfile_suffix
    os.chdir(tmp_path)
    settings.tempfile_suffix = "; %s; #" % _snippet("intervals.bed", "marker")
    try:
        _call_liftover(BedTool(str(source)), chainfile="missing.chain")
    finally:
        settings.tempfile_suffix = previous_suffix
        os.chdir(previous_cwd)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_generated_output_prefix_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    previous_cwd = os.getcwd()
    previous_prefix = settings.tempfile_prefix
    os.chdir(tmp_path)
    settings.tempfile_prefix = "p; %s; #" % _snippet("intervals.bed", "marker")
    try:
        _call_liftover(BedTool(str(source)), chainfile="missing.chain")
    finally:
        settings.tempfile_prefix = previous_prefix
        os.chdir(previous_cwd)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])


def test_tempdir_path_is_not_executed(tmp_path):
    source = _write_interval(tmp_path / "intervals.bed")
    marker = tmp_path / "marker"
    scratch = tmp_path / ("scratch; %s; #" % _snippet("intervals.bed", "marker"))
    scratch.mkdir()
    previous_cwd = os.getcwd()
    previous_tempdir = tempfile.tempdir
    os.chdir(tmp_path)
    try:
        pybedtools.set_tempdir(str(scratch))
        _call_liftover(BedTool(str(source)), chainfile="missing.chain")
    finally:
        tempfile.tempdir = previous_tempdir
        os.chdir(previous_cwd)
    _assert_shell_did_not_run(marker, [_watch(source, INTERVAL)])
