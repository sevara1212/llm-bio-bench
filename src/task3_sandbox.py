"""Task 3 sandbox: run agent-written Python in an isolated folder that holds only the raw PBMC3k counts.

Isolation (macOS sandbox-exec / Seatbelt profile, deny-by-default):
  - read access only to the Python installation, the project venv (packages) and the run's own folder;
    the rest of the project (answer key, results, .env), the home folder and all other files are denied
  - write access only inside the run folder (plus /dev/null)
  - all network access denied
  - a fresh process per execution, started with a minimal environment (no API keys, HOME = run folder)
  - 60-second timeout per execution; the whole process group is killed on timeout
Each execution's output (stdout + stderr) is truncated to ~50 lines / 4,000 characters before the agent
sees it; the untruncated output (capped at 20,000 characters) is kept for the trajectory.
"""
import os
import shutil
import signal
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

PROJECT = Path(__file__).resolve().parents[1]
VENV = PROJECT / "venv"
PYTHON = VENV / "bin" / "python"
PY_BASE = Path(sys.base_prefix).resolve()  # the Homebrew Python framework the venv is built on
TIMEOUT_S = 60
MAX_LINES, MAX_CHARS, KEEP_CHARS = 50, 4000, 20000
RAW = PROJECT / "data" / "pbmc3k_raw.h5ad"  # scanpy's cached copy of 10x filtered_gene_bc_matrices


def profile(work: Path) -> str:
    """Seatbelt profile: deny everything, then allow only what Python needs plus the run folder."""
    # /opt/homebrew holds only installed software (Python and the C libraries it links: OpenSSL, xz, ...).
    reads = [str(work), str(VENV), str(PY_BASE), "/opt/homebrew", "/usr/lib", "/usr/share", "/System",
             "/Library/Apple", "/private/var/db/dyld", "/dev"]
    read_rules = "\n".join(f'  (subpath "{p}")' for p in reads)
    return f"""(version 1)
(deny default)
(allow process-fork)
(allow process-exec (subpath "{PY_BASE}") (subpath "{VENV}") (subpath "/opt/homebrew/Cellar/python@3.14"))
(allow sysctl-read)
(allow mach-lookup)
(allow ipc-posix-shm)
(allow signal (target self))
(allow file-read-metadata)
(allow file-read*
{read_rules}
  (literal "/")
  (literal "/private")
  (literal "/private/etc/localtime"))
(allow file-write* (subpath "{work}") (literal "/dev/null"))
(deny network*)
"""


def write_raw_10x(dest: Path):
    """Write the raw PBMC3k counts as a 10x filtered_gene_bc_matrices/hg19 folder (matrix.mtx,
    genes.tsv, barcodes.tsv) - the same input scanpy's pbmc3k() tutorial data starts from."""
    import anndata as ad
    a = ad.read_h5ad(RAW)
    out = dest / "filtered_gene_bc_matrices" / "hg19"
    out.mkdir(parents=True)
    scipy.io.mmwrite(out / "matrix.mtx", sp.csr_matrix(a.X).T.astype(np.int32), field="integer")
    pd.DataFrame({"id": a.var.gene_ids.values, "symbol": a.var_names}).to_csv(out / "genes.tsv", sep="\t", header=False, index=False)
    pd.Series(a.obs_names).to_csv(out / "barcodes.tsv", header=False, index=False)


# Work folders live OUTSIDE the project: some libraries (anndata's settings via pydantic-settings) search
# parent folders for a .env file, and inside the project they would find the API-key file. The sandbox
# blocks that read anyway, but it crashes the import, so the folders are kept where no parent has a .env.
WORK_ROOT = Path("/private/tmp/llmbio_task3_runs")


def new_run_folder(run_id: str) -> Path:
    """Fresh work folder (outside the project) containing only the raw counts."""
    work = WORK_ROOT / run_id / "work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_raw_10x(work)
    return work


def truncate(text: str) -> str:
    lines = text.splitlines()
    if len(lines) > MAX_LINES:
        head, tail = lines[:20], lines[-(MAX_LINES - 20):]
        lines = head + [f"... [{len(lines) - MAX_LINES} lines truncated] ..."] + tail
    out = "\n".join(lines)
    if len(out) > MAX_CHARS:
        out = out[:1500] + f"\n... [{len(out) - MAX_CHARS} characters truncated] ...\n" + out[-(MAX_CHARS - 1500):]
    return out


def run_python(work: Path, code: str, n: int) -> dict:
    """Run one code cell in the sandbox. Returns output (full + truncated), exit code, timing."""
    work = work.resolve()
    cell = work / f"_cell_{n:02d}.py"
    cell.write_text(code)
    tmp = work / ".tmp"
    tmp.mkdir(exist_ok=True)
    env = {"HOME": str(work), "TMPDIR": str(tmp), "MPLCONFIGDIR": str(tmp), "NUMBA_CACHE_DIR": str(tmp),
           "MPLBACKEND": "Agg", "PATH": f"{VENV}/bin:/usr/bin:/bin", "PYTHONDONTWRITEBYTECODE": "1",
           "OMP_NUM_THREADS": "4", "LANG": "en_US.UTF-8"}
    cmd = ["/usr/bin/sandbox-exec", "-p", profile(work), str(PYTHON), cell.name]
    proc = subprocess.Popen(cmd, cwd=work, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, start_new_session=True)
    try:
        out, _ = proc.communicate(timeout=TIMEOUT_S)
        timed_out = False
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        out, _ = proc.communicate()
        timed_out = True
    out = (out or "") + (f"\n[TIMEOUT: execution stopped after {TIMEOUT_S} s]" if timed_out else "")
    error = timed_out or proc.returncode != 0 or "Traceback (most recent call last)" in out
    return {"code": code, "output_full": out[:KEEP_CHARS], "output_shown": truncate(out) or "(no output)",
            "returncode": proc.returncode, "timed_out": timed_out, "error": error}
