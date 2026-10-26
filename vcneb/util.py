"""Small shared helpers: json io, timing, subprocess, formatting."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime


class VcnebError(RuntimeError):
    """Error with an actionable message for the user."""


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def _json_default(obj):
    try:
        import numpy as np

        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
    except Exception:
        pass
    return str(obj)


def jdump(obj, path: str, indent: int = 1) -> None:
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=indent, default=_json_default)
    os.replace(tmp, path)


def jload(path: str, default=None):
    if not os.path.isfile(path):
        return default
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def now_str() -> str:
    return datetime.now().strftime("%H:%M:%S")


def fmt_duration(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h:d}h{m:02d}m"
    if m:
        return f"{m:d}m{s:02d}s"
    return f"{s:d}s"


class Timer:
    def __init__(self):
        self.t0 = time.time()

    def __call__(self) -> float:
        return time.time() - self.t0


def tail(path: str, n: int = 25) -> str:
    if not os.path.isfile(path):
        return ""
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()
    except OSError:
        return ""
    return "".join(lines[-n:])


def run_cmd(cmd, cwd: str, log_path: str | None = None, timeout: float | None = None,
            env: dict | None = None, echo: bool = False) -> dict:
    """Run a command; return dict(rc, stdout, stderr, seconds, cmd)."""
    if isinstance(cmd, str):
        cmd = cmd.split()
    cmd = [str(c) for c in cmd]
    exe = cmd[0]
    if shutil.which(exe) is None and not os.path.isabs(exe):
        raise VcnebError(
            f"executable '{exe}' not found in PATH.\n"
            f"  command was: {' '.join(cmd)}\n"
            f"  fix: set calc.<engine>.command in the config, or add {exe} to PATH."
        )
    env = dict(os.environ if env is None else env)
    env.setdefault("OMP_NUM_THREADS", "1")
    t0 = time.time()
    if echo:
        print(f"[cmd] {' '.join(cmd)}  (cwd={cwd})")
    try:
        proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                              timeout=timeout, env=env)
        out, err, rc = proc.stdout, proc.stderr, proc.returncode
    except subprocess.TimeoutExpired as exc:
        out = (exc.stdout or "") if isinstance(exc.stdout, str) else ""
        err = (exc.stderr or "") if isinstance(exc.stderr, str) else ""
        rc = -9
        err += f"\n[vcneb] TIMEOUT after {timeout:.0f}s\n"
    seconds = time.time() - t0
    if log_path:
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(f"# {' '.join(cmd)}\n# cwd: {cwd}\n# rc: {rc}   time: {seconds:.1f}s\n")
            fh.write("=" * 70 + "\n[stdout]\n")
            fh.write(out)
            fh.write("\n[stderr]\n")
            fh.write(err)
    return {"rc": rc, "stdout": out, "stderr": err, "seconds": seconds, "cmd": cmd}


def print_exc_brief(exc: BaseException) -> str:
    return "".join(traceback.format_exception_only(type(exc), exc)).strip()


class Logger:
    """Tiny logger that mirrors everything to a file."""

    def __init__(self, path: str | None = None):
        self.path = path
        self.t0 = time.time()
        if path:
            ensure_dir(os.path.dirname(path) or ".")
            open(path, "a").close()

    def __call__(self, msg: str = "", level: str = "info") -> None:
        prefix = {"info": "", "warn": "WARNING: ", "error": "ERROR: ", "ok": ""}.get(level, "")
        line = f"{prefix}{msg}"
        print(line, flush=True)
        if self.path:
            with open(self.path, "a", encoding="utf-8") as fh:
                fh.write(f"[{now_str()} +{self.t0 and time.time() - self.t0:6.1f}s] {line}\n")

    def warn(self, msg: str) -> None:
        self(msg, "warn")

    def error(self, msg: str) -> None:
        self(msg, "error")

    def rule(self, char: str = "-", width: int = 78) -> None:
        self(char * width)


def is_windows() -> bool:
    return sys.platform.startswith("win")
