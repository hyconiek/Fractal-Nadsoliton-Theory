#!/usr/bin/env python3
"""Portable Colab/Kaggle runner for FIN P486/P485/P487 and optional P475.

Each stage writes a log, status JSON, master state, and a checkpoint ZIP.  A
keyboard interrupt is forwarded to the active CAS child and recorded.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import zipfile


STAGES = {
    "P486": "fin_program_486.py",
    "P485": "fin_program_485_unlimited.py",
    "P487": "fin_program_487_unlimited.py",
    "P475": "fin_program_475_unlimited.py",
}


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def checkpoint(root: Path, label: str) -> Path:
    destination = root / f"FIN_Remote_Checkpoint_{label}.zip"
    patterns = ("FIN_Program_*.json", "FIN_Program_*.log", "*.txt.gz",
                "FIN_Remote_Campaign_State.json")
    files = sorted({path for pattern in patterns for path in root.glob(pattern)})
    temporary = destination.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            if path != temporary and path != destination:
                archive.write(path, path.name)
    temporary.replace(destination)
    return destination


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--programs", default="P486,P485,P487")
    parser.add_argument("--root", default=".")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    requested = [value.strip().upper() for value in args.programs.split(",") if value.strip()]
    unknown = [value for value in requested if value not in STAGES]
    if unknown:
        raise SystemExit(f"unknown programs: {unknown}")
    state_path = root / "FIN_Remote_Campaign_State.json"
    state = {
        "campaign": "FIN portable exact algebra v2",
        "programs": requested,
        "status": "running",
        "started_unix_time": time.time(),
        "platform": sys.platform,
        "python": sys.version,
        "stages": {},
    }
    atomic_json(state_path, state)
    environment = dict(os.environ)
    environment["PYTHONUNBUFFERED"] = "1"
    environment.setdefault("MPLCONFIGDIR", str(root / ".matplotlib"))
    for name in requested:
        script = root / STAGES[name]
        log_path = root / f"FIN_Program_{name}_Remote.log"
        stage = {"script": script.name, "status": "running", "started_unix_time": time.time()}
        state["current_stage"] = name
        state["stages"][name] = stage
        atomic_json(state_path, state)
        print(f"\n=== {name} started; log: {log_path.name} ===", flush=True)
        with log_path.open("a", encoding="utf-8") as log:
            process = subprocess.Popen(
                [sys.executable, "-u", str(script)], cwd=root, env=environment,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                bufsize=1, start_new_session=True,
            )
            try:
                assert process.stdout is not None
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    print(line, end="", flush=True)
                return_code = process.wait()
                stage["status"] = "completed" if return_code == 0 else "failed"
                stage["return_code"] = return_code
            except KeyboardInterrupt:
                os.killpg(process.pid, signal.SIGINT)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait(timeout=15)
                stage["status"] = "stopped_by_user"
                stage["return_code"] = process.returncode
                state["status"] = "stopped_by_user"
                state["current_stage"] = None
                stage["finished_unix_time"] = time.time()
                atomic_json(state_path, state)
                path = checkpoint(root, f"after_{name}_STOPPED")
                print(f"\nStopped safely. Checkpoint: {path}")
                return 130
        stage["finished_unix_time"] = time.time()
        state["current_stage"] = None
        atomic_json(state_path, state)
        path = checkpoint(root, f"after_{name}")
        print(f"=== {name}: {stage['status']}; checkpoint: {path.name} ===", flush=True)
    state["status"] = "completed"
    state["finished_unix_time"] = time.time()
    state["current_stage"] = None
    atomic_json(state_path, state)
    checkpoint(root, "FINAL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
