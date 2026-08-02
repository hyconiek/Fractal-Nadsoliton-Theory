#!/usr/bin/env python3
"""Sequential durable runner for FIN P486 -> P485 -> P487 -> P475.

No subprocess timeout or resource limit is applied.  Sequential execution is
intentional on the declared i3/16-GB host to avoid simultaneous CAS memory
pressure.  A failed stage is recorded and the next independent stage starts.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parent
STATE_PATH = ROOT / "FIN_Unlimited_Research_Pipeline_State.json"
DONE_PATH = ROOT / "FIN_Unlimited_Research_Pipeline_DONE.json"
PID_PATH = ROOT / "FIN_Unlimited_Research_Pipeline.pid"

STAGES = (
    ("P486", "fin_program_486.py", "FIN_Program_486_Pipeline.log"),
    ("P485", "fin_program_485_unlimited.py", "FIN_Program_485_Unlimited.log"),
    ("P487", "fin_program_487_unlimited.py", "FIN_Program_487_Unlimited.log"),
    ("P475-unlimited", "fin_program_475_unlimited.py", "FIN_Program_475_Unlimited.log"),
)


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value, indent=2)+"\n", encoding="utf-8")
    temporary.replace(path)


def persist(state: dict[str, object]) -> None:
    state["updated_unix_time"] = time.time()
    write_json(STATE_PATH, state)


def main() -> None:
    PID_PATH.write_text(f"{os.getpid()}\n", encoding="ascii")
    state: dict[str, object] = {
        "pipeline": "FIN-P486-P485-P487-P475-unlimited-v1",
        "pid": os.getpid(),
        "started_unix_time": time.time(),
        "status": "running",
        "execution_order": [stage[0] for stage in STAGES],
        "reason_for_sequential_execution": (
            "Prevent concurrent exact-CAS jobs from competing for 16 GB RAM."
        ),
        "no_application_time_limit": True,
        "no_application_memory_limit": True,
        "network_used": False,
        "stages": {},
    }
    persist(state)
    environment = dict(os.environ)
    environment.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    environment["PYTHONUNBUFFERED"] = "1"
    for name, script, log_name in STAGES:
        stage_state = {
            "script": script,
            "log": log_name,
            "status": "running",
            "started_unix_time": time.time(),
        }
        state["current_stage"] = name
        state["stages"][name] = stage_state  # type: ignore[index]
        persist(state)
        with (ROOT / log_name).open("a", encoding="utf-8") as log:
            log.write(
                f"\n=== {name} started at unix {stage_state['started_unix_time']} "
                f"by pipeline pid {os.getpid()} ===\n"
            )
            log.flush()
            completed = subprocess.run(
                [sys.executable, "-u", str(ROOT / script)],
                cwd=ROOT,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
            log.write(
                f"\n=== {name} exited with code {completed.returncode} "
                f"at unix {time.time()} ===\n"
            )
        stage_state["finished_unix_time"] = time.time()
        stage_state["return_code"] = completed.returncode
        stage_state["status"] = "completed" if completed.returncode == 0 else "failed"
        persist(state)
    state["status"] = "completed"
    state["current_stage"] = None
    state["finished_unix_time"] = time.time()
    persist(state)
    write_json(DONE_PATH, state)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        failure = {
            "pipeline": "FIN-P486-P485-P487-P475-unlimited-v1",
            "pid": os.getpid(),
            "status": "pipeline_failed",
            "unix_time": time.time(),
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
        }
        write_json(STATE_PATH, failure)
        raise
