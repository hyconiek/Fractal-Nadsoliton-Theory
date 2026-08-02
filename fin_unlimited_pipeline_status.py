#!/usr/bin/env python3
"""Print FIN unlimited-pipeline state and short log tails without mutation."""

from __future__ import annotations

import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PID_PATH = ROOT / "FIN_Unlimited_Research_Pipeline.pid"
STATE_PATH = ROOT / "FIN_Unlimited_Research_Pipeline_State.json"


def process_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, ValueError):
        return False
    except PermissionError:
        return True
    return True


def tail(path: Path, lines: int = 12) -> list[str]:
    if not path.exists():
        return []
    return path.read_text(encoding="utf-8", errors="replace").splitlines()[-lines:]


def main() -> None:
    pid = int(PID_PATH.read_text(encoding="ascii").strip()) if PID_PATH.exists() else -1
    state = json.loads(STATE_PATH.read_text(encoding="utf-8")) if STATE_PATH.exists() else {}
    logs = {}
    for name in (
        "FIN_Program_486_Pipeline.log",
        "FIN_Program_485_Unlimited.log",
        "FIN_Program_487_Unlimited.log",
        "FIN_Program_475_Unlimited.log",
    ):
        if (ROOT / name).exists():
            logs[name] = tail(ROOT / name)
    print(json.dumps({
        "pid": pid,
        "process_alive_same_pid_namespace_only": process_alive(pid) if pid > 0 else False,
        "pid_namespace_warning": (
            "Codex tool commands may use different PID namespaces. A false PID "
            "probe does not override a running durable state or advancing logs."
        ),
        "state": state,
        "log_tails": logs,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
