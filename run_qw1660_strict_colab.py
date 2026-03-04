#!/usr/bin/env python3
"""
Colab-oriented strict runner for QW-1660 methodology repair branch.

Runs:
1) phase61_full_null_model_strict.py
2) phase63_whitening_test_strict.py
3) phase65_micro_timeshift_strict.py
4) QW_2116_GW1660_METHOD_REPAIR_GATE.py
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PY = os.environ.get("PYTHON_BIN", "python3")


def run(cmd: list[str], env: dict[str, str]) -> None:
    print(f"[RUN] {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=ROOT, env=env, check=False)
    if proc.returncode != 0:
        raise SystemExit(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


def main() -> None:
    n_samples = os.environ.get("QW1660_STRICT_N_SAMPLES", "524288")
    n_trials = os.environ.get("QW1660_STRICT_N_TRIALS", "400")
    n_perm = os.environ.get("QW1660_STRICT_N_PERM", "400")

    env = os.environ.copy()
    env.update(
        {
            "PHASE61_STRICT_N_SAMPLES": n_samples,
            "PHASE61_STRICT_N_TRIALS": n_trials,
            "PHASE63_STRICT_N_SAMPLES": n_samples,
            "PHASE63_STRICT_N_TRIALS": n_trials,
            "PHASE65_STRICT_N_SAMPLES": n_samples,
            "PHASE65_STRICT_N_PERM": n_perm,
        }
    )

    run([PY, "phase61_full_null_model_strict.py"], env)
    run([PY, "phase63_whitening_test_strict.py"], env)
    run([PY, "phase65_micro_timeshift_strict.py"], env)
    run([PY, "QW_2116_GW1660_METHOD_REPAIR_GATE.py"], env)

    gate = ROOT / "report_qw2116_gw1660_method_repair_gate.json"
    if gate.exists():
        obj = json.loads(gate.read_text(encoding="utf-8"))
        print(json.dumps(
            {
                "verdict": obj.get("verdict"),
                "strict_protocol": obj.get("strict_protocol"),
                "post_repair_branch_status": obj.get("post_repair_branch_status"),
            },
            ensure_ascii=False,
            indent=2,
        ))
    else:
        print("Missing report_qw2116_gw1660_method_repair_gate.json")


if __name__ == "__main__":
    main()

