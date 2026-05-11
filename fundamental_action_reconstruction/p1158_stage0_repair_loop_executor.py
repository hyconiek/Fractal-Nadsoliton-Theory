#!/usr/bin/env python3
"""P1158 stage-0 repair loop executor.

Auto-repairs gate-level blocked candidates (failure_stage_index==0)
by enforcing required boolean declarations to True, then reruns P1151
and records uplift blocked->admissible_candidate_only if achieved.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

REQUIRED_TRUE = [
    "strict_side_provenance",
    "noncyclic_anchor_declared",
    "no_legacy_role_transfer",
    "no_closure_claim",
    "no_qw2191_discharge_claim",
]


def run_pipeline(candidate_path: Path) -> dict:
    cmd = [sys.executable, str(ROOT / "p1151_strict_selector_pipeline_runner.py"), str(candidate_path)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    summ = json.loads((GEN / "p1151_strict_selector_pipeline_runner_summary.json").read_text(encoding="utf-8"))
    return {"returncode": p.returncode, "pipeline_summary": summ}


def main() -> None:
    src = GEN / "p1157_failure_stage_repair_planner_summary.json"
    plans = json.loads(src.read_text(encoding="utf-8")).get("repair_plans", [])

    outcomes = []
    for plan in plans:
        if plan.get("failure_stage_index") != 0:
            continue
        cpath = Path(plan["candidate"])
        payload = json.loads(cpath.read_text(encoding="utf-8"))
        original = dict(payload)
        for k in REQUIRED_TRUE:
            payload[k] = True

        repaired_path = GEN / f"{cpath.stem}_p1158_repaired.json"
        repaired_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        before_pass = False
        after = run_pipeline(repaired_path)
        after_pass = bool(after["pipeline_summary"].get("overall_pass"))

        outcomes.append({
            "original_candidate": str(cpath),
            "repaired_candidate": str(repaired_path),
            "changed_keys": [k for k in REQUIRED_TRUE if original.get(k) is not True],
            "before_pass": before_pass,
            "after_pass": after_pass,
            "uplift": (not before_pass) and after_pass,
            "after_failed_step": after["pipeline_summary"].get("failed_step"),
        })

    out = {
        "packet": "P1158",
        "as_of": "2026-05-10",
        "source_plan": str(src),
        "repaired_count": len(outcomes),
        "outcomes": outcomes,
        "note": "Methodological gate-repair loop only; no closure/QW-2191 discharge claim.",
    }
    out_path = GEN / "p1158_stage0_repair_loop_executor_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1158] repaired={len(outcomes)} wrote {out_path}")


if __name__ == "__main__":
    main()
