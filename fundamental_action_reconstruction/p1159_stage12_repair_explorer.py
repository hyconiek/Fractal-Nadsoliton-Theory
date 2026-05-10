#!/usr/bin/env python3
"""P1159 stage1/2 repair explorer.

Framework for conservative repair attempts for candidates failing at probe stages
1 or 2. If no such candidates are present, emits an explicit skipped report.
No closure/discharge claims.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_pipeline(candidate: Path) -> dict:
    cmd = [sys.executable, str(ROOT / "p1151_strict_selector_pipeline_runner.py"), str(candidate)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    summ = json.loads((GEN / "p1151_strict_selector_pipeline_runner_summary.json").read_text(encoding="utf-8"))
    return {"returncode": p.returncode, "summary": summ}


def main() -> None:
    src = GEN / "p1152_strict_candidate_registry_runner_summary.json"
    registry = json.loads(src.read_text(encoding="utf-8"))

    targets = []
    for r in registry.get("results", []):
        fs = r.get("failure_stage")
        if isinstance(fs, dict) and fs.get("index") in (1, 2):
            targets.append(r)

    attempts = []
    for t in targets:
        cpath = Path(t["candidate"])
        payload = json.loads(cpath.read_text(encoding="utf-8"))

        # conservative tweak scaffold (only if numeric hint fields exist)
        variant = dict(payload)
        changed = []
        for k in ("omega_hint", "phi_hint", "beta_hint", "eta_hint"):
            if isinstance(variant.get(k), (int, float)):
                variant[k] = float(variant[k]) * 0.98
                changed.append(k)

        vpath = GEN / f"{cpath.stem}_p1159_stage12_variant.json"
        vpath.write_text(json.dumps(variant, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        res = run_pipeline(vpath)
        attempts.append({
            "original_candidate": str(cpath),
            "variant_candidate": str(vpath),
            "changed_hint_keys": changed,
            "variant_overall_pass": bool(res["summary"].get("overall_pass")),
            "variant_failed_step": res["summary"].get("failed_step"),
        })

    out = {
        "packet": "P1159",
        "as_of": "2026-05-10",
        "source_registry_summary": str(src),
        "target_count": len(targets),
        "attempt_count": len(attempts),
        "attempts": attempts,
        "skipped_reason": "no stage1/stage2 failures in current registry" if not targets else None,
        "note": "Exploratory stage1/2 repair framework only; no closure/QW-2191 discharge claim.",
    }

    out_path = GEN / "p1159_stage12_repair_explorer_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1159] targets={len(targets)} attempts={len(attempts)} wrote {out_path}")


if __name__ == "__main__":
    main()
