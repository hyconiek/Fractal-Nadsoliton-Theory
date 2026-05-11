#!/usr/bin/env python3
"""P1153 strict candidate quality ranking.

Builds a transparent ranking from P1152 registry-run results.
Scoring (higher is better):
- quality_score: +1.0 gate pass, +1.0 pipeline pass, +0.5 audit pass
- repair_priority_score: for failed candidates, later failure stage => higher repair priority
This is methodological ranking only, not physical closure evidence.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"



def safe_region_margin(candidate_path: str) -> float | None:
    try:
        c = json.loads(Path(candidate_path).read_text(encoding="utf-8"))
        b = json.loads((GEN / "p1172_safe_operating_region_summary.json").read_text(encoding="utf-8"))["safe_bounds"]
        vals = {
            "omega": c.get("omega_hint"),
            "phi": c.get("phi_hint"),
            "sigma": c.get("sigma_hint"),
            "kappa": c.get("kappa_hint"),
        }
        margins = []
        for k,v in vals.items():
            if v is None:
                return None
            lo,hi = b[k]["min"], b[k]["max"]
            if not (lo <= v <= hi):
                return -1.0
            margins.append(min(v-lo, hi-v))
        return min(margins)
    except Exception:
        return None

def main() -> None:
    args = sys.argv[1:]
    src = Path(args[0]).resolve() if args else (GEN / "p1152_strict_candidate_registry_runner_summary.json").resolve()
    reg = json.loads(src.read_text(encoding="utf-8"))

    ranking = []
    for r in reg.get("results", []):
        passed = bool(r.get("pass"))
        failure_stage = r.get("failure_stage")
        stage_index = failure_stage.get("index") if isinstance(failure_stage, dict) else None
        gate_pass = passed or (stage_index is not None and stage_index > 0)
        audit_pass = passed
        score = (1.0 if gate_pass else 0.0) + (1.0 if passed else 0.0) + (0.5 if audit_pass else 0.0)
        repair_priority_score = float(stage_index) if (not passed and stage_index is not None) else (-1.0 if passed else 0.0)
        margin = safe_region_margin(r.get("candidate"))
        ranking.append(
            {
                "candidate": r.get("candidate"),
                "gate_pass": gate_pass,
                "pipeline_pass": passed,
                "audit_pass": audit_pass,
                "quality_score": score,
                "repair_priority_score": repair_priority_score,
                "status": "admissible_candidate_only" if passed else "blocked",
                "failure_stage": failure_stage,
                "safe_region_margin": margin,
            }
        )

    ranking.sort(key=lambda x: (-x["quality_score"], -(x["safe_region_margin"] if isinstance(x["safe_region_margin"], (int,float)) else -999), -x["repair_priority_score"], x["candidate"] or ""))
    out = {
        "packet": "P1153",
        "as_of": "2026-05-10",
        "source_summary": str(src),
        "ranking": ranking,
        "top_recommendation": ranking[0]["candidate"] if ranking else None,
        "note": "Methodological ranking only; no strict-core closure or QW-2191 discharge claim.",
    }

    out_path = GEN / "p1153_strict_candidate_quality_ranking_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1153] ranked {len(ranking)} candidates, wrote {out_path}")


if __name__ == "__main__":
    main()
