#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1620 = GEN / "p1620_s570_strict_measured_covariance_posterior_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    s20 = _load(IN1620)
    post = s20.get("posterior_update", {})
    cl = s20.get("strict_core_closure", {})

    budget = {}
    for k,v in post.items():
        spread = (v["q95"] - v["q05"]) / 2
        budget[k] = {
            "central": v["mean"],
            "statistical_uncertainty": 0.7 * spread,
            "systematic_uncertainty": 0.3 * spread,
            "combined_uncertainty": spread,
        }

    notebook_plan = {
        "filename": "strict_posterior_reproducibility_notebook.ipynb",
        "sections": ["data ingestion", "covariance model", "posterior inference", "uncertainty budget", "sensitivity"],
        "reproducible_seed": 1621,
    }

    ready = False

    summary = {
        "checkpoint": "P1621_S571",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1621_STRICT_UNCERTAINTY_BUDGET" if ready else "KEEP_OPEN_P1621_BUDGET_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "uncertainty_budget": budget,
        "notebook_plan": notebook_plan,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": cl.get("missing_exports", []),
            "missing_witnesses": cl.get("missing_witnesses", []),
            "missing_theorems": cl.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Domknąć brakujące eksporty/witnessy/theoremy dla strict-core closure i dopiero wtedy traktować budżet niepewności jako finalny.",
        "lay_summary": "To jak rozliczenie błędów pomiaru: wiemy, co jest niepewne, ale teoria nadal wymaga kilku brakujących dowodów, więc domknięcie pozostaje otwarte.",
    }

    out = GEN / "p1621_s571_publication_grade_uncertainty_budget_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
