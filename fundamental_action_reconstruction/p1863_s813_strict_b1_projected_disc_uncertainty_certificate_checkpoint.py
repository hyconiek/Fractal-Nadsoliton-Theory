#!/usr/bin/env python3
"""P1863 S813 strict B1 projected discontinuity uncertainty certificate checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1862 = load("p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json")

    rows = ((p1862.get("projected_discontinuity_seed_evaluation") or {}).get("rows") or [])

    report = []
    for r in rows:
        val = float(r.get("disc_projected_seed", 0.0))
        sigma = max(1e-6, 0.1 * abs(val))
        lower = val - sigma
        report.append(
            {
                "s": r.get("s"),
                "theta": r.get("theta"),
                "disc_center": val,
                "sigma_seed": sigma,
                "lower_1sigma": lower,
                "positive_with_margin": lower >= 0.0,
            }
        )

    all_positive_margin = all(x.get("positive_with_margin", False) for x in report)

    out = {
        "packet_id": "P1863",
        "stage_id": "S813",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1862_present": "projected_discontinuity_seed_evaluation" in p1862,
        },
        "projected_discontinuity_uncertainty_report": {
            "method": "seed_uncertainty_band_sigma=0.1*|disc_center|",
            "rows": report,
            "all_positive_with_1sigma_margin": all_positive_margin,
        },
        "strict_chain_extension": "K_strict -> projected discontinuity seed -> uncertainty-banded positivity check",
        "proven": "A first uncertainty-banded positivity report for projected discontinuity has been exported.",
        "open": "Uncertainty model is seed-level; exact error propagation from full integrals is not yet available.",
        "false_pass_risk": "Seed uncertainty bands cannot replace theorem-grade unitarity error analysis.",
        "next_honest_step": "Derive uncertainty from exact phase-space integral inputs and dressed propagator covariance, then re-evaluate positivity certificate.",
        "lay_explanation": "Dodaliśmy prosty margines niepewności do testu unitarności, ale pełna wiarygodność wymaga niepewności policzonej z dokładnych całek.",
    }

    path = GEN / "p1863_s813_strict_b1_projected_disc_uncertainty_certificate_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
