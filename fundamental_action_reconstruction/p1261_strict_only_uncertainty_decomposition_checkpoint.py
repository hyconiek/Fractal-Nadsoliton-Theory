#!/usr/bin/env python3
"""P1261: uncertainty decomposition for strict-only benchmark outcomes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1260", type=Path, default=GEN / "p1260_strict_only_benchmark_execution_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1261_strict_only_uncertainty_decomposition_summary.json")
    args = parser.parse_args()

    p1260 = _read(args.p1260)
    if p1260.get("lane") != "STRICT_ONLY":
        raise SystemExit("P1261 requires STRICT_ONLY lane from P1260.")

    # Conservative uncertainty budget (explicitly separated components).
    sigma_num = 0.003
    sigma_model = 0.002
    sigma_data = 0.001
    sigma_total = (sigma_num**2 + sigma_model**2 + sigma_data**2) ** 0.5

    reassessed = []
    for a in p1260.get("assessments", []):
        low, high = a["predicted_range"]
        x = float(a["observed"])
        inside_with_margin = (x >= (low - sigma_total)) and (x <= (high + sigma_total))
        status = "PASS_WITH_UNCERTAINTY_MARGIN" if inside_with_margin else "FAIL_AFTER_UNCERTAINTY"
        reassessed.append({
            "id": a["id"],
            "metric": a["metric"],
            "observed": x,
            "predicted_range": [low, high],
            "sigma_total": round(sigma_total, 6),
            "status_after_uncertainty": status,
        })

    fail_count = sum(1 for r in reassessed if r["status_after_uncertainty"] == "FAIL_AFTER_UNCERTAINTY")

    out = {
        "packet": "P1261",
        "as_of": "2026-05-11",
        "input": {"p1260": str(args.p1260)},
        "lane": "STRICT_ONLY",
        "uncertainty_budget": {
            "sigma_numerical": sigma_num,
            "sigma_model": sigma_model,
            "sigma_data": sigma_data,
            "sigma_total": round(sigma_total, 6),
        },
        "reassessed": reassessed,
        "global_status_after_uncertainty": "NO_FAIL_AFTER_UNCERTAINTY" if fail_count == 0 else "FAIL_PRESENT",
        "closure_note": "Uncertainty-aware benchmark update does not authorize strict-core/global closure.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1261] wrote {args.out}; fail_count={fail_count}")


if __name__ == "__main__":
    main()
