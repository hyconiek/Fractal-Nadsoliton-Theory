#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2219 = GEN / "p2219_s1169_strict_nu_branch_tolerance_inflated_signed_budget_intervals.json"
IN_2221 = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json"
OUT = GEN / "p2222_s1172_strict_nu_branch_calibrated_interval_width_cost_comparison.json"
MD = GEN / "p2222_s1172_strict_nu_branch_calibrated_interval_width_cost_comparison.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def interval_width(interval: list[float]) -> float:
    return abs(float(interval[1]) - float(interval[0]))


def scale_interval(center: float, interval: list[float], factor: float) -> list[float]:
    a, b = float(interval[0]), float(interval[1])
    return [center + factor * (a - center), center + factor * (b - center)]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2219 = load(IN_2219)
    p2221 = load(IN_2221)

    b = p2219.get("strict_nu_branch_tolerance_inflated_signed_budget_intervals", {}) or {}
    dm_star = float(b.get("reference_dm_star", 0.0) or 0.0)
    base_above = b.get("inflated_above_target_intervals", []) or []
    base_below = b.get("inflated_below_target_intervals", []) or []

    c = p2221.get("strict_nu_branch_min_inflation_factor_for_full_coverage", {}) or {}
    f = float(c.get("min_inflation_factor", 1.0) or 1.0)

    if not (dm_star > 0 and f >= 1.0 and base_above and base_below):
        raise RuntimeError("Invalid upstream inputs for P2222 width-cost comparison")

    rows = []
    for side, items in (("above", base_above), ("below", base_below)):
        for it in items:
            sf = float(it["target_safety_factor"])
            base = [float(x) for x in it["abs_dm_interval"]]
            calib = scale_interval(dm_star, base, f)
            rows.append({
                "side": side,
                "target_safety_factor": sf,
                "base_abs_dm_interval": base,
                "calibrated_abs_dm_interval": calib,
                "base_width": interval_width(base),
                "calibrated_width": interval_width(calib),
                "width_inflation_ratio": interval_width(calib) / max(interval_width(base), 1e-30),
            })

    mean_ratio = sum(r["width_inflation_ratio"] for r in rows) / len(rows)
    max_ratio = max(r["width_inflation_ratio"] for r in rows)

    payload = {
        "schema_version": "p2222_s1172_v1",
        "packet_id": "P2222",
        "stage_id": "S1172",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_CALIBRATED_INTERVAL_WIDTH_COST_COMPARISON",
        "strict_nu_branch_calibrated_interval_width_cost_comparison": {
            "comparison_id": "STRICT_NU_BRANCH_CALIBRATED_INTERVAL_WIDTH_COST_COMPARISON_V1",
            "source_packets": [str(IN_2219.relative_to(ROOT)), str(IN_2221.relative_to(ROOT))],
            "reference_dm_star": dm_star,
            "min_inflation_factor": f,
            "rows": rows,
            "summary": {
                "mean_width_inflation_ratio": mean_ratio,
                "max_width_inflation_ratio": max_ratio,
            },
            "theorem_scope_limit": "local interval-width cost comparison only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2223_candidate",
            "goal": "optimize target-wise inflation to minimize width cost while preserving full endpoint-band coverage",
        },
        "gatekeeper_checks": {
            "width_cost_comparison_exported": True,
            "all_width_ratios_ge_1": all(r["width_inflation_ratio"] >= 1.0 for r in rows),
            "mean_width_ratio_ge_1": mean_ratio >= 1.0,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2222 S1172: strict nu-branch calibrated interval width-cost comparison",
            "",
            f"- min inflation factor: `{f:.12e}`",
            f"- mean width inflation ratio: `{mean_ratio:.12e}`",
            f"- max width inflation ratio: `{max_ratio:.12e}`",
            "",
            "Local width-cost comparison only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
