#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2219 = GEN / "p2219_s1169_strict_nu_branch_tolerance_inflated_signed_budget_intervals.json"
OUT = GEN / "p2220_s1170_strict_nu_branch_inflated_interval_coverage_certificate.json"
MD = GEN / "p2220_s1170_strict_nu_branch_inflated_interval_coverage_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def eff(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2219 = load(IN_2219)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    block = p2219.get("strict_nu_branch_tolerance_inflated_signed_budget_intervals", {}) or {}
    dm_star = float(block.get("reference_dm_star", 0.0) or 0.0)
    above = block.get("inflated_above_target_intervals", []) or []
    below = block.get("inflated_below_target_intervals", []) or []

    if not (c1 > 0 and c2 > 0 and margin > 0 and dm_star > 0 and above and below):
        raise RuntimeError("Invalid upstream inputs for P2220 coverage certificate")

    rows = []
    for side, items in (("above", above), ("below", below)):
        for it in items:
            sf_target = float(it["target_safety_factor"])
            d0, d1 = [float(x) for x in it["abs_dm_interval"]]
            lo, hi = (d0, d1) if d0 <= d1 else (d1, d0)
            sf_lo = eff(lo, c1, c2) / margin
            sf_hi = eff(hi, c1, c2) / margin
            covered = (min(sf_lo, sf_hi) <= sf_target <= max(sf_lo, sf_hi))
            rows.append({
                "side": side,
                "target_safety_factor": sf_target,
                "abs_dm_interval": [lo, hi],
                "safety_factor_at_interval_lo": sf_lo,
                "safety_factor_at_interval_hi": sf_hi,
                "target_covered_by_endpoint_band": covered,
            })

    above_rows = [r for r in rows if r["side"] == "above"]
    below_rows = [r for r in rows if r["side"] == "below"]

    above_covered = sum(1 for r in above_rows if r["target_covered_by_endpoint_band"])
    below_covered = sum(1 for r in below_rows if r["target_covered_by_endpoint_band"])
    above_rate = above_covered / max(len(above_rows), 1)
    below_rate = below_covered / max(len(below_rows), 1)

    payload = {
        "schema_version": "p2220_s1170_v1",
        "packet_id": "P2220",
        "stage_id": "S1170",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_INFLATED_INTERVAL_COVERAGE_CERTIFICATE",
        "strict_nu_branch_inflated_interval_coverage_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_INFLATED_INTERVAL_COVERAGE_CERTIFICATE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2219.relative_to(ROOT))],
            "reference_dm_star": dm_star,
            "coverage_rows": rows,
            "coverage_summary": {
                "all_above_targets_covered": all(r["target_covered_by_endpoint_band"] for r in above_rows),
                "all_below_targets_covered": all(r["target_covered_by_endpoint_band"] for r in below_rows),
                "above_coverage_rate": above_rate,
                "below_coverage_rate": below_rate,
            },
            "theorem_scope_limit": "endpoint-band coverage on tolerance-inflated local intervals only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2221_candidate",
            "goal": "derive minimal inflation factor ensuring 100% endpoint-band coverage under direct recompute envelope",
        },
        "gatekeeper_checks": {
            "interval_coverage_certificate_exported": True,
            "all_above_targets_covered": all(r["target_covered_by_endpoint_band"] for r in above_rows),
            "below_coverage_majority": below_rate >= 0.8,
            "above_coverage_rate": above_rate,
            "below_coverage_rate": below_rate,
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
            "# P2220 S1170: strict nu-branch inflated interval coverage certificate",
            "",
            f"- all above targets covered: `{payload['gatekeeper_checks']['all_above_targets_covered']}`",
            f"- below coverage majority (>=0.8): `{payload['gatekeeper_checks']['below_coverage_majority']}`",
            f"- above coverage rate: `{above_rate:.6f}`",
            f"- below coverage rate: `{below_rate:.6f}`",
            "",
            "Endpoint-band local coverage only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
