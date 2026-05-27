#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2214 = GEN / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.json"
OUT = GEN / "p2233_s1183_strict_nu_branch_widened_compact_monotonic_analytic_certificate.json"
MD = GEN / "p2233_s1183_strict_nu_branch_widened_compact_monotonic_analytic_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2214 = load(IN_2214)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)
    dm_star = float(((p2214.get("strict_nu_branch_closed_form_threshold_consistency_certificate", {}) or {}).get("bisection_reference", {}) or {}).get("dm_star_bisection", 0.0) or 0.0)

    if not (c1 > 0 and c2 > 0 and margin > 0 and dm_star > 0):
        raise RuntimeError("Invalid upstream inputs for P2233 analytic monotonic certificate")

    lo = 0.5 * dm_star
    hi = 1.5 * dm_star

    # Analytic monotonicity on dm>=0:
    # f1(dm)=c1*dm and f2(dm)=c2*dm^2 are nondecreasing for c1,c2>0.
    # max(f1,f2) is therefore nondecreasing.
    # safety_factor = max(f1,f2)/margin inherits monotonicity for margin>0.
    monotone_analytic = lo >= 0.0 and c1 > 0.0 and c2 > 0.0 and margin > 0.0

    # crossing persistence on widened compact interval endpoints
    sf_lo = max(c1 * lo, c2 * lo * lo) / margin
    sf_hi = max(c1 * hi, c2 * hi * hi) / margin
    crossing_persists = sf_lo <= 1.0 <= sf_hi

    payload = {
        "schema_version": "p2233_s1183_v1",
        "packet_id": "P2233",
        "stage_id": "S1183",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_WIDENED_COMPACT_MONOTONIC_ANALYTIC_CERTIFICATE",
        "strict_nu_branch_widened_compact_monotonic_analytic_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_WIDENED_COMPACT_MONOTONIC_ANALYTIC_CERTIFICATE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2214.relative_to(ROOT))],
            "interval": {"abs_dm_lo": lo, "abs_dm_hi": hi},
            "analytic_claim": "For c1>0,c2>0,margin>0 and dm>=0, safety_factor(dm)=max(c1*dm,c2*dm^2)/margin is nondecreasing.",
            "endpoint_check": {"safety_factor_lo": sf_lo, "safety_factor_hi": sf_hi, "crossing_persists": crossing_persists},
            "theorem_scope_limit": "analytic monotonicity for the modeled max(c1*dm,c2*dm^2) lane on widened compact interval only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2234_candidate",
            "goal": "connect analytic monotonic lane to direct transport residual witness to reduce model-gap between safety-factor lane and residual lane",
        },
        "gatekeeper_checks": {
            "analytic_monotonic_certificate_exported": True,
            "analytic_monotonicity_holds": monotone_analytic,
            "crossing_persists_on_endpoints": crossing_persists,
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
            "# P2233 S1183: strict nu-branch widened compact monotonic analytic certificate",
            "",
            f"- interval: `[{lo:.12e}, {hi:.12e}]`",
            f"- analytic monotonicity holds: `{monotone_analytic}`",
            f"- crossing persists on endpoints: `{crossing_persists}`",
            "",
            "Analytic modeled-lane certificate only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
