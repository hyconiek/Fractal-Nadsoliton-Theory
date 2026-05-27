#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2233 = GEN / "p2233_s1183_strict_nu_branch_widened_compact_monotonic_analytic_certificate.json"
OUT = GEN / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.json"
MD = GEN / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2203 = load(IN_2203)
    p2233 = load(IN_2233)

    rblock = p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}
    rows = rblock.get("residual_map_rows", []) or []
    summary = rblock.get("summary", {}) or {}
    max_res = float(summary.get("max_residual", 0.0) or 0.0)
    mean_res = float(summary.get("mean_residual", 0.0) or 0.0)

    mblock = p2233.get("strict_nu_branch_widened_compact_monotonic_analytic_certificate", {}) or {}
    e = mblock.get("endpoint_check", {}) or {}
    sf_lo = float(e.get("safety_factor_lo", 0.0) or 0.0)
    sf_hi = float(e.get("safety_factor_hi", 0.0) or 0.0)

    if not (rows and max_res > 0.0 and sf_hi >= sf_lo >= 0.0):
        raise RuntimeError("Invalid upstream inputs for P2234 lane-gap probe")

    # Honest quantification only: no direct theorem-level mapping between residual-lane and safety-factor lane.
    residual_dynamic_ratio = max_res / max(mean_res, 1e-30)
    safety_dynamic_span = sf_hi - sf_lo

    payload = {
        "schema_version": "p2234_s1184_v1",
        "packet_id": "P2234",
        "stage_id": "S1184",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_MODELED_TO_RESIDUAL_LANE_GAP_QUANTIFICATION_PROBE",
        "strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe": {
            "probe_id": "STRICT_NU_BRANCH_MODELED_TO_RESIDUAL_LANE_GAP_QUANTIFICATION_PROBE_V1",
            "source_packets": [str(IN_2203.relative_to(ROOT)), str(IN_2233.relative_to(ROOT))],
            "residual_lane_observables": {
                "max_residual": max_res,
                "mean_residual": mean_res,
                "max_over_mean_ratio": residual_dynamic_ratio,
            },
            "modeled_lane_observables": {
                "safety_factor_lo": sf_lo,
                "safety_factor_hi": sf_hi,
                "safety_factor_span": safety_dynamic_span,
            },
            "bridge_gap_statement": "No theorem-grade operator bridge exported between residual lane and modeled safety-factor lane; current packet exports only side-by-side quantitative diagnostics.",
            "theorem_scope_limit": "diagnostic gap quantification only; not a bridge theorem and not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2235_candidate",
            "goal": "construct explicit residual->safety surrogate map or prove non-existence on current strict lane",
        },
        "gatekeeper_checks": {
            "gap_quantification_exported": True,
            "no_bridge_theorem_claimed": True,
            "diagnostic_observables_finite": residual_dynamic_ratio > 0.0 and safety_dynamic_span >= 0.0,
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
            "# P2234 S1184: strict nu-branch modeled-lane vs residual-lane gap quantification probe",
            "",
            f"- residual max/mean ratio: `{residual_dynamic_ratio:.12e}`",
            f"- modeled safety span: `{safety_dynamic_span:.12e}`",
            "- bridge theorem exported: `False`",
            "",
            "Diagnostic-only gap quantification; no bridge theorem and no global closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
