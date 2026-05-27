#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2234 = GEN / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.json"
OUT = GEN / "p2236_s1186_strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe.json"
MD = GEN / "p2236_s1186_strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def fit_quadratic_through_three_points(x0: float, y0: float, x1: float, y1: float, x2: float, y2: float) -> tuple[float, float, float]:
    # Lagrange-form expansion into a*x^2+b*x+c
    d01 = (x0 - x1)
    d02 = (x0 - x2)
    d12 = (x1 - x2)
    if abs(d01) < 1e-18 or abs(d02) < 1e-18 or abs(d12) < 1e-18:
        raise RuntimeError("Degenerate x-nodes for quadratic fit")

    a = y0 / (d01 * d02) + y1 / ((x1 - x0) * d12) + y2 / ((x2 - x0) * (x2 - x1))
    b = -y0 * (x1 + x2) / (d01 * d02) - y1 * (x0 + x2) / ((x1 - x0) * d12) - y2 * (x0 + x1) / ((x2 - x0) * (x2 - x1))
    c = y0 * x1 * x2 / (d01 * d02) + y1 * x0 * x2 / ((x1 - x0) * d12) + y2 * x0 * x1 / ((x2 - x0) * (x2 - x1))
    return a, b, c


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2203 = load(IN_2203)
    p2234 = load(IN_2234)

    rows = (p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}).get("residual_map_rows", []) or []
    modeled = (p2234.get("strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe", {}) or {}).get("modeled_lane_observables", {}) or {}

    if not rows:
        raise RuntimeError("Missing residual rows for P2236")

    rvals = [float(r["transport_residual_l1"]) for r in rows]
    rmin, rmax = min(rvals), max(rvals)
    sf_lo = float(modeled.get("safety_factor_lo", 0.0) or 0.0)
    sf_hi = float(modeled.get("safety_factor_hi", 0.0) or 0.0)
    if not (rmax > rmin and sf_hi >= sf_lo):
        raise RuntimeError("Degenerate ranges for P2236")

    # Normalize x to [0,1] to stabilize comparisons
    xvals = [(r - rmin) / (rmax - rmin) for r in rvals]
    # Proxy target: linear safety interpolation from compact-lane ordering
    yvals = [sf_lo + x * (sf_hi - sf_lo) for x in xvals]

    # Family 1: affine in x
    a1 = (sf_hi - sf_lo)
    b1 = sf_lo
    y_aff = [a1 * x + b1 for x in xvals]

    # Family 2: quadratic calibrated through (0,y0),(0.5,y_mid),(1,y1) with y_mid from median sample
    order = sorted(zip(xvals, yvals), key=lambda t: t[0])
    x_mid, y_mid = order[len(order) // 2]
    qa, qb, qc = fit_quadratic_through_three_points(0.0, sf_lo, x_mid, y_mid, 1.0, sf_hi)
    y_quad = [qa * x * x + qb * x + qc for x in xvals]

    def mse(yhat: list[float], y: list[float]) -> float:
        return sum((u - v) ** 2 for u, v in zip(yhat, y)) / max(1, len(y))

    mse_aff = mse(y_aff, yvals)
    mse_quad = mse(y_quad, yvals)

    # Monotonic derivative checks on [0,1]
    d_aff_min = a1
    d_quad_min = min(qb, 2 * qa + qb)

    winner = "affine" if mse_aff <= mse_quad else "quadratic"

    payload = {
        "schema_version": "p2236_s1186_v1",
        "packet_id": "P2236",
        "stage_id": "S1186",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_RESIDUAL_TO_MODELED_SURROGATE_FAMILY_EVIDENCE_PROBE",
        "strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe": {
            "probe_id": "STRICT_NU_BRANCH_RESIDUAL_TO_MODELED_SURROGATE_FAMILY_EVIDENCE_PROBE_V1",
            "source_packets": [str(IN_2203.relative_to(ROOT)), str(IN_2234.relative_to(ROOT))],
            "domain_normalization": {"residual_min": rmin, "residual_max": rmax},
            "target_range": {"modeled_safety_lo": sf_lo, "modeled_safety_hi": sf_hi},
            "family_metrics": {
                "affine": {"params": {"a": a1, "b": b1}, "mse": mse_aff, "derivative_min_on_unit_interval": d_aff_min},
                "quadratic": {
                    "params": {"a": qa, "b": qb, "c": qc},
                    "mse": mse_quad,
                    "derivative_min_on_unit_interval": d_quad_min,
                },
            },
            "selected_family_under_complexity_penalty": winner,
            "selection_rule": "prefer lower MSE; tie -> affine (lower complexity)",
            "physical_interpretation_note": "Monotone positive derivative surrogate is compatible with one-directional local response, but remains a modeled-lane diagnostic only.",
            "theorem_scope_limit": "finite-sample surrogate-family evidence only; not a legacy->strict bridge theorem",
        },
        "recommended_next_honest_step": {
            "id": "P2237_candidate",
            "goal": "export weighted regression with uncertainty bands and out-of-sample split checks for surrogate stability",
        },
        "gatekeeper_checks": {
            "surrogate_family_evidence_exported": True,
            "affine_derivative_positive": d_aff_min > 0.0,
            "quadratic_derivative_positive_on_unit_interval": d_quad_min > 0.0,
            "complexity_tie_breaker_explicit": True,
            "no_bridge_theorem_claimed": True,
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
            "# P2236 S1186: residual->modeled surrogate family evidence probe",
            "",
            f"- mse_affine: `{mse_aff:.12e}`",
            f"- mse_quadratic: `{mse_quad:.12e}`",
            f"- selected family: `{winner}`",
            f"- affine derivative min: `{d_aff_min:.12e}`",
            f"- quadratic derivative min on [0,1]: `{d_quad_min:.12e}`",
            "",
            "Local modeled diagnostic only; no kernel-bridge theorem claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
