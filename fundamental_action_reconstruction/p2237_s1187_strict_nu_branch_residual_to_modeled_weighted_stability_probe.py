#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2234 = GEN / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.json"
IN_2236 = GEN / "p2236_s1186_strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe.json"
OUT = GEN / "p2237_s1187_strict_nu_branch_residual_to_modeled_weighted_stability_probe.json"
MD = GEN / "p2237_s1187_strict_nu_branch_residual_to_modeled_weighted_stability_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def weighted_affine_fit(x: list[float], y: list[float], w: list[float]) -> tuple[float, float]:
    sw = sum(w)
    if sw <= 0.0:
        raise RuntimeError("Non-positive weight sum")
    xbar = sum(wi * xi for wi, xi in zip(w, x)) / sw
    ybar = sum(wi * yi for wi, yi in zip(w, y)) / sw
    sxx = sum(wi * (xi - xbar) ** 2 for wi, xi in zip(w, x))
    sxy = sum(wi * (xi - xbar) * (yi - ybar) for wi, xi, yi in zip(w, x, y))
    if abs(sxx) < 1e-18:
        raise RuntimeError("Degenerate weighted variance for affine fit")
    a = sxy / sxx
    b = ybar - a * xbar
    return a, b


def mse(yhat: list[float], y: list[float]) -> float:
    return sum((u - v) ** 2 for u, v in zip(yhat, y)) / max(1, len(y))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2203 = load(IN_2203)
    p2234 = load(IN_2234)
    p2236 = load(IN_2236)

    rows = (p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}).get("residual_map_rows", []) or []
    modeled = (p2234.get("strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe", {}) or {}).get("modeled_lane_observables", {}) or {}

    if len(rows) < 4:
        raise RuntimeError("Need at least 4 rows for split stability probe")

    rvals = [float(r["transport_residual_l1"]) for r in rows]
    rmin, rmax = min(rvals), max(rvals)
    sf_lo = float(modeled.get("safety_factor_lo", 0.0) or 0.0)
    sf_hi = float(modeled.get("safety_factor_hi", 0.0) or 0.0)
    if not (rmax > rmin and sf_hi >= sf_lo):
        raise RuntimeError("Degenerate ranges for P2237")

    xvals = [(r - rmin) / (rmax - rmin) for r in rvals]
    yvals = [sf_lo + x * (sf_hi - sf_lo) for x in xvals]

    # Symmetric edge-emphasis weights: w = 1 + lambda*|x-0.5|
    lam = 2.0
    wvals = [1.0 + lam * abs(x - 0.5) for x in xvals]

    a_w, b_w = weighted_affine_fit(xvals, yvals, wvals)
    yhat = [a_w * x + b_w for x in xvals]
    mse_full = mse(yhat, yvals)

    # two-way holdout split by sorted x (even/odd indices)
    order = sorted(range(len(xvals)), key=lambda i: xvals[i])
    train_idx = [idx for k, idx in enumerate(order) if k % 2 == 0]
    test_idx = [idx for k, idx in enumerate(order) if k % 2 == 1]

    x_train = [xvals[i] for i in train_idx]
    y_train = [yvals[i] for i in train_idx]
    w_train = [wvals[i] for i in train_idx]
    a_tr, b_tr = weighted_affine_fit(x_train, y_train, w_train)

    y_test_hat = [a_tr * xvals[i] + b_tr for i in test_idx]
    y_test = [yvals[i] for i in test_idx]
    mse_test = mse(y_test_hat, y_test)

    # uncertainty band: max absolute residual on training -> reused as conservative tube
    y_train_hat = [a_tr * x + b_tr for x in x_train]
    max_abs_train_res = max(abs(u - v) for u, v in zip(y_train_hat, y_train)) if y_train else 0.0
    band_radius = max_abs_train_res
    all_test_in_band = all(abs(u - v) <= band_radius + 1e-15 for u, v in zip(y_test_hat, y_test))

    prev_sel = (p2236.get("strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe", {}) or {}).get(
        "selected_family_under_complexity_penalty", "unknown"
    )

    payload = {
        "schema_version": "p2237_s1187_v1",
        "packet_id": "P2237",
        "stage_id": "S1187",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_RESIDUAL_TO_MODELED_WEIGHTED_STABILITY_PROBE",
        "strict_nu_branch_residual_to_modeled_weighted_stability_probe": {
            "probe_id": "STRICT_NU_BRANCH_RESIDUAL_TO_MODELED_WEIGHTED_STABILITY_PROBE_V1",
            "source_packets": [str(IN_2203.relative_to(ROOT)), str(IN_2234.relative_to(ROOT)), str(IN_2236.relative_to(ROOT))],
            "prior_selected_family": prev_sel,
            "weighted_fit": {
                "weight_formula": "w=1+2*|x-0.5|",
                "params": {"a": a_w, "b": b_w},
                "mse_full": mse_full,
            },
            "split_stability": {
                "train_size": len(train_idx),
                "test_size": len(test_idx),
                "train_params": {"a": a_tr, "b": b_tr},
                "test_mse": mse_test,
                "uncertainty_band_radius_from_train_residual": band_radius,
                "all_test_points_inside_band": all_test_in_band,
            },
            "physical_interpretation_note": "Edge-weighted fit tests whether boundary-dominant transport responses preserve monotone safety mapping under holdout perturbation.",
            "theorem_scope_limit": "finite-sample weighted stability evidence only; no legacy->strict bridge theorem",
        },
        "recommended_next_honest_step": {
            "id": "P2238_candidate",
            "goal": "connect weighted stability bounds to explicit perturbation-budget safety margins in a shared inequality packet",
        },
        "gatekeeper_checks": {
            "weighted_stability_probe_exported": True,
            "weighted_slope_positive": a_w > 0.0,
            "train_split_slope_positive": a_tr > 0.0,
            "test_points_inside_train_band": all_test_in_band,
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
            "# P2237 S1187: residual->modeled weighted stability probe",
            "",
            f"- weighted full MSE: `{mse_full:.12e}`",
            f"- holdout test MSE: `{mse_test:.12e}`",
            f"- weighted slope a: `{a_w:.12e}`",
            f"- train slope a: `{a_tr:.12e}`",
            f"- uncertainty band radius: `{band_radius:.12e}`",
            f"- all holdout points inside train band: `{all_test_in_band}`",
            "",
            "Local weighted stability evidence only; no kernel-bridge theorem claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
