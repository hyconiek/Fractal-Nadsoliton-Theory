#!/usr/bin/env python3
"""
QW-1710: Flavor closure test for FIN using a topological overlap operator.

What this script does:
1) Tests a "derived-like" setting with frozen FIN parameters only.
2) Tests a "calibrated" setting where operator parameters are optimized.
3) Compares predicted |CKM| and |PMNS| against reference matrices.
4) Reports whether flavor sector can be considered "closed".
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1710_flavor_topological_operator.json"
OUT_MD = ROOT / "RAPORT_QW1710_FLAVOR_TOPOLOGICAL_OPERATOR.md"


# Reference matrices (central values, absolute elements; compact benchmark set)
CKM_REF = np.array(
    [
        [0.97401, 0.22650, 0.00361],
        [0.22636, 0.97320, 0.04053],
        [0.00854, 0.03978, 0.99917],
    ],
    dtype=float,
)

PMNS_REF = np.array(
    [
        [0.821, 0.550, 0.150],
        [0.432, 0.582, 0.693],
        [0.378, 0.598, 0.707],
    ],
    dtype=float,
)


def unitarize(m: np.ndarray) -> np.ndarray:
    u, _, vh = np.linalg.svd(m)
    return u @ vh


def build_overlap_matrix(
    q_left: np.ndarray,
    q_right: np.ndarray,
    eta: float,
    phase: float,
    skew: float,
) -> np.ndarray:
    n = len(q_left)
    m = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            dq = abs(q_left[i] - q_right[j])
            gen_gap = abs(i - j)
            amp = math.exp(-eta * dq) * (1.0 + skew * (i - j))
            m[i, j] = amp * np.exp(1j * phase * gen_gap)
    return unitarize(m)


def matrix_error(pred_abs: np.ndarray, ref: np.ndarray) -> dict:
    abs_err = np.abs(pred_abs - ref)
    rel_err = abs_err / np.clip(ref, 1e-12, None)
    return {
        "mae": float(np.mean(abs_err)),
        "max_abs": float(np.max(abs_err)),
        "mean_rel_pct": float(np.mean(rel_err) * 100.0),
        "max_rel_pct": float(np.max(rel_err) * 100.0),
    }


def search_best(
    q_left: np.ndarray,
    q_right: np.ndarray,
    ref: np.ndarray,
) -> dict:
    best = None
    # Compact brute-force grid: robust, transparent, deterministic.
    etas = np.linspace(0.01, 1.20, 32)
    phases = np.linspace(0.0, math.pi, 32)
    skews = np.linspace(-0.40, 0.40, 21)
    for eta in etas:
        for phase in phases:
            for skew in skews:
                pred = np.abs(build_overlap_matrix(q_left, q_right, float(eta), float(phase), float(skew)))
                e = matrix_error(pred, ref)
                score = e["mean_rel_pct"] + 0.25 * e["max_rel_pct"]
                if best is None or score < best["score"]:
                    best = {
                        "score": float(score),
                        "eta": float(eta),
                        "phase": float(phase),
                        "skew": float(skew),
                        "pred_abs": pred,
                        "error": e,
                    }
    return best


def main() -> None:
    # FIN frozen parameters
    alpha_geo = 4.0 * math.log(2.0)
    beta_tors = 0.01
    phi = math.pi / 6.0

    # Topological charge assignment used in this closure probe
    # Up sector, down sector, neutrino sector, charged-lepton sector
    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)

    # A) Derived-like test (strictly frozen)
    ckm_derived = np.abs(build_overlap_matrix(q_up, q_down, beta_tors, phi, 0.0))
    pmns_derived = np.abs(build_overlap_matrix(q_nu, q_lep, beta_tors, phi, 0.0))
    ckm_derived_err = matrix_error(ckm_derived, CKM_REF)
    pmns_derived_err = matrix_error(pmns_derived, PMNS_REF)

    # B) Calibrated best-fit test (operator closure potential)
    ckm_best = search_best(q_up, q_down, CKM_REF)
    pmns_best = search_best(q_nu, q_lep, PMNS_REF)

    # C) Shared-parameter closure test (single operator for both sectors)
    shared_best = None
    for eta in np.linspace(0.01, 1.20, 32):
        for phase in np.linspace(0.0, math.pi, 32):
            for skew in np.linspace(-0.40, 0.40, 21):
                ckm_p = np.abs(build_overlap_matrix(q_up, q_down, float(eta), float(phase), float(skew)))
                pmns_p = np.abs(build_overlap_matrix(q_nu, q_lep, float(eta), float(phase), float(skew)))
                e_ckm = matrix_error(ckm_p, CKM_REF)
                e_pmns = matrix_error(pmns_p, PMNS_REF)
                # balanced score, both sectors must be good
                score = (
                    e_ckm["mean_rel_pct"]
                    + e_pmns["mean_rel_pct"]
                    + 0.20 * (e_ckm["max_rel_pct"] + e_pmns["max_rel_pct"])
                )
                if shared_best is None or score < shared_best["score"]:
                    shared_best = {
                        "score": float(score),
                        "eta": float(eta),
                        "phase": float(phase),
                        "skew": float(skew),
                        "ckm_error": e_ckm,
                        "pmns_error": e_pmns,
                        "ckm_pred_abs": ckm_p,
                        "pmns_pred_abs": pmns_p,
                    }

    # Verdict logic
    closure_threshold_mean_rel_pct = 12.0
    derived_ok = (
        ckm_derived_err["mean_rel_pct"] < closure_threshold_mean_rel_pct
        and pmns_derived_err["mean_rel_pct"] < closure_threshold_mean_rel_pct
    )
    calibrated_shared_ok = (
        shared_best["ckm_error"]["mean_rel_pct"] < closure_threshold_mean_rel_pct
        and shared_best["pmns_error"]["mean_rel_pct"] < closure_threshold_mean_rel_pct
    )

    if derived_ok:
        verdict = "FLAVOR_CLOSURE_DERIVED_PLAUSIBLE"
    elif calibrated_shared_ok:
        verdict = "FLAVOR_CLOSURE_REQUIRES_CALIBRATED_OPERATOR"
    else:
        verdict = "FLAVOR_NOT_CLOSED_YET"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "frozen_params": {
                "alpha_geo": alpha_geo,
                "beta_tors": beta_tors,
                "phi": phi,
            },
            "q_assignments": {
                "q_up": q_up.tolist(),
                "q_down": q_down.tolist(),
                "q_nu": q_nu.tolist(),
                "q_lep": q_lep.tolist(),
            },
        },
        "derived_like_test": {
            "ckm_pred_abs": ckm_derived.tolist(),
            "pmns_pred_abs": pmns_derived.tolist(),
            "ckm_error": ckm_derived_err,
            "pmns_error": pmns_derived_err,
        },
        "calibrated_independent_best": {
            "ckm_best": {
                "eta": ckm_best["eta"],
                "phase": ckm_best["phase"],
                "skew": ckm_best["skew"],
                "error": ckm_best["error"],
            },
            "pmns_best": {
                "eta": pmns_best["eta"],
                "phase": pmns_best["phase"],
                "skew": pmns_best["skew"],
                "error": pmns_best["error"],
            },
        },
        "shared_operator_best": {
            "eta": shared_best["eta"],
            "phase": shared_best["phase"],
            "skew": shared_best["skew"],
            "ckm_error": shared_best["ckm_error"],
            "pmns_error": shared_best["pmns_error"],
            "ckm_pred_abs": shared_best["ckm_pred_abs"].tolist(),
            "pmns_pred_abs": shared_best["pmns_pred_abs"].tolist(),
        },
        "verdict": verdict,
        "thresholds": {
            "closure_threshold_mean_rel_pct": closure_threshold_mean_rel_pct,
        },
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1710: FLAVOR TOPOLOGICAL OPERATOR",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## 1) Test derived-like (tylko zamrożone FIN)",
        f"- CKM mean rel. error: {ckm_derived_err['mean_rel_pct']:.2f}%",
        f"- PMNS mean rel. error: {pmns_derived_err['mean_rel_pct']:.2f}%",
        "",
        "## 2) Najlepszy operator współdzielony (kalibracja)",
        f"- eta={shared_best['eta']:.4f}, phase={shared_best['phase']:.4f}, skew={shared_best['skew']:.4f}",
        f"- CKM mean rel. error: {shared_best['ckm_error']['mean_rel_pct']:.2f}%",
        f"- PMNS mean rel. error: {shared_best['pmns_error']['mean_rel_pct']:.2f}%",
        "",
        "## 3) Interpretacja",
        "- Jeśli test derived-like nie przechodzi, flavor nie jest domknięty z samego jądra.",
        "- Jeśli operator współdzielony przechodzi po kalibracji, to istnieje spójny kierunek domknięcia, ale klasy `C/M`, nie `D`.",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1710] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1710] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
