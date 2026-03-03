#!/usr/bin/env python3
"""
QW-1719: Flavor model with locality and group constraints.

Model:
1) Unitary matrix built via standard sequence U = R23 * U13(delta) * R12.
2) Angles derived from topological locality distances on generation manifold.
3) Shared parameter set for CKM and PMNS (single FIN-style mechanism).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1719_flavor_locality_group_constraint.json"
OUT_MD = ROOT / "RAPORT_QW1719_FLAVOR_LOCALITY_GROUP_CONSTRAINT.md"

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


def unitary_from_angles(th12: float, th23: float, th13: float, delta: float) -> np.ndarray:
    c12, s12 = np.cos(th12), np.sin(th12)
    c23, s23 = np.cos(th23), np.sin(th23)
    c13, s13 = np.cos(th13), np.sin(th13)
    e = np.exp(-1j * delta)
    r12 = np.array([[c12, s12, 0.0], [-s12, c12, 0.0], [0.0, 0.0, 1.0]], dtype=complex)
    r23 = np.array([[1.0, 0.0, 0.0], [0.0, c23, s23], [0.0, -s23, c23]], dtype=complex)
    u13 = np.array(
        [[c13, 0.0, s13 * e], [0.0, 1.0, 0.0], [-s13 * np.conj(e), 0.0, c13]],
        dtype=complex,
    )
    return r23 @ u13 @ r12


def matrix_error(pred_abs: np.ndarray, ref: np.ndarray) -> dict:
    abs_err = np.abs(pred_abs - ref)
    rel = abs_err / np.clip(ref, 1e-12, None)
    return {
        "mae": float(np.mean(abs_err)),
        "max_abs": float(np.max(abs_err)),
        "mean_rel_pct": float(np.mean(rel) * 100.0),
        "max_rel_pct": float(np.max(rel) * 100.0),
    }


def angles_from_locality(p: np.ndarray, kappa: float, lam: float, rho23: float, eta13: float):
    d12 = abs(p[0] - p[1])
    d23 = abs(p[1] - p[2])
    d13 = abs(p[0] - p[2])
    th12 = kappa / (1.0 + lam * d12)
    th23 = rho23 * kappa / (1.0 + lam * d23)
    th13 = eta13 * (kappa / (1.0 + lam * d13)) ** 2
    # keep physical range
    th12 = float(np.clip(th12, 0.0, 1.4))
    th23 = float(np.clip(th23, 0.0, 1.4))
    th13 = float(np.clip(th13, 0.0, 0.7))
    return th12, th23, th13


def main() -> None:
    # generation topology proxies (shared logic, different sectors)
    p_ckm = np.array([7.0, 9.0, 14.0], dtype=float)   # quark-like ladder proxy
    p_pmns = np.array([0.0, 1.0, 2.0], dtype=float)   # neutrino-like ladder proxy

    best = None
    for kappa in np.linspace(0.05, 1.30, 31):
        for lam in np.linspace(0.01, 2.50, 31):
            for rho23 in np.linspace(0.15, 1.20, 21):
                for eta13 in np.linspace(0.2, 2.5, 21):
                    delta = math.pi / 3.0

                    t12, t23, t13 = angles_from_locality(p_ckm, float(kappa), float(lam), float(rho23), float(eta13))
                    u_ckm = np.abs(unitary_from_angles(t12, t23, t13, delta))
                    e_ckm = matrix_error(u_ckm, CKM_REF)

                    n12, n23, n13 = angles_from_locality(p_pmns, float(kappa), float(lam), float(rho23), float(eta13))
                    u_pmns = np.abs(unitary_from_angles(n12, n23, n13, delta))
                    e_pmns = matrix_error(u_pmns, PMNS_REF)

                    reg = 0.03 * (abs(kappa) + abs(lam) + abs(rho23) + abs(eta13))
                    score = e_ckm["mean_rel_pct"] + e_pmns["mean_rel_pct"] + 0.10 * (e_ckm["max_rel_pct"] + e_pmns["max_rel_pct"]) + reg

                    if best is None or score < best["score"]:
                        best = {
                            "score": float(score),
                            "kappa": float(kappa),
                            "lam": float(lam),
                            "rho23": float(rho23),
                            "eta13": float(eta13),
                            "delta": float(delta),
                            "ckm_error": e_ckm,
                            "pmns_error": e_pmns,
                            "ckm_pred_abs": u_ckm,
                            "pmns_pred_abs": u_pmns,
                            "angles_ckm": {"th12": t12, "th23": t23, "th13": t13},
                            "angles_pmns": {"th12": n12, "th23": n23, "th13": n13},
                        }

    threshold = 15.0
    pass_shared = (best["ckm_error"]["mean_rel_pct"] < threshold) and (best["pmns_error"]["mean_rel_pct"] < threshold)
    verdict = "FLAVOR_LOCALITY_GROUP_CLOSED" if pass_shared else "FLAVOR_LOCALITY_GROUP_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "best_shared_params": {
            "kappa": best["kappa"],
            "lam": best["lam"],
            "rho23": best["rho23"],
            "eta13": best["eta13"],
            "delta": best["delta"],
        },
        "angles_ckm": best["angles_ckm"],
        "angles_pmns": best["angles_pmns"],
        "ckm_error": best["ckm_error"],
        "pmns_error": best["pmns_error"],
        "ckm_pred_abs": best["ckm_pred_abs"].tolist(),
        "pmns_pred_abs": best["pmns_pred_abs"].tolist(),
        "threshold_mean_rel_pct": threshold,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1719: FLAVOR LOCALITY + GROUP CONSTRAINT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Best shared parameter set",
        f"- kappa={best['kappa']:.4f}, lam={best['lam']:.4f}, rho23={best['rho23']:.4f}, eta13={best['eta13']:.4f}, delta={best['delta']:.4f}",
        f"- CKM mean rel. error: {best['ckm_error']['mean_rel_pct']:.2f}%",
        f"- PMNS mean rel. error: {best['pmns_error']['mean_rel_pct']:.2f}%",
        "",
        "## Interpretacja",
        "- Model narzuca unitarność (group constraint) i lokalność topologiczną, ale musi przejść jednocześnie CKM i PMNS na wspólnym zestawie parametrów.",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1719] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1719] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
