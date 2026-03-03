#!/usr/bin/env python3
"""
QW-1943: Topological Q-assignment audit for mass+flavor under single frozen kernel.

Purpose:
- test whether current blockers can be explained mainly by Q assignments,
- keep one strict operator mapping (no per-sector fitting),
- scan a constrained, physically local family of integer Q assignments.
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import polar


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1943_topological_q_assignment_audit.json"
OUT_MD = ROOT / "RAPORT_QW1943_TOPOLOGICAL_Q_ASSIGNMENT_AUDIT.md"


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

MASS_EXP = {
    "Top": 173_000.0,
    "Bottom": 4_180.0,
    "Tau": 1_776.9,
    "Charm": 1_270.0,
    "Muon": 105.7,
    "Electron": 0.511,
}


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def cyclic_distance_matrix(q_left: np.ndarray, q_right: np.ndarray, modulus: int = 24) -> np.ndarray:
    dq = np.abs(q_left[:, None] - q_right[None, :]) % float(modulus)
    return np.minimum(dq, float(modulus) - dq)


def derive_shared_params(omega: float, phi: float, beta: float, eta: float) -> Dict[str, float]:
    # Same strict mapping as QW-1940 to isolate assignment effect.
    d = np.arange(1.0, 13.0, dtype=float)
    k = np.abs(kernel_fn(d, omega, phi, beta, eta))
    w = k / max(np.sum(k), 1e-15)
    mean_d = float(np.sum(w * d))
    var_d = float(np.sum(w * (d - mean_d) ** 2))
    decay_ratio = float(k[3] / max(k[0], 1e-15))

    p_amp = float(np.clip(1.0 + 0.60 * np.tanh((mean_d - 2.0) / 2.0), 0.30, 2.20))
    r_dist = float(np.clip(np.tanh((var_d - 1.0) / 2.5), -1.20, 1.80))
    phase_scale = float(np.clip(1.0 + 0.70 * np.tanh(abs(phi)) + 0.25 * np.tanh(1.0 - decay_ratio), 0.00, 3.00))
    return {"p_amp": p_amp, "r_dist": r_dist, "phase_scale": phase_scale}


def flavor_prediction_abs(
    q_left: np.ndarray,
    q_right: np.ndarray,
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
) -> np.ndarray:
    n = len(q_left)
    d = 1.0 + cyclic_distance_matrix(q_left, q_right, modulus=24)
    k = kernel_fn(d, omega=omega, phi=phi, beta=beta, eta=eta)

    amp = np.sign(k) * ((np.abs(k) ** p_amp) * (d**r_dist))
    idx = np.arange(n, dtype=float)
    gap = np.abs(idx[:, None] - idx[None, :])
    sign = np.where((idx[:, None] - idx[None, :]) < 0.0, 1.0, -1.0)
    phase = np.exp(1j * (phi + phase_scale * omega * sign * gap))

    m = amp * phase
    u = polar(m)[0]
    return np.abs(u)


def flavor_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, float]:
    ckm = flavor_prediction_abs(
        q_up,
        q_down,
        params["p_amp"],
        params["r_dist"],
        params["phase_scale"],
        kernel["omega"],
        kernel["phi"],
        kernel["beta"],
        kernel["eta"],
    )
    pmns = flavor_prediction_abs(
        q_nu,
        q_lep,
        params["p_amp"],
        params["r_dist"],
        params["phase_scale"],
        kernel["omega"],
        kernel["phi"],
        kernel["beta"],
        kernel["eta"],
    )

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)
    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
    }


def mass_metrics_from_assignment(mass_q: Dict[str, float], kernel: Dict[str, float]) -> Dict[str, float]:
    # Exact hard mass formula with gamma from d=1..4 (as QW-1939 primary).
    k1 = float(abs(kernel_fn(np.array([1.0]), kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    k4 = float(abs(kernel_fn(np.array([4.0]), kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    gamma = float(-4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0)

    errs = []
    for name, m_exp in MASS_EXP.items():
        q = float(mass_q[name])
        pred = 173_000.0 * (4.0 ** (-(gamma * q / 4.0)))
        err = float(abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0)
        errs.append(err)
    return {
        "gamma": gamma,
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
    }


def build_assignment_space() -> Dict[str, List[Tuple[int, int, int]]]:
    q_up = []
    for q2 in [8, 9, 10]:
        for q1 in [13, 14, 15]:
            if len({0, q2, q1}) == 3:
                q_up.append((0, q2, q1))

    q_down = []
    for q3 in [6, 7, 8]:
        for q2 in [8, 9, 10]:
            for q1 in [13, 14, 15]:
                if len({q3, q2, q1}) == 3:
                    q_down.append((q3, q2, q1))

    q_nu = list(itertools.permutations([0, 1, 2], 3))

    q_lep = []
    for qe in [23, 24]:
        for qmu in [13, 14, 15]:
            for qtau in [8, 9, 10]:
                if len({qe, qmu, qtau}) == 3:
                    q_lep.append((qe, qmu, qtau))

    return {"q_up": q_up, "q_down": q_down, "q_nu": q_nu, "q_lep": q_lep}


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }
    params = derive_shared_params(**kernel)

    thresholds = {
        "mass_mean_rel_pct_max": 15.0,
        "mass_max_rel_pct_max": 75.0,
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
    }

    space = build_assignment_space()

    rows = []
    best_joint = None
    best_mass = None
    best_flavor = None

    mass_pass_count = 0
    flavor_pass_count = 0
    joint_pass_count = 0

    for qu in space["q_up"]:
        q_up = np.array(qu, dtype=float)
        for qd in space["q_down"]:
            q_down = np.array(qd, dtype=float)
            for qn in space["q_nu"]:
                q_nu = np.array(qn, dtype=float)
                for ql in space["q_lep"]:
                    q_lep = np.array(ql, dtype=float)

                    f = flavor_metrics(q_up, q_down, q_nu, q_lep, params=params, kernel=kernel)

                    mass_q = {
                        "Top": float(q_up[0]),
                        "Charm": float(q_up[1]),
                        "Bottom": float(q_down[0]),
                        "Electron": float(q_lep[0]),
                        "Muon": float(q_lep[1]),
                        "Tau": float(q_lep[2]),
                    }
                    m = mass_metrics_from_assignment(mass_q=mass_q, kernel=kernel)

                    mass_pass = bool(
                        m["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]
                        and m["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]
                    )
                    flavor_pass = bool(
                        f["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]
                        and f["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]
                    )
                    joint_pass = bool(mass_pass and flavor_pass)

                    if mass_pass:
                        mass_pass_count += 1
                    if flavor_pass:
                        flavor_pass_count += 1
                    if joint_pass:
                        joint_pass_count += 1

                    mass_loss = m["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"] + m["max_rel_err_pct"] / thresholds["mass_max_rel_pct_max"]
                    flavor_loss = f["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"] + f["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
                    joint_loss = float(mass_loss + flavor_loss)

                    row = {
                        "q_up": [int(x) for x in q_up],
                        "q_down": [int(x) for x in q_down],
                        "q_nu": [int(x) for x in q_nu],
                        "q_lep": [int(x) for x in q_lep],
                        "mass_q": {k: int(v) for k, v in mass_q.items()},
                        "metrics": {
                            "mass_mean_rel_pct": m["mean_rel_err_pct"],
                            "mass_max_rel_pct": m["max_rel_err_pct"],
                            "ckm_mean_rel_pct": f["ckm_mean_rel_pct"],
                            "pmns_mean_rel_pct": f["pmns_mean_rel_pct"],
                        },
                        "flags": {
                            "mass_pass": mass_pass,
                            "flavor_pass": flavor_pass,
                            "joint_pass": joint_pass,
                        },
                        "loss": {
                            "mass_loss": mass_loss,
                            "flavor_loss": flavor_loss,
                            "joint_loss": joint_loss,
                        },
                    }
                    rows.append(row)

                    if best_joint is None or joint_loss < best_joint["loss"]["joint_loss"]:
                        best_joint = row
                    if best_mass is None or mass_loss < best_mass["loss"]["mass_loss"]:
                        best_mass = row
                    if best_flavor is None or flavor_loss < best_flavor["loss"]["flavor_loss"]:
                        best_flavor = row

    rows_sorted_joint = sorted(rows, key=lambda r: r["loss"]["joint_loss"])
    top20 = rows_sorted_joint[:20]

    baseline_signature = {
        "q_up": [0, 9, 14],
        "q_down": [7, 9, 14],
        "q_nu": [0, 1, 2],
        "q_lep": [24, 14, 9],
    }
    baseline_row = None
    for r in rows:
        if (
            r["q_up"] == baseline_signature["q_up"]
            and r["q_down"] == baseline_signature["q_down"]
            and r["q_nu"] == baseline_signature["q_nu"]
            and r["q_lep"] == baseline_signature["q_lep"]
        ):
            baseline_row = r
            break

    verdict = (
        "Q_ASSIGNMENT_AUDIT_HAS_JOINT_PASS_CANDIDATE"
        if joint_pass_count > 0
        else "Q_ASSIGNMENT_AUDIT_NO_JOINT_PASS_IN_LOCAL_PHYSICS_CONSTRAINED_SPACE"
    )
    required_next = (
        "PROMOTE_AUDITED_ASSIGNMENT_AND_RETEST_TRIPLE_SECTOR_GATE"
        if joint_pass_count > 0
        else "REWORK_FLAVOR_OPERATOR_OR_MASS_LINK_BEYOND_ASSIGNMENT_LOCALITY"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "shared_operator_params": params,
        "space_sizes": {k: int(len(v)) for k, v in space.items()},
        "n_total_candidates": int(len(rows)),
        "thresholds": thresholds,
        "counts": {
            "mass_pass_count": int(mass_pass_count),
            "flavor_pass_count": int(flavor_pass_count),
            "joint_pass_count": int(joint_pass_count),
        },
        "baseline_row": baseline_row,
        "best_mass_row": best_mass,
        "best_flavor_row": best_flavor,
        "best_joint_row": best_joint,
        "top20_joint": top20,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    b = baseline_row
    lines = [
        "# RAPORT QW-1943: TOPOLOGICAL Q-ASSIGNMENT AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- n_total_candidates: {len(rows)}",
        f"- space sizes (up/down/nu/lep): {len(space['q_up'])}/{len(space['q_down'])}/{len(space['q_nu'])}/{len(space['q_lep'])}",
        "",
        "## Pass Counts",
        f"- mass_pass_count: {mass_pass_count}",
        f"- flavor_pass_count: {flavor_pass_count}",
        f"- joint_pass_count: {joint_pass_count}",
        "",
        "## Baseline Assignment",
    ]
    if b is not None:
        lines.extend(
            [
                f"- q_up/down/nu/lep: {b['q_up']} / {b['q_down']} / {b['q_nu']} / {b['q_lep']}",
                (
                    "- mass mean/max rel% | CKM/PMNS mean rel%: "
                    f"{b['metrics']['mass_mean_rel_pct']:.3f}/{b['metrics']['mass_max_rel_pct']:.3f} | "
                    f"{b['metrics']['ckm_mean_rel_pct']:.3f}/{b['metrics']['pmns_mean_rel_pct']:.3f}"
                ),
            ]
        )
    else:
        lines.append("- baseline not found in scanned space.")

    bj = best_joint
    lines.extend(
        [
            "",
            "## Best Joint Candidate",
            f"- q_up/down/nu/lep: {bj['q_up']} / {bj['q_down']} / {bj['q_nu']} / {bj['q_lep']}",
            (
                "- mass mean/max rel% | CKM/PMNS mean rel%: "
                f"{bj['metrics']['mass_mean_rel_pct']:.3f}/{bj['metrics']['mass_max_rel_pct']:.3f} | "
                f"{bj['metrics']['ckm_mean_rel_pct']:.3f}/{bj['metrics']['pmns_mean_rel_pct']:.3f}"
            ),
            f"- joint_loss: {bj['loss']['joint_loss']:.4f}",
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1943] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1943] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1943] verdict={verdict} n={len(rows)} joint_pass_count={joint_pass_count}")


if __name__ == "__main__":
    main()

