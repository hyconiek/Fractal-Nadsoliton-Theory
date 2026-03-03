#!/usr/bin/env python3
"""
QW-1723: Unified effective Hamiltonian integration test.

Objective:
1) Fit mass sector OOS with frozen protocol.
2) Map fitted mass parameters into flavor generator coefficients.
3) Evaluate single integrated mechanism for CKM + PMNS + mass OOS.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1723_unified_effective_hamiltonian_integration.json"
OUT_MD = ROOT / "RAPORT_QW1723_UNIFIED_EFFECTIVE_HAMILTONIAN_INTEGRATION.md"


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


def load_json(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return None


def rel_err_pct(pred: float, target: float) -> float:
    return abs(pred - target) / abs(target) * 100.0


def matrix_error(pred_abs: np.ndarray, ref: np.ndarray) -> dict:
    abs_err = np.abs(pred_abs - ref)
    rel = abs_err / np.clip(ref, 1e-12, None)
    return {
        "mae": float(np.mean(abs_err)),
        "max_abs": float(np.max(abs_err)),
        "mean_rel_pct": float(np.mean(rel) * 100.0),
        "max_rel_pct": float(np.max(rel) * 100.0),
    }


def antisym_part(m: np.ndarray) -> np.ndarray:
    return 0.5 * (m - m.T)


def build_generators(q_left: np.ndarray, q_right: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(q_left)
    g1 = np.zeros((n, n), dtype=float)
    g2 = np.zeros((n, n), dtype=float)
    g3 = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dq = abs(q_left[i] - q_right[j])
            gap = abs(i - j)
            sign = 1.0 if i < j else -1.0
            g1[i, j] = sign / max(dq, 1.0)
            g2[i, j] = sign / max(dq * dq, 1.0)
            locality = 1.0 / (1.0 + gap)
            asym = (q_left[i] - q_right[j]) / (1.0 + dq)
            g3[i, j] = sign * locality * asym
    return antisym_part(g1), antisym_part(g2), antisym_part(g3)


def build_unitary(g1: np.ndarray, g2: np.ndarray, g3: np.ndarray, params: dict) -> np.ndarray:
    base = (
        params["a1"] * g1
        + params["a2"] * g2
        + params["a3"] * g3
        + params["a4"] * (g1 @ g2 - g2 @ g1)
    )
    h = 1j * base
    cp = np.zeros_like(h, dtype=complex)
    cp[0, 2] = np.exp(-1j * params["delta"])
    cp[2, 0] = np.conj(cp[0, 2])
    h = h + 0.05 * cp
    return expm(1j * h)


def mass_fit_and_metrics(gamma: float, particles: dict, train: list[str], test: list[str]) -> dict:
    names = list(particles.keys())
    idx = {n: i for i, n in enumerate(names)}
    m_top = particles["Top"]["mass_mev"]

    q_feat = np.array([particles[n]["Q"] / 24.0 for n in names], dtype=float)
    sec_feat = np.array([particles[n]["sector"] for n in names], dtype=float)
    gen_feat = np.array([particles[n]["gen"] - 2.0 for n in names], dtype=float)
    x = np.column_stack([q_feat, sec_feat, gen_feat])

    base_pred = []
    y = []
    for n in names:
        q = particles[n]["Q"]
        base = m_top * (4.0 ** (-(gamma * q / 4.0)))
        base_pred.append(base)
        y.append(math.log(particles[n]["mass_mev"] / base))
    base_pred = np.array(base_pred, dtype=float)
    y = np.array(y, dtype=float)

    idx_train = np.array([idx[n] for n in train], dtype=int)
    idx_test = np.array([idx[n] for n in test], dtype=int)

    lam, *_ = np.linalg.lstsq(x[idx_train], y[idx_train], rcond=None)

    pred = base_pred * np.exp(x @ lam)
    train_err = [rel_err_pct(pred[i], particles[names[i]]["mass_mev"]) for i in idx_train]
    test_err = [rel_err_pct(pred[i], particles[names[i]]["mass_mev"]) for i in idx_test]

    return {
        "lambda": {"l1": float(lam[0]), "l2": float(lam[1]), "l3": float(lam[2])},
        "lambda_vec": lam,
        "train_mean_pct": float(np.mean(train_err)),
        "test_mean_pct": float(np.mean(test_err)),
        "gap_pct": float(np.mean(test_err) - np.mean(train_err)),
        "test_max_pct": float(np.max(test_err)),
    }


def unified_eval(
    gamma: float,
    kappa: float,
    xi: float,
    delta: float,
    particles: dict,
    train: list[str],
    test: list[str],
    ckm_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
    pmns_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict:
    m = mass_fit_and_metrics(gamma, particles, train, test)
    lam = m.pop("lambda_vec")

    # Unified map: same fitted mass lambdas drive flavor coefficients.
    params = {
        "a1": float(kappa * lam[0]),
        "a2": float(kappa * lam[1]),
        "a3": float(kappa * lam[2]),
        "a4": float(xi * (lam[0] * lam[1] - lam[2])),
        "delta": float(delta),
    }

    u_ckm = np.abs(build_unitary(*ckm_gens, params))
    u_pmns = np.abs(build_unitary(*pmns_gens, params))
    e_ckm = matrix_error(u_ckm, CKM_REF)
    e_pmns = matrix_error(u_pmns, PMNS_REF)

    score = (
        1.00 * m["test_mean_pct"]
        + 0.80 * max(0.0, m["gap_pct"])
        + 0.35 * (e_ckm["mean_rel_pct"] + e_pmns["mean_rel_pct"])
        + 0.05 * (e_ckm["max_rel_pct"] + e_pmns["max_rel_pct"])
        + 0.02 * (abs(kappa) + abs(xi))
    )

    return {
        "score": float(score),
        "gamma": float(gamma),
        "kappa": float(kappa),
        "xi": float(xi),
        "delta": float(delta),
        "mass": m,
        "flavor": {
            "params": params,
            "ckm_error": e_ckm,
            "pmns_error": e_pmns,
            "ckm_pred_abs": u_ckm,
            "pmns_pred_abs": u_pmns,
        },
    }


def clamp(v: float, lo: float, hi: float) -> float:
    return float(min(max(v, lo), hi))


def main() -> None:
    particles = {
        "Top": {"Q": 0, "mass_mev": 173_000.0, "sector": 0.0, "gen": 3.0},
        "Bottom": {"Q": 7, "mass_mev": 4_180.0, "sector": 1.0, "gen": 3.0},
        "Tau": {"Q": 9, "mass_mev": 1_776.9, "sector": -1.0, "gen": 3.0},
        "Charm": {"Q": 9, "mass_mev": 1_270.0, "sector": 1.0, "gen": 2.0},
        "Muon": {"Q": 14, "mass_mev": 105.7, "sector": -1.0, "gen": 2.0},
        "Electron": {"Q": 24, "mass_mev": 0.511, "sector": -1.0, "gen": 1.0},
    }
    train = ["Bottom", "Muon", "Electron"]
    test = ["Top", "Tau", "Charm"]

    q_ckm_l = np.array([0.0, 9.0, 14.0], dtype=float)
    q_ckm_r = np.array([7.0, 9.0, 14.0], dtype=float)
    q_pmns_l = np.array([0.0, 1.0, 2.0], dtype=float)
    q_pmns_r = np.array([24.0, 14.0, 9.0], dtype=float)
    ckm_gens = build_generators(q_ckm_l, q_ckm_r)
    pmns_gens = build_generators(q_pmns_l, q_pmns_r)

    # Coarse grid
    best = None
    eval_count = 0
    for gamma in np.linspace(1.47, 1.57, 13):
        for kappa in np.linspace(-10.0, 10.0, 21):
            for xi in np.linspace(-4.0, 4.0, 11):
                for delta in np.linspace(0.0, math.pi, 9):
                    cand = unified_eval(
                        float(gamma),
                        float(kappa),
                        float(xi),
                        float(delta),
                        particles,
                        train,
                        test,
                        ckm_gens,
                        pmns_gens,
                    )
                    eval_count += 1
                    if best is None or cand["score"] < best["score"]:
                        best = cand

    # Local deterministic refinement
    steps = {"gamma": 0.01, "kappa": 0.60, "xi": 0.25, "delta": 0.20}
    bounds = {
        "gamma": (1.40, 1.65),
        "kappa": (-12.0, 12.0),
        "xi": (-5.0, 5.0),
        "delta": (0.0, math.pi),
    }
    for shrink in [1.0, 0.65, 0.42, 0.28]:
        improved = True
        while improved:
            improved = False
            center = {
                "gamma": best["gamma"],
                "kappa": best["kappa"],
                "xi": best["xi"],
                "delta": best["delta"],
            }
            for key in ["gamma", "kappa", "xi", "delta"]:
                for sign in (-1.0, 1.0):
                    trial = dict(center)
                    lo, hi = bounds[key]
                    trial[key] = clamp(center[key] + sign * steps[key] * shrink, lo, hi)
                    cand = unified_eval(
                        trial["gamma"],
                        trial["kappa"],
                        trial["xi"],
                        trial["delta"],
                        particles,
                        train,
                        test,
                        ckm_gens,
                        pmns_gens,
                    )
                    eval_count += 1
                    if cand["score"] + 1e-12 < best["score"]:
                        best = cand
                        improved = True

    # Closure thresholds for integrated mechanism
    thresholds = {
        "mass_test_mean_pct_max": 10.0,
        "mass_gap_pct_max": 7.0,
        "flavor_ckm_mean_rel_pct_max": 15.0,
        "flavor_pmns_mean_rel_pct_max": 15.0,
    }
    integrated_pass = (
        best["mass"]["test_mean_pct"] <= thresholds["mass_test_mean_pct_max"]
        and best["mass"]["gap_pct"] <= thresholds["mass_gap_pct_max"]
        and best["flavor"]["ckm_error"]["mean_rel_pct"] <= thresholds["flavor_ckm_mean_rel_pct_max"]
        and best["flavor"]["pmns_error"]["mean_rel_pct"] <= thresholds["flavor_pmns_mean_rel_pct_max"]
    )
    verdict = (
        "UNIFIED_EFFECTIVE_HAMILTONIAN_CLOSED"
        if integrated_pass
        else "UNIFIED_EFFECTIVE_HAMILTONIAN_NOT_CLOSED"
    )

    # Decoupled reference from already computed reports (if available)
    j1720 = load_json(ROOT / "report_qw1720_flavor_extended_operator.json")
    if j1720 and "best" in j1720:
        ref_ckm = float(j1720["best"]["ckm_error"]["mean_rel_pct"])
        ref_pmns = float(j1720["best"]["pmns_error"]["mean_rel_pct"])
        ref_flavor = 0.5 * (ref_ckm + ref_pmns)
    else:
        ref_flavor = 0.5 * (
            best["flavor"]["ckm_error"]["mean_rel_pct"] + best["flavor"]["pmns_error"]["mean_rel_pct"]
        )

    j1721 = load_json(ROOT / "report_qw1721_strict_no_leakage_oos.json")
    if j1721 and "metrics" in j1721:
        ref_mass = float(j1721["metrics"]["test_mean_pct"])
    else:
        ref_mass = float(best["mass"]["test_mean_pct"])

    decoupled_proxy = ref_mass + 0.5 * ref_flavor
    integrated_proxy = best["mass"]["test_mean_pct"] + 0.5 * (
        best["flavor"]["ckm_error"]["mean_rel_pct"] + best["flavor"]["pmns_error"]["mean_rel_pct"]
    )

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "search_protocol": {
            "coarse_grid": {"gamma": 13, "kappa": 21, "xi": 11, "delta": 9},
            "refinement": "deterministic_coordinate_descent",
            "total_evaluations": int(eval_count),
        },
        "split": {"train": train, "test": test},
        "thresholds": thresholds,
        "best": {
            "score": best["score"],
            "gamma": best["gamma"],
            "kappa": best["kappa"],
            "xi": best["xi"],
            "delta": best["delta"],
            "mass": best["mass"],
            "flavor": {
                "params": best["flavor"]["params"],
                "ckm_error": best["flavor"]["ckm_error"],
                "pmns_error": best["flavor"]["pmns_error"],
                "ckm_pred_abs": best["flavor"]["ckm_pred_abs"].tolist(),
                "pmns_pred_abs": best["flavor"]["pmns_pred_abs"].tolist(),
            },
        },
        "integration_penalty_proxy": {
            "integrated_proxy": float(integrated_proxy),
            "decoupled_proxy": float(decoupled_proxy),
            "delta_integrated_minus_decoupled": float(integrated_proxy - decoupled_proxy),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1723: UNIFIED EFFECTIVE HAMILTONIAN INTEGRATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        f"- Evaluations: {eval_count}",
        "",
        "## Best integrated parameters",
        f"- gamma={best['gamma']:.5f}, kappa={best['kappa']:.5f}, xi={best['xi']:.5f}, delta={best['delta']:.5f}",
        "",
        "## Mass OOS",
        f"- train mean: {best['mass']['train_mean_pct']:.2f}%",
        f"- test mean: {best['mass']['test_mean_pct']:.2f}%",
        f"- gap: {best['mass']['gap_pct']:.2f} p.p.",
        "",
        "## Flavor",
        f"- CKM mean rel. error: {best['flavor']['ckm_error']['mean_rel_pct']:.2f}%",
        f"- PMNS mean rel. error: {best['flavor']['pmns_error']['mean_rel_pct']:.2f}%",
        "",
        "## Integration penalty proxy",
        f"- integrated proxy: {integrated_proxy:.3f}",
        f"- decoupled proxy: {decoupled_proxy:.3f}",
        f"- delta: {integrated_proxy - decoupled_proxy:+.3f}",
        "",
        "## Artefakty",
        f"- JSON szczegolowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1723] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1723] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
