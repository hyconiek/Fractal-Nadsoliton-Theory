#!/usr/bin/env python3
"""
QW-1722: Parameter stability under perturbations.

Scope:
1) Mass sector: sensitivity of (l1,l2,l3) and OOS metrics to data + gamma perturbations.
2) Flavor sector: sensitivity of shared operator parameters (a1,a2,a3,a4,delta)
   to perturbed CKM/PMNS targets.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1722_parameter_stability_perturbation.json"
OUT_MD = ROOT / "RAPORT_QW1722_PARAMETER_STABILITY_PERTURBATION.md"


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

PARAM_BOUNDS = {
    "a1": (-2.0, 2.0),
    "a2": (-1.2, 1.2),
    "a3": (-1.2, 1.2),
    "a4": (-0.5, 0.5),
    "delta": (0.0, math.pi),
}
PARAM_KEYS = ["a1", "a2", "a3", "a4", "delta"]


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


def flavor_score(
    params: dict,
    ckm_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
    pmns_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
    ckm_ref: np.ndarray,
    pmns_ref: np.ndarray,
) -> dict:
    u_ckm = np.abs(build_unitary(*ckm_gens, params))
    u_pmns = np.abs(build_unitary(*pmns_gens, params))
    e_ckm = matrix_error(u_ckm, ckm_ref)
    e_pmns = matrix_error(u_pmns, pmns_ref)
    reg = 0.03 * (
        abs(params["a1"])
        + abs(params["a2"])
        + abs(params["a3"])
        + abs(params["a4"])
    )
    score = (
        e_ckm["mean_rel_pct"]
        + e_pmns["mean_rel_pct"]
        + 0.10 * (e_ckm["max_rel_pct"] + e_pmns["max_rel_pct"])
        + reg
    )
    return {
        "score": float(score),
        "ckm_error": e_ckm,
        "pmns_error": e_pmns,
    }


def clamp_param(key: str, value: float) -> float:
    lo, hi = PARAM_BOUNDS[key]
    return float(min(max(value, lo), hi))


def perturb_matrix(ref: np.ndarray, rel_sigma: float, rng: np.random.Generator) -> np.ndarray:
    out = ref * (1.0 + rng.normal(0.0, rel_sigma, size=ref.shape))
    return np.clip(out, 1e-6, 1.0)


def local_refit(
    init_params: dict,
    ckm_target: np.ndarray,
    pmns_target: np.ndarray,
    ckm_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
    pmns_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict:
    best_params = dict(init_params)
    best_eval = flavor_score(best_params, ckm_gens, pmns_gens, ckm_target, pmns_target)

    step0 = {"a1": 0.35, "a2": 0.22, "a3": 0.22, "a4": 0.10, "delta": 0.25}
    for shrink in [1.0, 0.65, 0.45, 0.30]:
        improved = True
        while improved:
            improved = False
            center = dict(best_params)
            for key in PARAM_KEYS:
                for sign in (-1.0, 1.0):
                    cand = dict(center)
                    cand[key] = clamp_param(key, center[key] + sign * step0[key] * shrink)
                    ev = flavor_score(cand, ckm_gens, pmns_gens, ckm_target, pmns_target)
                    if ev["score"] + 1e-12 < best_eval["score"]:
                        best_params = cand
                        best_eval = ev
                        improved = True
    return {"params": best_params, "eval": best_eval}


def robust_scale_ratio(samples: np.ndarray, denom: float, eps: float = 1e-12) -> float:
    q25, q75 = np.quantile(samples, [0.25, 0.75])
    iqr = float(q75 - q25)
    return iqr / (abs(denom) + eps)


def param_summary(values: np.ndarray, keys: list[str], bounds: dict) -> dict:
    out = {}
    for i, k in enumerate(keys):
        arr = values[:, i]
        mean = float(np.mean(arr))
        med = float(np.median(arr))
        std = float(np.std(arr))
        q25, q75 = np.quantile(arr, [0.25, 0.75])
        iqr = float(q75 - q25)
        lo, hi = bounds[k]
        norm_iqr = float(iqr / (hi - lo))
        sign_ref = 1.0 if med >= 0.0 else -1.0
        sign_flip = float(np.mean(np.sign(arr) != sign_ref))
        out[k] = {
            "mean": mean,
            "median": med,
            "std": std,
            "iqr": iqr,
            "normalized_iqr_to_range": norm_iqr,
            "sign_flip_rate": sign_flip,
        }
    return out


def main() -> None:
    rng = np.random.default_rng(1722)

    # ---------- Mass perturbation block ----------
    particles = {
        "Top": {"Q": 0, "mass_mev": 173_000.0, "sigma_mev": 300.0, "sector": 0.0, "gen": 3.0},
        "Bottom": {"Q": 7, "mass_mev": 4_180.0, "sigma_mev": 30.0, "sector": 1.0, "gen": 3.0},
        "Tau": {"Q": 9, "mass_mev": 1_776.9, "sigma_mev": 0.2, "sector": -1.0, "gen": 3.0},
        "Charm": {"Q": 9, "mass_mev": 1_270.0, "sigma_mev": 20.0, "sector": 1.0, "gen": 2.0},
        "Muon": {"Q": 14, "mass_mev": 105.7, "sigma_mev": 0.0002, "sector": -1.0, "gen": 2.0},
        "Electron": {"Q": 24, "mass_mev": 0.511, "sigma_mev": 1e-6, "sector": -1.0, "gen": 1.0},
    }
    names = list(particles.keys())
    idx = {n: i for i, n in enumerate(names)}

    train_names = ["Bottom", "Muon", "Electron"]
    test_names = ["Top", "Tau", "Charm"]
    idx_train = np.array([idx[n] for n in train_names], dtype=int)
    idx_test = np.array([idx[n] for n in test_names], dtype=int)

    q_feat = np.array([particles[n]["Q"] / 24.0 for n in names], dtype=float)
    sec_feat = np.array([particles[n]["sector"] for n in names], dtype=float)
    gen_feat = np.array([particles[n]["gen"] - 2.0 for n in names], dtype=float)
    x_all = np.column_stack([q_feat, sec_feat, gen_feat])

    m_top = particles["Top"]["mass_mev"]
    n_mass_trials = 800
    mass_samples = []

    for _ in range(n_mass_trials):
        gamma = float(rng.uniform(1.49, 1.55))

        obs_mass = {}
        for n in names:
            mu = particles[n]["mass_mev"]
            sigma = particles[n]["sigma_mev"]
            draw = float(rng.normal(mu, sigma))
            obs_mass[n] = max(draw, 1e-12)

        base_pred = {}
        y = []
        for n in names:
            q = particles[n]["Q"]
            base = m_top * (4.0 ** (-(gamma * q / 4.0)))
            base_pred[n] = base
            y.append(math.log(obs_mass[n] / base))
        y = np.array(y, dtype=float)

        x_train = x_all[idx_train]
        y_train = y[idx_train]
        lam, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        errs_test = []
        for i in idx_test:
            n = names[i]
            delta = float(x_all[i] @ lam)
            pred_true = base_pred[n] * math.exp(delta)
            errs_test.append(rel_err_pct(pred_true, particles[n]["mass_mev"]))

        mass_samples.append(
            {
                "gamma": gamma,
                "l1": float(lam[0]),
                "l2": float(lam[1]),
                "l3": float(lam[2]),
                "test_mean_pct": float(np.mean(errs_test)),
                "test_max_pct": float(np.max(errs_test)),
            }
        )

    arr_mass = np.array([[r["l1"], r["l2"], r["l3"]] for r in mass_samples], dtype=float)
    arr_test_mean = np.array([r["test_mean_pct"] for r in mass_samples], dtype=float)
    arr_test_max = np.array([r["test_max_pct"] for r in mass_samples], dtype=float)
    arr_gamma = np.array([r["gamma"] for r in mass_samples], dtype=float)

    mass_param_stats = param_summary(
        arr_mass,
        keys=["l1", "l2", "l3"],
        bounds={"l1": (-2.0, 2.0), "l2": (-2.0, 2.0), "l3": (-2.0, 2.0)},
    )
    mass_metrics = {
        "n_trials": n_mass_trials,
        "gamma_range_tested": [float(np.min(arr_gamma)), float(np.max(arr_gamma))],
        "test_mean_median_pct": float(np.median(arr_test_mean)),
        "test_mean_p90_pct": float(np.quantile(arr_test_mean, 0.90)),
        "test_max_p90_pct": float(np.quantile(arr_test_max, 0.90)),
    }

    mass_thresholds = {
        "test_mean_median_pct_max": 10.0,
        "test_mean_p90_pct_max": 20.0,
        "test_max_p90_pct_max": 40.0,
        "max_sign_flip_rate": 0.15,
        "max_normalized_iqr": 0.25,
    }

    mass_pass = (
        mass_metrics["test_mean_median_pct"] <= mass_thresholds["test_mean_median_pct_max"]
        and mass_metrics["test_mean_p90_pct"] <= mass_thresholds["test_mean_p90_pct_max"]
        and mass_metrics["test_max_p90_pct"] <= mass_thresholds["test_max_p90_pct_max"]
        and max(v["sign_flip_rate"] for v in mass_param_stats.values()) <= mass_thresholds["max_sign_flip_rate"]
        and max(v["normalized_iqr_to_range"] for v in mass_param_stats.values()) <= mass_thresholds["max_normalized_iqr"]
    )

    # ---------- Flavor perturbation block ----------
    q_ckm_l = np.array([0.0, 9.0, 14.0], dtype=float)
    q_ckm_r = np.array([7.0, 9.0, 14.0], dtype=float)
    q_pmns_l = np.array([0.0, 1.0, 2.0], dtype=float)
    q_pmns_r = np.array([24.0, 14.0, 9.0], dtype=float)
    ckm_gens = build_generators(q_ckm_l, q_ckm_r)
    pmns_gens = build_generators(q_pmns_l, q_pmns_r)

    j1720 = load_json(ROOT / "report_qw1720_flavor_extended_operator.json")
    if j1720 and "best" in j1720 and "params" in j1720["best"]:
        init = {k: float(j1720["best"]["params"][k]) for k in PARAM_KEYS}
    else:
        init = {"a1": 0.0, "a2": 0.0, "a3": 0.0, "a4": 0.0, "delta": math.pi / 3.0}

    n_flavor_trials = 120
    flavor_samples = []
    for _ in range(n_flavor_trials):
        target_ckm = perturb_matrix(CKM_REF, rel_sigma=0.004, rng=rng)
        target_pmns = perturb_matrix(PMNS_REF, rel_sigma=0.020, rng=rng)

        start = dict(init)
        for k in PARAM_KEYS:
            span = PARAM_BOUNDS[k][1] - PARAM_BOUNDS[k][0]
            start[k] = clamp_param(k, start[k] + 0.03 * span * float(rng.normal(0.0, 1.0)))

        fit = local_refit(start, target_ckm, target_pmns, ckm_gens, pmns_gens)
        p = fit["params"]

        canonical_eval = flavor_score(p, ckm_gens, pmns_gens, CKM_REF, PMNS_REF)
        flavor_samples.append(
            {
                "a1": p["a1"],
                "a2": p["a2"],
                "a3": p["a3"],
                "a4": p["a4"],
                "delta": p["delta"],
                "fit_score_perturbed": fit["eval"]["score"],
                "canonical_ckm_mean_rel_pct": canonical_eval["ckm_error"]["mean_rel_pct"],
                "canonical_pmns_mean_rel_pct": canonical_eval["pmns_error"]["mean_rel_pct"],
            }
        )

    arr_flavor = np.array([[r[k] for k in PARAM_KEYS] for r in flavor_samples], dtype=float)
    arr_fit = np.array([r["fit_score_perturbed"] for r in flavor_samples], dtype=float)
    arr_ckm_can = np.array([r["canonical_ckm_mean_rel_pct"] for r in flavor_samples], dtype=float)
    arr_pmns_can = np.array([r["canonical_pmns_mean_rel_pct"] for r in flavor_samples], dtype=float)

    flavor_param_stats = param_summary(arr_flavor, PARAM_KEYS, PARAM_BOUNDS)
    flavor_metrics = {
        "n_trials": n_flavor_trials,
        "perturbed_fit_score_median": float(np.median(arr_fit)),
        "canonical_ckm_median_pct": float(np.median(arr_ckm_can)),
        "canonical_pmns_median_pct": float(np.median(arr_pmns_can)),
    }
    flavor_thresholds = {
        "max_sign_flip_rate": 0.10,
        "max_normalized_iqr": 0.18,
        "canonical_ckm_median_pct_max": 25.0,
        "canonical_pmns_median_pct_max": 35.0,
    }
    flavor_pass = (
        max(v["sign_flip_rate"] for v in flavor_param_stats.values()) <= flavor_thresholds["max_sign_flip_rate"]
        and max(v["normalized_iqr_to_range"] for v in flavor_param_stats.values()) <= flavor_thresholds["max_normalized_iqr"]
        and flavor_metrics["canonical_ckm_median_pct"] <= flavor_thresholds["canonical_ckm_median_pct_max"]
        and flavor_metrics["canonical_pmns_median_pct"] <= flavor_thresholds["canonical_pmns_median_pct_max"]
    )

    verdict = (
        "PARAMETERS_STABLE_UNDER_PERTURBATION"
        if (mass_pass and flavor_pass)
        else "PARAMETERS_UNSTABLE_UNDER_PERTURBATION"
    )

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "seed": 1722,
        "mass_block": {
            "thresholds": mass_thresholds,
            "metrics": mass_metrics,
            "param_stats": mass_param_stats,
            "pass": mass_pass,
            "sample_rows_head_40": mass_samples[:40],
        },
        "flavor_block": {
            "thresholds": flavor_thresholds,
            "metrics": flavor_metrics,
            "param_stats": flavor_param_stats,
            "pass": flavor_pass,
            "sample_rows_head_40": flavor_samples[:40],
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1722: PARAMETER STABILITY PERTURBATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Seed: {output['seed']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## 1) Mass block",
        f"- pass: {mass_pass}",
        f"- test_mean median/p90: {mass_metrics['test_mean_median_pct']:.2f}% / {mass_metrics['test_mean_p90_pct']:.2f}%",
        f"- test_max p90: {mass_metrics['test_max_p90_pct']:.2f}%",
        "- param stats (normalized IQR, sign flip):",
    ]
    for k in ["l1", "l2", "l3"]:
        s = mass_param_stats[k]
        lines.append(
            f"  - {k}: norm_iqr={s['normalized_iqr_to_range']:.3f}, sign_flip={s['sign_flip_rate']:.3f}"
        )

    lines.extend(
        [
            "",
            "## 2) Flavor block",
            f"- pass: {flavor_pass}",
            f"- canonical median errors CKM/PMNS: {flavor_metrics['canonical_ckm_median_pct']:.2f}% / {flavor_metrics['canonical_pmns_median_pct']:.2f}%",
            f"- perturbed fit score median: {flavor_metrics['perturbed_fit_score_median']:.2f}",
            "- param stats (normalized IQR, sign flip):",
        ]
    )
    for k in PARAM_KEYS:
        s = flavor_param_stats[k]
        lines.append(
            f"  - {k}: norm_iqr={s['normalized_iqr_to_range']:.3f}, sign_flip={s['sign_flip_rate']:.3f}"
        )

    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegolowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1722] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1722] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
