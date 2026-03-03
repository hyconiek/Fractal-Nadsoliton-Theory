#!/usr/bin/env python3
"""
QW-1920: High-power identifiability experiment for triad interior stability.

Goal:
- answer Stage-B next step from QW-1919:
  RUN_HIGH_POWER_IDENTIFIABILITY_EXPERIMENT_FOR_TRIAD_INTERIOR_STABILITY

Protocol (two-arm):
A) Extended-profile real-data fit
   - regenerate signed-dynamic micromodel profiles with longer distance horizon (dmax=24),
   - estimate shared (omega, phi, beta) without ansatz,
   - quantify boundary pressure / multimodality / conditioning.

B) Synthetic power audit of the same estimator
   - generate interior triads under controlled noise,
   - recover triad with exactly the same estimator,
   - test whether estimator itself can recover interior solutions.

Interpretation:
- if Arm A remains boundary-pushed while Arm B has strong recovery,
  then limitation is data/information geometry, not estimator mechanics.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1920_high_power_identifiability_interior_stability.json"
OUT_MD = ROOT / "RAPORT_QW1920_HIGH_POWER_IDENTIFIABILITY_INTERIOR_STABILITY.md"


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_1920", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_1920"] = mod
    spec.loader.exec_module(mod)
    return mod


def build_extended_profiles(
    n_grid: List[int],
    seeds_per_n: int,
    dmax: int,
) -> Dict[str, object]:
    mod = load_qw1739_module()

    rows = []
    for n in n_grid:
        for sidx in range(seeds_per_n):
            seed = 192000 + 1000 * n + sidx
            cfg = mod.MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)
            theta, q = mod.build_micro_state(cfg, rng)
            w = mod.build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = mod.effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = mod.profile_from_matrix(g, dmax=min(int(dmax), n // 2))
            rows.append(
                {
                    "n_nodes": int(n),
                    "seed": int(seed),
                    "profile": y.tolist(),
                }
            )

    Y = np.array([np.array(r["profile"], dtype=float) for r in rows], dtype=float)
    d = np.arange(1, Y.shape[1] + 1, dtype=float)
    w = np.array([1.0 / max(float(np.var(y)), 1e-5) for y in Y], dtype=float)

    return {
        "rows": rows,
        "d": d,
        "Y": Y,
        "weights": w,
    }


def objective(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    omega, phi, beta = theta
    if omega < 0.02 or omega > 1.90 or beta <= 1e-4 or beta > 2.0:
        return float("inf")

    b = np.cos(omega * d + phi) / (1.0 + beta * d)
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")

    a = (y @ b) / bb
    pred = a[:, None] * b[None, :]
    res = y - pred
    sse = np.sum(weights[:, None] * (res ** 2))
    return float(sse)


def safe_objective(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    f = objective(theta, d, y, weights)
    if not np.isfinite(f):
        return 1e30
    return float(f)


def coordinate_refine(
    start: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    quick: bool = False,
) -> Tuple[Tuple[float, float, float], float]:
    omega, phi, beta = start
    omega = float(np.clip(omega, 0.02, 1.90))
    phi = wrap_phi(phi)
    beta = float(np.clip(beta, 1e-4, 2.0))

    cur = (omega, phi, beta)
    fcur = safe_objective(cur, d, y, w)

    if quick:
        step_plan = [
            (0.20, 0.70, 0.40),
            (0.08, 0.28, 0.15),
            (0.030, 0.10, 0.05),
        ]
    else:
        step_plan = [
            (0.30, 1.00, 0.50),
            (0.12, 0.45, 0.22),
            (0.045, 0.16, 0.08),
            (0.015, 0.07, 0.03),
        ]

    for so, sp, sb in step_plan:
        improved = True
        while improved:
            improved = False

            cand_omega = np.linspace(max(0.02, cur[0] - so), min(1.90, cur[0] + so), 9)
            best_local = (cur, fcur)
            for om in cand_omega:
                th = (float(om), cur[1], cur[2])
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

            cand_phi = np.linspace(cur[1] - sp, cur[1] + sp, 9)
            best_local = (cur, fcur)
            for ph in cand_phi:
                th = (cur[0], wrap_phi(float(ph)), cur[2])
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

            cand_beta = np.linspace(max(1e-4, cur[2] - sb), min(2.0, cur[2] + sb), 9)
            best_local = (cur, fcur)
            for be in cand_beta:
                th = (cur[0], cur[1], float(be))
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

    return cur, fcur


def unique_modes(params: List[Tuple[float, float, float]], tol: Tuple[float, float, float]) -> List[Tuple[float, float, float]]:
    modes: List[Tuple[float, float, float]] = []
    for p in params:
        add = True
        for m in modes:
            d0 = abs(p[0] - m[0]) < tol[0]
            d1 = abs(((p[1] - m[1] + math.pi) % (2.0 * math.pi)) - math.pi) < tol[1]
            d2 = abs(p[2] - m[2]) < tol[2]
            if d0 and d1 and d2:
                add = False
                break
        if add:
            modes.append(p)
    return modes


def bootstrap_boundary(
    theta0: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    n_boot: int,
    seed: int,
) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    n = y.shape[0]

    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        wb = w[idx]
        th, _ = coordinate_refine(theta0, d, yb, wb, quick=True)
        vals.append(th)

    arr = np.array(vals, dtype=float)
    om = arr[:, 0]
    ph = arr[:, 1]
    be = arr[:, 2]

    out = {
        "n_boot": int(n_boot),
        "omega": {
            "median": float(np.median(om)),
            "q05": float(np.quantile(om, 0.05)),
            "q95": float(np.quantile(om, 0.95)),
        },
        "phi": {
            "median": float(np.median(np.unwrap(ph))),
            "q05": float(np.quantile(np.unwrap(ph), 0.05)),
            "q95": float(np.quantile(np.unwrap(ph), 0.95)),
        },
        "beta": {
            "median": float(np.median(be)),
            "q05": float(np.quantile(be, 0.05)),
            "q95": float(np.quantile(be, 0.95)),
        },
        "boundary_fraction": {
            "omega_low": float(np.mean(om <= 0.03)),
            "omega_high": float(np.mean(om >= 1.88)),
            "beta_low": float(np.mean(be <= 0.005)),
            "beta_high": float(np.mean(be >= 1.80)),
        },
    }
    return out


def hessian_numeric(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, w: np.ndarray) -> Dict[str, object]:
    x = np.array(theta, dtype=float)
    h = np.array([0.006, 0.030, 0.018], dtype=float)
    f0 = safe_objective(tuple(x), d, y, w)
    hess = np.zeros((3, 3), dtype=float)

    for i in range(3):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h[i]
        xm[i] -= h[i]
        if i == 1:
            xp[i] = wrap_phi(float(xp[i]))
            xm[i] = wrap_phi(float(xm[i]))
        fp = safe_objective(tuple(xp), d, y, w)
        fm = safe_objective(tuple(xm), d, y, w)
        hess[i, i] = (fp - 2.0 * f0 + fm) / (h[i] ** 2)

    for i in range(3):
        for j in range(i + 1, 3):
            xpp = x.copy()
            xpm = x.copy()
            xmp = x.copy()
            xmm = x.copy()
            xpp[i] += h[i]
            xpp[j] += h[j]
            xpm[i] += h[i]
            xpm[j] -= h[j]
            xmp[i] -= h[i]
            xmp[j] += h[j]
            xmm[i] -= h[i]
            xmm[j] -= h[j]
            if i == 1:
                xpp[i] = wrap_phi(float(xpp[i]))
                xpm[i] = wrap_phi(float(xpm[i]))
                xmp[i] = wrap_phi(float(xmp[i]))
                xmm[i] = wrap_phi(float(xmm[i]))
            if j == 1:
                xpp[j] = wrap_phi(float(xpp[j]))
                xpm[j] = wrap_phi(float(xpm[j]))
                xmp[j] = wrap_phi(float(xmp[j]))
                xmm[j] = wrap_phi(float(xmm[j]))
            fpp = safe_objective(tuple(xpp), d, y, w)
            fpm = safe_objective(tuple(xpm), d, y, w)
            fmp = safe_objective(tuple(xmp), d, y, w)
            fmm = safe_objective(tuple(xmm), d, y, w)
            val = (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j])
            hess[i, j] = val
            hess[j, i] = val

    hsym = 0.5 * (hess + hess.T)
    hsym = np.nan_to_num(hsym, nan=1e6, posinf=1e6, neginf=-1e6)
    eig = np.linalg.eigvalsh(hsym)
    min_eig = float(np.min(eig))
    reg = 0.0
    if min_eig <= 1e-8:
        reg = (1e-8 - min_eig) + 1e-6
    hreg = hsym + reg * np.eye(3)
    cond = float(np.linalg.cond(hreg))

    return {
        "matrix": hsym.tolist(),
        "regularization_added": float(reg),
        "eigenvalues": [float(v) for v in eig],
        "condition_number": cond,
    }


def estimate_global_triad(d: np.ndarray, y: np.ndarray, w: np.ndarray) -> Dict[str, object]:
    starts: List[Tuple[float, float, float]] = [
        (0.118, 0.475, 0.300),
        (0.088, 0.890, 0.400),
        (0.75, 0.50, 0.02),
        (1.10, -0.30, 0.20),
        (0.20, -0.90, 0.80),
    ]

    rng = np.random.default_rng(192001)
    for _ in range(22):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
            )
        )

    sols: List[Tuple[Tuple[float, float, float], float]] = []
    for st in starts:
        th, f = coordinate_refine(st, d, y, w, quick=False)
        sols.append((th, f))

    sols = sorted(sols, key=lambda x: x[1])
    theta_best = sols[0][0]
    f_best = float(sols[0][1])

    params = [s[0] for s in sols]
    modes = unique_modes(params, tol=(0.05, 0.18, 0.10))

    hdiag = hessian_numeric(theta_best, d, y, w)
    bstrap = bootstrap_boundary(theta_best, d, y, w, n_boot=100, seed=192002)

    return {
        "optimum": {
            "omega": float(theta_best[0]),
            "phi": float(theta_best[1]),
            "beta": float(theta_best[2]),
            "objective_sse_weighted": f_best,
        },
        "multistart": {
            "n_starts": len(starts),
            "n_unique_modes": len(modes),
            "unique_modes": [[float(x) for x in m] for m in modes],
            "top_solutions": [
                {"theta": [float(x) for x in th], "objective": float(f)}
                for th, f in sols[:12]
            ],
        },
        "hessian": hdiag,
        "bootstrap": bstrap,
    }


def estimate_one_profile(y: np.ndarray, d: np.ndarray, rng: np.random.Generator) -> Tuple[float, float, float]:
    y2 = y.reshape(1, -1)
    w2 = np.array([1.0], dtype=float)

    starts: List[Tuple[float, float, float]] = [
        (0.118, 0.475, 0.300),
        (0.088, 0.890, 0.400),
        (0.75, 0.50, 0.02),
    ]
    for _ in range(6):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
            )
        )

    sols = []
    for st in starts:
        th, f = coordinate_refine(st, d, y2, w2, quick=False)
        sols.append((th, f))
    sols = sorted(sols, key=lambda x: x[1])

    th = sols[0][0]
    return float(th[0]), float(th[1]), float(th[2])


def circular_abs_diff_mod_pi(a: float, b: float) -> float:
    # In y = a*cos(omega*d + phi), nuisance amplitude a can flip sign.
    # Therefore phase is identifiable modulo pi, not modulo 2*pi.
    d = (a - b + 0.5 * math.pi) % math.pi - 0.5 * math.pi
    return float(abs(d))


def synthetic_power_audit(d: np.ndarray, sigma_noise: float, n_rep: int) -> Dict[str, object]:
    rng = np.random.default_rng(192010)

    rows = []
    for _ in range(n_rep):
        omega_true = float(np.clip(rng.normal(math.pi / 4.0, 0.08), 0.35, 1.10))
        if rng.random() < 0.5:
            phi_center = math.pi / 6.0
        else:
            phi_center = -math.pi / 6.0
        phi_true = float(np.clip(rng.normal(phi_center, 0.18), -1.4, 1.4))
        beta_true = float(np.clip(rng.normal(0.020, 0.010), 0.006, 0.060))

        amp_true = float(np.clip(rng.normal(1.0, 0.15), 0.4, 1.8))
        y_clean = amp_true * np.cos(omega_true * d + phi_true) / (1.0 + beta_true * d)
        y_obs = y_clean + rng.normal(0.0, sigma_noise, size=len(d))

        omega_hat, phi_hat, beta_hat = estimate_one_profile(y_obs, d, rng)

        eom = float(abs(omega_hat - omega_true))
        eph = circular_abs_diff_mod_pi(phi_hat, phi_true)
        ebe = float(abs(beta_hat - beta_true))

        rows.append(
            {
                "true": {"omega": omega_true, "phi": phi_true, "beta": beta_true},
                "hat": {"omega": omega_hat, "phi": phi_hat, "beta": beta_hat},
                "errors": {"omega_abs": eom, "phi_abs": eph, "beta_abs": ebe},
                "hat_boundary": bool(beta_hat >= 1.80 or beta_hat <= 0.005),
                "joint_hit": bool(eom <= 0.10 and eph <= 0.35 and ebe <= 0.020),
            }
        )

    eom = np.array([r["errors"]["omega_abs"] for r in rows], dtype=float)
    eph = np.array([r["errors"]["phi_abs"] for r in rows], dtype=float)
    ebe = np.array([r["errors"]["beta_abs"] for r in rows], dtype=float)
    hit = np.array([1.0 if r["joint_hit"] else 0.0 for r in rows], dtype=float)
    bnd = np.array([1.0 if r["hat_boundary"] else 0.0 for r in rows], dtype=float)

    summary = {
        "n_rep": int(n_rep),
        "sigma_noise": float(sigma_noise),
        "median_abs_error": {
            "omega": float(np.median(eom)),
            "phi": float(np.median(eph)),
            "beta": float(np.median(ebe)),
        },
        "joint_hit_rate": float(np.mean(hit)),
        "boundary_estimate_rate": float(np.mean(bnd)),
        "omega_tol_rate_0p10": float(np.mean(eom <= 0.10)),
        "phi_tol_rate_0p35": float(np.mean(eph <= 0.35)),
        "beta_tol_rate_0p02": float(np.mean(ebe <= 0.020)),
    }

    return {
        "summary": summary,
        "rows_head": rows[:80],
    }


def main() -> None:
    cfg = {
        "n_grid": [64, 96, 128],
        "seeds_per_n": 14,
        "dmax": 24,
        "synthetic_n_rep": 180,
    }

    built = build_extended_profiles(
        n_grid=list(cfg["n_grid"]),
        seeds_per_n=int(cfg["seeds_per_n"]),
        dmax=int(cfg["dmax"]),
    )

    d = built["d"]
    Y = built["Y"]
    W = built["weights"]

    real_est = estimate_global_triad(d, Y, W)

    # Noise scale for synthetic audit from residual spread of real rows.
    om = float(real_est["optimum"]["omega"])
    ph = float(real_est["optimum"]["phi"])
    be = float(real_est["optimum"]["beta"])
    b = np.cos(om * d + ph) / (1.0 + be * d)
    bb = float(np.dot(b, b))
    a = (Y @ b) / max(bb, 1e-12)
    pred = a[:, None] * b[None, :]
    residual = Y - pred
    sigma_noise = float(np.median(np.std(residual, axis=1)))
    sigma_noise = float(np.clip(sigma_noise, 0.008, 0.050))

    synth = synthetic_power_audit(d, sigma_noise=sigma_noise, n_rep=int(cfg["synthetic_n_rep"]))

    bf = real_est["bootstrap"]["boundary_fraction"]
    real_beta_boundary = float(bf["beta_high"])
    real_cond = float(real_est["hessian"]["condition_number"])
    real_modes = int(real_est["multistart"]["n_unique_modes"])
    beta_hat = float(real_est["optimum"]["beta"])

    synth_joint = float(synth["summary"]["joint_hit_rate"])
    synth_boundary = float(synth["summary"]["boundary_estimate_rate"])

    estimator_power_strong = bool(synth_joint >= 0.65 and synth_boundary <= 0.15)
    interior_stable_real = bool(beta_hat < 1.20 and real_beta_boundary <= 0.40 and real_cond <= 1e10 and real_modes <= 4)

    if estimator_power_strong and interior_stable_real:
        verdict = "HIGH_POWER_IDENTIFIABILITY_INTERIOR_PASS"
        required_next = "UPDATE_STAGE_B_GATE_WITH_INTERIOR_STABILITY_PASS"
    elif estimator_power_strong and not interior_stable_real:
        verdict = "HIGH_POWER_IDENTIFIABILITY_DATA_LIMITED_INTERIOR_NOT_ESTABLISHED"
        required_next = "DESIGN_ORTHOGONAL_BETA_OBSERVABLE_ON_BLIND_EXTERNAL_WITH_INTERVENTION_PROTOCOL"
    else:
        verdict = "HIGH_POWER_IDENTIFIABILITY_METHOD_OR_DATA_WEAK"
        required_next = "STRENGTHEN_ESTIMATOR_AND_REPEAT_POWER_AUDIT"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": cfg,
        "real_arm": {
            "n_profiles": int(Y.shape[0]),
            "n_points_per_profile": int(Y.shape[1]),
            "estimate": real_est,
        },
        "synthetic_power_arm": synth,
        "derived_flags": {
            "estimator_power_strong": bool(estimator_power_strong),
            "interior_stable_real": bool(interior_stable_real),
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1920: HIGH-POWER IDENTIFIABILITY INTERIOR STABILITY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Required next step: `{required_next}`",
        "",
        "## Real Arm (extended profiles)",
        f"- n_profiles: {Y.shape[0]}",
        f"- n_points_per_profile: {Y.shape[1]}",
        f"- optimum: omega={om:.6f}, phi={ph:.6f}, beta={be:.6f}",
        f"- unique_modes: {real_modes}",
        f"- hessian_cond: {real_cond:.3e}",
        f"- boundary beta_high fraction: {real_beta_boundary:.3f}",
        "",
        "## Synthetic Power Arm",
        f"- n_rep: {synth['summary']['n_rep']}",
        f"- sigma_noise: {synth['summary']['sigma_noise']:.4f}",
        f"- joint_hit_rate: {synth_joint:.3f}",
        f"- boundary_estimate_rate: {synth_boundary:.3f}",
        f"- median_abs_error omega/phi/beta: "
        f"{synth['summary']['median_abs_error']['omega']:.4f} / "
        f"{synth['summary']['median_abs_error']['phi']:.4f} / "
        f"{synth['summary']['median_abs_error']['beta']:.4f}",
        "",
        "## Flags",
        f"- estimator_power_strong: {estimator_power_strong}",
        f"- interior_stable_real: {interior_stable_real}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1920] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1920] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1920] verdict={verdict}")


if __name__ == "__main__":
    main()
