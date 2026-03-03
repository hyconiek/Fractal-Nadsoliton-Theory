#!/usr/bin/env python3
"""
QW-1959: Numerical test of minimal EFT term without label fitting.

Protocol:
- freeze kernel from QW-1932,
- freeze baseline operator map from QW-1955,
- sweep EFT coefficients (lambda_I, mu_I) in natural range,
- evaluate primary+stress cohorts, no optimization to labels.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.signal import fftconvolve


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1959_eft_term_numerical_test_without_label_fit.json"
OUT_MD = ROOT / "RAPORT_QW1959_EFT_TERM_NUMERICAL_TEST_WITHOUT_LABEL_FIT.md"


THRESHOLDS = {
    "dual_accuracy_min": 0.70,
    "dual_mean_corr_min": 0.45,
    "dual_info_bits_min": 0.20,
    "acc_gain_vs_control_min": 0.10,
    "info_gain_vs_control_min": 0.10,
    "channel_complementarity_min": 0.02,
}


CLASS_NAMES = [
    "horizontal_stripes",
    "vertical_stripes",
    "checker",
    "ring",
    "diagonal",
    "spots",
]


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def normalize01(img: np.ndarray) -> np.ndarray:
    mn = float(np.min(img))
    mx = float(np.max(img))
    return (img - mn) / max(mx - mn, 1e-15)


def corr2(a: np.ndarray, b: np.ndarray) -> float:
    av = a.reshape(-1).astype(float)
    bv = b.reshape(-1).astype(float)
    av -= float(np.mean(av))
    bv -= float(np.mean(bv))
    den = math.sqrt(float(np.sum(av * av) * np.sum(bv * bv)))
    if den <= 1e-15:
        return 0.0
    return float(np.sum(av * bv) / den)


def build_surface(class_id: int, idx: int, cohort_id: int, n: int = 56) -> np.ndarray:
    rng = np.random.default_rng(201_000 + 2_000 * cohort_id + 1_000 * class_id + idx)
    x = np.linspace(-1.0, 1.0, n, dtype=float)
    y = np.linspace(-1.0, 1.0, n, dtype=float)
    xx, yy = np.meshgrid(x, y)

    freq = int(2 + (idx % 4))
    phase = float(rng.uniform(0.0, 2.0 * np.pi))

    if class_id == 0:
        base = 0.5 + 0.5 * np.cos(np.pi * freq * yy + phase)
    elif class_id == 1:
        base = 0.5 + 0.5 * np.cos(np.pi * freq * xx + phase)
    elif class_id == 2:
        base = np.sign(np.sin(np.pi * freq * xx + phase) * np.sin(np.pi * freq * yy + phase))
        base = 0.5 + 0.5 * base
    elif class_id == 3:
        rr = np.sqrt(xx * xx + yy * yy)
        base = 0.5 + 0.5 * np.cos(2.0 * np.pi * freq * rr + phase)
    elif class_id == 4:
        base = 0.5 + 0.5 * np.cos(np.pi * freq * (xx + yy) / np.sqrt(2.0) + phase)
    else:
        base = np.zeros_like(xx)
        n_spots = 4 + cohort_id
        for _ in range(n_spots):
            cx = float(rng.uniform(-0.7, 0.7))
            cy = float(rng.uniform(-0.7, 0.7))
            s = float(rng.uniform(0.06, 0.18))
            base += np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * s * s))
        base = base / max(float(np.max(base)), 1e-15)

    sx = int(rng.integers(-4, 5))
    sy = int(rng.integers(-4, 5))
    base = np.roll(base, shift=sx, axis=1)
    base = np.roll(base, shift=sy, axis=0)

    contrast = float(rng.uniform(0.82, 1.18))
    bias = float(rng.uniform(-0.10, 0.10))
    texture = (0.03 + 0.01 * cohort_id) * rng.standard_normal(size=base.shape)
    return np.clip(contrast * base + bias + texture, 0.0, 1.0)


def build_repaired_psf(
    radial_w: np.ndarray,
    a2: float,
    b1: float,
    b3: float,
    psi0: float,
    phase_shift: float,
    sign_odd: float,
    size: int = 21,
) -> np.ndarray:
    c = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    x = xx - c
    y = yy - c
    rr = np.sqrt(x * x + y * y)
    th = np.arctan2(y, x)
    ridx = np.clip(np.rint(rr).astype(int), 0, len(radial_w) - 1)

    even_mode = a2 * np.cos(2.0 * (th - psi0))
    odd_mode = sign_odd * (b1 * np.sin((th - psi0) + phase_shift) + b3 * np.sin(3.0 * (th - psi0) + phase_shift))
    ang = even_mode + odd_mode
    psf = radial_w[ridx] * (1.0 + ang)
    psf = np.clip(psf, 1e-15, None)
    psf = psf / max(float(np.sum(psf)), 1e-15)
    return psf


def gaussian_psf_from_radial(radial_w: np.ndarray, size: int = 21) -> np.ndarray:
    r = np.arange(len(radial_w), dtype=float)
    r2_mean = float(np.sum(radial_w * (r**2)) / max(np.sum(radial_w), 1e-15))
    sigma = max(math.sqrt(max(r2_mean, 1e-15) / 2.0), 0.25)
    c = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    rr2 = (yy - c) ** 2 + (xx - c) ** 2
    g = np.exp(-rr2 / (2.0 * sigma * sigma))
    g = g / max(float(np.sum(g)), 1e-15)
    return g


def image_features(img: np.ndarray) -> np.ndarray:
    gx = np.diff(img, axis=1)
    gy = np.diff(img, axis=0)
    lap = (
        -4.0 * img
        + np.roll(img, 1, axis=0)
        + np.roll(img, -1, axis=0)
        + np.roll(img, 1, axis=1)
        + np.roll(img, -1, axis=1)
    )
    hist, _ = np.histogram(img, bins=16, range=(0.0, 1.0), density=True)
    p = hist / max(np.sum(hist), 1e-15)
    ent = float(-np.sum(p * np.log2(np.clip(p, 1e-15, None))))
    return np.array(
        [
            float(np.mean(img)),
            float(np.std(img)),
            float(np.mean(np.abs(gx))),
            float(np.mean(np.abs(gy))),
            float(np.mean(np.abs(lap))),
            float(np.mean(gx * gx)),
            float(np.mean(gy * gy)),
            ent,
        ],
        dtype=float,
    )


def centroid_accuracy(features: np.ndarray, labels: np.ndarray, train_mask: np.ndarray) -> float:
    x_train = features[train_mask]
    y_train = labels[train_mask]
    x_test = features[~train_mask]
    y_test = labels[~train_mask]
    cls = sorted(set(int(v) for v in y_train))
    centroids = {c: np.mean(x_train[y_train == c], axis=0) for c in cls}
    ok = 0
    for x, y in zip(x_test, y_test):
        best_c = None
        best_d = float("inf")
        for c in cls:
            d = float(np.sum((x - centroids[c]) ** 2))
            if d < best_d:
                best_d = d
                best_c = c
        ok += int(best_c == int(y))
    return float(ok / max(len(y_test), 1))


def fisher_ratio(features: np.ndarray, labels: np.ndarray) -> float:
    mu = np.mean(features, axis=0)
    classes = sorted(set(int(v) for v in labels))
    between = 0.0
    within = 0.0
    for c in classes:
        x = features[labels == c]
        muc = np.mean(x, axis=0)
        between += float(len(x)) * float(np.sum((muc - mu) ** 2))
        within += float(np.sum((x - muc) ** 2))
    between /= max(len(classes), 1)
    within /= max(len(features), 1)
    return float(between / max(within, 1e-15))


def info_bits_from_fisher(fr: float) -> float:
    return float(0.5 * math.log2(1.0 + max(fr, 0.0)))


def simulate_cohort(
    psf_base: np.ndarray,
    psf_plus: np.ndarray,
    psf_minus: np.ndarray,
    p: Dict[str, float],
    cohort_id: int,
    mode: str,
    n_per_class: int = 24,
    n: int = 56,
) -> Dict[str, float]:
    labels: List[int] = []
    feats: List[np.ndarray] = []
    corrs: List[float] = []
    pm_corrs: List[float] = []

    x = np.linspace(0.0, 1.0, n, dtype=float)
    xx, yy = np.meshgrid(x, x)
    ux = float(np.cos(p["incident_amp"] * 5.0 + 0.5 + 0.07 * cohort_id))
    uy = float(np.sin(p["incident_amp"] * 5.0 - 0.3 - 0.05 * cohort_id))

    noise_scale = 1.0 + 0.20 * cohort_id

    for c in range(len(CLASS_NAMES)):
        for i in range(n_per_class):
            surf = build_surface(c, i, cohort_id=cohort_id, n=n)
            rng = np.random.default_rng(211_000 + 2_000 * cohort_id + 1_000 * c + i)

            phase = 2.0 * np.pi * (0.37 * ux * xx + 0.37 * uy * yy)
            incident = 1.0 + p["incident_amp"] * np.cos(phase)
            pre = np.clip((1.0 - p["absorption"]) * surf * incident, 0.0, None)

            if mode == "control":
                base = fftconvolve(pre, psf_base, mode="same")
                obs = base
                plus = base
                minus = base
            else:
                plus = fftconvolve(pre, psf_plus, mode="same")
                minus = fftconvolve(pre, psf_minus, mode="same")
                if mode == "single":
                    obs = 0.5 * (plus + minus)
                elif mode == "dual":
                    obs = 0.5 * (plus + minus) + 0.20 * (plus - minus)
                else:
                    raise ValueError(f"unknown mode={mode}")

            noise = noise_scale * p["noise_sigma"] * (0.4 + np.sqrt(np.clip(obs, 0.0, None))) * rng.standard_normal(obs.shape)
            obs_noisy = normalize01(np.clip(obs + noise, 0.0, None))

            labels.append(c)
            feats.append(image_features(obs_noisy))
            corrs.append(corr2(obs_noisy, surf))
            pm_corrs.append(corr2(normalize01(np.clip(plus, 0.0, None)), normalize01(np.clip(minus, 0.0, None))))

    labels_np = np.array(labels, dtype=int)
    feats_np = np.vstack(feats).astype(float)
    train_mask = np.zeros_like(labels_np, dtype=bool)
    for c in range(len(CLASS_NAMES)):
        idx = np.where(labels_np == c)[0]
        train_mask[idx[: int(0.7 * len(idx))]] = True

    acc = centroid_accuracy(feats_np, labels_np, train_mask)
    fr = fisher_ratio(feats_np, labels_np)
    info_bits = info_bits_from_fisher(fr)

    return {
        "accuracy": float(acc),
        "mean_reconstruction_corr": float(np.mean(corrs)),
        "fisher_ratio": float(fr),
        "info_bits": float(info_bits),
        "plus_minus_corr_mean": float(np.mean(pm_corrs)),
        "channel_complementarity": float(1.0 - np.mean(pm_corrs)),
        "n_samples": int(len(labels)),
    }


def score_strict(row: Dict[str, object]) -> float:
    # lower is better; sum misses across primary and stress
    misses = []
    for cohort in ("primary", "stress"):
        m = row["cohorts"][cohort]["dual"]
        c = row["cohorts"][cohort]
        misses.extend(
            [
                max(0.0, THRESHOLDS["dual_accuracy_min"] - m["accuracy"]) / THRESHOLDS["dual_accuracy_min"],
                max(0.0, THRESHOLDS["dual_mean_corr_min"] - m["mean_reconstruction_corr"]) / THRESHOLDS["dual_mean_corr_min"],
                max(0.0, THRESHOLDS["dual_info_bits_min"] - m["info_bits"]) / THRESHOLDS["dual_info_bits_min"],
                max(0.0, THRESHOLDS["acc_gain_vs_control_min"] - c["acc_gain_dual_vs_control"]) / THRESHOLDS["acc_gain_vs_control_min"],
                max(0.0, THRESHOLDS["info_gain_vs_control_min"] - c["info_gain_dual_vs_control"]) / THRESHOLDS["info_gain_vs_control_min"],
                max(0.0, THRESHOLDS["channel_complementarity_min"] - m["channel_complementarity"]) / THRESHOLDS["channel_complementarity_min"],
            ]
        )
    return float(np.sum(misses))


def main() -> None:
    d1932 = load("report_qw1932_physical_reparameterization_eta_scan.json")
    d1955 = load("report_qw1955_nogo_and_minimal_operator_repair.json")
    d1958 = load("report_qw1958_formal_nogo_and_minimal_lagrangian_term.json")

    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }
    p_base = d1955["minimal_repair_params"]
    lam0 = float(d1958["minimal_eft_repair_term"]["deterministic_coefficients"]["lambda_I"])
    mu0 = float(d1958["minimal_eft_repair_term"]["deterministic_coefficients"]["mu_I"])

    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0
    w_base = ka / max(float(np.sum(ka)), 1e-15)
    w_plus_base = np.clip(ka * (1.0 + float(p_base["eps_pol"]) * s), 1e-15, None)
    w_minus_base = np.clip(ka * (1.0 - float(p_base["eps_pol"]) * s), 1e-15, None)
    w_plus_base /= max(float(np.sum(w_plus_base)), 1e-15)
    w_minus_base /= max(float(np.sum(w_minus_base)), 1e-15)

    lambda_grid = np.linspace(max(0.05, 0.35 * lam0), min(1.20, 1.60 * lam0), 7)
    mu_grid = np.linspace(max(0.03, 0.35 * mu0), min(1.20, 1.60 * mu0), 7)

    rows: List[Dict[str, object]] = []

    for lam in lambda_grid:
        for mu in mu_grid:
            # Deterministic EFT->operator map (no label fit).
            r_l = float(lam / max(lam0, 1e-15))
            r_m = float(mu / max(mu0, 1e-15))
            a2 = float(np.clip(float(p_base["a2_even_mode"]) * (0.90 + 0.10 * r_m), 0.05, 0.95))
            b1 = float(np.clip(float(p_base["b1_odd_mode"]) * r_l, 0.01, 1.20))
            b3 = float(np.clip(float(p_base["b3_odd_mode"]) * (0.60 * r_l + 0.40 * r_m), 0.01, 1.20))
            phase = float(np.clip(float(p_base["retard_phase"]) * (0.90 + 0.15 * r_l), 0.005, 1.20))
            psi0 = float(p_base["orientation_psi0"])
            p = {
                "incident_amp": float(p_base["incident_amp"]),
                "absorption": float(p_base["absorption"]),
                "noise_sigma": float(p_base["noise_sigma"]),
                "a2": a2,
                "b1": b1,
                "b3": b3,
                "phase": phase,
                "psi0": psi0,
            }

            psf_plus = build_repaired_psf(w_plus_base, a2, b1, b3, psi0, +phase, +1.0, size=21)
            psf_minus = build_repaired_psf(w_minus_base, a2, b1, b3, psi0, -phase, -1.0, size=21)
            psf_single = build_repaired_psf(w_base, 0.0, 0.0, 0.0, psi0, 0.0, 0.0, size=21)
            psf_control = gaussian_psf_from_radial(w_base, size=21)

            cohorts = {}
            for cname, cid in [("primary", 0), ("stress", 1)]:
                dual = simulate_cohort(psf_single, psf_plus, psf_minus, p, cohort_id=cid, mode="dual")
                control = simulate_cohort(psf_control, psf_control, psf_control, p, cohort_id=cid, mode="control")
                single = simulate_cohort(psf_single, psf_plus, psf_minus, p, cohort_id=cid, mode="single")
                cohorts[cname] = {
                    "dual": dual,
                    "single": single,
                    "control": control,
                    "acc_gain_dual_vs_control": float(dual["accuracy"] - control["accuracy"]),
                    "info_gain_dual_vs_control": float((dual["info_bits"] - control["info_bits"]) / max(control["info_bits"], 1e-15)),
                }

            flags = {}
            for cname in ("primary", "stress"):
                m = cohorts[cname]["dual"]
                c = cohorts[cname]
                flags[cname] = {
                    "dual_accuracy_ge_min": bool(m["accuracy"] >= THRESHOLDS["dual_accuracy_min"]),
                    "dual_mean_corr_ge_min": bool(m["mean_reconstruction_corr"] >= THRESHOLDS["dual_mean_corr_min"]),
                    "dual_info_bits_ge_min": bool(m["info_bits"] >= THRESHOLDS["dual_info_bits_min"]),
                    "acc_gain_vs_control_ge_min": bool(c["acc_gain_dual_vs_control"] >= THRESHOLDS["acc_gain_vs_control_min"]),
                    "info_gain_vs_control_ge_min": bool(c["info_gain_dual_vs_control"] >= THRESHOLDS["info_gain_vs_control_min"]),
                    "channel_complementarity_ge_min": bool(m["channel_complementarity"] >= THRESHOLDS["channel_complementarity_min"]),
                }
                flags[cname]["all_pass"] = bool(all(flags[cname].values()))

            strict_both_pass = bool(flags["primary"]["all_pass"] and flags["stress"]["all_pass"])
            row = {
                "lambda_I": float(lam),
                "mu_I": float(mu),
                "mapped_operator": {"a2_even_mode": a2, "b1_odd_mode": b1, "b3_odd_mode": b3, "retard_phase": phase},
                "cohorts": cohorts,
                "flags": flags,
                "strict_both_pass": strict_both_pass,
            }
            row["strict_loss"] = score_strict(row)
            rows.append(row)

    best = min(rows, key=lambda r: (float(r["strict_loss"]), 0 if r["strict_both_pass"] else 1))
    pass_count = int(sum(1 for r in rows if bool(r["strict_both_pass"])))

    verdict = "EFT_NUMERICAL_TEST_PASS_STRICT" if pass_count > 0 else "EFT_NUMERICAL_TEST_FAIL_STRICT"
    required_next = (
        "PROMOTE_BEST_EFT_POINT_TO_GLOBAL_STAGE_C_PLUS_INFO_GATE"
        if pass_count > 0
        else "REJECT_CURRENT_EFT_TERM_CLASS_OR_ADD_NEXT_ORDER_STRUCTURE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "baseline_source": {
            "q1955_minimal_repair_params": {
                "a2_even_mode": float(p_base["a2_even_mode"]),
                "b1_odd_mode": float(p_base["b1_odd_mode"]),
                "b3_odd_mode": float(p_base["b3_odd_mode"]),
                "retard_phase": float(p_base["retard_phase"]),
            },
            "q1958_deterministic_coeffs": {"lambda_I": lam0, "mu_I": mu0},
        },
        "thresholds": THRESHOLDS,
        "grid": {
            "lambda_I": [float(x) for x in lambda_grid],
            "mu_I": [float(x) for x in mu_grid],
            "n_total": int(len(rows)),
        },
        "summary": {
            "strict_both_pass_count": pass_count,
            "best_point": best,
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    bp = best
    lines = [
        "# RAPORT QW-1959: EFT TERM NUMERICAL TEST (NO LABEL FIT)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- strict_both_pass_count: {pass_count}/{len(rows)}",
        "",
        "## Best Point",
        f"- lambda_I/mu_I: {bp['lambda_I']:.6f}/{bp['mu_I']:.6f}",
        f"- strict_loss: {bp['strict_loss']:.6f}",
        f"- strict_both_pass: {bp['strict_both_pass']}",
        "",
        "## Primary Cohort (best)",
        f"- dual acc/info/complementarity: {bp['cohorts']['primary']['dual']['accuracy']:.4f}/"
        f"{bp['cohorts']['primary']['dual']['info_bits']:.4f}/"
        f"{bp['cohorts']['primary']['dual']['channel_complementarity']:.6f}",
        f"- acc_gain/info_gain vs control: {bp['cohorts']['primary']['acc_gain_dual_vs_control']:.4f}/"
        f"{bp['cohorts']['primary']['info_gain_dual_vs_control']:.4f}",
        "",
        "## Stress Cohort (best)",
        f"- dual acc/info/complementarity: {bp['cohorts']['stress']['dual']['accuracy']:.4f}/"
        f"{bp['cohorts']['stress']['dual']['info_bits']:.4f}/"
        f"{bp['cohorts']['stress']['dual']['channel_complementarity']:.6f}",
        f"- acc_gain/info_gain vs control: {bp['cohorts']['stress']['acc_gain_dual_vs_control']:.4f}/"
        f"{bp['cohorts']['stress']['info_gain_dual_vs_control']:.4f}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1959] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1959] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1959] verdict={verdict} strict_pass={pass_count}/{len(rows)}")


if __name__ == "__main__":
    main()

