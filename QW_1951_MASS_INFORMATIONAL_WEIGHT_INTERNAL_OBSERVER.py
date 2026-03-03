#!/usr/bin/env python3
"""
QW-1951: Internal observer closed loop with mass as informational weight.

Extension over QW-1950:
- observer remains internal/emergent from nadsoliton,
- mass hierarchy contributes explicit informational weighting in state update.
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
OUT_JSON = ROOT / "report_qw1951_mass_informational_weight_internal_observer.json"
OUT_MD = ROOT / "RAPORT_QW1951_MASS_INFORMATIONAL_WEIGHT_INTERNAL_OBSERVER.md"


THRESHOLDS = {
    "closed_accuracy_min": 0.70,
    "closed_info_bits_min": 0.20,
    "closed_mean_corr_min": 0.45,
    "closed_acc_gain_vs_open_min": 0.03,
    "closed_info_gain_vs_control_min": 0.10,
    "observer_loop_gain_min": 0.02,
    "observer_state_stability_max": 0.40,
    "mass_info_coupling_min": 0.02,
}


CLASS_NAMES = [
    "horizontal_stripes",
    "vertical_stripes",
    "checker",
    "ring",
    "diagonal",
    "spots",
]

PARTICLE_ORDER = ["Top", "Bottom", "Tau", "Charm", "Muon", "Electron"]


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


def build_radial_psf(weights: np.ndarray, size: int = 21) -> np.ndarray:
    c = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    rr = np.sqrt((yy - c) ** 2 + (xx - c) ** 2)
    ridx = np.clip(np.rint(rr).astype(int), 0, len(weights) - 1)
    psf = weights[ridx]
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


def derive_internal_params(kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0

    flips = float(np.mean((s[1:] != s[:-1]).astype(float)))
    even = float(np.sum(ka[(d.astype(int) % 2) == 0]))
    odd = float(np.sum(ka[(d.astype(int) % 2) == 1]))
    parity_imb = float(abs(even - odd) / max(even + odd, 1e-15))
    short_decay = float(np.mean(np.diff(np.log(np.clip(ka[:8], 1e-15, None))) * -1.0))
    mem_short = float(np.sum(ka[:6]) / max(np.sum(ka), 1e-15))

    eps_pol = float(np.clip(0.12 * flips + 0.08 * parity_imb + 0.03 * abs(math.sin(kernel["phi"])), 0.005, 0.20))
    noise_sigma = float(np.clip(0.015 + 0.05 * flips + 0.03 * parity_imb, 0.01, 0.08))
    absorption = float(np.clip(0.08 + 0.12 * short_decay, 0.03, 0.35))
    incident_amp = float(np.clip(0.04 + 0.08 * abs(math.sin(kernel["omega"])) + 0.04 * abs(ka[0] - ka[1]), 0.03, 0.24))

    tau = float(np.clip(2.0 + 8.0 * mem_short, 2.5, 11.0))
    fb_gain = float(np.clip(0.05 + 0.25 * mem_short + 0.12 * flips, 0.05, 0.55))
    fb_theta = float(np.clip(0.30 + 0.25 * parity_imb, 0.20, 0.60))

    return {
        "eps_pol": eps_pol,
        "noise_sigma": noise_sigma,
        "absorption": absorption,
        "incident_amp": incident_amp,
        "observer_tau": tau,
        "observer_feedback_gain": fb_gain,
        "observer_feedback_theta": fb_theta,
        "flip_rate": flips,
        "parity_imbalance": parity_imb,
        "short_memory_fraction": mem_short,
    }


def load_mass_informational_weights() -> Dict[str, object]:
    d1943 = json.loads((ROOT / "report_qw1943_topological_q_assignment_audit.json").read_text(encoding="utf-8"))
    d1945 = json.loads((ROOT / "report_qw1945_physical_gamma_extraction_hard_mass.json").read_text(encoding="utf-8"))

    mass_q = d1943["baseline_row"]["mass_q"]
    primary_label = d1945["primary_variant"]
    gamma = None
    for v in d1945["variants"]:
        if v["label"] == primary_label:
            gamma = float(v["gamma"])
            break
    if gamma is None:
        raise RuntimeError("Cannot find primary gamma in QW-1945 report")

    pred = {}
    for p in PARTICLE_ORDER:
        q = float(mass_q[p])
        pred[p] = float(173_000.0 * (4.0 ** (-(gamma * q / 4.0))))

    lv = np.array([math.log1p(max(pred[p], 0.0)) for p in PARTICLE_ORDER], dtype=float)
    w = lv / max(float(np.sum(lv)), 1e-15)

    heavy = float(np.sum(w[:3]))
    light = float(np.sum(w[3:]))
    entropy = float(-np.sum(w * np.log2(np.clip(w, 1e-15, None))))
    mass_gain = float(np.clip(0.05 + 0.25 * (heavy - light + 0.5) + 0.02 * entropy, 0.05, 0.55))

    return {
        "mass_q": {k: int(v) for k, v in mass_q.items()},
        "gamma": gamma,
        "pred_mev": pred,
        "informational_weights": {p: float(wi) for p, wi in zip(PARTICLE_ORDER, w)},
        "heavy_weight_sum": heavy,
        "light_weight_sum": light,
        "weight_entropy_bits": entropy,
        "mass_gain": mass_gain,
    }


def build_surface(class_id: int, idx: int, n: int = 64) -> np.ndarray:
    rng = np.random.default_rng(43_000 + 700 * class_id + idx)
    x = np.linspace(-1.0, 1.0, n, dtype=float)
    y = np.linspace(-1.0, 1.0, n, dtype=float)
    xx, yy = np.meshgrid(x, y)
    freq = int(2 + (idx % 3))
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
        for _ in range(4):
            cx = float(rng.uniform(-0.7, 0.7))
            cy = float(rng.uniform(-0.7, 0.7))
            s = float(rng.uniform(0.06, 0.18))
            base += np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * s * s))
        base = base / max(float(np.max(base)), 1e-15)

    sx = int(rng.integers(-4, 5))
    sy = int(rng.integers(-4, 5))
    base = np.roll(base, shift=sx, axis=1)
    base = np.roll(base, shift=sy, axis=0)

    contrast = float(rng.uniform(0.85, 1.15))
    bias = float(rng.uniform(-0.08, 0.08))
    texture = 0.03 * rng.standard_normal(size=base.shape)
    return np.clip(contrast * base + bias + texture, 0.0, 1.0)


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


def mass_info_scalar(img: np.ndarray, weights: np.ndarray) -> float:
    gx = np.diff(img, axis=1)
    gy = np.diff(img, axis=0)
    lap = (
        -4.0 * img
        + np.roll(img, 1, axis=0)
        + np.roll(img, -1, axis=0)
        + np.roll(img, 1, axis=1)
        + np.roll(img, -1, axis=1)
    )
    val = np.array(
        [
            float(np.mean(img)),
            float(np.std(img)),
            float(np.mean(np.abs(gx))),
            float(np.mean(np.abs(gy))),
            float(np.mean(np.abs(lap))),
            float(np.mean(np.abs(gx[:-1, :] * gy[:, :-1]))),
        ],
        dtype=float,
    )
    val = np.clip(val, 0.0, None)
    val = val / max(float(np.sum(val)), 1e-15)
    return float(np.dot(weights, val))


def simulate(
    psf_base: np.ndarray,
    psf_plus: np.ndarray,
    psf_minus: np.ndarray,
    p: Dict[str, float],
    mass_info: Dict[str, object],
    mode: str,
    n_per_class: int = 40,
    n: int = 64,
    t_steps: int = 4,
) -> Dict[str, object]:
    labels: List[int] = []
    feats: List[np.ndarray] = []
    corrs: List[float] = []
    pm_corrs: List[float] = []
    loop_gains: List[float] = []
    z_stabilities: List[float] = []
    mass_couplings: List[float] = []

    w_mass = np.array([mass_info["informational_weights"][k] for k in PARTICLE_ORDER], dtype=float)
    mass_gain = float(mass_info["mass_gain"])
    heavy_bias = float(mass_info["heavy_weight_sum"] - mass_info["light_weight_sum"])

    x = np.linspace(0.0, 1.0, n, dtype=float)
    xx, yy = np.meshgrid(x, x)
    ux = float(np.cos(p["incident_amp"] * 5.0 + 0.5))
    uy = float(np.sin(p["incident_amp"] * 5.0 - 0.3))

    tau = p["observer_tau"]
    fb_gain = p["observer_feedback_gain"]
    fb_theta = p["observer_feedback_theta"]

    for c in range(len(CLASS_NAMES)):
        for i in range(n_per_class):
            surf = build_surface(c, i, n=n)
            rng = np.random.default_rng(91_000 + 5_000 * c + i)
            z = 0.0
            z_hist = []
            y_hist = []
            m_hist = []

            obs_final = None
            plus_final = None
            minus_final = None

            for _ in range(t_steps):
                fb = fb_gain * math.tanh(z - fb_theta) if mode == "closed" else 0.0
                phase = 2.0 * np.pi * (0.37 * ux * xx + 0.37 * uy * yy)
                incident = 1.0 + p["incident_amp"] * (1.0 + fb) * np.cos(phase)
                pre = np.clip((1.0 - p["absorption"]) * surf * incident, 0.0, None)

                if mode == "control":
                    base = fftconvolve(pre, psf_base, mode="same")
                    obs = base
                    plus = base
                    minus = base
                else:
                    plus = fftconvolve(pre, psf_plus, mode="same")
                    minus = fftconvolve(pre, psf_minus, mode="same")
                    if mode == "closed":
                        w_dyn = 0.5 + 0.35 * math.tanh(z - 0.5 * fb_theta) + 0.10 * heavy_bias
                    else:
                        w_dyn = 0.5 + 0.10 * heavy_bias
                    w_dyn = float(np.clip(w_dyn, 0.05, 0.95))
                    obs = w_dyn * plus + (1.0 - w_dyn) * minus

                noise = p["noise_sigma"] * (0.4 + np.sqrt(np.clip(obs, 0.0, None))) * rng.standard_normal(obs.shape)
                obs_noisy = normalize01(np.clip(obs + noise, 0.0, None))

                y_proxy = float(np.mean(obs_noisy) + 0.5 * np.std(obs_noisy))
                m_proxy = mass_info_scalar(obs_noisy, w_mass)
                z = (1.0 - 1.0 / tau) * z + (1.0 / tau) * (y_proxy + mass_gain * m_proxy)

                z_hist.append(z)
                y_hist.append(y_proxy)
                m_hist.append(m_proxy)

                obs_final = obs_noisy
                plus_final = normalize01(np.clip(plus, 0.0, None))
                minus_final = normalize01(np.clip(minus, 0.0, None))

            assert obs_final is not None
            assert plus_final is not None
            assert minus_final is not None

            labels.append(c)
            feats.append(image_features(obs_final))
            corrs.append(corr2(obs_final, surf))
            pm_corrs.append(corr2(plus_final, minus_final))

            z_np = np.array(z_hist, dtype=float)
            y_np = np.array(y_hist, dtype=float)
            m_np = np.array(m_hist, dtype=float)
            loop_gains.append(float(np.mean(np.abs(np.diff(z_np))) / max(np.mean(np.abs(y_np)), 1e-15)))
            z_stabilities.append(float(np.std(np.diff(z_np))))
            mass_couplings.append(float(np.corrcoef(z_np, m_np)[0, 1]) if np.std(z_np) > 0 and np.std(m_np) > 0 else 0.0)

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
        "observer_loop_gain": float(np.mean(loop_gains)),
        "observer_state_stability": float(np.mean(z_stabilities)),
        "mass_info_coupling": float(np.mean(np.abs(mass_couplings))),
        "n_samples": int(len(labels)),
    }


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }
    params = derive_internal_params(kernel)
    mass_info = load_mass_informational_weights()

    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0

    w_base = ka / max(float(np.sum(ka)), 1e-15)
    w_plus = np.clip(ka * (1.0 + params["eps_pol"] * s), 1e-15, None)
    w_minus = np.clip(ka * (1.0 - params["eps_pol"] * s), 1e-15, None)
    w_plus /= max(float(np.sum(w_plus)), 1e-15)
    w_minus /= max(float(np.sum(w_minus)), 1e-15)

    psf_base = build_radial_psf(w_base, size=21)
    psf_plus = build_radial_psf(w_plus, size=21)
    psf_minus = build_radial_psf(w_minus, size=21)
    psf_control = gaussian_psf_from_radial(w_base, size=21)

    res_closed = simulate(psf_base, psf_plus, psf_minus, params, mass_info, mode="closed")
    res_open = simulate(psf_base, psf_plus, psf_minus, params, mass_info, mode="open")
    res_control = simulate(psf_control, psf_control, psf_control, params, mass_info, mode="control")

    metrics = {
        "closed": res_closed,
        "open": res_open,
        "control": res_control,
        "closed_acc_gain_vs_open": float(res_closed["accuracy"] - res_open["accuracy"]),
        "closed_info_gain_vs_control": float(
            (res_closed["info_bits"] - res_control["info_bits"]) / max(res_control["info_bits"], 1e-15)
        ),
    }

    flags = {
        "closed_accuracy_ge_min": bool(res_closed["accuracy"] >= THRESHOLDS["closed_accuracy_min"]),
        "closed_info_bits_ge_min": bool(res_closed["info_bits"] >= THRESHOLDS["closed_info_bits_min"]),
        "closed_mean_corr_ge_min": bool(res_closed["mean_reconstruction_corr"] >= THRESHOLDS["closed_mean_corr_min"]),
        "closed_acc_gain_vs_open_ge_min": bool(metrics["closed_acc_gain_vs_open"] >= THRESHOLDS["closed_acc_gain_vs_open_min"]),
        "closed_info_gain_vs_control_ge_min": bool(
            metrics["closed_info_gain_vs_control"] >= THRESHOLDS["closed_info_gain_vs_control_min"]
        ),
        "observer_loop_gain_ge_min": bool(res_closed["observer_loop_gain"] >= THRESHOLDS["observer_loop_gain_min"]),
        "observer_state_stability_le_max": bool(
            res_closed["observer_state_stability"] <= THRESHOLDS["observer_state_stability_max"]
        ),
        "mass_info_coupling_ge_min": bool(res_closed["mass_info_coupling"] >= THRESHOLDS["mass_info_coupling_min"]),
    }
    pass_all = bool(all(flags.values()))

    if pass_all:
        verdict = "MASS_INFORMATIONAL_INTERNAL_OBSERVER_PASS"
        next_step = "MERGE_MASS_INFORMATIONAL_INTERNAL_OBSERVER_IN_STAGE_C"
    else:
        verdict = "MASS_INFORMATIONAL_INTERNAL_OBSERVER_FAIL"
        next_step = "REWORK_MASS_TO_INFORMATION_MAPPING_AND_CHANNEL_DEGENERACY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "thresholds": THRESHOLDS,
        "derived_internal_params": params,
        "mass_informational_weights": mass_info,
        "metrics": metrics,
        "flags": flags,
        "pass": pass_all,
        "verdict": verdict,
        "required_next_step": next_step,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1951: MASS INFORMATIONAL INTERNAL OBSERVER",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Closed Loop Metrics",
        f"- accuracy: {res_closed['accuracy']:.4f}",
        f"- mean reconstruction corr: {res_closed['mean_reconstruction_corr']:.4f}",
        f"- info bits: {res_closed['info_bits']:.4f}",
        f"- mass-info coupling: {res_closed['mass_info_coupling']:.4f}",
        f"- observer loop gain: {res_closed['observer_loop_gain']:.4f}",
        "",
        "## Relative Gains",
        f"- accuracy gain (closed-open): {metrics['closed_acc_gain_vs_open']:.4f}",
        f"- info gain (closed-control): {metrics['closed_info_gain_vs_control']:.4f}",
        "",
        "## Flags",
        f"- pass: {pass_all}",
        f"- mass_info_coupling_ge_min: {flags['mass_info_coupling_ge_min']}",
        "",
        "## Required Next Step",
        f"- {next_step}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1951] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1951] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1951] verdict={verdict}")


if __name__ == "__main__":
    main()
