#!/usr/bin/env python3
"""
QW-1952: Information-channel de-degeneracy operator (L->M).

Goal:
- break near-degeneracy of +/- channels under frozen kernel,
- keep deterministic maps only (no label fitting),
- quantify whether informational channel gains become non-trivial.
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
OUT_JSON = ROOT / "report_qw1952_information_channel_dedegeneracy_operator.json"
OUT_MD = ROOT / "RAPORT_QW1952_INFORMATION_CHANNEL_DEDEGENERACY_OPERATOR.md"


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


def derive_params(kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0

    flips = float(np.mean((s[1:] != s[:-1]).astype(float)))
    even = float(np.sum(ka[(d.astype(int) % 2) == 0]))
    odd = float(np.sum(ka[(d.astype(int) % 2) == 1]))
    parity_imb = float(abs(even - odd) / max(even + odd, 1e-15))
    decay = float(np.mean(np.diff(np.log(np.clip(ka[:8], 1e-15, None))) * -1.0))

    eps_pol = float(np.clip(0.11 * flips + 0.07 * parity_imb + 0.04 * abs(math.sin(kernel["phi"])), 0.01, 0.25))
    noise_sigma = float(np.clip(0.012 + 0.05 * flips + 0.02 * parity_imb, 0.01, 0.08))
    incident_amp = float(np.clip(0.05 + 0.10 * abs(math.sin(kernel["omega"])) + 0.03 * abs(ka[0] - ka[1]), 0.03, 0.25))
    absorption = float(np.clip(0.06 + 0.14 * decay, 0.03, 0.35))

    # De-degeneracy-specific terms (all deterministic from kernel invariants).
    anis = float(np.clip(0.18 + 0.60 * parity_imb + 0.40 * flips, 0.12, 0.78))
    retard = float(np.clip(0.03 + 0.18 * abs(kernel["phi"]) / math.pi + 0.10 * flips, 0.02, 0.35))
    psi0 = float(np.mod(0.5 * kernel["phi"] + 0.8 * kernel["omega"], 2.0 * math.pi))

    return {
        "eps_pol": eps_pol,
        "noise_sigma": noise_sigma,
        "incident_amp": incident_amp,
        "absorption": absorption,
        "anisotropy_strength": anis,
        "retard_phase": retard,
        "orientation_psi0": psi0,
        "flip_rate": flips,
        "parity_imbalance": parity_imb,
    }


def build_surface(class_id: int, idx: int, n: int = 64) -> np.ndarray:
    rng = np.random.default_rng(111_000 + 1_000 * class_id + idx)
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


def build_anisotropic_psf(radial_w: np.ndarray, anis: float, psi0: float, phase: float, size: int = 21) -> np.ndarray:
    c = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    x = xx - c
    y = yy - c
    rr = np.sqrt(x * x + y * y)
    th = np.arctan2(y, x)
    ridx = np.clip(np.rint(rr).astype(int), 0, len(radial_w) - 1)

    angular = np.cos(2.0 * (th - psi0) + phase)
    psf = radial_w[ridx] * (1.0 + anis * angular)
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


def simulate(
    psf_base: np.ndarray,
    psf_plus: np.ndarray,
    psf_minus: np.ndarray,
    p: Dict[str, float],
    mode: str,
    n_per_class: int = 40,
    n: int = 64,
) -> Dict[str, object]:
    labels: List[int] = []
    feats: List[np.ndarray] = []
    corrs: List[float] = []
    pm_corrs: List[float] = []

    x = np.linspace(0.0, 1.0, n, dtype=float)
    xx, yy = np.meshgrid(x, x)
    ux = float(np.cos(p["incident_amp"] * 5.0 + 0.5))
    uy = float(np.sin(p["incident_amp"] * 5.0 - 0.3))

    for c in range(len(CLASS_NAMES)):
        for i in range(n_per_class):
            surf = build_surface(c, i, n=n)
            rng = np.random.default_rng(131_000 + 1_000 * c + i)

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
                    # Mild nonlinear mixing to preserve both channels.
                    obs = 0.5 * (plus + minus) + 0.25 * np.sqrt(np.clip(plus * minus, 0.0, None))
                else:
                    raise ValueError(f"unknown mode={mode}")

            noise = p["noise_sigma"] * (0.4 + np.sqrt(np.clip(obs, 0.0, None))) * rng.standard_normal(obs.shape)
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


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }
    p = derive_params(kernel)

    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0

    w_base = ka / max(float(np.sum(ka)), 1e-15)
    w_plus = np.clip(ka * (1.0 + p["eps_pol"] * s), 1e-15, None)
    w_minus = np.clip(ka * (1.0 - p["eps_pol"] * s), 1e-15, None)
    w_plus /= max(float(np.sum(w_plus)), 1e-15)
    w_minus /= max(float(np.sum(w_minus)), 1e-15)

    psf_plus = build_anisotropic_psf(w_plus, p["anisotropy_strength"], p["orientation_psi0"], +p["retard_phase"], size=21)
    psf_minus = build_anisotropic_psf(w_minus, p["anisotropy_strength"], p["orientation_psi0"], -p["retard_phase"], size=21)
    psf_single = build_radial_single = build_anisotropic_psf(w_base, 0.0, p["orientation_psi0"], 0.0, size=21)
    psf_control = gaussian_psf_from_radial(w_base, size=21)

    res_dual = simulate(psf_single, psf_plus, psf_minus, p, mode="dual")
    res_single = simulate(psf_single, psf_plus, psf_minus, p, mode="single")
    res_control = simulate(psf_control, psf_control, psf_control, p, mode="control")

    metrics = {
        "dual": res_dual,
        "single": res_single,
        "control": res_control,
        "acc_gain_dual_vs_control": float(res_dual["accuracy"] - res_control["accuracy"]),
        "info_gain_dual_vs_control": float((res_dual["info_bits"] - res_control["info_bits"]) / max(res_control["info_bits"], 1e-15)),
    }

    flags = {
        "dual_accuracy_ge_min": bool(res_dual["accuracy"] >= THRESHOLDS["dual_accuracy_min"]),
        "dual_mean_corr_ge_min": bool(res_dual["mean_reconstruction_corr"] >= THRESHOLDS["dual_mean_corr_min"]),
        "dual_info_bits_ge_min": bool(res_dual["info_bits"] >= THRESHOLDS["dual_info_bits_min"]),
        "acc_gain_vs_control_ge_min": bool(metrics["acc_gain_dual_vs_control"] >= THRESHOLDS["acc_gain_vs_control_min"]),
        "info_gain_vs_control_ge_min": bool(metrics["info_gain_dual_vs_control"] >= THRESHOLDS["info_gain_vs_control_min"]),
        "channel_complementarity_ge_min": bool(res_dual["channel_complementarity"] >= THRESHOLDS["channel_complementarity_min"]),
    }
    dedeg_pass = bool(all(flags.values()))

    verdict = "INFORMATION_CHANNEL_DEDEGENERACY_PASS" if dedeg_pass else "INFORMATION_CHANNEL_DEDEGENERACY_FAIL"
    required_next = (
        "PROMOTE_DEDEGENERACY_OPERATOR_TO_INTERNAL_OBSERVER_TEST"
        if dedeg_pass
        else "REWORK_ANGULAR_SCATTERING_AND_RETRADATION_MAPPING"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "thresholds": THRESHOLDS,
        "derived_params": p,
        "metrics": metrics,
        "flags": flags,
        "dedeg_pass": dedeg_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1952: INFORMATION CHANNEL DE-DEGENERACY OPERATOR",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Dual Metrics",
        f"- accuracy: {res_dual['accuracy']:.4f}",
        f"- mean reconstruction corr: {res_dual['mean_reconstruction_corr']:.4f}",
        f"- info bits: {res_dual['info_bits']:.4f}",
        f"- channel complementarity: {res_dual['channel_complementarity']:.6f}",
        "",
        "## Gains vs Control",
        f"- acc gain: {metrics['acc_gain_dual_vs_control']:.4f}",
        f"- info gain: {metrics['info_gain_dual_vs_control']:.4f}",
        "",
        "## Flags",
        f"- dedeg_pass: {dedeg_pass}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1952] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1952] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1952] verdict={verdict}")


if __name__ == "__main__":
    main()

