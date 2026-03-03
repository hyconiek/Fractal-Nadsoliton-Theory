#!/usr/bin/env python3
"""
QW-1950: Internal emergent observer closed-loop informational channel.

Core hypothesis (as requested):
- observer is emergent from nadsoliton and therefore INTERNAL,
- informational transfer must be tested as a closed loop:
    L -> M -> O -> L
  with observer state derived from the same frozen kernel.

Strict rules:
- frozen kernel from QW-1932,
- deterministic maps only,
- no supervised parameter fit to labels.
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
OUT_JSON = ROOT / "report_qw1950_internal_emergent_observer_closed_loop.json"
OUT_MD = ROOT / "RAPORT_QW1950_INTERNAL_EMERGENT_OBSERVER_CLOSED_LOOP.md"


THRESHOLDS = {
    "closed_accuracy_min": 0.70,
    "closed_info_bits_min": 0.20,
    "closed_mean_corr_min": 0.45,
    "closed_acc_gain_vs_open_min": 0.03,
    "closed_info_gain_vs_control_min": 0.10,
    "observer_loop_gain_min": 0.02,
    "observer_state_stability_max": 0.40,
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


def derive_internal_observer_params(kernel: Dict[str, float]) -> Dict[str, float]:
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

    obs_gain_plus = float(0.50 * (1.0 + 0.50 * eps_pol))
    obs_gain_minus = float(0.50 * (1.0 - 0.50 * eps_pol))

    # Internal observer loop from kernel moments.
    tau = float(np.clip(2.0 + 8.0 * mem_short, 2.5, 11.0))
    fb_gain = float(np.clip(0.05 + 0.25 * mem_short + 0.12 * flips, 0.05, 0.55))
    fb_theta = float(np.clip(0.30 + 0.25 * parity_imb, 0.20, 0.60))

    return {
        "eps_pol": eps_pol,
        "noise_sigma": noise_sigma,
        "absorption": absorption,
        "incident_amp": incident_amp,
        "observer_gain_plus": obs_gain_plus,
        "observer_gain_minus": obs_gain_minus,
        "observer_tau": tau,
        "observer_feedback_gain": fb_gain,
        "observer_feedback_theta": fb_theta,
        "flip_rate": flips,
        "parity_imbalance": parity_imb,
        "short_memory_fraction": mem_short,
    }


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


def normalize01(img: np.ndarray) -> np.ndarray:
    mn = float(np.min(img))
    mx = float(np.max(img))
    return (img - mn) / max(mx - mn, 1e-15)


def corr2(a: np.ndarray, b: np.ndarray) -> float:
    av = a.reshape(-1).astype(float)
    bv = b.reshape(-1).astype(float)
    av = av - float(np.mean(av))
    bv = bv - float(np.mean(bv))
    den = math.sqrt(float(np.sum(av * av) * np.sum(bv * bv)))
    if den <= 1e-15:
        return 0.0
    return float(np.sum(av * bv) / den)


def build_surface(class_id: int, idx: int, n: int = 64) -> np.ndarray:
    rng = np.random.default_rng(21_000 + 700 * class_id + idx)
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
    out = np.clip(contrast * base + bias + texture, 0.0, 1.0)
    return out


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


def simulate_closed_loop(
    psf_base: np.ndarray,
    psf_plus: np.ndarray,
    psf_minus: np.ndarray,
    params: Dict[str, float],
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

    x = np.linspace(0.0, 1.0, n, dtype=float)
    xx, yy = np.meshgrid(x, x)
    ux = float(np.cos(params["incident_amp"] * 5.0 + 0.5))
    uy = float(np.sin(params["incident_amp"] * 5.0 - 0.3))

    tau = params["observer_tau"]
    fb_gain = params["observer_feedback_gain"]
    fb_theta = params["observer_feedback_theta"]

    for c in range(len(CLASS_NAMES)):
        for i in range(n_per_class):
            surf = build_surface(c, i, n=n)
            rng = np.random.default_rng(81_000 + 5_000 * c + i)

            z = 0.0
            z_hist = []
            y_hist = []

            obs_final = None
            plus_final = None
            minus_final = None

            for _ in range(t_steps):
                if mode == "closed":
                    fb = fb_gain * math.tanh(z - fb_theta)
                else:
                    fb = 0.0

                phase = 2.0 * np.pi * (0.37 * ux * xx + 0.37 * uy * yy)
                incident = 1.0 + params["incident_amp"] * (1.0 + fb) * np.cos(phase)
                pre = np.clip((1.0 - params["absorption"]) * surf * incident, 0.0, None)

                if mode == "control":
                    base = fftconvolve(pre, psf_base, mode="same")
                    obs = base
                    plus = base
                    minus = base
                else:
                    plus = fftconvolve(pre, psf_plus, mode="same")
                    minus = fftconvolve(pre, psf_minus, mode="same")
                    if mode == "open":
                        obs = params["observer_gain_plus"] * plus + params["observer_gain_minus"] * minus
                    elif mode == "closed":
                        # Internal observer state modulates readout weights.
                        w_dyn = 0.5 + 0.35 * math.tanh(z - 0.5 * fb_theta)
                        obs = w_dyn * plus + (1.0 - w_dyn) * minus
                    else:
                        raise ValueError(f"unknown mode={mode}")

                noise = params["noise_sigma"] * (0.4 + np.sqrt(np.clip(obs, 0.0, None))) * rng.standard_normal(obs.shape)
                obs_noisy = normalize01(np.clip(obs + noise, 0.0, None))

                y_proxy = float(np.mean(obs_noisy) + 0.5 * np.std(obs_noisy))
                z = (1.0 - 1.0 / tau) * z + (1.0 / tau) * y_proxy

                z_hist.append(z)
                y_hist.append(y_proxy)
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

            z_hist_np = np.array(z_hist, dtype=float)
            y_hist_np = np.array(y_hist, dtype=float)
            loop_gain = float(np.mean(np.abs(np.diff(z_hist_np))) / max(np.mean(np.abs(y_hist_np)), 1e-15))
            z_stability = float(np.std(np.diff(z_hist_np)))
            loop_gains.append(loop_gain)
            z_stabilities.append(z_stability)

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
    params = derive_internal_observer_params(kernel)

    d = np.arange(1.0, 25.0, dtype=float)
    k_signed = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    k_abs = np.abs(k_signed)
    sign = np.sign(k_signed)
    sign[sign == 0.0] = 1.0

    w_base = k_abs / max(float(np.sum(k_abs)), 1e-15)
    w_plus = np.clip(k_abs * (1.0 + params["eps_pol"] * sign), 1e-15, None)
    w_minus = np.clip(k_abs * (1.0 - params["eps_pol"] * sign), 1e-15, None)
    w_plus /= max(float(np.sum(w_plus)), 1e-15)
    w_minus /= max(float(np.sum(w_minus)), 1e-15)

    psf_kernel = build_radial_psf(w_base, size=21)
    psf_plus = build_radial_psf(w_plus, size=21)
    psf_minus = build_radial_psf(w_minus, size=21)
    psf_control = gaussian_psf_from_radial(w_base, size=21)

    res_closed = simulate_closed_loop(psf_kernel, psf_plus, psf_minus, params, mode="closed")
    res_open = simulate_closed_loop(psf_kernel, psf_plus, psf_minus, params, mode="open")
    res_control = simulate_closed_loop(psf_control, psf_control, psf_control, params, mode="control")

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
        "closed_acc_gain_vs_open_ge_min": bool(
            metrics["closed_acc_gain_vs_open"] >= THRESHOLDS["closed_acc_gain_vs_open_min"]
        ),
        "closed_info_gain_vs_control_ge_min": bool(
            metrics["closed_info_gain_vs_control"] >= THRESHOLDS["closed_info_gain_vs_control_min"]
        ),
        "observer_loop_gain_ge_min": bool(res_closed["observer_loop_gain"] >= THRESHOLDS["observer_loop_gain_min"]),
        "observer_state_stability_le_max": bool(
            res_closed["observer_state_stability"] <= THRESHOLDS["observer_state_stability_max"]
        ),
    }
    closed_loop_pass = bool(all(flags.values()))

    p1947 = ROOT / "report_qw1947_coupling_completeness_and_mass_operator_frontier.json"
    mass_status = None
    if p1947.exists():
        d1947 = json.loads(p1947.read_text(encoding="utf-8"))
        mass_status = bool(d1947.get("summary", {}).get("strict_pass_exists", False))

    if closed_loop_pass and mass_status is True:
        verdict = "INTERNAL_OBSERVER_CHANNEL_PASS_WITH_MASS_READY"
        required_next = "INTEGRATE_INTERNAL_OBSERVER_GATE_IN_STAGE_C"
    elif closed_loop_pass:
        verdict = "INTERNAL_OBSERVER_CHANNEL_PASS_MASS_STILL_BLOCKED"
        required_next = "KEEP_INTERNAL_OBSERVER_CHANNEL_LOCKED_AND_REPAIR_MASS"
    else:
        verdict = "INTERNAL_OBSERVER_CHANNEL_FAIL"
        required_next = "REWORK_INTERNAL_OBSERVER_FEEDBACK_OPERATOR_AT_KERNEL_LEVEL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "class_names": CLASS_NAMES,
        "thresholds": THRESHOLDS,
        "derived_internal_observer_params": params,
        "metrics": metrics,
        "flags": flags,
        "closed_loop_pass": closed_loop_pass,
        "mass_strict_pass_from_qw1947": mass_status,
        "verdict": verdict,
        "required_next_step": required_next,
        "notes": {
            "closed_loop_model": "L->M->O->L with observer state internal to nadsoliton dynamics",
            "no_fit_statement": "All maps deterministic from frozen kernel; no label fitting.",
        },
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1950: INTERNAL EMERGENT OBSERVER CLOSED LOOP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Closed Loop Metrics",
        f"- accuracy: {res_closed['accuracy']:.4f}",
        f"- mean reconstruction corr: {res_closed['mean_reconstruction_corr']:.4f}",
        f"- info bits: {res_closed['info_bits']:.4f}",
        f"- observer loop gain: {res_closed['observer_loop_gain']:.4f}",
        f"- observer state stability: {res_closed['observer_state_stability']:.4f}",
        "",
        "## Relative Gains",
        f"- accuracy gain (closed-open): {metrics['closed_acc_gain_vs_open']:.4f}",
        f"- info gain (closed-control): {metrics['closed_info_gain_vs_control']:.4f}",
        "",
        "## Flags",
        f"- closed_accuracy_ge_min: {flags['closed_accuracy_ge_min']}",
        f"- closed_info_bits_ge_min: {flags['closed_info_bits_ge_min']}",
        f"- closed_mean_corr_ge_min: {flags['closed_mean_corr_ge_min']}",
        f"- closed_acc_gain_vs_open_ge_min: {flags['closed_acc_gain_vs_open_ge_min']}",
        f"- closed_info_gain_vs_control_ge_min: {flags['closed_info_gain_vs_control_ge_min']}",
        f"- observer_loop_gain_ge_min: {flags['observer_loop_gain_ge_min']}",
        f"- observer_state_stability_le_max: {flags['observer_state_stability_le_max']}",
        f"- closed_loop_pass: {closed_loop_pass}",
        "",
        "## Link to Mass",
        f"- mass_strict_pass_from_qw1947: {mass_status}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1950] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1950] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1950] verdict={verdict}")


if __name__ == "__main__":
    main()

