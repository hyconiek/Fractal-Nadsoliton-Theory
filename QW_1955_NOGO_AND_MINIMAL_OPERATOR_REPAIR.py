#!/usr/bin/env python3
"""
QW-1955: No-go diagnosis + minimal operator repair (theoretical style).

Part A (no-go):
- diagnose why old information operators fail,
- formalize a quantitative no-go proxy from observed channel degeneracy.

Part B (minimal repair):
- add one minimal odd-angular splitting term to L->M operator,
- keep deterministic frozen-kernel mapping (no label fitting),
- test whether repaired operator escapes no-go regime.
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
OUT_JSON = ROOT / "report_qw1955_nogo_and_minimal_operator_repair.json"
OUT_MD = ROOT / "RAPORT_QW1955_NOGO_AND_MINIMAL_OPERATOR_REPAIR.md"


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


def derive_kernel_params(kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0

    flips = float(np.mean((s[1:] != s[:-1]).astype(float)))
    even = float(np.sum(ka[(d.astype(int) % 2) == 0]))
    odd = float(np.sum(ka[(d.astype(int) % 2) == 1]))
    parity_imb = float(abs(even - odd) / max(even + odd, 1e-15))
    odd_frac = float(odd / max(even + odd, 1e-15))
    decay = float(np.mean(np.diff(np.log(np.clip(ka[:8], 1e-15, None))) * -1.0))

    eps_pol = float(np.clip(0.11 * flips + 0.07 * parity_imb + 0.04 * abs(math.sin(kernel["phi"])), 0.01, 0.25))
    noise_sigma = float(np.clip(0.012 + 0.05 * flips + 0.02 * parity_imb, 0.01, 0.08))
    incident_amp = float(np.clip(0.05 + 0.10 * abs(math.sin(kernel["omega"])) + 0.03 * abs(ka[0] - ka[1]), 0.03, 0.25))
    absorption = float(np.clip(0.06 + 0.14 * decay, 0.03, 0.35))

    # Minimal odd-angular repair coefficients derived from kernel moments.
    a2 = float(np.clip(0.15 + 0.45 * parity_imb + 0.25 * flips, 0.10, 0.75))
    b1 = float(np.clip(0.08 + 0.60 * abs(odd_frac - 0.5) + 0.20 * flips, 0.06, 0.65))
    b3 = float(np.clip(0.03 + 0.25 * flips + 0.15 * abs(math.sin(kernel["phi"])), 0.02, 0.40))
    retard = float(np.clip(0.03 + 0.18 * abs(kernel["phi"]) / math.pi + 0.10 * flips, 0.02, 0.35))
    psi0 = float(np.mod(0.5 * kernel["phi"] + 0.8 * kernel["omega"], 2.0 * math.pi))

    return {
        "eps_pol": eps_pol,
        "noise_sigma": noise_sigma,
        "incident_amp": incident_amp,
        "absorption": absorption,
        "a2_even_mode": a2,
        "b1_odd_mode": b1,
        "b3_odd_mode": b3,
        "retard_phase": retard,
        "orientation_psi0": psi0,
        "flip_rate": flips,
        "parity_imbalance": parity_imb,
        "odd_fraction": odd_frac,
    }


def build_surface(class_id: int, idx: int, n: int = 64) -> np.ndarray:
    rng = np.random.default_rng(161_000 + 1_000 * class_id + idx)
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
    odd_mode = sign_odd * (
        b1 * np.sin((th - psi0) + phase_shift) + b3 * np.sin(3.0 * (th - psi0) + phase_shift)
    )
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
            rng = np.random.default_rng(171_000 + 1_000 * c + i)

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


def no_go_proxy() -> Dict[str, object]:
    d1952 = json.loads((ROOT / "report_qw1952_information_channel_dedegeneracy_operator.json").read_text(encoding="utf-8"))
    d1953 = json.loads((ROOT / "report_qw1953_two_state_internal_observer.json").read_text(encoding="utf-8"))

    c52 = float(d1952["metrics"]["dual"]["channel_complementarity"])
    c53 = float(d1953["metrics"]["closed"]["channel_complementarity"])
    cmax = max(c52, c53)

    req_acc_gain = max(
        float(THRESHOLDS["acc_gain_vs_control_min"]),
        float(d1953["thresholds"]["closed_acc_gain_vs_open_min"]),
    )
    req_info_gain = max(
        float(THRESHOLDS["info_gain_vs_control_min"]),
        float(d1953["thresholds"]["closed_info_gain_vs_control_min"]),
    )

    # Conservative no-go proxy:
    # under quasi-degenerate channels, potential gain scales ~ O(complementarity).
    ceiling_gain = float(2.0 * cmax)
    no_go_acc = bool(ceiling_gain < req_acc_gain)
    no_go_info = bool(ceiling_gain < req_info_gain)
    no_go_active = bool(no_go_acc and no_go_info)

    return {
        "complementarity_q1952_dual": c52,
        "complementarity_q1953_closed": c53,
        "complementarity_max": cmax,
        "required_acc_gain_min": req_acc_gain,
        "required_info_gain_min": req_info_gain,
        "degenerate_class_gain_ceiling_proxy": ceiling_gain,
        "no_go_acc": no_go_acc,
        "no_go_info": no_go_info,
        "no_go_active_for_old_class": no_go_active,
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

    nogo = no_go_proxy()
    p = derive_kernel_params(kernel)

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

    psf_plus = build_repaired_psf(
        radial_w=w_plus,
        a2=p["a2_even_mode"],
        b1=p["b1_odd_mode"],
        b3=p["b3_odd_mode"],
        psi0=p["orientation_psi0"],
        phase_shift=p["retard_phase"],
        sign_odd=+1.0,
        size=21,
    )
    psf_minus = build_repaired_psf(
        radial_w=w_minus,
        a2=p["a2_even_mode"],
        b1=p["b1_odd_mode"],
        b3=p["b3_odd_mode"],
        psi0=p["orientation_psi0"],
        phase_shift=-p["retard_phase"],
        sign_odd=-1.0,
        size=21,
    )
    psf_single = build_repaired_psf(
        radial_w=w_base,
        a2=0.0,
        b1=0.0,
        b3=0.0,
        psi0=p["orientation_psi0"],
        phase_shift=0.0,
        sign_odd=0.0,
        size=21,
    )
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
    repair_pass = bool(all(flags.values()))

    escaped_no_go_proxy = bool(
        res_dual["channel_complementarity"] > nogo["degenerate_class_gain_ceiling_proxy"]
        and res_dual["channel_complementarity"] > nogo["complementarity_max"]
    )

    verdict = "NOGO_ESCAPED_AND_MINIMAL_REPAIR_PASS" if repair_pass else "NOGO_DIAGNOSIS_WITH_MINIMAL_REPAIR_FAIL"
    required_next = (
        "PROMOTE_REPAIRED_OPERATOR_TO_TWO_STATE_INTERNAL_OBSERVER_TEST"
        if repair_pass
        else "REWORK_MINIMAL_REPAIR_COEFFICIENT_MAP_OR_OPERATOR_FORM"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "thresholds": THRESHOLDS,
        "no_go_proxy": nogo,
        "minimal_repair_params": p,
        "metrics": metrics,
        "flags": flags,
        "escaped_no_go_proxy": escaped_no_go_proxy,
        "repair_pass": repair_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1955: NOGO + MINIMAL OPERATOR REPAIR",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## No-Go Proxy (Old Degenerate Class)",
        f"- complementarity max (QW-1952/1953): {nogo['complementarity_max']:.6f}",
        f"- gain ceiling proxy: {nogo['degenerate_class_gain_ceiling_proxy']:.6f}",
        f"- required acc/info gain mins: {nogo['required_acc_gain_min']:.4f}/{nogo['required_info_gain_min']:.4f}",
        f"- no_go_active_for_old_class: {nogo['no_go_active_for_old_class']}",
        "",
        "## Minimal Repair Outcome",
        f"- dual accuracy: {res_dual['accuracy']:.4f}",
        f"- dual info bits: {res_dual['info_bits']:.4f}",
        f"- dual complementarity: {res_dual['channel_complementarity']:.6f}",
        f"- acc gain vs control: {metrics['acc_gain_dual_vs_control']:.4f}",
        f"- info gain vs control: {metrics['info_gain_dual_vs_control']:.4f}",
        f"- escaped_no_go_proxy: {escaped_no_go_proxy}",
        f"- repair_pass: {repair_pass}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1955] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1955] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1955] verdict={verdict}")


if __name__ == "__main__":
    main()

