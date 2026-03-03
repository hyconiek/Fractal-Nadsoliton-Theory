#!/usr/bin/env python3
"""
QW-1726: Retest of FIN projection hypothesis in GW channels.

Hypothesis under test:
"A true shared background with H=0.23 appears as H_cross~0.31 after detector response."

Protocol:
1) Build synthetic channels with explicit shared/non-shared latent background.
2) Apply empirical detector amplitude envelopes (H1/L1) in Fourier domain.
3) Measure Cross-MF-DFA q=0 under fixed settings and SNR grid.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import detrend


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_strain_unfiltered"
OUT_JSON = ROOT / "report_qw1726_gw_fin_projection_retest.json"
OUT_MD = ROOT / "RAPORT_QW1726_GW_FIN_PROJECTION_RETEST.md"


def read_strain(path: Path, n: int) -> np.ndarray:
    with h5py.File(path, "r") as f:
        x = f["strain"][:n] if "strain" in f else f[list(f.keys())[0]][:n]
    return detrend(np.asarray(x, dtype=float))


def scales_for_n(n: int) -> np.ndarray:
    return np.unique(np.logspace(3, np.log10(n // 4), 15).astype(int))


def cross_mfdfa_q0(x: np.ndarray, y: np.ndarray, scales: np.ndarray) -> float:
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    n = len(z)
    vals = []
    used = []
    for s in scales:
        ns = n // s
        if ns == 0:
            continue
        t = np.arange(s)
        rms = []
        for i in range(ns):
            seg = z[i * s : (i + 1) * s]
            p = np.polyfit(t, seg, 1)
            trend = np.polyval(p, t)
            var = np.mean((seg - trend) ** 2)
            if var > 0:
                rms.append(var)
        if rms:
            r = np.array(rms, dtype=float)
            vals.append(np.exp(0.5 * np.mean(np.log(r + 1e-300))))
            used.append(s)
    if len(vals) < 3:
        return float("nan")
    return float(np.polyfit(np.log(used), np.log(vals), 1)[0])


def generate_fgn(n: int, h: float, rng: np.random.Generator) -> np.ndarray:
    # Spectral synthesis (consistent with previous project scripts).
    f = np.fft.fftfreq(n)
    f[0] = 1e-10
    power = np.abs(f) ** (-(2.0 * h - 1.0))
    power[0] = 0.0
    ph = rng.uniform(0.0, 2.0 * np.pi, n)
    spec = np.sqrt(power) * np.exp(1j * ph)
    if n % 2 == 0:
        spec[n // 2] = np.abs(spec[n // 2])
    for i in range(1, n // 2):
        spec[n - i] = np.conj(spec[i])
    spec[0] = 0.0
    x = np.real(np.fft.ifft(spec))
    x = (x - np.mean(x)) / (np.std(x) + 1e-30)
    return x


def apply_empirical_envelope(x: np.ndarray, amp: np.ndarray) -> np.ndarray:
    xf = np.fft.rfft(x)
    out = np.fft.irfft(amp * np.exp(1j * np.angle(xf)), n=len(x))
    return detrend(out)


def summarize(arr: np.ndarray) -> dict:
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "median": float(np.median(arr)),
        "p10": float(np.quantile(arr, 0.10)),
        "p90": float(np.quantile(arr, 0.90)),
        "prob_near_031_pm_002": float(np.mean(np.abs(arr - 0.31) <= 0.02)),
    }


def main() -> None:
    rng = np.random.default_rng(1726)
    gps = 1266965117
    n = 131072
    n_trials = 24
    h_shared = 0.23
    h_local = 0.04
    snr_grid = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

    p_h1 = RAW / f"H1_unfiltered_{gps}.h5"
    p_l1 = RAW / f"L1_unfiltered_{gps}.h5"
    if not (p_h1.exists() and p_l1.exists()):
        raise FileNotFoundError("Brak danych H1/L1 dla testu 1726.")

    x_h1 = read_strain(p_h1, n)
    x_l1 = read_strain(p_l1, n)
    a_h1 = np.abs(np.fft.rfft(x_h1))
    a_l1 = np.abs(np.fft.rfft(x_l1))
    scales = scales_for_n(n)

    rows = []
    for snr in snr_grid:
        vals_shared = []
        vals_unshared = []
        for _ in range(n_trials):
            bg = generate_fgn(n, h_shared, rng)
            n1 = generate_fgn(n, h_local, rng)
            n2 = generate_fgn(n, h_local, rng)
            s1 = bg + snr * n1
            s2 = bg + snr * n2
            y1 = apply_empirical_envelope(s1, a_h1)
            y2 = apply_empirical_envelope(s2, a_l1)
            vals_shared.append(cross_mfdfa_q0(y1, y2, scales))

            # Control: independent backgrounds (same marginal statistics, no shared field)
            bg1 = generate_fgn(n, h_shared, rng)
            bg2 = generate_fgn(n, h_shared, rng)
            u1 = bg1 + snr * n1
            u2 = bg2 + snr * n2
            z1 = apply_empirical_envelope(u1, a_h1)
            z2 = apply_empirical_envelope(u2, a_l1)
            vals_unshared.append(cross_mfdfa_q0(z1, z2, scales))

        a = np.array(vals_shared, dtype=float)
        b = np.array(vals_unshared, dtype=float)
        rows.append(
            {
                "snr_local_over_shared": float(snr),
                "shared_background": summarize(a),
                "unshared_background_control": summarize(b),
                "delta_mean_shared_minus_unshared": float(np.mean(a) - np.mean(b)),
            }
        )

    # Decide projection support:
    # require at least one SNR with shared mean near 0.31 and clear separation from unshared.
    support = False
    best_row = None
    best_dist = None
    for r in rows:
        mu = r["shared_background"]["mean"]
        dist = abs(mu - 0.31)
        sep = r["delta_mean_shared_minus_unshared"]
        if best_row is None or dist < best_dist:
            best_row = r
            best_dist = dist
        if dist <= 0.02 and sep >= 0.05:
            support = True

    verdict = "FIN_023_TO_031_PROJECTION_SUPPORTED" if support else "FIN_023_TO_031_PROJECTION_NOT_SUPPORTED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": {
            "gps": gps,
            "n_samples": n,
            "n_trials_per_snr": n_trials,
            "h_shared": h_shared,
            "h_local": h_local,
            "snr_grid": snr_grid,
            "seed": 1726,
        },
        "rows": rows,
        "best_row_closest_to_031": best_row,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1726: GW FIN PROJECTION RETEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Wyniki po SNR",
    ]
    for r in rows:
        snr = r["snr_local_over_shared"]
        m_s = r["shared_background"]["mean"]
        sd_s = r["shared_background"]["std"]
        m_u = r["unshared_background_control"]["mean"]
        sd_u = r["unshared_background_control"]["std"]
        d = r["delta_mean_shared_minus_unshared"]
        lines.append(
            f"- SNR={snr:.1f}: shared={m_s:.3f}±{sd_s:.3f}, unshared={m_u:.3f}±{sd_u:.3f}, delta={d:+.3f}"
        )
    if best_row is not None:
        lines.extend(
            [
                "",
                "## Najblizszy punkt do 0.31",
                (
                    f"- SNR={best_row['snr_local_over_shared']:.1f}, "
                    f"shared mean={best_row['shared_background']['mean']:.3f}, "
                    f"delta vs unshared={best_row['delta_mean_shared_minus_unshared']:+.3f}"
                ),
            ]
        )
    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegolowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1726] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1726] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
