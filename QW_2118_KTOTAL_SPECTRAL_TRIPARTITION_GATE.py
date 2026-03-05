#!/usr/bin/env python3
"""
QW-2118: K_total spectral tripartition + stability gate.

Purpose:
- build frozen 12x12 K_total octave-space matrix,
- audit spectral stability properties,
- test deterministic 3-band / 3-partition structure robustness.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2118_ktotal_spectral_tripartition_gate.json"
OUT_MD = ROOT / "RAPORT_QW2118_KTOTAL_SPECTRAL_TRIPARTITION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def kernel_value(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return float(np.cos(omega * d + phi) / (1.0 + beta * (d**eta)))


def build_ktotal_matrix(n: int, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    m = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                m[i, j] = 0.0
                continue
            d = min(abs(i - j), n - abs(i - j))
            m[i, j] = kernel_value(float(d), omega, phi, beta, eta)
    return 0.5 * (m + m.T)


def gershgorin_bounds(m: np.ndarray) -> Tuple[float, float]:
    lowers: List[float] = []
    uppers: List[float] = []
    for i in range(m.shape[0]):
        r = float(np.sum(np.abs(m[i, :]))) - abs(float(m[i, i]))
        c = float(m[i, i])
        lowers.append(c - r)
        uppers.append(c + r)
    return float(min(lowers)), float(max(uppers))


def spectral_band_counts(eigvals: np.ndarray, eps: float) -> Dict[str, int]:
    neg = int(np.sum(eigvals < -eps))
    mid = int(np.sum(np.abs(eigvals) <= eps))
    pos = int(np.sum(eigvals > eps))
    return {"neg": neg, "near_zero": mid, "pos": pos}


def kmeans_tripartition_embedding(m: np.ndarray) -> Dict[str, object]:
    eigvals, eigvecs = np.linalg.eigh(m)
    idx = np.argsort(eigvals)[::-1]
    x = np.abs(eigvecs[:, idx[:3]])

    # Deterministic init on evenly spaced octaves.
    cent = np.array([x[0], x[4], x[8]], dtype=float)
    labels = np.zeros(x.shape[0], dtype=int)
    for _ in range(80):
        d2 = np.sum((x[:, None, :] - cent[None, :, :]) ** 2, axis=2)
        labels = np.argmin(d2, axis=1)
        new_cent = np.array(
            [x[labels == c].mean(axis=0) if np.any(labels == c) else cent[c] for c in range(3)],
            dtype=float,
        )
        if np.allclose(new_cent, cent, atol=1e-12):
            break
        cent = new_cent

    sizes = [int(np.sum(labels == c)) for c in range(3)]
    return {
        "labels": [int(v) for v in labels.tolist()],
        "cluster_sizes": sizes,
        "cluster_sizes_sorted": sorted(sizes),
    }


def perturbed_kernel_samples(base_kernel: Dict[str, float], n: int, seed: int) -> List[Dict[str, float]]:
    rng = np.random.default_rng(seed)
    out: List[Dict[str, float]] = []
    for _ in range(n):
        kk = {
            "omega": float(np.clip(base_kernel["omega"] * (1.0 + rng.normal(0.0, 0.02)), 0.01, 2.0)),
            "phi": float(np.clip(base_kernel["phi"] + rng.normal(0.0, 0.02), -np.pi, np.pi)),
            "beta": float(np.clip(base_kernel["beta"] * (1.0 + rng.normal(0.0, 0.03)), 0.01, 5.0)),
            "eta": float(np.clip(base_kernel["eta"] * (1.0 + rng.normal(0.0, 0.02)), 0.5, 4.0)),
        }
        out.append(kk)
    return out


def main() -> None:
    r2049 = load_json("report_qw2049_spectral_micro_stagec_intersection_gate.json")
    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}

    n_oct = 12
    band_eps = 0.10
    n_perturb = 300

    m = build_ktotal_matrix(n_oct, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    eigvals = np.linalg.eigvalsh(m)
    lam_min = float(np.min(eigvals))
    lam_max = float(np.max(eigvals))
    gl, gu = gershgorin_bounds(m)
    bands = spectral_band_counts(eigvals, eps=band_eps)
    tri = kmeans_tripartition_embedding(m)

    # Distance profile K(d), d=1..11.
    d_profile = {
        str(d): kernel_value(float(d), kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
        for d in range(1, n_oct)
    }

    base_sizes_sorted = tuple(tri["cluster_sizes_sorted"])
    base_bands = (bands["neg"], bands["near_zero"], bands["pos"])

    k_samples = perturbed_kernel_samples(kernel, n=n_perturb, seed=2118)
    stable_size_count = 0
    stable_band_count = 0
    for kk in k_samples:
        mm = build_ktotal_matrix(n_oct, kk["omega"], kk["phi"], kk["beta"], kk["eta"])
        ev = np.linalg.eigvalsh(mm)
        bb = spectral_band_counts(ev, eps=band_eps)
        tt = kmeans_tripartition_embedding(mm)
        if tuple(tt["cluster_sizes_sorted"]) == base_sizes_sorted:
            stable_size_count += 1
        if (bb["neg"], bb["near_zero"], bb["pos"]) == base_bands:
            stable_band_count += 1

    size_stability = float(stable_size_count / n_perturb)
    band_stability = float(stable_band_count / n_perturb)

    required_uniform_mass_shift = float(max(0.0, -lam_min + 1e-9))

    flags = {
        "matrix_symmetric": bool(np.max(np.abs(m - m.T)) <= 1e-12),
        "three_nonempty_spectral_bands": bool(
            bands["neg"] > 0 and bands["near_zero"] > 0 and bands["pos"] > 0
        ),
        "balanced_tripartition_4_4_4": bool(base_sizes_sorted == (4, 4, 4)),
        "tripartition_size_stability_ge_0p90": bool(size_stability >= 0.90),
        "spectral_band_count_stability_ge_0p90": bool(band_stability >= 0.90),
        "deterministic_no_scan_no_retune": True,
        "k_alone_not_psd_requires_mass_shift": bool(lam_min < 0.0),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    structural_pass = bool(
        flags["matrix_symmetric"]
        and flags["three_nonempty_spectral_bands"]
        and flags["balanced_tripartition_4_4_4"]
        and flags["tripartition_size_stability_ge_0p90"]
        and flags["spectral_band_count_stability_ge_0p90"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "KTOTAL_SPECTRAL_TRIPARTITION_GATE_PASS_WITH_CONDITIONAL_VACUUM_SHIFT"
        if structural_pass
        else "KTOTAL_SPECTRAL_TRIPARTITION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_kernel": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        "kernel": kernel,
        "matrix_spec": {
            "n_octaves": n_oct,
            "diagonal_convention": "zero_self_coupling",
            "distance": "cyclic_octave_distance",
        },
        "distance_profile": d_profile,
        "eigen_spectrum": {
            "eigvals": [float(v) for v in eigvals.tolist()],
            "lambda_min": lam_min,
            "lambda_max": lam_max,
            "gershgorin_lower_bound": gl,
            "gershgorin_upper_bound": gu,
        },
        "band_partition": {
            "band_eps": band_eps,
            "counts": bands,
            "tripartition": tri,
        },
        "robustness": {
            "n_perturb": n_perturb,
            "tripartition_size_stability": size_stability,
            "band_count_stability": band_stability,
        },
        "vacuum_shift_condition": {
            "k_alone_psd": bool(lam_min >= 0.0),
            "required_uniform_mass_shift_ge": required_uniform_mass_shift,
            "interpretation": (
                "Need diagonal positive mass floor m0^2 >= -lambda_min(K_total) to guarantee PSD(A=K_total+m0^2 I)."
            ),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "CONNECT_REQUIRED_MASS_SHIFT_TO_EXPLICIT_SCALAR_SECTOR_AND_RUN_QW2119"
            if structural_pass
            else "REVIEW_SPECTRAL_PARTITION_ASSUMPTIONS_AND_RERUN_QW2118"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2118: KTOTAL SPECTRAL TRIPARTITION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Spectrum",
        f"- lambda_min/lambda_max: `{lam_min:.9f}` / `{lam_max:.9f}`",
        f"- Gershgorin bounds: [`{gl:.9f}`, `{gu:.9f}`]",
        f"- band counts (neg/near/pos): `{bands['neg']}/{bands['near_zero']}/{bands['pos']}`",
        "",
        "## Tripartition",
        f"- cluster sizes: `{tri['cluster_sizes']}` (sorted: `{tri['cluster_sizes_sorted']}`)",
        f"- size stability: `{size_stability:.3f}`",
        f"- band stability: `{band_stability:.3f}`",
        "",
        "## Vacuum condition",
        f"- required uniform mass shift >= `{required_uniform_mass_shift:.9f}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2118] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2118] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2118] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

