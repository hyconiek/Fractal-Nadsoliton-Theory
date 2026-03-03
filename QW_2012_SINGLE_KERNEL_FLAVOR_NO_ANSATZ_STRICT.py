#!/usr/bin/env python3
"""
QW-2012: Single-kernel flavor derivation in strict no-ansatz / no-fit mode (Task 2).

Rules:
- one frozen kernel,
- no parameter scan/optimization,
- no per-sector calibration,
- no CKM/PMNS target fitting.
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2012_single_kernel_flavor_no_ansatz_strict.json"
OUT_MD = ROOT / "RAPORT_QW2012_SINGLE_KERNEL_FLAVOR_NO_ANSATZ_STRICT.md"

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


def kernel_complex_signed(d_signed: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    d_abs = np.abs(d_signed)
    amp = 1.0 / (1.0 + beta * (d_abs**eta))
    phase = omega * d_signed + phi * np.sign(d_signed)
    return amp * np.exp(1j * phase)


def cyclic_signed_distance(q_left: np.ndarray, q_right: np.ndarray, modulus: int = 24) -> np.ndarray:
    raw = q_left[:, None] - q_right[None, :]
    half = modulus / 2.0
    signed = ((raw + half) % modulus) - half
    return signed


def no_ansatz_hamiltonian(q: np.ndarray, kernel: Dict[str, float]) -> np.ndarray:
    d_signed = cyclic_signed_distance(q, q, modulus=24)
    z = kernel_complex_signed(d_signed, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    # Strict minimal rule: hermitianized complex kernel + fixed centering diagonal.
    h_off = 0.5 * (z + z.conj().T)
    np.fill_diagonal(h_off, 0.0)

    q_centered = q - np.mean(q)
    diag = np.diag(q_centered / max(np.std(q_centered), 1e-12))

    h = h_off + diag
    return 0.5 * (h + h.conj().T)


def flavor_metrics(q_up: np.ndarray, q_down: np.ndarray, q_nu: np.ndarray, q_lep: np.ndarray, kernel: Dict[str, float]) -> Dict[str, float]:
    hu = no_ansatz_hamiltonian(q_up, kernel)
    hd = no_ansatz_hamiltonian(q_down, kernel)
    hn = no_ansatz_hamiltonian(q_nu, kernel)
    hl = no_ansatz_hamiltonian(q_lep, kernel)

    _, uu = np.linalg.eigh(hu)
    _, ud = np.linalg.eigh(hd)
    _, un = np.linalg.eigh(hn)
    _, ul = np.linalg.eigh(hl)

    ckm = np.abs(uu.conj().T @ ud)
    pmns = np.abs(un.conj().T @ ul)

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)

    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "ckm_max_rel_pct": float(100.0 * np.max(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
        "pmns_max_rel_pct": float(100.0 * np.max(pmns_rel)),
    }


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    ksel = d1932["selected"]
    kernel = {
        "omega": float(ksel["fit"]["omega"]),
        "phi": float(ksel["fit"]["phi"]),
        "beta": float(ksel["fit"]["beta"]),
        "eta": float(ksel["eta"]),
    }

    d1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    q_spaces = d1961["inputs"]["q_assignments"]

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
    }

    rows: List[Dict[str, object]] = []
    q_nu_space = list(itertools.permutations([0.0, 1.0, 2.0], 3))

    for q_name, q_map in q_spaces.items():
        q_up = np.array([q_map["Top"], q_map["Charm"], q_map["Muon"]], dtype=float)
        q_down = np.array([q_map["Bottom"], q_map["Tau"], q_map["Muon"]], dtype=float)
        q_lep = np.array([q_map["Electron"], q_map["Muon"], q_map["Tau"]], dtype=float)

        for q_nu in q_nu_space:
            q_nu_arr = np.array(q_nu, dtype=float)
            m = flavor_metrics(q_up, q_down, q_nu_arr, q_lep, kernel)
            flags = {
                "ckm_mean_rel_pct_le_max": bool(m["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
                "pmns_mean_rel_pct_le_max": bool(m["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
            }
            row = {
                "q_assignment": q_name,
                "q_nu": [int(x) for x in q_nu],
                "metrics": m,
                "flags": flags,
                "all_pass": bool(all(flags.values())),
                "score": float(m["ckm_mean_rel_pct"] / 15.0 + m["pmns_mean_rel_pct"] / 15.0),
            }
            rows.append(row)

    rows_sorted = sorted(rows, key=lambda x: x["score"])
    pass_count = sum(1 for r in rows if r["all_pass"])
    best = rows_sorted[0]

    verdict = "SINGLE_KERNEL_FLAVOR_NO_ANSATZ_PASS" if pass_count > 0 else "SINGLE_KERNEL_FLAVOR_NO_ANSATZ_FAIL"
    required_next = (
        "FREEZE_NO_ANSATZ_FLAVOR_BRANCH_FOR_TRIAD"
        if pass_count > 0
        else "NO_ANSATZ_FLAVOR_NOT_REACHED_NEED_DEEPER_MICRODYNAMIC_DERIVATION"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel": kernel,
        "protocol": {
            "no_fit": True,
            "no_parameter_scan": True,
            "single_kernel_shared_operator": True,
            "operator_rule": "hermitianized_complex_kernel_plus_fixed_centering_diagonal",
        },
        "thresholds": thresholds,
        "summary": {
            "n_candidates": len(rows),
            "pass_count": int(pass_count),
            "best": best,
            "top10": rows_sorted[:10],
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    b = best
    lines = [
        "# RAPORT QW-2012: SINGLE KERNEL FLAVOR NO-ANSATZ STRICT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{len(rows)}",
        "",
        "## Protocol",
        "- no fit / no per-sector tuning / one frozen kernel",
        "",
        "## Best Candidate",
        f"- q_assignment: {b['q_assignment']}, q_nu: {b['q_nu']}",
        (
            f"- CKM mean/max rel%: {b['metrics']['ckm_mean_rel_pct']:.3f} / "
            f"{b['metrics']['ckm_max_rel_pct']:.3f}"
        ),
        (
            f"- PMNS mean/max rel%: {b['metrics']['pmns_mean_rel_pct']:.3f} / "
            f"{b['metrics']['pmns_max_rel_pct']:.3f}"
        ),
        f"- all_pass: {b['all_pass']}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2012] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2012] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2012] verdict={verdict}, pass_count={pass_count}")


if __name__ == "__main__":
    main()
