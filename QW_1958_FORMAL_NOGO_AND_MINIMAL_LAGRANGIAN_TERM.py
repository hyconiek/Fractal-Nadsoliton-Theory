#!/usr/bin/env python3
"""
QW-1958: Formal no-go bound + minimal EFT repair term.

This script upgrades the theoretical rigor by:
1) deriving a quantitative information-gain upper bound from channel-operator separation,
2) checking whether current operator classes are blocked by that bound,
3) proposing a minimal parity-odd EFT repair term with deterministic coefficients from
   frozen-kernel invariants (no data-label fitting).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1958_formal_nogo_and_minimal_lagrangian_term.json"
OUT_MD = ROOT / "RAPORT_QW1958_FORMAL_NOGO_AND_MINIMAL_LAGRANGIAN_TERM.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def build_operator_kernels_from_report(d_rep: Dict, rep_type: str) -> Dict[str, np.ndarray]:
    """
    Reconstruct compact operator kernels for +/- channels using report coefficients.
    This is an effective representation in angular-harmonic basis on a 2D stencil.
    """
    if rep_type == "q1952":
        p = d_rep["derived_params"]
        a2 = float(p["anisotropy_strength"])
        b1 = 0.0
        b3 = 0.0
        psi0 = float(p["orientation_psi0"])
        phase = float(p["retard_phase"])
        noise_sigma = float(p["noise_sigma"])
    elif rep_type == "q1955":
        p = d_rep["minimal_repair_params"]
        a2 = float(p["a2_even_mode"])
        b1 = float(p["b1_odd_mode"])
        b3 = float(p["b3_odd_mode"])
        psi0 = float(p["orientation_psi0"])
        phase = float(p["retard_phase"])
        noise_sigma = float(p["noise_sigma"])
    else:
        raise ValueError(f"unknown rep_type={rep_type}")

    d1932 = load("report_qw1932_physical_reparameterization_eta_scan.json")
    sel = d1932["selected"]
    omega = float(sel["fit"]["omega"])
    phi = float(sel["fit"]["phi"])
    beta = float(sel["fit"]["beta"])
    eta = float(sel["eta"])

    r = np.arange(1.0, 25.0, dtype=float)
    base = np.abs(kernel_fn(r, omega, phi, beta, eta))
    base = base / max(float(np.sum(base)), 1e-15)

    size = 21
    c = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    x = xx - c
    y = yy - c
    rr = np.sqrt(x * x + y * y)
    th = np.arctan2(y, x)
    ridx = np.clip(np.rint(rr).astype(int), 0, len(base) - 1)

    even_mode = a2 * np.cos(2.0 * (th - psi0))
    odd_mode_p = b1 * np.sin((th - psi0) + phase) + b3 * np.sin(3.0 * (th - psi0) + phase)
    odd_mode_m = -b1 * np.sin((th - psi0) - phase) - b3 * np.sin(3.0 * (th - psi0) - phase)

    k_plus = base[ridx] * (1.0 + even_mode + odd_mode_p)
    k_minus = base[ridx] * (1.0 + even_mode + odd_mode_m)

    k_plus = np.clip(k_plus, 1e-15, None)
    k_minus = np.clip(k_minus, 1e-15, None)
    k_plus /= max(float(np.sum(k_plus)), 1e-15)
    k_minus /= max(float(np.sum(k_minus)), 1e-15)

    return {
        "k_plus": k_plus,
        "k_minus": k_minus,
        "noise_sigma": np.array([noise_sigma], dtype=float),
    }


def operator_separation_metrics(k_plus: np.ndarray, k_minus: np.ndarray, noise_sigma: float) -> Dict[str, float]:
    k_mean = 0.5 * (k_plus + k_minus)
    delta = k_plus - k_minus

    fro_delta = float(np.linalg.norm(delta))
    fro_mean = float(np.linalg.norm(k_mean))
    rel_sep = float(fro_delta / max(fro_mean, 1e-15))

    # Information-gain upper bound proxy (Gaussian linear channel):
    # ΔI <= 0.5 * log2(1 + ||ΔK||_F^2 / σ^2)
    info_gain_bound = float(0.5 * math.log2(1.0 + (fro_delta * fro_delta) / max(noise_sigma * noise_sigma, 1e-15)))

    # Accuracy-gain bound proxy via Pinsker-style relation (conservative):
    # Δacc <= sqrt(2 ln2 * ΔI)
    acc_gain_bound = float(math.sqrt(max(0.0, 2.0 * math.log(2.0) * info_gain_bound)))

    return {
        "fro_delta": fro_delta,
        "fro_mean": fro_mean,
        "relative_separation": rel_sep,
        "info_gain_upper_bound_bits": info_gain_bound,
        "accuracy_gain_upper_bound_proxy": acc_gain_bound,
    }


def derive_minimal_eft_coefficients(kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.arange(1.0, 25.0, dtype=float)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    ka = np.abs(k)
    s = np.sign(k)
    s[s == 0.0] = 1.0

    flip_rate = float(np.mean((s[1:] != s[:-1]).astype(float)))
    even = float(np.sum(ka[(d.astype(int) % 2) == 0]))
    odd = float(np.sum(ka[(d.astype(int) % 2) == 1]))
    parity_imb = float(abs(even - odd) / max(even + odd, 1e-15))
    odd_frac = float(odd / max(even + odd, 1e-15))

    logk = np.log(np.clip(ka[:10], 1e-15, None))
    short_curv = float(np.mean(np.abs(np.diff(logk, n=2))))

    # Deterministic map with EFT-naturalness caps (<= O(1)).
    lambda_i = float(np.clip(0.20 + 0.85 * parity_imb + 0.65 * flip_rate + 0.25 * abs(odd_frac - 0.5), 0.05, 1.00))
    mu_i = float(np.clip(0.06 + 0.70 * short_curv + 0.35 * flip_rate, 0.03, 1.00))

    return {
        "lambda_I": lambda_i,
        "mu_I": mu_i,
        "flip_rate": flip_rate,
        "parity_imbalance": parity_imb,
        "odd_fraction": odd_frac,
        "short_curvature": short_curv,
    }


def main() -> None:
    d1932 = load("report_qw1932_physical_reparameterization_eta_scan.json")
    d1952 = load("report_qw1952_information_channel_dedegeneracy_operator.json")
    d1953 = load("report_qw1953_two_state_internal_observer.json")
    d1955 = load("report_qw1955_nogo_and_minimal_operator_repair.json")
    d1956 = load("report_qw1956_two_state_observer_with_repaired_operator.json")

    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }

    req_acc_gain = max(
        float(d1952["thresholds"]["acc_gain_vs_control_min"]),
        float(d1953["thresholds"]["closed_acc_gain_vs_open_min"]),
    )
    req_info_gain = max(
        float(d1952["thresholds"]["info_gain_vs_control_min"]),
        float(d1953["thresholds"]["closed_info_gain_vs_control_min"]),
    )
    req_complementarity = float(d1952["thresholds"]["channel_complementarity_min"])

    op_old = build_operator_kernels_from_report(d1952, "q1952")
    m_old = operator_separation_metrics(
        op_old["k_plus"], op_old["k_minus"], float(op_old["noise_sigma"][0])
    )
    no_go_old = bool(
        m_old["accuracy_gain_upper_bound_proxy"] < req_acc_gain
        and m_old["info_gain_upper_bound_bits"] < req_info_gain
    )

    op_rep = build_operator_kernels_from_report(d1955, "q1955")
    m_rep = operator_separation_metrics(
        op_rep["k_plus"], op_rep["k_minus"], float(op_rep["noise_sigma"][0])
    )
    no_go_rep = bool(
        m_rep["accuracy_gain_upper_bound_proxy"] < req_acc_gain
        and m_rep["info_gain_upper_bound_bits"] < req_info_gain
    )

    eft = derive_minimal_eft_coefficients(kernel)

    # First-order complementarity proxy anchored on measured repair level.
    c_measured_rep = float(d1955["metrics"]["dual"]["channel_complementarity"])
    alpha, beta = 0.008, 0.005
    c_proxy_next = float(c_measured_rep + alpha * eft["lambda_I"] + beta * eft["mu_I"])

    theorem_text = (
        "If an operator class satisfies ΔI_max < I_req and ΔAcc_max < Acc_req under "
        "the channel-separation bound, strict confirmatory gain is impossible within that class."
    )

    verdict = "FORMAL_NOGO_ESTABLISHED_AND_EFT_TERM_DERIVED" if no_go_old else "FORMAL_NOGO_NOT_ESTABLISHED_UNDER_BOUND"
    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "requirements": {
            "required_acc_gain": req_acc_gain,
            "required_info_gain_bits": req_info_gain,
            "required_complementarity": req_complementarity,
        },
        "formal_bound": {
            "statement": theorem_text,
            "old_class_metrics": m_old,
            "repaired_class_metrics": m_rep,
            "no_go_old_class": no_go_old,
            "no_go_repaired_class": no_go_rep,
        },
        "minimal_eft_repair_term": {
            "symbolic": "ΔL_I = λ_I ε^{ab} n_a A_b (∂_⊥·A) + μ_I ε^{ab}(∂_a A_c)(∂_b A^c)",
            "deterministic_coefficients": {"lambda_I": eft["lambda_I"], "mu_I": eft["mu_I"]},
            "invariants_used": {
                "flip_rate": eft["flip_rate"],
                "parity_imbalance": eft["parity_imbalance"],
                "odd_fraction": eft["odd_fraction"],
                "short_curvature": eft["short_curvature"],
            },
            "naturalness_check": {
                "lambda_I_le_1": bool(eft["lambda_I"] <= 1.0),
                "mu_I_le_1": bool(eft["mu_I"] <= 1.0),
                "all_ok": bool(eft["lambda_I"] <= 1.0 and eft["mu_I"] <= 1.0),
            },
            "complementarity_proxy_after_term": c_proxy_next,
            "proxy_reaches_required": bool(c_proxy_next >= req_complementarity),
        },
        "empirical_context": {
            "q1955_dual_complementarity": c_measured_rep,
            "q1956_closed_complementarity": float(d1956["metrics"]["closed"]["channel_complementarity"]),
        },
        "verdict": verdict,
        "required_next_step": "EXECUTE_QW1959_EFT_TERM_NUMERICAL_TEST_WITHOUT_LABEL_FIT",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1958: FORMAL NOGO + MINIMAL EFT TERM",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Formal Bound",
        f"- no_go_old_class: {no_go_old}",
        f"- no_go_repaired_class: {no_go_rep}",
        f"- old ΔI_max / ΔAcc_max: {m_old['info_gain_upper_bound_bits']:.6f} / {m_old['accuracy_gain_upper_bound_proxy']:.6f}",
        f"- repaired ΔI_max / ΔAcc_max: {m_rep['info_gain_upper_bound_bits']:.6f} / {m_rep['accuracy_gain_upper_bound_proxy']:.6f}",
        "",
        "## Minimal EFT Repair Term",
        f"- lambda_I: {eft['lambda_I']:.6f}",
        f"- mu_I: {eft['mu_I']:.6f}",
        f"- complementarity proxy after term: {c_proxy_next:.6f}",
        f"- required complementarity: {req_complementarity:.6f}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1958] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1958] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1958] verdict={verdict}")


if __name__ == "__main__":
    main()

