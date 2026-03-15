#!/usr/bin/env python3
from __future__ import annotations

import cmath
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "p436_current_strict_vacuum_eom_restricted_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p436_current_strict_vacuum_eom_restricted_canonical_local_diagonal_mode2_defect_underdetermination_audit_probe_summary.json"
)


def max_abs(values: list[float]) -> float:
    return float(max(abs(float(v)) for v in values)) if values else 0.0


def complex_to_json(z: complex) -> dict[str, float]:
    return {"re": float(z.real), "im": float(z.imag), "abs": float(abs(z))}


def build_ksym_nearest_neighbor_circulant(*, n: int, weight: float) -> list[list[float]]:
    ksym: list[list[float]] = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        j = (i + 1) % n
        ksym[i][j] = float(weight)
        ksym[j][i] = float(weight)
    return ksym


def mixing_row_dot(*, i: int, ksym: list[list[float]], vpsi: list[float]) -> float:
    n = len(vpsi)
    return float(sum(ksym[i][j] * vpsi[j] for j in range(n) if j != i))


def solve_m2_from_constant_vacuum_eom(
    *,
    ksym: list[list[float]],
    vpsi: list[float],
    vphi: float,
    g4: list[float],
    g6: list[float],
    gY: list[float],
) -> tuple[list[float], list[float]]:
    n = len(vpsi)
    m2: list[float] = []
    residuals: list[float] = []
    for i in range(n):
        mix = mixing_row_dot(i=i, ksym=ksym, vpsi=vpsi)
        m2_i = float(
            -(mix / vpsi[i])
            - g4[i] * (vpsi[i] ** 2)
            - g6[i] * (vpsi[i] ** 4)
            - 2.0 * gY[i] * (vphi**2)
        )
        m2.append(m2_i)
        res = float(
            mix
            + g4[i] * (vpsi[i] ** 3)
            + g6[i] * (vpsi[i] ** 5)
            + 2.0 * gY[i] * (vphi**2) * vpsi[i]
            + m2_i * vpsi[i]
        )
        residuals.append(res)
    return m2, residuals


def canonical_diag_entry(
    *,
    vpsi_i: float,
    vphi: float,
    g4_i: float,
    g6_i: float,
    gY_i: float,
    m2_i: float,
) -> float:
    return float(
        3.0 * g4_i * (vpsi_i**2) + 5.0 * g6_i * (vpsi_i**4) + 2.0 * gY_i * (vphi**2) + m2_i
    )


def compute_d_profile(
    *,
    vpsi: list[float],
    vphi: float,
    g4: list[float],
    g6: list[float],
    gY: list[float],
    m2: list[float],
    m0_sq: float,
) -> list[float]:
    n = len(vpsi)
    d: list[float] = []
    for i in range(n):
        h_i = canonical_diag_entry(
            vpsi_i=vpsi[i], vphi=vphi, g4_i=g4[i], g6_i=g6[i], gY_i=gY[i], m2_i=m2[i]
        )
        d.append(float(h_i - m0_sq))
    return d


def compute_sigma_opposite_pair_sums(*, d: list[float]) -> dict[str, float]:
    if len(d) != 12:
        raise ValueError("This probe is scoped to n=12.")
    out: dict[str, float] = {}
    for k in range(6):
        out[f"Sigma_psi{k}_psi{k+6}"] = float(d[k] + d[k + 6])
    return out


def compute_f2_from_d(*, d: list[float]) -> complex:
    n = len(d)
    # Mode-2 coefficient: sum_k d_k * exp(i * 4πk/n).
    return sum(complex(dk) * cmath.exp(1j * (4.0 * cmath.pi * k / n)) for k, dk in enumerate(d))


def compute_f2_from_sigmas(*, sigmas: dict[str, float]) -> complex:
    # For n=12, reduction from N467: F2(d) = sum_{k=0..5} Sigma_k * exp(i π k/3).
    total: complex = 0.0 + 0.0j
    for k in range(6):
        total += complex(sigmas[f"Sigma_psi{k}_psi{k+6}"]) * cmath.exp(1j * (cmath.pi * k / 3.0))
    return total


def run_case(*, case_id: str, ksym: list[list[float]], vpsi: list[float], vphi: float, g4: list[float], g6: list[float], gY: list[float], m0_sq: float) -> dict[str, Any]:
    m2, eom_residuals = solve_m2_from_constant_vacuum_eom(ksym=ksym, vpsi=vpsi, vphi=vphi, g4=g4, g6=g6, gY=gY)
    d = compute_d_profile(vpsi=vpsi, vphi=vphi, g4=g4, g6=g6, gY=gY, m2=m2, m0_sq=m0_sq)
    sigmas = compute_sigma_opposite_pair_sums(d=d)
    f2_d = compute_f2_from_d(d=d)
    f2_sig = compute_f2_from_sigmas(sigmas=sigmas)
    f2_consistency_error = complex(f2_d - f2_sig)
    return {
        "case_id": case_id,
        "inputs": {
            "vpsi": vpsi,
            "vphi": vphi,
            "g4": g4,
            "g6": g6,
            "gY": gY,
            "m0_sq": m0_sq,
            "ksym": ksym,
        },
        "derived": {
            "m2_solved_from_constant_vacuum_eom": m2,
            "d_profile": d,
            "sigma_opposite_pair_sums": sigmas,
            "F2_from_d": complex_to_json(f2_d),
            "F2_from_sigmas": complex_to_json(f2_sig),
            "F2_consistency_error": complex_to_json(f2_consistency_error),
        },
        "checks": {
            "max_abs_constant_vacuum_eom_residual": max_abs(eom_residuals),
        },
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n = 12
    tol = 1e-12

    # Use one fixed translation-invariant symmetric mixing matrix (circulant nearest-neighbor toy kernel).
    ksym = build_ksym_nearest_neighbor_circulant(n=n, weight=1.0)

    # Keep a constant nonzero vacuum vector to isolate the remaining local-family freedom.
    vpsi = [1.0 for _ in range(n)]
    vphi = 1.0
    g6 = [0.0 for _ in range(n)]
    gY = [0.0 for _ in range(n)]
    m0_sq = 0.0

    # Case A: uniform local couplings => constant diagonal profile => F2(d)=0.
    g4_uniform = [0.0 for _ in range(n)]
    case_f2_zero = run_case(
        case_id="CASE_F2_ZERO_UNIFORM_LOCAL_COUPLINGS",
        ksym=ksym,
        vpsi=vpsi,
        vphi=vphi,
        g4=g4_uniform,
        g6=g6,
        gY=gY,
        m0_sq=m0_sq,
    )

    # Case B: one-site local defect (g4_psi0) => non-translation-invariant diagonal profile => F2(d)≠0.
    g4_defect = [0.0 for _ in range(n)]
    g4_defect[0] = 1.0
    case_f2_nonzero = run_case(
        case_id="CASE_F2_NONZERO_SINGLE_SITE_G4_DEFECT",
        ksym=ksym,
        vpsi=vpsi,
        vphi=vphi,
        g4=g4_defect,
        g6=g6,
        gY=gY,
        m0_sq=m0_sq,
    )

    def f2_abs(case: dict[str, Any]) -> float:
        return float(case["derived"]["F2_from_d"]["abs"])

    artifact: dict[str, Any] = {
        "probe_id": "P436_CURRENT_STRICT_VACUUM_EOM_RESTRICTED_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_UNDERDETERMINATION_AUDIT_PROBE",
        "as_of": "2026-03-13",
        "no_false_pass": True,
        "n_sites": n,
        "tolerance": tol,
        "toy_kernel_note": "Nearest-neighbor symmetric circulant Ksym. This is a toy instantiation, not a strict-derived kernel claim.",
        "cases": {
            "case_f2_zero": case_f2_zero,
            "case_f2_nonzero": case_f2_nonzero,
        },
        "checks": {
            "case_f2_zero_passes_at_tolerance": bool(
                f2_abs(case_f2_zero) <= tol and case_f2_zero["checks"]["max_abs_constant_vacuum_eom_residual"] <= tol
            ),
            "case_f2_nonzero_passes_eom_at_tolerance": bool(
                case_f2_nonzero["checks"]["max_abs_constant_vacuum_eom_residual"] <= tol
            ),
            "case_f2_nonzero_is_nonzero": bool(f2_abs(case_f2_nonzero) > tol),
        },
        "theorem_pointers": ["N474", "N475", "N476"],
        "target_pointer": "T166",
        "notes": [
            "Toy audit only: this does not claim any strict-derived or physically admissible instantiation.",
            "Both cases satisfy the constant-vacuum EoM by construction (m2_i solved per-site).",
            "The two cases differ only in a single local coupling entry (g4_psi0).",
            "This demonstrates underdetermination of F2(d) even under constant-vacuum stationarity + vpsi_i≠0 restrictions.",
        ],
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_TOY_AUDIT_NO_FALSE_PASS",
        "n_sites": n,
        "tolerance": tol,
        "case_f2_zero_abs": artifact["cases"]["case_f2_zero"]["derived"]["F2_from_d"]["abs"],
        "case_f2_nonzero_abs": artifact["cases"]["case_f2_nonzero"]["derived"]["F2_from_d"]["abs"],
        "case_f2_zero_passes_at_tolerance": artifact["checks"]["case_f2_zero_passes_at_tolerance"],
        "case_f2_nonzero_passes_eom_at_tolerance": artifact["checks"]["case_f2_nonzero_passes_eom_at_tolerance"],
        "case_f2_nonzero_is_nonzero": artifact["checks"]["case_f2_nonzero_is_nonzero"],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

