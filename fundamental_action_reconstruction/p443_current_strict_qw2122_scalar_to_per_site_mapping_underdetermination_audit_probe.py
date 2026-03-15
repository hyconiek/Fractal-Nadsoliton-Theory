#!/usr/bin/env python3
from __future__ import annotations

import cmath
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

R14_JSON = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
R15_JSON = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
QW2122_JSON = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"

OUT_JSON = (
    GENERATED
    / "p443_current_strict_qw2122_scalar_to_per_site_mapping_underdetermination_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p443_current_strict_qw2122_scalar_to_per_site_mapping_underdetermination_audit_probe_summary.json"
)


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def complex_to_json(z: complex) -> dict[str, float]:
    return {"re": float(z.real), "im": float(z.imag), "abs": float(abs(z))}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def mixing_row_dot(*, i: int, ksym: list[list[float]], vpsi: list[float]) -> float:
    n = len(vpsi)
    return float(sum(float(ksym[i][j]) * float(vpsi[j]) for j in range(n) if j != i))


def solve_m2_from_constant_vacuum_eom(
    *,
    ksym: list[list[float]],
    vpsi: list[float],
    g4: list[float],
    g6: list[float],
    vphi: float = 0.0,
    gY: list[float] | None = None,
) -> tuple[list[float], list[float]]:
    # This is a toy check: we set gY=0 by default and solve m2_i per-site from the stationarity equation.
    n = len(vpsi)
    gY_eff = gY if gY is not None else [0.0 for _ in range(n)]
    m2: list[float] = []
    residuals: list[float] = []
    for i in range(n):
        mix = mixing_row_dot(i=i, ksym=ksym, vpsi=vpsi)
        m2_i = float(
            -(mix / vpsi[i])
            - float(g4[i]) * (vpsi[i] ** 2)
            - float(g6[i]) * (vpsi[i] ** 4)
            - 2.0 * float(gY_eff[i]) * (vphi**2)
        )
        m2.append(m2_i)
        res = float(
            mix
            + float(g4[i]) * (vpsi[i] ** 3)
            + float(g6[i]) * (vpsi[i] ** 5)
            + 2.0 * float(gY_eff[i]) * (vphi**2) * vpsi[i]
            + m2_i * vpsi[i]
        )
        residuals.append(res)
    return m2, residuals


def compute_d_profile_n477(
    *,
    ksym: list[list[float]],
    vpsi: list[float],
    g4: list[float],
    g6: list[float],
    m0_sq: float,
) -> list[float]:
    # N477: Yukawa-free, K_total-numeric diagonal residual (conditional on stationarity + vpsi_k != 0).
    n = len(vpsi)
    d: list[float] = []
    for k in range(n):
        denom = float(vpsi[k])
        mix_ratio = 0.0
        for j in range(n):
            if j == k:
                continue
            mix_ratio += float(ksym[k][j]) * (float(vpsi[j]) / denom)
        val = -mix_ratio + 2.0 * float(g4[k]) * (denom**2) + 4.0 * float(g6[k]) * (denom**4) - float(m0_sq)
        d.append(float(val))
    return d


def compute_sigma_opposite_pair_sums(*, d: list[float]) -> dict[str, float]:
    if len(d) != 12:
        raise ValueError("This probe is scoped to n=12.")
    return {f"Sigma_psi{k}_psi{k+6}": float(d[k] + d[k + 6]) for k in range(6)}


def compute_f2_from_sigmas(*, sigmas: dict[str, float]) -> complex:
    # N467 reduction: F2(d) = sum_{k=0..5} Sigma_k * exp(i π k/3).
    total: complex = 0.0 + 0.0j
    for k in range(6):
        total += complex(sigmas[f"Sigma_psi{k}_psi{k+6}"]) * cmath.exp(1j * (math.pi * k / 3.0))
    return total


def build_uniform_vpsi_from_rho_star_sq(*, rho_star_sq: float) -> list[float]:
    v = math.sqrt(float(rho_star_sq) / 12.0)
    return [float(v) for _ in range(12)]


def build_two_site_perturbation_same_norm(*, v_base: float, eps: float) -> list[float]:
    # Keep v0^2 + v6^2 = 2 * v_base^2 with v0 = v_base + eps (eps small enough), others unchanged.
    a = float(v_base + eps)
    b_sq = 2.0 * (v_base**2) - (a**2)
    if not (b_sq > 0.0):
        raise ValueError("eps too large: cannot keep two-site norm constraint with real v6")
    b = math.sqrt(b_sq)
    vpsi = [float(v_base) for _ in range(12)]
    vpsi[0] = float(a)
    vpsi[6] = float(b)
    return vpsi


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in [R14_JSON, R15_JSON, QW2122_JSON]:
        if not p.exists():
            missing.append(str(p.relative_to(REPO)) if p.is_absolute() else str(p))
    if missing:
        summary = {
            "stage": "P443",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing,
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    r14 = load_json(R14_JSON)
    r15 = load_json(R15_JSON)
    q2122 = load_json(QW2122_JSON)

    K_total = r14.get("host_kernel_rows")
    m0_sq = (r15.get("diagonal_decomposition") or {}).get("host_diagonal_floor_value")
    rho_star_sq = ((q2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("rho_star_sq")
    lambda_psi_strict = (q2122.get("inputs") or {}).get("lambda_psi_strict")

    if not (
        isinstance(K_total, list)
        and len(K_total) == 12
        and all(isinstance(r, list) and len(r) == 12 for r in K_total)
        and all(is_number(x) for row in K_total for x in row)
    ):
        missing.append("R14.host_kernel_rows (12x12 numeric K_total matrix)")
    if not is_number(m0_sq):
        missing.append("R15.diagonal_decomposition.host_diagonal_floor_value (numeric m0^2)")
    if not is_number(rho_star_sq):
        missing.append("QW-2122.branch_results.broken_branch_strict.rho_star_sq (numeric rho_*^2)")
    if not is_number(lambda_psi_strict):
        missing.append("QW-2122.inputs.lambda_psi_strict (numeric lambda_psi_strict)")

    if missing:
        summary = {
            "stage": "P443",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_VALUES",
            "missing_or_invalid": sorted(set(missing)),
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    assert isinstance(K_total, list)
    m0_sq_f = float(m0_sq)
    rho_star_sq_f = float(rho_star_sq)
    lambda_psi_f = float(lambda_psi_strict)

    # Explicit premises for the audit (not strict-derived):
    # - uniform per-site quartic matching along the uniform ray (AX24-style): g4_i := 12*lambda_psi_strict, g6_i := 0
    g4 = [12.0 * lambda_psi_f for _ in range(12)]
    g6 = [0.0 for _ in range(12)]

    vpsi_uniform = build_uniform_vpsi_from_rho_star_sq(rho_star_sq=rho_star_sq_f)
    v_base = float(vpsi_uniform[0])
    eps = 0.05
    vpsi_perturbed = build_two_site_perturbation_same_norm(v_base=v_base, eps=eps)

    def check_norm(vpsi: list[float]) -> float:
        return float(sum(float(x) ** 2 for x in vpsi))

    def case(case_id: str, vpsi: list[float]) -> dict[str, Any]:
        if any(float(x) == 0.0 for x in vpsi):
            raise ValueError("vpsi contains zero entry; N477 division requires vpsi_k != 0")
        d = compute_d_profile_n477(ksym=K_total, vpsi=vpsi, g4=g4, g6=g6, m0_sq=m0_sq_f)
        sigmas = compute_sigma_opposite_pair_sums(d=d)
        F2 = compute_f2_from_sigmas(sigmas=sigmas)
        m2, eom_residuals = solve_m2_from_constant_vacuum_eom(ksym=K_total, vpsi=vpsi, g4=g4, g6=g6)
        return {
            "case_id": case_id,
            "inputs": {
                "vpsi": vpsi,
                "rho_star_sq_target": rho_star_sq_f,
                "rho_star_sq_actual": check_norm(vpsi),
                "g4": g4,
                "g6": g6,
                "m0_sq": m0_sq_f,
                "eps_used": eps if case_id.endswith("PERTURBED_VPSI") else 0.0,
            },
            "derived": {
                "m2_solved_from_constant_vacuum_eom": m2,
                "max_abs_constant_vacuum_eom_residual": float(max(abs(float(r)) for r in eom_residuals)),
                "d_local_residual_profile": d,
                "sigma_opposite_pair_sums": sigmas,
                "F2_from_sigmas": complex_to_json(F2),
            },
        }

    case_uniform = case("CASE_A_UNIFORM_VPSI_SAME_RHO_STAR_SQ", vpsi_uniform)
    case_nonuniform = case("CASE_B_TWO_SITE_PERTURBED_VPSI", vpsi_perturbed)

    artifact: dict[str, Any] = {
        "stage": "P443",
        "as_of": "2026-03-14",
        "no_false_pass": True,
        "goal": "toy_audit_that_qw2122_scalar_outputs_do_not_canonically_fix_the_canonical_per_site_arrays_needed_by_T168/P437",
        "strict_sources_used_as_inputs": {
            "QW-2122": str(QW2122_JSON.name),
            "R14": str(R14_JSON.relative_to(REPO)),
            "R15": str(R15_JSON.relative_to(REPO)),
        },
        "inputs_from_strict_sources": {
            "rho_star_sq": rho_star_sq_f,
            "lambda_psi_strict": lambda_psi_f,
            "m0_sq": m0_sq_f,
        },
        "explicit_mapping_premises": {
            "uniform_self_couplings": "g4_i := 12*lambda_psi_strict, g6_i := 0 for i=0..11 (premise; not strict-derived)",
            "rho_to_vpsi_only_fixes_norm": "rho_star_sq fixes only Σ_i vpsi_i^2; per-site distribution/orientation remains free without an additional strict mapping/selector.",
        },
        "cases": {"case_uniform": case_uniform, "case_nonuniform": case_nonuniform},
        "conclusion": {
            "case_uniform_F2_abs": case_uniform["derived"]["F2_from_sigmas"]["abs"],
            "case_nonuniform_F2_abs": case_nonuniform["derived"]["F2_from_sigmas"]["abs"],
            "mapping_underdetermined": True,
            "recommended_next_strict_target": "T168",
        },
    }

    summary = {
        "stage": "P443",
        "status": "EXECUTED_TOY_AUDIT_NO_FALSE_PASS",
        "rho_star_sq": rho_star_sq_f,
        "lambda_psi_strict": lambda_psi_f,
        "m0_sq": m0_sq_f,
        "case_uniform_F2_abs": artifact["conclusion"]["case_uniform_F2_abs"],
        "case_nonuniform_F2_abs": artifact["conclusion"]["case_nonuniform_F2_abs"],
        "recommended_next_strict_target": "T168",
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

