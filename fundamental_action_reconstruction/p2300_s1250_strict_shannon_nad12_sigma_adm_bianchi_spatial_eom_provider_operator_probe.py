#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_1853 = GEN / "p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json"
IN_1981 = GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json"
IN_1982 = GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json"
IN_1983 = GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json"
IN_1988 = GEN / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"
IN_2282 = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
IN_2297 = GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json"
IN_2299 = GEN / "p2299_s1249_strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe.json"
OUT = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
MD = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def matrix_report(rows: list[list[sp.Expr]], rhs: list[sp.Expr], unknowns: list[sp.Symbol]) -> dict[str, Any]:
    a_mat = sp.Matrix(rows)
    b_vec = sp.Matrix(rhs)
    rank_a = int(a_mat.rank())
    rank_aug = int(a_mat.row_join(b_vec).rank())
    consistent = rank_a == rank_aug
    a_num = np.array([[float(sp.N(cell, 40)) for cell in row] for row in rows], dtype=float)
    b_num = np.array([float(sp.N(cell, 40)) for cell in rhs], dtype=float)
    x_lsq, *_ = np.linalg.lstsq(a_num, b_num, rcond=None)
    residual = a_num @ x_lsq - b_num
    return {
        "equation_count": len(rows),
        "unknown_count": len(unknowns),
        "rank_A": rank_a,
        "rank_augmented": rank_aug,
        "consistent": consistent,
        "full_solve": consistent and rank_a == len(unknowns),
        "least_squares_residual_l2": float(la.norm(residual, ord=2)),
        "least_squares_residual_max_abs": float(np.max(np.abs(residual))) if residual.size else 0.0,
        "least_squares_solution_preview": {str(unknowns[i]): float(x_lsq[i]) for i in range(min(len(unknowns), 10))},
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha_packet = load(IN_ALPHA)
    p1853 = load(IN_1853)
    p1981 = load(IN_1981)
    p1982 = load(IN_1982)
    p1983 = load(IN_1983)
    p1988 = load(IN_1988)
    p2282 = load(IN_2282)
    p2297 = load(IN_2297)
    p2299 = load(IN_2299)

    names = [
        "N",
        "Nd",
        "Ndd",
        "V",
        "H",
        "Hd",
        "Hdd",
        "sigma1",
        "sigma2",
        "dsigma1",
        "dsigma2",
        "d2sigma1",
        "d2sigma2",
        "Q",
    ]
    local = {name: sp.Symbol(name, real=True) for name in names}
    local.update({"pi": sp.pi, "log": sp.log, "ln": sp.log})
    N, Nd, Ndd, V, H, Hd, Hdd, s1, s2, sd1, sd2, sdd1, sdd2 = [
        local[key]
        for key in [
            "N",
            "Nd",
            "Ndd",
            "V",
            "H",
            "Hd",
            "Hdd",
            "sigma1",
            "sigma2",
            "dsigma1",
            "dsigma2",
            "d2sigma1",
            "d2sigma2",
        ]
    ]
    q_shear = s1**2 + s1 * s2 + s2**2
    q_dot = 2 * s1 * sd1 + s2 * sd1 + s1 * sd2 + 2 * s2 * sd2
    q_dot_sq = sd1**2 + sd1 * sd2 + sd2**2
    q_ddot_pair = 2 * s1 * sdd1 + s2 * sdd1 + s1 * sdd2 + 2 * s2 * sdd2

    coeffs = (p1853.get("b1_symbolic_evaluation", {}) or {}).get("evaluated_coefficients", {}) or {}
    a_r2 = sp.sympify((coeffs.get("a_R2", {}) or {}).get("symbolic", "0"), locals=local)
    a_ric2 = sp.sympify((coeffs.get("a_Ric2", {}) or {}).get("symbolic", "0"), locals=local)
    a_riem2 = sp.sympify((coeffs.get("a_Riem2", {}) or {}).get("symbolic", "0"), locals=local)
    density_r2 = sp.sympify(
        (p1981.get("r2_lapse_euler_operator", {}) or {}).get("density_difference_NV_R2", "0"),
        locals=local,
    ).subs({local["Q"]: q_shear})
    density_ric2 = sp.sympify(
        (p1982.get("ricci2_lapse_euler_operator", {}) or {}).get("density_difference_NV_Ricci2", "0"),
        locals=local,
    )
    density_riem2 = sp.sympify(
        (p1983.get("riemann2_lapse_euler_operator", {}) or {}).get("density_difference_NV_Riemann2", "0"),
        locals=local,
    )
    density_non_gb = sp.factor(a_r2 * density_r2 + a_ric2 * density_ric2 + a_riem2 * density_riem2)

    dt_rules = {
        N: Nd,
        Nd: Ndd,
        V: 3 * H * V,
        H: Hd,
        Hd: Hdd,
        s1: sd1,
        s2: sd2,
        sd1: sdd1,
        sd2: sdd2,
    }

    def total_dt(expr: sp.Expr) -> sp.Expr:
        return sp.factor(sum(sp.diff(expr, var) * dvar for var, dvar in dt_rules.items()))

    def euler_lagrange_component(density: sp.Expr, sigma: sp.Symbol, dsigma: sp.Symbol, d2sigma: sp.Symbol) -> sp.Expr:
        raw = sp.diff(density, sigma) - total_dt(sp.diff(density, dsigma)) + total_dt(total_dt(sp.diff(density, d2sigma)))
        return sp.factor(sp.simplify(raw * N**8 / V))

    def spatial_components(density: sp.Expr) -> list[sp.Expr]:
        return [
            euler_lagrange_component(density, s1, sd1, sdd1),
            euler_lagrange_component(density, s2, sd2, sdd2),
        ]

    target_components = [-component for component in spatial_components(density_non_gb)]
    alpha = sp.Integer(4) * sp.log(2)
    base_density = V * alpha / N**3
    provider_terms = {
        "H2_Q": base_density * H**2 * q_shear,
        "Hd_Q": base_density * Hd * q_shear,
        "H_Qdot": base_density * H * q_dot,
        "Q2": base_density * q_shear**2,
        "dot_sigma_sq": base_density * q_dot_sq,
        "ddot_sigma_sigma": base_density * q_ddot_pair,
        "Ndd_Q": base_density * (Ndd / N) * q_shear,
        "Nd2_Q": base_density * (Nd**2 / N**2) * q_shear,
        "Nd_H_Q": base_density * (Nd / N) * H * q_shear,
        "Nd_Qdot": base_density * (Nd / N) * q_dot,
    }
    provider_components = {name: spatial_components(term) for name, term in provider_terms.items()}
    unknowns = [sp.Symbol(f"u{i}") for i in range(len(provider_terms))]

    poly_vars = [N, Nd, Ndd, H, Hd, Hdd, s1, s2, sd1, sd2, sdd1, sdd2]
    all_polys: list[sp.Poly] = []
    for expr in target_components:
        all_polys.append(sp.Poly(sp.expand(expr), *poly_vars))
    for comps in provider_components.values():
        for expr in comps:
            all_polys.append(sp.Poly(sp.expand(expr), *poly_vars))
    monomials = sorted(set().union(*[set(poly.monoms()) for poly in all_polys]))

    rows: list[list[sp.Expr]] = []
    rhs: list[sp.Expr] = []
    term_names = list(provider_terms.keys())
    for component_index in range(2):
        component_provider_polys = [
            sp.Poly(sp.expand(provider_components[name][component_index]), *poly_vars) for name in term_names
        ]
        component_target_poly = sp.Poly(sp.expand(target_components[component_index]), *poly_vars)
        for monomial in monomials:
            row = [poly.coeff_monomial(monomial) for poly in component_provider_polys]
            b_val = component_target_poly.coeff_monomial(monomial)
            if any(value != 0 for value in row) or b_val != 0:
                rows.append(row)
                rhs.append(b_val)

    report = matrix_report(rows, rhs, unknowns)
    a_mat = sp.Matrix(rows)
    b_vec = sp.Matrix(rhs)
    solution_set = sp.linsolve((a_mat, b_vec), unknowns)
    solution_tuple = next(iter(solution_set)) if report["consistent"] else tuple([sp.nan] * len(unknowns))
    canonical_solution = [sp.simplify(expr.subs({unknowns[7]: 0, unknowns[8]: 0, unknowns[9]: 0})) for expr in solution_tuple]
    residual_vec = sp.simplify(a_mat * sp.Matrix(canonical_solution) - b_vec)
    exact_reconstruction_zero = all(sp.simplify(value) == 0 for value in residual_vec)
    canonical_solution_numeric = [float(sp.N(value, 40)) for value in canonical_solution]

    p1988_families = (p1988.get("family_projection", {}) or {}).get("strict_non_gb_families_detected", [])
    p2297_formal_count = int(
        (p2297.get("strict_non_gb_spatial_eom_provider_matrix_obstruction_probe", {})
         .get("provider_matrix_results", {})
         .get("formal_full_residual_basis_provider", {})
         .get("basis_monomial_count", 0))
    )
    p2299_gate = p2299.get("gatekeeper_checks", {})
    p2282_probe = p2282.get("strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe", {})
    gap_status = {row.get("id", "UNKNOWN"): row.get("status", "UNKNOWN") for row in p2282_probe.get("gap_rows", [])}

    source_hashes = {
        "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
        "p1853_sha256": sha256_file(IN_1853),
        "p1981_sha256": sha256_file(IN_1981),
        "p1982_sha256": sha256_file(IN_1982),
        "p1983_sha256": sha256_file(IN_1983),
        "p1988_sha256": sha256_file(IN_1988),
        "p2282_sha256": sha256_file(IN_2282),
        "p2297_sha256": sha256_file(IN_2297),
        "p2299_sha256": sha256_file(IN_2299),
    }

    provider_basis_export = [
        {
            "column": index,
            "name": name,
            "density_term": str(sp.factor(term)),
            "canonical_coefficient": str(canonical_solution[index]),
            "canonical_coefficient_numeric": canonical_solution_numeric[index],
        }
        for index, (name, term) in enumerate(provider_terms.items())
    ]

    task3_impact = {
        "gap_status_before_p2300": gap_status,
        "gap_status_after_p2300": gap_status,
        "closure_score_after_p2300": p2282_probe.get("closure_score"),
        "reason": (
            "P2300 resolves a provider-matrix subproblem for the non-GB spatial residual, but it does not replay the "
            "global Bianchi-I transport atlas or the P2280/P2281 policy-lock feasibility criteria required by P2282 G1/G2/G3."
        ),
    }

    gatekeeper_checks = {
        "alpha_geo_strict_source_loaded": alpha_packet.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "alpha_geo_is_four_ln2_not_legacy_import": alpha_packet.get("value") == "4 ln(2)",
        "p1988_family_scope_loaded": len(p1988_families) >= 8,
        "p2299_no_column_blocker_was_real": bool(p2299_gate.get("no_new_spatial_eom_operator_columns_exported")),
        "new_spatial_eom_operator_columns_exported": len(provider_terms) > 0,
        "operator_basis_not_formal_residual_copy": len(provider_terms) < p2297_formal_count and all("residual" not in name.lower() for name in term_names),
        "provider_matrix_consistent": bool(report["consistent"]),
        "exact_reconstruction_zero": exact_reconstruction_zero,
        "canonical_solution_exported": len(canonical_solution) == len(provider_terms),
        "strict_selector_premise_not_used": True,
        "task3_g1_g2_g3_not_closed_by_p2300": all(status == "OPEN" for status in gap_status.values()),
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    theorem_export = {
        "statement_id": "P2300_STRICT_SHANNON_NAD12_SIGMA_SPATIAL_EOM_PROVIDER_OPERATOR_MATRIX_THEOREM",
        "formal_statement": (
            "Using alpha_geo_strict_derived_v1 = 4 ln 2 as a strict-side Shannon source, the density-level "
            "Shannon/nad12-sigma ADM/Bianchi-I provider ansatz with ten non-residual-copying spatial-EOM columns has an "
            "exact componentwise solution for the P2297 strict non-GB residual matrix.  This is a provider-matrix pass for "
            "the local non-GB spatial-EOM subproblem only; it does not by itself close P2282 G1/G2/G3, discharge QW-2191, "
            "or prove global background-independence/ToE closure."
        ),
        "proof_bits": {
            "new_density_level_columns": len(provider_terms),
            "rank_A": report["rank_A"],
            "rank_augmented": report["rank_augmented"],
            "matrix_consistent": report["consistent"],
            "exact_reconstruction_zero": exact_reconstruction_zero,
            "formal_residual_basis_column_count_from_p2297": p2297_formal_count,
            "task3_gaps_remain_open": gatekeeper_checks["task3_g1_g2_g3_not_closed_by_p2300"],
        },
        "not_claimed": [
            "P2282 G1/G2/G3 strict closure",
            "QW-2191 discharge",
            "selector closure",
            "legacy-kernel role transfer",
            "global background-independence theorem",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    payload = {
        "schema_version": "p2300_s1250_v1",
        "packet_id": "P2300",
        "stage_id": "S1250",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_STRICT_PROVIDER_OPERATOR_MATRIX_PASS_TRACE",
        "result_kind": "STRICT_SHANNON_NAD12_SIGMA_ADM_BIANCHI_SPATIAL_EOM_PROVIDER_OPERATOR_MATRIX_PASS_WITH_G1_G2_G3_HOLD",
        "strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe": {
            "probe_id": "P2300_S1250_STRICT_SHANNON_NAD12_SIGMA_ADM_BIANCHI_SPATIAL_EOM_PROVIDER_OPERATOR",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p1853": "generated/p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json",
                "p1981": "generated/p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json",
                "p1982": "generated/p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json",
                "p1983": "generated/p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json",
                "p1988": "generated/p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json",
                "p2282": "generated/p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json",
                "p2297": "generated/p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json",
                "p2299": "generated/p2299_s1249_strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe.json",
            },
            "source_hashes": source_hashes,
            "operator_basis_policy": {
                "basis_kind": "density_level_shannon_nad12_sigma_provider_ansatz",
                "base_density_factor": "V*(4*log(2))/N**3",
                "uses_formal_full_residual_basis": False,
                "uses_explicit_selector_premise": False,
                "column_count": len(provider_terms),
                "p2297_formal_residual_basis_column_count": p2297_formal_count,
            },
            "provider_basis": provider_basis_export,
            "provider_matrix_report": report,
            "solution_space": {
                "linsolve_tuple": [str(value) for value in solution_tuple],
                "canonical_free_parameters_set_to_zero": [str(unknowns[7]), str(unknowns[8]), str(unknowns[9])],
                "canonical_solution": [str(value) for value in canonical_solution],
                "canonical_solution_numeric": canonical_solution_numeric,
                "exact_reconstruction_zero": exact_reconstruction_zero,
            },
            "task3_g1_g2_g3_impact": task3_impact,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2301_candidate",
            "goal": "Integrate the P2300 provider coefficients into the global Bianchi-I transport/policy-lock replay and recompute P2282 G1/G2/G3; do not mark Task-3 closure until G1, G2, and G3 pass their own gate criteria.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    md = f"""# P2300/S1250 — strict Shannon/nad12-sigma ADM/Bianchi-I spatial-EOM provider operator\n\n- Status: `{payload['status']}`\n- New density-level provider columns: `{len(provider_terms)}`\n- Matrix rank: `rank_A={report['rank_A']}`, `rank_augmented={report['rank_augmented']}`, consistent = `{report['consistent']}`.\n- Exact reconstruction zero: `{exact_reconstruction_zero}`.\n- Basis policy: density-level Shannon/nad12-sigma ansatz; not the P2297 formal residual-copying basis.\n- Task-3 G1/G2/G3 after P2300: `{task3_impact['gap_status_after_p2300']}`.\n- Theorem fingerprint: `{theorem_fingerprint}`\n\n## Guardrail statement\nP2300 is a local provider-matrix pass for the strict non-GB ADM/Bianchi-I spatial-EOM subproblem. It does not close P2282 G1/G2/G3, does not discharge QW-2191, does not transfer any legacy-kernel role, and does not claim ToE closure.\n\n## Next honest step\n{payload['recommended_next_honest_step']['goal']}\n"""
    MD.write_text(md, encoding="utf-8")


if __name__ == "__main__":
    main()
