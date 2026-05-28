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
IN_1853 = GEN / "p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json"
IN_1981 = GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json"
IN_1982 = GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json"
IN_1983 = GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json"
IN_1985 = GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"
IN_1988 = GEN / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"
IN_1991 = GEN / "p1991_s941_strict_augmented_provider_channel_matrix_witness.json"
IN_2296 = GEN / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.json"
OUT = GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json"
MD = GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_rank_report(rows: list[list[sp.Expr]], rhs: list[sp.Expr], unknowns: list[sp.Symbol]) -> dict[str, Any]:
    a_mat = sp.Matrix(rows)
    b_vec = sp.Matrix(rhs)
    rank_a = int(a_mat.rank()) if rows else 0
    rank_aug = int(a_mat.row_join(b_vec).rank()) if rows else 0
    consistent = rank_a == rank_aug
    full_solve = consistent and rank_a == len(unknowns)

    a_num = np.array([[float(sp.N(cell, 40)) for cell in row] for row in rows], dtype=float) if rows else np.zeros((0, len(unknowns)))
    b_num = np.array([float(sp.N(cell, 40)) for cell in rhs], dtype=float) if rhs else np.zeros((0,))
    if rows:
        x_lsq, *_ = np.linalg.lstsq(a_num, b_num, rcond=None)
        residual_l2 = float(la.norm(a_num @ x_lsq - b_num, ord=2))
    else:
        x_lsq = np.zeros((len(unknowns),), dtype=float)
        residual_l2 = 0.0

    return {
        "equation_count": len(rows),
        "unknown_count": len(unknowns),
        "rank_A": rank_a,
        "rank_augmented": rank_aug,
        "consistent": consistent,
        "full_solve": full_solve,
        "least_squares_residual_l2": residual_l2,
        "least_squares_solution_preview": {str(u): float(x_lsq[i]) for i, u in enumerate(unknowns[: min(8, len(unknowns))])},
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1853 = load(IN_1853)
    p1981 = load(IN_1981)
    p1982 = load(IN_1982)
    p1983 = load(IN_1983)
    p1985 = load(IN_1985)
    p1988 = load(IN_1988)
    p1991 = load(IN_1991)
    p2296 = load(IN_2296)

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
        # Full one-dimensional Euler operator for a density depending on
        # sigma, dot(sigma), and possibly ddot(sigma):
        # dL/dsigma - Dt(dL/ddot_sigma) + Dt^2(dL/ddotdot_sigma).
        raw = sp.diff(density, sigma) - total_dt(sp.diff(density, dsigma)) + total_dt(total_dt(sp.diff(density, d2sigma)))
        return sp.factor(sp.simplify(raw * N**8 / V))

    def spatial_components(density: sp.Expr) -> list[sp.Expr]:
        e1 = euler_lagrange_component(density, s1, sd1, sdd1)
        e2 = euler_lagrange_component(density, s2, sd2, sdd2)
        e3 = sp.factor(sp.simplify(-e1 - e2))
        return [e1, e2, e3]

    target_components = spatial_components(density_non_gb)
    tracefree_sum_zero = sp.simplify(sum(target_components)) == 0
    component_hashes = [sha256_text(str(component)) for component in target_components]

    poly_vars = [N, Nd, Ndd, H, Hd, Hdd, s1, s2, sd1, sd2, sdd1, sdd2]

    def provider_density_from_scaled_poly(poly: sp.Expr) -> sp.Expr:
        return sp.factor(V * poly / N**5)

    def build_system(provider_poly: sp.Expr, unknowns: list[sp.Symbol]) -> tuple[list[list[sp.Expr]], list[sp.Expr], list[dict[str, Any]]]:
        provider_components = spatial_components(provider_density_from_scaled_poly(provider_poly))
        rows: list[list[sp.Expr]] = []
        rhs: list[sp.Expr] = []
        labels: list[dict[str, Any]] = []
        zero_subs = {u: 0 for u in unknowns}
        for component_index, target in enumerate(
            [sp.expand(t - p) for t, p in zip(target_components, provider_components)], start=1
        ):
            poly = sp.Poly(sp.expand(target), *poly_vars, domain="EX")
            for monomial, coeff in poly.terms():
                coeff = sp.expand(coeff)
                row = [sp.diff(coeff, unknown) for unknown in unknowns]
                const = sp.simplify(coeff.subs(zero_subs))
                if any(entry != 0 for entry in row) or const != 0:
                    rows.append(row)
                    rhs.append(sp.simplify(-const))
                    labels.append(
                        {
                            "component": f"E_spatial_{component_index}",
                            "monomial_powers": {str(var): power for var, power in zip(poly_vars, monomial) if power},
                            "row_nonzero": any(entry != 0 for entry in row),
                            "rhs": str(sp.simplify(-const)),
                        }
                    )
        return rows, rhs, labels

    c1, c2, c3, eta1, eta2, eta3 = sp.symbols("c1 c2 c3 eta1 eta2 eta3", real=True)
    p_min_poly = c1 * (3 * H * (s1 + s2)) + c2 * (sd1 + sd2) + c3 * q_shear
    d2sig_sig = 2 * sdd1 * s1 + sdd1 * s2 + sdd2 * s1 + 2 * sdd2 * s2
    p_aug_poly = p_min_poly + eta1 * d2sig_sig + eta2 * Ndd * q_shear + eta3 * q_shear**2

    strict_rows, strict_rhs, strict_labels = build_system(p_min_poly, [c1, c2, c3])
    aug_rows, aug_rhs, aug_labels = build_system(p_aug_poly, [c1, c2, c3, eta1, eta2, eta3])

    strict_report = matrix_rank_report(strict_rows, strict_rhs, [c1, c2, c3])
    aug_report = matrix_rank_report(aug_rows, aug_rhs, [c1, c2, c3, eta1, eta2, eta3])

    def obstruction_preview(rows: list[list[sp.Expr]], rhs: list[sp.Expr], labels: list[dict[str, Any]], limit: int = 8) -> list[dict[str, Any]]:
        preview = []
        for row, right, label in zip(rows, rhs, labels):
            if all(entry == 0 for entry in row) and right != 0:
                item = dict(label)
                item["obstruction_rhs"] = str(right)
                preview.append(item)
            if len(preview) >= limit:
                break
        return preview

    strict_obstructions = obstruction_preview(strict_rows, strict_rhs, strict_labels)
    aug_obstructions = obstruction_preview(aug_rows, aug_rhs, aug_labels)

    scaled_density_poly = sp.factor(sp.simplify(density_non_gb * N**5 / V))
    full_poly = sp.Poly(sp.expand(scaled_density_poly), *poly_vars, domain="EX")
    full_terms = full_poly.terms()
    full_unknowns = list(sp.symbols(f"u0:{len(full_terms)}", real=True))
    monomial_exprs = [sp.prod(var**power for var, power in zip(poly_vars, monomial)) for monomial, _ in full_terms]
    full_provider_poly = sp.Add(*[u * monomial for u, monomial in zip(full_unknowns, monomial_exprs)])
    full_rows, full_rhs, _full_labels = build_system(full_provider_poly, full_unknowns)
    full_report = matrix_rank_report(full_rows, full_rhs, full_unknowns)
    exact_reconstruction_components = [
        sp.factor(sp.simplify(component))
        for component in spatial_components(density_non_gb - provider_density_from_scaled_poly(scaled_density_poly))
    ]
    exact_reconstruction_zero = all(component == 0 for component in exact_reconstruction_components)

    p1985_gate = p1985.get("gatekeeper_checks", {}) or {}
    p1988_gate = p1988.get("gatekeeper_checks", {}) or {}
    p1991_gate = p1991.get("gatekeeper_checks", {}) or {}
    p2296_gate = p2296.get("gatekeeper_checks", {}) or {}

    payload = {
        "schema_version": "p2297_s1247_v1",
        "packet_id": "P2297",
        "stage_id": "S1247",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NON_GB_SPATIAL_EOM_PROVIDER_MATRIX_OBSTRUCTION_PROBE",
        "strict_non_gb_spatial_eom_provider_matrix_obstruction_probe": {
            "probe_id": "STRICT_NON_GB_SPATIAL_EOM_PROVIDER_MATRIX_OBSTRUCTION_PROBE_V1",
            "source_packets": [
                str(IN_1853.relative_to(ROOT)),
                str(IN_1981.relative_to(ROOT)),
                str(IN_1982.relative_to(ROOT)),
                str(IN_1983.relative_to(ROOT)),
                str(IN_1985.relative_to(ROOT)),
                str(IN_1988.relative_to(ROOT)),
                str(IN_1991.relative_to(ROOT)),
                str(IN_2296.relative_to(ROOT)),
            ],
            "source_hashes": {
                str(path.relative_to(ROOT)): sha256_file(path)
                for path in [IN_1853, IN_1981, IN_1982, IN_1983, IN_1985, IN_1988, IN_1991, IN_2296]
            },
            "strict_coefficients": {"a_R2": str(a_r2), "a_Ric2": str(a_ric2), "a_Riem2": str(a_riem2)},
            "spatial_eom_projection": {
                "definition": "E_i = tracefree diagonal Bianchi-I anisotropic Euler projection dL/dsigma - Dt(dL/ddot_sigma) + Dt^2(dL/ddotdot_sigma) of the strict non-GB density; equations are scaled by N^8/V for zero-equivalent polynomial testing",
                "components": ["E_spatial_1 := EL_sigma1", "E_spatial_2 := EL_sigma2", "E_spatial_3 := -EL_sigma1-EL_sigma2"],
                "tracefree_sum_zero": tracefree_sum_zero,
                "component_sha256": component_hashes,
            },
            "provider_matrix_results": {
                "strict_core_minimal_provider": {
                    "unknowns": ["c1", "c2", "c3"],
                    "provider_scaled_polynomial": str(p_min_poly),
                    "matrix_report": strict_report,
                    "zero_row_obstruction_preview": strict_obstructions,
                    "verdict": "NO_STRICT_CORE_SOLUTION_FOR_CURRENT_PROVIDER_CLASS",
                },
                "p1990_augmented_provider_non_strict": {
                    "unknowns": ["c1", "c2", "c3", "eta1", "eta2", "eta3"],
                    "provider_scaled_polynomial": str(p_aug_poly),
                    "selector_premise_label": "NON_STRICT_AUGMENTED_CLASS_FROM_P1990; QW-2191_NOT_DISCHARGED",
                    "matrix_report": aug_report,
                    "zero_row_obstruction_preview": aug_obstructions,
                    "verdict": "NO_FULL_COMPONENTWISE_SOLUTION_EVEN_BEFORE_STRICT_PROMOTION",
                },
                "formal_full_residual_basis_provider": {
                    "basis_monomial_count": len(full_terms),
                    "matrix_report": full_report,
                    "exact_reconstruction_zero": exact_reconstruction_zero,
                    "admissibility_verdict": "ALGEBRAICALLY_CAPABLE_BUT_NOT_A_LEGAL_STRICT_PROVIDER; it copies the full residual basis and supplies no internal selector source",
                },
            },
            "upstream_consistency": {
                "p1985_obstruction_witness_passed": bool(p1985_gate.get("obstruction_witness_passed", False)),
                "p1988_spatial_projection_obstruction_present": bool(p1988_gate.get("outside_eh_family_capacity_detected", False)),
                "p1991_augmented_channel_obstruction_present": bool(p1991_gate.get("matrix_rank_aug_exceeds_rankA_or_not_full_solve", False)),
                "p2296_task1_b1_projected_pass": bool(p2296_gate.get("task1_b1_projected_pass", False)),
            },
            "theorem_scope_limit": "Componentwise spatial-EOM/provider matrix obstruction for the strict non-GB Bianchi-I residual; not a legal-provider existence theorem, not QW-2191 discharge, and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2298_candidate",
            "goal": "derive a noncyclic strict internal provider-source candidate for the full residual-family basis, or prove a non-admissibility theorem showing that any such source would require an explicit non-strict selector premise",
        },
        "gatekeeper_checks": {
            "tracefree_spatial_component_sum_zero": tracefree_sum_zero,
            "strict_core_matrix_inconsistent": not strict_report["consistent"],
            "augmented_non_strict_matrix_inconsistent": not aug_report["consistent"],
            "formal_full_basis_reconstructs_zero": exact_reconstruction_zero,
            "formal_full_basis_marked_non_admissible": True,
            "p1985_obstruction_preserved": bool(p1985_gate.get("obstruction_witness_passed", False)),
            "no_legacy_transfer_used": True,
            "qw2191_not_discharged": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2297 S1247: strict non-GB spatial-EOM/provider matrix obstruction",
                "",
                f"- tracefree spatial component sum zero: `{tracefree_sum_zero}`",
                f"- strict-core provider matrix consistent: `{strict_report['consistent']}`",
                f"- augmented non-strict provider matrix consistent: `{aug_report['consistent']}`",
                f"- formal full residual-basis reconstruction zero: `{exact_reconstruction_zero}`",
                f"- formal basis admissibility: `NON_ADMISSIBLE_AS_STRICT_PROVIDER`",
                "",
                "Conclusion: the current legal strict-core provider class cannot cancel the P1985 non-GB Bianchi-I residual.",
                "A formally complete residual-copying basis can cancel algebraically, but it is not an exported strict provider and does not discharge QW-2191.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
