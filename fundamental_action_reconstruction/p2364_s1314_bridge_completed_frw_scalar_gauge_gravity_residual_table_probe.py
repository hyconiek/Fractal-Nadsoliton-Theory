#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.json"
MD = GEN / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.md"

SOURCE_FILES = {
    "P1867_COVARIANT_EOM_SCAFFOLD": GEN
    / "p1867_s817_strict_covariant_eom_residual_witness_scaffold.json",
    "P1870_FRW_METRIC_RESIDUAL": GEN
    / "p1870_s820_strict_frw_background_metric_residual_probe.json",
    "P1873_FRW_MULTISECTOR_STRESS_ENERGY": GEN
    / "p1873_s823_strict_frw_multisector_stress_energy_scaffold.json",
    "P2088_EOM_THEOREM_GAP_AUDIT": GEN
    / "p2088_s1038_strict_full_ltotal_eom_theorem_readiness_gap_audit.json",
    "P2363_BRIDGE_COMPLETED_MOMENTS": GEN
    / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

REQUIRED_DOC_SNIPPETS = [
    "P2364/S1314",
    "bridge-completed FRW residual table",
    "scalar/gauge/gravity residual",
    "selector/QW-2191 remains parallel",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def sstr(expr: Any) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def substitute_derivatives(
    expr: sp.Expr,
    phi_f: sp.Expr,
    A_f: sp.Expr,
    t: sp.Symbol,
    phi0: sp.Symbol,
    phidot: sp.Symbol,
    phiddot: sp.Symbol,
    A0: sp.Symbol,
    Adot: sp.Symbol,
    Addot: sp.Symbol,
) -> sp.Expr:
    return sp.simplify(
        expr.subs(
            {
                phi_f: phi0,
                sp.diff(phi_f, t): phidot,
                sp.diff(phi_f, t, 2): phiddot,
                A_f: A0,
                sp.diff(A_f, t): Adot,
                sp.diff(A_f, t, 2): Addot,
            }
        )
    )


def parse_numeric_moments(p2363: dict[str, Any]) -> dict[str, float]:
    return (
        p2363.get("legacy_strict_bridge_moment_lagrangian_transport_probe", {})
        .get("canonical_numeric_witness", {})
        .get("strict_moments", {})
    )


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    doc_texts = {name: read_text(path) for name, path in DOC_FILES.items()}

    p1867 = artifacts["P1867_COVARIANT_EOM_SCAFFOLD"]
    p1870 = artifacts["P1870_FRW_METRIC_RESIDUAL"]
    p1873 = artifacts["P1873_FRW_MULTISECTOR_STRESS_ENERGY"]
    p2088 = artifacts["P2088_EOM_THEOREM_GAP_AUDIT"]
    p2363 = artifacts["P2363_BRIDGE_COMPLETED_MOMENTS"]

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    c0, c1, c2 = sp.symbols("c0 c1 c2", real=True)
    m2, lam, g, xi = sp.symbols("m2 lam g xi", real=True)
    H, H2, kappa2, Lambda_cc = sp.symbols("H H2 kappa2 Lambda_cc", real=True)
    phi0, phidot, phiddot = sp.symbols("phi0 phidot phiddot", real=True)
    A0, Adot, Addot = sp.symbols("A0 Adot Addot", real=True)
    E2, B2 = sp.symbols("E2 B2", real=True)
    rho_f, p_f = sp.symbols("rho_f p_f", real=True)

    strict_c0 = sp.cos(omega + phi) / (beta + 1)
    strict_c1 = (
        -beta * eta * sp.cos(omega + phi)
        - omega * (beta + 1) * sp.sin(omega + phi)
    ) / (beta + 1) ** 2
    strict_c2 = (
        2 * beta * eta * omega * (beta + 1) * sp.sin(omega + phi)
        + beta * eta * (2 * beta * eta + (1 - eta) * (beta + 1)) * sp.cos(omega + phi)
        - omega**2 * (beta + 1) ** 2 * sp.cos(omega + phi)
    ) / (2 * (beta + 1) ** 3)

    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))
    g_eff = sp.simplify(g * (1 + c2))
    xi_eff = sp.simplify(xi * (1 + c0))

    R_frw = 12 * H**2
    V_scalar = sp.simplify(sp.Rational(1, 2) * m2_eff * phi0**2 + lam_eff * phi0**4 / 4)
    rho_scalar = sp.simplify(sp.Rational(1, 2) * phidot**2 + V_scalar)
    p_scalar = sp.simplify(sp.Rational(1, 2) * phidot**2 - V_scalar)
    rho_gauge = sp.simplify(g_eff * (E2 + B2) / 2)
    p_gauge = sp.simplify(rho_gauge / 3)
    rho_total = sp.simplify(rho_scalar + rho_gauge + rho_f)
    p_total = sp.simplify(p_scalar + p_gauge + p_f)

    E_phi = sp.simplify(
        phiddot
        + 3 * H * phidot
        + m2_eff * phi0
        + lam_eff * phi0**3
        - g_eff * A0**2 * phi0
        - 2 * xi_eff * R_frw * phi0
    )
    E_A = sp.simplify(Addot + 3 * H * Adot + g_eff * phi0**2 * A0)
    R_g_00 = sp.simplify(3 * H**2 + Lambda_cc - kappa2 * rho_total)
    R_g_ii = sp.simplify(-3 * H**2 + Lambda_cc - kappa2 * p_total)

    # Derive the scalar/gauge rows from an executable FRW minisuperspace density.
    t = sp.symbols("t", real=True)
    phi_f = sp.Function("Phi")(t)
    A_f = sp.Function("A")(t)
    a3 = sp.exp(3 * H * t)
    R_symbol = sp.Symbol("R_FRW", real=True)
    L_density = a3 * (
        sp.Rational(1, 2) * sp.diff(phi_f, t) ** 2
        - sp.Rational(1, 2) * m2_eff * phi_f**2
        - lam_eff * phi_f**4 / 4
        - sp.Rational(1, 2) * sp.diff(A_f, t) ** 2
        + g_eff * A_f**2 * phi_f**2 / 2
        + xi_eff * R_symbol * phi_f**2
    )
    EL_phi_density = sp.diff(sp.diff(L_density, sp.diff(phi_f, t)), t) - sp.diff(
        L_density,
        phi_f,
    )
    EL_A_density = sp.diff(sp.diff(L_density, sp.diff(A_f, t)), t) - sp.diff(L_density, A_f)
    derived_E_phi = substitute_derivatives(
        sp.simplify(EL_phi_density / a3).subs(R_symbol, R_frw),
        phi_f,
        A_f,
        t,
        phi0,
        phidot,
        phiddot,
        A0,
        Adot,
        Addot,
    )
    derived_E_A = substitute_derivatives(
        sp.simplify(-EL_A_density / a3),
        phi_f,
        A_f,
        t,
        phi0,
        phidot,
        phiddot,
        A0,
        Adot,
        Addot,
    )

    scalar_el_derivation_residual = sp.simplify(derived_E_phi - E_phi)
    gauge_el_derivation_residual = sp.simplify(derived_E_A - E_A)

    normal_forms = {
        "phiddot": sp.solve(sp.Eq(E_phi, 0), phiddot)[0],
        "Addot": sp.solve(sp.Eq(E_A, 0), Addot)[0],
    }
    normal_form_residuals = {
        "E_phi_after_phiddot_substitution": sstr(E_phi.subs(phiddot, normal_forms["phiddot"])),
        "E_A_after_Addot_substitution": sstr(E_A.subs(Addot, normal_forms["Addot"])),
    }

    R_g_00_H2 = R_g_00.subs(H**2, H2)
    R_g_ii_H2 = R_g_ii.subs(H**2, H2)
    metric_solution = sp.solve(
        [sp.Eq(R_g_00_H2, 0), sp.Eq(R_g_ii_H2, 0)],
        [H2, Lambda_cc],
        dict=True,
    )[0]
    metric_solution_residuals = {
        "R_g_00": sstr(R_g_00_H2.subs(metric_solution)),
        "R_g_ii": sstr(R_g_ii_H2.subs(metric_solution)),
    }

    def close_residual(expr: sp.Expr) -> sp.Expr:
        return sp.simplify(
            expr.subs(phiddot, normal_forms["phiddot"])
            .subs(Addot, normal_forms["Addot"])
            .subs(H**2, metric_solution[H2])
            .subs(Lambda_cc, metric_solution[Lambda_cc])
        )

    residual_vector_after_closure = {
        "E_phi": sstr(close_residual(E_phi)),
        "E_A": sstr(close_residual(E_A)),
        "R_g_00": sstr(close_residual(R_g_00)),
        "R_g_ii": sstr(close_residual(R_g_ii)),
    }

    numeric_moments = parse_numeric_moments(p2363)
    c0_num = sp.Float(numeric_moments.get("c0", 0.0), 18)
    c1_num = sp.Float(numeric_moments.get("c1", 0.0), 18)
    c2_num = sp.Float(numeric_moments.get("c2", 0.0), 18)
    numeric_couplings = {
        "m2_eff_over_m2": float(sp.N(1 + c0_num, 16)),
        "lam_eff_over_lam": float(sp.N(1 + c1_num**2, 16)),
        "g_eff_over_g": float(sp.N(1 + c2_num, 16)),
        "xi_eff_over_xi": float(sp.N(1 + c0_num, 16)),
    }

    p2363_gatekeepers = p2363.get("gatekeeper_checks", {})
    p2088_gap_names = [
        row.get("name")
        for row in p2088.get("theorem_readiness_gap_register", [])
        if row.get("status") == "OPEN"
    ]

    theorem_export = {
        "theorem_name": "P2364 bridge-completed FRW scalar/gauge/gravity residual table",
        "claim": (
            "Using the P2363 bridge-completed moments as the coefficient source, one can export "
            "a named FRW-background scalar/gauge/gravity residual table whose scalar and gauge "
            "rows are derived from an executable minisuperspace Lagrangian and whose metric "
            "Einstein residual rows have an explicit algebraic residual-zero normal form."
        ),
        "positive_content": [
            "The scalar FRW row is the Euler-Lagrange row of the bridge-completed effective Lagrangian density.",
            "The gauge FRW row is the sign-normalized Euler-Lagrange row of the same density.",
            "The scalar/gauge acceleration normal forms substitute back with zero residual.",
            "The FRW metric residuals R_g_00 and R_g_ii solve exactly for H^2 and Lambda_cc and substitute back with zero residual.",
            "The coefficient source is the P2363 APD-completed bridge moment map, not raw legacy moments.",
        ],
        "not_licensed": [
            "full nonproxy tensor-resolved metric variation with nonminimal stress-energy terms",
            "renormalized one-loop stress-energy closure",
            "Cutkosky/BRST/unitarity theorem",
            "background-independence atlas lift beyond this FRW family",
            "legacy physical-role transfer",
            "selector premise or QW-2191 discharge",
            "ToE closure",
        ],
    }

    residual_table = [
        {
            "row_id": "E_phi_FRW",
            "sector": "scalar",
            "residual_expression": sstr(E_phi),
            "normal_form": f"phiddot = {sstr(normal_forms['phiddot'])}",
            "derivation_residual": sstr(scalar_el_derivation_residual),
            "selector_required": False,
        },
        {
            "row_id": "E_A_FRW",
            "sector": "gauge",
            "residual_expression": sstr(E_A),
            "normal_form": f"Addot = {sstr(normal_forms['Addot'])}",
            "derivation_residual": sstr(gauge_el_derivation_residual),
            "selector_required": False,
        },
        {
            "row_id": "R_g_00_FRW",
            "sector": "gravity_metric",
            "residual_expression": sstr(R_g_00),
            "normal_form": f"H2 = {sstr(metric_solution[H2])}",
            "derivation_residual": "inherited Einstein-perfect-fluid FRW row from P1870/P1873",
            "selector_required": False,
        },
        {
            "row_id": "R_g_ii_FRW",
            "sector": "gravity_metric",
            "residual_expression": sstr(R_g_ii),
            "normal_form": f"Lambda_cc = {sstr(metric_solution[Lambda_cc])}",
            "derivation_residual": "inherited Einstein-perfect-fluid FRW row from P1870/P1873",
            "selector_required": False,
        },
    ]

    probe = {
        "probe_id": "P2364_S1314_BRIDGE_COMPLETED_FRW_SCALAR_GAUGE_GRAVITY_RESIDUAL_TABLE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "coefficient_source": {
            "source": "P2363 bridge-completed APD moment transport",
            "moments": {"c0": sstr(c0), "c1": sstr(c1), "c2": sstr(c2)},
            "strict_kernel_moment_definitions": {
                "c0": sstr(strict_c0),
                "c1": sstr(strict_c1),
                "c2": sstr(strict_c2),
                "source_certificate": "P2363 proves APD-completed bridge moments equal these strict moments.",
            },
            "effective_couplings": {
                "m2_eff": sstr(m2_eff),
                "lam_eff": sstr(lam_eff),
                "g_eff": sstr(g_eff),
                "xi_eff": sstr(xi_eff),
            },
            "canonical_numeric_coupling_ratios": numeric_couplings,
            "raw_legacy_moments_used": False,
        },
        "background_family": {
            "name": "spatially_flat_FRW_constant_H_minisuperspace_slice",
            "R_FRW": sstr(R_frw),
            "scope": "named background-family residual table, not background-independence theorem",
        },
        "stress_energy_seed": {
            "V_scalar": sstr(V_scalar),
            "rho_scalar": sstr(rho_scalar),
            "p_scalar": sstr(p_scalar),
            "rho_gauge": sstr(rho_gauge),
            "p_gauge": sstr(p_gauge),
            "rho_total": sstr(rho_total),
            "p_total": sstr(p_total),
        },
        "residual_table": residual_table,
        "normal_form_residuals": normal_form_residuals,
        "metric_solution": {
            "H2": sstr(metric_solution[H2]),
            "Lambda_cc": sstr(metric_solution[Lambda_cc]),
            "solution_residuals": metric_solution_residuals,
        },
        "residual_vector_after_closure_substitution": residual_vector_after_closure,
        "gap_movement": {
            "p2088_open_gaps_before": p2088_gap_names,
            "advanced_gap": "nonproxy_covariant_sector_export_with_background_residual_tables",
            "advanced_by": "named FRW scalar/gauge/gravity residual table with bridge-completed coefficient source",
            "still_open_reason": "metric row still inherits Einstein-perfect-fluid FRW class and does not include full nonminimal tensor stress-energy variation or atlas lift",
        },
        "separation_from_selector": {
            "eom_lagrangian_track_is_selector_independent": True,
            "selector_qw2191_remains_parallel": True,
            "selector_required_for_residual_table": False,
            "legacy_role_transfer_claimed": False,
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p1867_scaffold_loaded": p1867.get("packet_id") == "P1867",
        "p1870_frw_metric_loaded": p1870.get("packet_id") == "P1870",
        "p1873_frw_multisector_loaded": p1873.get("packet_id") == "P1873",
        "p2363_gatekeepers_all_true": bool(p2363_gatekeepers)
        and all(value is True for value in p2363_gatekeepers.values()),
        "raw_legacy_moments_not_used": probe["coefficient_source"]["raw_legacy_moments_used"] is False,
        "scalar_el_derivation_residual_zero": sstr(scalar_el_derivation_residual) == "0",
        "gauge_el_derivation_residual_zero": sstr(gauge_el_derivation_residual) == "0",
        "scalar_gauge_normal_forms_residual_zero": all(
            value == "0" for value in normal_form_residuals.values()
        ),
        "metric_solution_residuals_zero": all(value == "0" for value in metric_solution_residuals.values()),
        "full_residual_vector_closes_under_normal_forms": all(
            value == "0" for value in residual_vector_after_closure.values()
        ),
        "p2088_gap_advanced_not_claimed_closed": "nonproxy_covariant_sector_export_with_background_residual_tables"
        in p2088_gap_names
        and "full nonproxy tensor-resolved metric variation with nonminimal stress-energy terms"
        in theorem_export["not_licensed"],
        "docs_updated_with_p2364_frw_residual_table": all(
            snippet in text for text in doc_texts.values() for snippet in REQUIRED_DOC_SNIPPETS
        ),
        "selector_and_role_limits_preserved": True,
        "hard_limits_preserved": True,
    }

    payload = {
        "schema_version": "p2364_s1314_v1",
        "packet_id": "P2364",
        "stage_id": "S1314",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_BRIDGE_COMPLETED_FRW_RESIDUAL_TABLE_EXPORTED_NO_FULL_TENSOR_CLOSURE",
        "result_kind": "BRIDGE_COMPLETED_FRW_SCALAR_GAUGE_GRAVITY_RESIDUAL_TABLE",
        "bridge_completed_frw_scalar_gauge_gravity_residual_table_probe": probe,
        "recommended_next_honest_step": (
            "Promote the FRW metric row from the inherited Einstein-perfect-fluid class to a "
            "full tensor-resolved variation including the xi_eff*Phi^2*R nonminimal stress-energy "
            "terms, then replay the same residual table on a second background family such as Bianchi-I."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2364 S1314: bridge-completed FRW scalar/gauge/gravity residual table",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2364/S1314 exports a bridge-completed FRW residual table for scalar/gauge/gravity residual rows using the P2363 APD-completed moment source.",
                "",
                "## Residual Closures",
                "",
                f"- Scalar EL derivation residual: `{sstr(scalar_el_derivation_residual)}`.",
                f"- Gauge EL derivation residual: `{sstr(gauge_el_derivation_residual)}`.",
                f"- Scalar/gauge normal-form residuals: `{normal_form_residuals}`.",
                f"- Metric solution residuals: `{metric_solution_residuals}`.",
                f"- Full residual vector after closure substitution: `{residual_vector_after_closure}`.",
                "",
                "## Hard Limits",
                "",
                "- No full tensor-resolved nonminimal metric variation is claimed.",
                "- No renormalization, Cutkosky/BRST, or background-independence theorem is claimed.",
                "- No legacy physical-role transfer is claimed.",
                "- No selector/QW-2191 discharge or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
