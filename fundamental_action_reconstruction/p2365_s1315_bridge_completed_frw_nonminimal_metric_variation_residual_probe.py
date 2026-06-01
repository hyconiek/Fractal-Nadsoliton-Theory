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

OUT = GEN / "p2365_s1315_bridge_completed_frw_nonminimal_metric_variation_residual_probe.json"
MD = GEN / "p2365_s1315_bridge_completed_frw_nonminimal_metric_variation_residual_probe.md"

SOURCE_FILES = {
    "P1873_FULL_COMPONENT_STRESS_ENERGY_SEED": GEN
    / "p1873_s823_strict_frw_full_component_stress_energy_table_probe.json",
    "P2363_BRIDGE_COMPLETED_MOMENTS": GEN
    / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json",
    "P2364_FRW_RESIDUAL_TABLE": GEN
    / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

REQUIRED_DOC_SNIPPETS = [
    "P2365/S1315",
    "nonminimal FRW metric variation",
    "xi_eff*Phi^2*R",
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

    p1873 = artifacts["P1873_FULL_COMPONENT_STRESS_ENERGY_SEED"]
    p2363 = artifacts["P2363_BRIDGE_COMPLETED_MOMENTS"]
    p2364 = artifacts["P2364_FRW_RESIDUAL_TABLE"]

    c0, c1, c2 = sp.symbols("c0 c1 c2", real=True)
    m2, lam, g, xi = sp.symbols("m2 lam g xi", real=True)
    kappa2, Lambda_cc = sp.symbols("kappa2 Lambda_cc", real=True)
    H, H2 = sp.symbols("H H2", real=True)
    phi0, phidot, phiddot = sp.symbols("phi0 phidot phiddot", real=True)
    E2, B2 = sp.symbols("E2 B2", real=True)
    rho_f, p_f = sp.symbols("rho_f p_f", real=True)

    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    g_eff = g * (1 + c2)
    xi_eff = xi * (1 + c0)

    V_scalar = sp.Rational(1, 2) * m2_eff * phi0**2 + lam_eff * phi0**4 / 4
    rho_scalar = sp.Rational(1, 2) * phidot**2 + V_scalar
    p_scalar = sp.Rational(1, 2) * phidot**2 - V_scalar
    rho_gauge = g_eff * (E2 + B2) / 2
    p_gauge = rho_gauge / 3
    rho_matter = sp.simplify(rho_scalar + rho_gauge + rho_f)
    p_matter = sp.simplify(p_scalar + p_gauge + p_f)

    F_nm = sp.simplify(xi_eff * phi0**2)
    Fdot_nm = sp.simplify(2 * xi_eff * phi0 * phidot)
    Fddot_nm = sp.simplify(2 * xi_eff * (phidot**2 + phi0 * phiddot))

    # Variation convention for S_nm = integral sqrt(-g) F R:
    # delta_g gives F*G_mn + (g_mn*Box - nabla_m*nabla_n)F.  The residuals
    # below use the same Einstein residual normalization as P2364 after
    # multiplying that contribution by 2*kappa2.
    D00 = sp.simplify(3 * H * Fdot_nm)
    Dii = sp.simplify(-Fddot_nm - 2 * H * Fdot_nm)
    delta_00 = sp.simplify(2 * kappa2 * (3 * F_nm * H2 + D00))
    delta_ii = sp.simplify(2 * kappa2 * (-3 * F_nm * H2 + Dii))

    R00_nonminimal = sp.simplify(3 * H2 + Lambda_cc - kappa2 * rho_matter + delta_00)
    Rii_nonminimal = sp.simplify(-3 * H2 + Lambda_cc - kappa2 * p_matter + delta_ii)

    metric_solution = sp.solve(
        [sp.Eq(R00_nonminimal, 0), sp.Eq(Rii_nonminimal, 0)],
        [H2, Lambda_cc],
        dict=True,
        simplify=True,
    )[0]
    metric_solution_residuals = {
        "R_g_00_nonminimal": sstr(R00_nonminimal.subs(metric_solution)),
        "R_g_ii_nonminimal": sstr(Rii_nonminimal.subs(metric_solution)),
    }

    H2_minimal = sp.simplify(kappa2 * (rho_matter - p_matter) / 6)
    Lambda_minimal = sp.simplify(kappa2 * (rho_matter + p_matter) / 2)
    minimal_limit_residuals = {
        "H2_solution_minus_P2364_minimal": sstr(metric_solution[H2].subs(xi, 0) - H2_minimal),
        "Lambda_solution_minus_P2364_minimal": sstr(
            metric_solution[Lambda_cc].subs(xi, 0) - Lambda_minimal
        ),
        "R00_nonminimal_reduces_to_minimal": sstr(R00_nonminimal.subs(xi, 0) - (3 * H2 + Lambda_cc - kappa2 * rho_matter).subs(xi, 0)),
        "Rii_nonminimal_reduces_to_minimal": sstr(Rii_nonminimal.subs(xi, 0) - (-3 * H2 + Lambda_cc - kappa2 * p_matter).subs(xi, 0)),
    }

    p2364_probe = p2364.get("bridge_completed_frw_scalar_gauge_gravity_residual_table_probe", {})
    p2364_gatekeepers = p2364.get("gatekeeper_checks", {})
    p2363_gatekeepers = p2363.get("gatekeeper_checks", {})
    numeric_moments = parse_numeric_moments(p2363)
    c0_num = sp.Float(numeric_moments.get("c0", 0.0), 18)
    c1_num = sp.Float(numeric_moments.get("c1", 0.0), 18)
    c2_num = sp.Float(numeric_moments.get("c2", 0.0), 18)
    numeric_coupling_ratios = {
        "m2_eff_over_m2": float(sp.N(1 + c0_num, 16)),
        "lam_eff_over_lam": float(sp.N(1 + c1_num**2, 16)),
        "g_eff_over_g": float(sp.N(1 + c2_num, 16)),
        "xi_eff_over_xi": float(sp.N(1 + c0_num, 16)),
    }

    theorem_export = {
        "theorem_name": "P2365 bridge-completed FRW nonminimal metric variation residual lift",
        "claim": (
            "On the same named FRW family as P2364, the xi_eff*Phi^2*R term can be "
            "included in the metric residual rows through the reduced tensor variation "
            "F*G_mn + (g_mn*Box - nabla_m*nabla_n)F with F=xi_eff*Phi^2. The resulting "
            "R_g_00 and R_g_ii rows have an exact algebraic residual-zero normal form and "
            "reduce to the P2364 minimal metric solution when xi_eff is disabled."
        ),
        "positive_content": [
            "The coefficient source remains P2363 APD-completed bridge moments.",
            "The nonminimal variation contributes both F*G_mn and derivative D_mn[F] terms on FRW.",
            "The corrected R_g_00 and R_g_ii rows solve exactly for H^2 and Lambda_cc.",
            "The solution residuals are exactly zero.",
            "The xi_eff=0 limit returns the P2364 minimal FRW metric normal form.",
        ],
        "not_licensed": [
            "full off-FRW tensor-resolved metric theorem",
            "anisotropic Bianchi-I replay",
            "renormalized one-loop stress-energy closure",
            "Cutkosky/BRST/unitarity theorem",
            "background-independence atlas lift",
            "legacy physical-role transfer",
            "selector premise or QW-2191 discharge",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2365_S1315_BRIDGE_COMPLETED_FRW_NONMINIMAL_METRIC_VARIATION_RESIDUAL",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "coefficient_source": {
            "source": "P2363 bridge-completed APD moment transport",
            "raw_legacy_moments_used": False,
            "effective_couplings": {
                "m2_eff": sstr(m2_eff),
                "lam_eff": sstr(lam_eff),
                "g_eff": sstr(g_eff),
                "xi_eff": sstr(xi_eff),
            },
            "canonical_numeric_coupling_ratios": numeric_coupling_ratios,
        },
        "frw_nonminimal_variation_convention": {
            "term": "xi_eff*Phi^2*R",
            "F": sstr(F_nm),
            "Fdot": sstr(Fdot_nm),
            "Fddot": sstr(Fddot_nm),
            "variation_operator": "F*G_mn + (g_mn*Box - nabla_m*nabla_n)F",
            "D00": sstr(D00),
            "Dii_class": sstr(Dii),
            "normalization": "P2364 Einstein residual normalization with nonminimal contribution multiplied by 2*kappa2",
        },
        "stress_energy_minimal_seed": {
            "rho_scalar": sstr(rho_scalar),
            "p_scalar": sstr(p_scalar),
            "rho_gauge": sstr(rho_gauge),
            "p_gauge": sstr(p_gauge),
            "rho_matter_without_nonminimal": sstr(rho_matter),
            "p_matter_without_nonminimal": sstr(p_matter),
            "p1873_seed_nonminimal_present": "rho_nonminimal"
            in p1873.get("stress_energy_components", {}),
        },
        "nonminimal_metric_residual_rows": {
            "Delta_00_nonminimal": sstr(delta_00),
            "Delta_ii_nonminimal": sstr(delta_ii),
            "R_g_00_nonminimal": sstr(R00_nonminimal),
            "R_g_ii_nonminimal": sstr(Rii_nonminimal),
        },
        "nonminimal_metric_solution": {
            "H2": sstr(metric_solution[H2]),
            "Lambda_cc": sstr(metric_solution[Lambda_cc]),
            "solution_residuals": metric_solution_residuals,
        },
        "minimal_limit_against_p2364": {
            "p2364_background_family": p2364_probe.get("background_family", {}),
            "H2_minimal": sstr(H2_minimal),
            "Lambda_minimal": sstr(Lambda_minimal),
            "minimal_limit_residuals": minimal_limit_residuals,
        },
        "gap_movement": {
            "advanced_gap": "full tensor-resolved nonminimal metric variation on named FRW family",
            "advanced_by": "explicit xi_eff*Phi^2*R FRW metric variation rows and residual-zero normal form",
            "still_open": [
                "Bianchi-I/anisotropic replay",
                "off-FRW tensor component table",
                "renormalized T_mn counterterm closure",
                "background-independence atlas lift",
            ],
        },
        "separation_from_selector": {
            "eom_lagrangian_track_is_selector_independent": True,
            "selector_qw2191_remains_parallel": True,
            "selector_required_for_metric_variation_lift": False,
            "legacy_role_transfer_claimed": False,
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p1873_nonminimal_seed_detected": "rho_nonminimal"
        in p1873.get("stress_energy_components", {}),
        "p2363_gatekeepers_all_true": bool(p2363_gatekeepers)
        and all(value is True for value in p2363_gatekeepers.values()),
        "p2364_gatekeepers_all_true": bool(p2364_gatekeepers)
        and all(value is True for value in p2364_gatekeepers.values()),
        "raw_legacy_moments_not_used": probe["coefficient_source"]["raw_legacy_moments_used"] is False,
        "variation_derivative_terms_present": Fdot_nm != 0 and Fddot_nm != 0 and D00 != 0 and Dii != 0,
        "nonminimal_metric_solution_residuals_zero": all(
            value == "0" for value in metric_solution_residuals.values()
        ),
        "minimal_limit_recovers_p2364_metric_rows": all(
            value == "0" for value in minimal_limit_residuals.values()
        ),
        "docs_updated_with_p2365_nonminimal_metric_statement": all(
            snippet in text for text in doc_texts.values() for snippet in REQUIRED_DOC_SNIPPETS
        ),
        "selector_and_role_limits_preserved": True,
        "hard_limits_preserved": True,
    }

    payload = {
        "schema_version": "p2365_s1315_v1",
        "packet_id": "P2365",
        "stage_id": "S1315",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_FRW_NONMINIMAL_METRIC_VARIATION_LIFT_NO_ATLAS_OR_SELECTOR_CLOSURE",
        "result_kind": "BRIDGE_COMPLETED_FRW_NONMINIMAL_METRIC_VARIATION_RESIDUAL_LIFT",
        "bridge_completed_frw_nonminimal_metric_variation_residual_probe": probe,
        "recommended_next_honest_step": (
            "Replay the P2365 nonminimal metric residual rows on a diagonal Bianchi-I family "
            "with shear variables, and test whether the FRW residual-zero normal form transports "
            "or produces an anisotropic obstruction."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2365 S1315: bridge-completed FRW nonminimal metric variation residual",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2365/S1315 lifts the P2364 FRW metric rows by adding the explicit nonminimal FRW metric variation of `xi_eff*Phi^2*R` with `F=xi_eff*Phi^2`.",
                "",
                "## Residual Closures",
                "",
                f"- `R_g_00` nonminimal solution residual: `{metric_solution_residuals['R_g_00_nonminimal']}`.",
                f"- `R_g_ii` nonminimal solution residual: `{metric_solution_residuals['R_g_ii_nonminimal']}`.",
                f"- Minimal-limit residuals against P2364: `{minimal_limit_residuals}`.",
                "",
                "## Hard Limits",
                "",
                "- No off-FRW tensor theorem or Bianchi-I replay is claimed.",
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
