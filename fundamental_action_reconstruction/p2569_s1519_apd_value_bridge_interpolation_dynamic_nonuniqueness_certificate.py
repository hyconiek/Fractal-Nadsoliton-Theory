#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate.json"
MD = GEN / "p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
    "P2568_PHASE_FREQUENCY_SELECTOR": GEN / "p2568_s1518_phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate.json",
}
DOMAIN = list(range(12))
LAMBDA_VALUES = ["0", "1e-9", "-1e-9", "1e-12"]
PROBE_POINTS = [sp.Rational(1, 2), sp.Rational(11, 2), sp.Rational(23, 2)]
NEGATIVE_EXPORT_FLAGS = [
    "apd_interpolation_dynamic_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2569|S1519|APD interpolation dynamic nonuniqueness|APD value bridge interpolation|A/P/D interpolation nonuniqueness",
        "intended_research_nonduplication": "APD.*interpolation|interpolation.*APD|APD.*nonuniqueness|finite APD.*dynamic|dynamic source.*APD|vanishing polynomial.*APD",
        "apd_precursors": "P2416|S1366|APD multiplicative bridge|APD value bridge|A/P/D|strict_dynamical_source_for_A_P_D",
        "frontier_precursors": "P2561|S1511|post-damping residual bridge|P2568|S1518|strict_phase_frequency_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def apd_rows(p2416: dict[str, Any]) -> list[dict[str, Any]]:
    cert = p2416.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {})
    rows = cert.get("pointwise_factor_rows", [])
    return [{"d": int(row["d"]), "apd_value": float(row["factor_product"])} for row in rows]


def interpolation_payload(rows: list[dict[str, Any]]) -> dict[str, Any]:
    x = sp.symbols("x")
    points = [(sp.Integer(row["d"]), sp.Rational(str(row["apd_value"]))) for row in rows]
    base_poly = sp.interpolate(points, x)
    vanish = sp.prod(x - d for d in DOMAIN)
    base_derivative = sp.diff(base_poly, x)
    vanish_derivative = sp.diff(vanish, x)
    family_rows = []
    for lambda_text in LAMBDA_VALUES:
        lam = sp.Rational(lambda_text)
        poly = base_poly + lam * vanish
        derivative = sp.diff(poly, x)
        node_residuals = [float(poly.subs(x, d) - sp.Rational(str(rows[d]["apd_value"]))) for d in DOMAIN]
        probe_values = []
        for point in PROBE_POINTS:
            probe_values.append({
                "x": str(point),
                "value": float(poly.subs(x, point)),
                "derivative": float(derivative.subs(x, point)),
                "base_value": float(base_poly.subs(x, point)),
                "base_derivative": float(base_derivative.subs(x, point)),
                "value_delta_from_base": float(lam * vanish.subs(x, point)),
                "derivative_delta_from_base": float(lam * vanish_derivative.subs(x, point)),
            })
        family_rows.append({
            "lambda": lambda_text,
            "max_abs_node_residual": max(abs(value) for value in node_residuals),
            "probe_values": probe_values,
            "same_finite_apd_values_as_base": max(abs(value) for value in node_residuals) == 0.0,
            "different_off_node_dynamics_from_base": any(abs(row["value_delta_from_base"]) > 0.0 or abs(row["derivative_delta_from_base"]) > 0.0 for row in probe_values),
        })
    nonzero = [row for row in family_rows if row["lambda"] != "0"]
    return {
        "sympy_version": sp.__version__,
        "domain": DOMAIN,
        "base_interpolation_degree": int(sp.degree(base_poly, gen=x)),
        "vanishing_polynomial_degree": int(sp.degree(vanish, gen=x)),
        "vanishing_polynomial_leading_coefficient": str(sp.LC(vanish, x)),
        "probe_points": [str(point) for point in PROBE_POINTS],
        "family_lambda_values": LAMBDA_VALUES,
        "family_rows": family_rows,
        "all_family_members_preserve_apd_nodes": all(row["same_finite_apd_values_as_base"] for row in family_rows),
        "nonzero_family_members_change_off_node_dynamics": all(row["different_off_node_dynamics_from_base"] for row in nonzero),
        "interpolation_family_cardinality_certified": "infinite one-parameter family q_lambda(x)=q_interp(x)+lambda*prod_{d=0}^{11}(x-d)",
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2568_payload = load_json(SOURCE_FILES["P2568_PHASE_FREQUENCY_SELECTOR"])
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    p2568 = theorem(p2568_payload, "phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate")
    rows = apd_rows(p2416_payload)
    interpolation = interpolation_payload(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2569_T1_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate",
        "audited_chain": ["P2416/S1366", "P2561/S1511", "P2568/S1518"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2561_apd_residual_atom_inherited": "strict_dynamical_source_for_A_P_D" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "p2568_phase_frequency_still_unsourced_inherited": p2568.get("strict_phase_frequency_source_exported") is False,
        "apd_node_rows": rows,
        "interpolation_nonuniqueness_audit": interpolation,
        "finite_apd_values_do_not_select_dynamic_law": interpolation["all_family_members_preserve_apd_nodes"] and interpolation["nonzero_family_members_change_off_node_dynamics"],
        "apd_value_bridge_remains_below_dynamic_source": True,
        "recommended_next_honest_step": (
            "Do not promote P2416 finite APD value exactness to strict A/P/D dynamics. The next honest step is to add a strict dynamical equation or regularity/minimal-action principle for A/P/D and test whether it rejects the vanishing-polynomial family; otherwise APD remains a value-level bridge component, not a source theorem."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2561_apd_atom_inherited": theorem_export["p2561_apd_residual_atom_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "base_degree_11": interpolation["base_interpolation_degree"] == 11,
        "vanishing_degree_12": interpolation["vanishing_polynomial_degree"] == 12,
        "all_family_members_preserve_nodes": interpolation["all_family_members_preserve_apd_nodes"],
        "nonzero_family_changes_off_node": interpolation["nonzero_family_members_change_off_node_dynamics"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2569",
        "stage_id": "S1519",
        "status": "P2569_APD_VALUE_BRIDGE_INTERPOLATION_DYNAMIC_NONUNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
                "P2568_PHASE_FREQUENCY_SELECTOR": sha256_json(p2568_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate"]["theorem_export"]
    interpolation = t["interpolation_nonuniqueness_audit"]
    lines = [
        "# P2569/S1519 APD value-bridge interpolation dynamic nonuniqueness certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2416 APD value bridge inherited: `{t['p2416_apd_value_bridge_inherited']}`.",
        f"- P2561 APD residual atom inherited: `{t['p2561_apd_residual_atom_inherited']}`.",
        f"- Base interpolation degree: `{interpolation['base_interpolation_degree']}`.",
        f"- Vanishing polynomial degree: `{interpolation['vanishing_polynomial_degree']}`.",
        f"- All family members preserve APD nodes: `{interpolation['all_family_members_preserve_apd_nodes']}`.",
        f"- Nonzero family members change off-node dynamics: `{interpolation['nonzero_family_members_change_off_node_dynamics']}`.",
        f"- Finite APD values select dynamic law: `{not t['finite_apd_values_do_not_select_dynamic_law']}`.", "",
        "## Interpretation", "",
        "The finite APD bridge values admit an infinite interpolation family `q_lambda(x)=q_interp(x)+lambda*prod_{d=0}^{11}(x-d)`.  Every member agrees with the audited APD values on the finite bridge domain, while nonzero lambda changes off-node values/derivatives.  Therefore finite APD exactness is not a strict dynamical source.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD interpolation dynamic source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2569/S1519` audits the residual `strict_dynamical_source_for_A_P_D` atom by separating P2416 finite APD value exactness from actual dynamics.  The `12` APD node values admit a degree-`11` interpolant, but adding `lambda*prod_{d=0}^{11}(x-d)` gives an infinite family that preserves every audited node value while changing off-node values and derivatives.  Thus APD value-level bridge bookkeeping is not a strict A/P/D dynamical source.
""".strip()
    lag_section = """
`P2569/S1519` blocks promotion of finite APD bridge values into role-bearing A/P/D dynamics in `L_total`.  A strict source must provide a dynamical equation, regularity principle, or minimal-action law that rejects the vanishing-polynomial family; otherwise APD remains value-level bridge evidence only.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2569/S1519 APD value-bridge interpolation dynamic nonuniqueness guard", "## P2569/S1519 APD value-bridge interpolation dynamic nonuniqueness guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2569/S1519 APD value-bridge interpolation dynamic nonuniqueness Ltotal guard", "## P2569/S1519 APD value-bridge interpolation dynamic nonuniqueness Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
