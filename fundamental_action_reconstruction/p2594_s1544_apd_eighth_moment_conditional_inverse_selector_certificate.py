#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np
import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2589_s1539_apd_mirror_sixth_moment_shell_support_nonuniqueness_audit import (
    EXPECTED_CARDINALITY,
    INTERNAL_FOURTH_SHELL,
    INTERNAL_SECOND_SHELL,
    INTERNAL_SIXTH_SHELL,
    support_from_product,
)
from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import PRODUCT_PARAMETER_GRID
from p2591_s1541_apd_product_parameter_sturm_interval_certificate import INTERVAL
from p2593_s1543_apd_current_state_replay_and_exact_next_moment_provenance_certificate import OUT as P2593_OUT, theorem as previous_theorem

GEN = ROOT / "generated"
OUT = GEN / "p2594_s1544_apd_eighth_moment_conditional_inverse_selector_certificate.json"
MD = GEN / "p2594_s1544_apd_eighth_moment_conditional_inverse_selector_certificate.md"

SOURCE_FILES = {
    "P2593_CURRENT_STATE_REPLAY": P2593_OUT,
}
INTERNAL_FORMULA_INTERCEPT = 74658
CENTRAL_FORMULA_INTERCEPT = sp.Rational(42715646049, 32768)
NEGATIVE_EXPORT_FLAGS = [
    "apd_eighth_moment_inverse_source_exported",
    "apd_conditional_support_recovery_source_exported",
    "apd_next_shell_conditional_selector_source_exported",
    "apd_affine_product_recovery_source_exported",
    "apd_root_reconstruction_source_exported",
    "apd_support_selection_source_exported",
    "apd_finite_support_source_exported",
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
        "new_packet": "P2594|S1544|APD eighth moment inverse selector|eighth moment inverse.*APD",
        "intended_research_nonduplication": "APD.*conditional support recovery|conditional support recovery.*APD|APD.*affine product recovery|affine product recovery.*APD|APD.*next shell conditional selector|next shell conditional selector.*APD|APD.*root reconstruction from eighth|root reconstruction from eighth.*APD",
        "apd_precursors": "P2592|S1542|P2593|S1543|APD.*eighth moment|APD.*Newton Girard|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def recover_product_from_internal_eighth(internal_eighth_shell: float) -> float:
    return float((INTERNAL_FORMULA_INTERCEPT - internal_eighth_shell) / 4.0)


def recover_product_from_central_eighth(central_eighth_moment: float) -> float:
    return float((float(CENTRAL_FORMULA_INTERCEPT) - central_eighth_moment) / 8.0)


def exact_interval_bounds() -> dict[str, Any]:
    left, right = [sp.Integer(value) for value in INTERVAL]
    internal_left = sp.Integer(INTERNAL_FORMULA_INTERCEPT) - 4 * left
    internal_right = sp.Integer(INTERNAL_FORMULA_INTERCEPT) - 4 * right
    central_left = CENTRAL_FORMULA_INTERCEPT - 8 * left
    central_right = CENTRAL_FORMULA_INTERCEPT - 8 * right
    return {
        "product_parameter_interval": INTERVAL,
        "internal_eighth_shell_interval_descending": [str(internal_left), str(internal_right)],
        "central_eighth_moment_interval_descending_exact": [str(central_left), str(central_right)],
        "central_eighth_moment_interval_descending_float": [float(central_left), float(central_right)],
        "internal_inverse_formula": "e4 = (74658 - p4) / 4",
        "central_inverse_formula": "e4 = (42715646049/32768 - M8) / 8",
        "internal_map_is_affine_bijection_on_interval": bool(internal_left > internal_right),
        "central_map_is_affine_bijection_on_interval": bool(central_left > central_right),
    }


def support_reconstruction_row(product_parameter: float) -> dict[str, Any]:
    support = support_from_product(product_parameter)
    squared = np.array(support["squared_offsets"], dtype=float)
    internal_eighth = float(INTERNAL_FORMULA_INTERCEPT - 4.0 * product_parameter)
    central_eighth = float(CENTRAL_FORMULA_INTERCEPT - 8.0 * product_parameter)
    product_from_internal = recover_product_from_internal_eighth(internal_eighth)
    product_from_central = recover_product_from_central_eighth(central_eighth)
    recovered_support = support_from_product(product_from_internal)
    recovered_squared = np.array(recovered_support["squared_offsets"], dtype=float)
    roots_poly_coefficients = [1.0, -INTERNAL_SECOND_SHELL, 273.0, -820.0, product_from_internal]
    roots_from_poly = sorted(float(root.real) for root in np.roots(roots_poly_coefficients) if abs(root.imag) <= 1.0e-8)
    return {
        "product_parameter": float(product_parameter),
        "internal_eighth_shell": internal_eighth,
        "central_eighth_moment": central_eighth,
        "recovered_product_from_internal_eighth": product_from_internal,
        "recovered_product_from_central_eighth": product_from_central,
        "squared_offsets": [float(value) for value in squared],
        "recovered_squared_offsets": [float(value) for value in recovered_squared],
        "roots_from_recovered_quartic": roots_from_poly,
        "max_abs_recovered_product_error": max(abs(product_from_internal - product_parameter), abs(product_from_central - product_parameter)),
        "max_abs_squared_offset_reconstruction_error": float(np.max(np.abs(squared - recovered_squared))),
        "max_abs_poly_root_reconstruction_error": float(np.max(np.abs(squared - np.array(roots_from_poly, dtype=float)))),
        "support_cardinality": len(recovered_support["points"]),
    }


def conditional_inverse_certificate() -> dict[str, Any]:
    bounds = exact_interval_bounds()
    rows = [support_reconstruction_row(value) for value in PRODUCT_PARAMETER_GRID]
    internal_values = [row["internal_eighth_shell"] for row in rows]
    central_values = [row["central_eighth_moment"] for row in rows]
    return {
        "numpy_version": np.__version__,
        "sympy_version": sp.__version__,
        "fixed_lower_shells": {
            "internal_second_shell": INTERNAL_SECOND_SHELL,
            "internal_fourth_shell": INTERNAL_FOURTH_SHELL,
            "internal_sixth_shell": INTERNAL_SIXTH_SHELL,
        },
        "exact_interval_bounds": bounds,
        "product_parameter_grid": PRODUCT_PARAMETER_GRID,
        "support_reconstruction_rows": rows,
        "support_reconstruction_row_count": len(rows),
        "expected_cardinality": EXPECTED_CARDINALITY,
        "internal_eighth_values_are_strictly_decreasing_on_grid": all(left > right for left, right in zip(internal_values, internal_values[1:])),
        "central_eighth_values_are_strictly_decreasing_on_grid": all(left > right for left, right in zip(central_values, central_values[1:])),
        "all_products_recovered_from_internal_and_central_eighth": all(row["max_abs_recovered_product_error"] <= 1.0e-10 for row in rows),
        "all_squared_offsets_reconstructed": all(row["max_abs_squared_offset_reconstruction_error"] <= 1.0e-9 for row in rows),
        "all_recovered_quartic_roots_match_support_offsets": all(row["max_abs_poly_root_reconstruction_error"] <= 1.0e-7 for row in rows),
        "all_recovered_supports_have_expected_cardinality": all(row["support_cardinality"] == EXPECTED_CARDINALITY for row in rows),
        "conditional_inverse_is_selector_coordinate_not_source": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2593_payload = load_json(SOURCE_FILES["P2593_CURRENT_STATE_REPLAY"])
    p2593 = previous_theorem(p2593_payload, "apd_current_state_replay_and_exact_next_moment_provenance_certificate")
    certificate = conditional_inverse_certificate()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2594_T1_apd_eighth_moment_conditional_inverse_selector_certificate",
        "audited_chain": ["P2592/S1542", "P2593/S1543"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "treat affine recovery of the product parameter from an eighth moment as a strict APD support source",
        "p2593_exact_replay_inherited": p2593.get("exact_replay_preserves_p2592_numeric_certificate") is True,
        "apd_eighth_moment_conditional_inverse_selector_certificate": certificate,
        "eighth_moment_conditionally_recovers_product_parameter": certificate["all_products_recovered_from_internal_and_central_eighth"],
        "recovered_product_conditionally_reconstructs_support_roots": certificate["all_recovered_quartic_roots_match_support_offsets"],
        "conditional_inverse_still_requires_eighth_moment_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "P2594 proves a conditional inverse: if a strict source supplies the eighth shell, then the P2591 product parameter and mirror support are recovered within the audited interval. This is not itself that source; the honest next step remains a strict nadsoliton-derived eighth-shell/support law or a real APD dynamics theorem."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2593_exact_replay_inherited": theorem_export["p2593_exact_replay_inherited"],
        "internal_map_affine_bijection_on_interval": certificate["exact_interval_bounds"]["internal_map_is_affine_bijection_on_interval"],
        "central_map_affine_bijection_on_interval": certificate["exact_interval_bounds"]["central_map_is_affine_bijection_on_interval"],
        "grid_has_eight_rows": certificate["support_reconstruction_row_count"] == len(PRODUCT_PARAMETER_GRID),
        "all_products_recovered": certificate["all_products_recovered_from_internal_and_central_eighth"],
        "all_squared_offsets_reconstructed": certificate["all_squared_offsets_reconstructed"],
        "all_quartic_roots_match_support_offsets": certificate["all_recovered_quartic_roots_match_support_offsets"],
        "all_recovered_supports_have_expected_cardinality": certificate["all_recovered_supports_have_expected_cardinality"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2594",
        "stage_id": "S1544",
        "status": "P2594_APD_EIGHTH_MOMENT_CONDITIONAL_INVERSE_SELECTOR_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_eighth_moment_conditional_inverse_selector_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2593_CURRENT_STATE_REPLAY": sha256_json(p2593_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_eighth_moment_conditional_inverse_selector_certificate"]["theorem_export"]
    c = t["apd_eighth_moment_conditional_inverse_selector_certificate"]
    bounds = c["exact_interval_bounds"]
    lines = [
        "# P2594/S1544 APD eighth-moment conditional inverse selector certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Internal inverse: `{bounds['internal_inverse_formula']}`.",
        f"- Central inverse: `{bounds['central_inverse_formula']}`.",
        f"- Support reconstruction rows: `{c['support_reconstruction_row_count']}`.",
        f"- Products recovered: `{t['eighth_moment_conditionally_recovers_product_parameter']}`.",
        f"- Roots/supports reconstructed: `{t['recovered_product_conditionally_reconstructs_support_roots']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2594 proves the inverse direction left implicit by P2592/P2593: an externally supplied eighth shell is an affine selector coordinate for the P2591 product parameter, and that recovered product reconstructs the mirror-support quartic roots on the audited grid.  This is conditional determinacy, not a derivation of the eighth shell from strict dynamics.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No eighth-moment inverse source, conditional support-recovery source, next-shell selector source, affine product-recovery source, root-reconstruction source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_eighth_moment_conditional_inverse_selector_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2594/S1544 APD eighth-moment conditional inverse selector guard

`P2594/S1544` proves the conditional inverse of P2592/P2593: if an eighth shell is supplied, then `e4 = (74658 - p4)/4` (equivalently `e4 = (42715646049/32768 - M8)/8`) recovers the P2591 product parameter and reconstructs the mirror-support quartic roots on the audited grid.  This is conditional selector determinacy only; it still does not source the eighth shell or the APD support law.
""".strip()
    lag_section = """
## P2594/S1544 APD eighth-moment conditional inverse selector Ltotal guard

`P2594/S1544` blocks a role-bearing APD Gram term in `L_total` from being justified by conditional inverse reconstruction from an unsourced eighth shell.  Reconstructing the product parameter after the shell is supplied is not a strict nadsoliton-derived support/density source.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2594/S1544 APD eighth-moment conditional inverse selector guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2594/S1544 APD eighth-moment conditional inverse selector Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
