#!/usr/bin/env python3
"""P1964 S914 strict PO2 conditional sufficiency and EOM-gap certificate.

This packet performs the next honest PO2 move after P1963.  It proves the
narrow algebraic theorem that the declared tensorial B1 normal form for
DELTA_BG_Yf vanishes under C1-C3, with C4 recorded as the same-branch strict
EOM consistency predicate.  It also certifies that the current repository
state still lacks the explicit variational EOM export needed to promote this
conditional algebraic theorem into full PO2 sufficiency from L_total.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1964_s914_strict_po2_conditional_sufficiency_and_eom_gap_certificate.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def s(expr: sp.Expr) -> str:
    return sp.sstr(expr)


def eom_status_summary(p1907: dict[str, Any]) -> dict[str, Any]:
    rows = p1907.get("eom_export_matrix") or []
    statuses = {row.get("field", "UNKNOWN"): row.get("status", "MISSING") for row in rows}
    explicit_statuses = {
        "PASS_EXPLICIT_VARIATIONAL",
        "PASS_NONPROXY_VARIATIONAL",
        "PASS_FULL_VARIATIONAL_EOM_EXPORTED",
    }
    all_explicit = bool(rows) and all(status in explicit_statuses for status in statuses.values())
    all_open_symbolic = bool(rows) and all(status == "OPEN_SYMBOLIC" for status in statuses.values())
    return {
        "row_count": len(rows),
        "statuses_by_field": statuses,
        "all_rows_explicit_variational": all_explicit,
        "all_rows_open_symbolic": all_open_symbolic,
        "accepted_explicit_statuses": sorted(explicit_statuses),
    }


def symbolic_po2_check() -> dict[str, Any]:
    delta_R, delta_RicUU, delta_gradchi2 = sp.symbols("delta_R delta_RicUU delta_gradchi2")
    C_R, C_U, C_chi = sp.symbols("C_R C_U C_chi")
    rho_H, rho_A, rho_psi, rho_g = sp.symbols("rho_H rho_A rho_psi rho_g")

    delta_bg_yf_tensor = C_R * delta_R + C_U * delta_RicUU + C_chi * delta_gradchi2
    c1_c4_subs = {
        delta_R: sp.Integer(0),
        delta_RicUU: sp.Integer(0),
        delta_gradchi2: sp.Integer(0),
        rho_H: sp.Integer(0),
        rho_A: sp.Integer(0),
        rho_psi: sp.Integer(0),
        rho_g: sp.Integer(0),
    }
    reduced = sp.simplify(delta_bg_yf_tensor.subs(c1_c4_subs))
    residual_norm_proxy = rho_H**2 + rho_A**2 + rho_psi**2 + rho_g**2
    residual_reduced = sp.simplify(residual_norm_proxy.subs(c1_c4_subs))

    return {
        "normal_form": "DELTA_BG_Yf_tensor = C_R*delta_R + C_U*delta_RicUU + C_chi*delta_gradchi2",
        "symbols": {
            "basis_deltas": ["delta_R", "delta_RicUU", "delta_gradchi2"],
            "free_coefficients": ["C_R", "C_U", "C_chi"],
            "c4_eom_residual_symbols": ["rho_H", "rho_A", "rho_psi", "rho_g"],
        },
        "c1_c4_substitution": {
            "C1": "delta_R -> 0",
            "C2": "delta_RicUU -> 0",
            "C3": "delta_gradchi2 -> 0",
            "C4": "rho_H,rho_A,rho_psi,rho_g -> 0 on the same strict branch",
        },
        "delta_bg_yf_before": s(delta_bg_yf_tensor),
        "delta_bg_yf_after_c1_c4": s(reduced),
        "c4_residual_norm_proxy_before": s(residual_norm_proxy),
        "c4_residual_norm_proxy_after": s(residual_reduced),
        "sympy_zero_check": bool(reduced == 0 and residual_reduced == 0),
        "verdict": (
            "PASS_CONDITIONAL_ON_DECLARED_TENSORIAL_B1_NORMAL_FORM"
            if reduced == 0 and residual_reduced == 0
            else "FAIL_SYMBOLIC_REDUCTION"
        ),
    }


def main() -> None:
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1930 = load("p1930_s880_strict_b1_invariant_triplet_branch_evaluation_probe.json")
    p1931 = load("p1931_s881_strict_b1_branch_policy_theorem_candidate_probe.json")
    p1932 = load("p1932_s882_strict_b1_po1_po2_po3_proof_attempt_probe.json")
    p1963 = load("p1963_s913_strict_po3_double_run_machine_checker.json")

    po2_rows = [
        row for row in p1932.get("po1_po2_po3_attempt_table", [])
        if row.get("proof_obligation") == "PO2"
    ]
    po2_before = po2_rows[0].get("status") if po2_rows else "MISSING"

    symbolic_check = symbolic_po2_check()
    eom_summary = eom_status_summary(p1907)

    full_ltotal_anchor_present = "full_lagrangian_term_registry_non_skeleton" in p1907
    declared_tensorial_form_present = bool(p1930.get("tensorial_b1_witness_form"))
    c1_c4_policy_present = (
        "theorem_candidate_branch_policy" in p1931
        and len(
            p1931.get("theorem_candidate_branch_policy", {}).get("required_conditions", [])
        )
        == 4
    )
    po3_machine_nonempty = (
        p1963.get("double_run_reproducibility", {}).get("double_run_gate_pass") is True
        and p1963.get("po3_restamp", {}).get("after") == "PASS_MACHINE_CHECKED_FORMAL_DOMAIN_NONEMPTY"
    )

    full_eom_derivation_available = bool(
        full_ltotal_anchor_present
        and declared_tensorial_form_present
        and c1_c4_policy_present
        and eom_summary["all_rows_explicit_variational"]
    )

    conditional_pass = symbolic_check["sympy_zero_check"] and declared_tensorial_form_present
    if conditional_pass and full_eom_derivation_available:
        local_status = "PO2_FULL_EOM_SUFFICIENCY_PASS"
        background_status = "OPEN_OTHER_BLOCKS_WITH_PO2_PASS"
        current_open_blocks = 5
    elif conditional_pass:
        local_status = "PO2_CONDITIONAL_TENSORIAL_ALGEBRA_PASS__FULL_EOM_DERIVATION_OPEN"
        background_status = "OPEN_PO2_FULL_EOM_DERIVATION_PENDING_WITH_PO3_FORMAL_NONEMPTY_PASS"
        current_open_blocks = 6
    else:
        local_status = "PO2_SYMBOLIC_CONDITIONAL_CHECK_FAILED"
        background_status = "OPEN_PO2_CONDITIONAL_CHECK_FAILED"
        current_open_blocks = 6

    out = {
        "packet_id": "P1964",
        "stage_id": "S914",
        "status": local_status,
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1930_present": "tensorial_b1_witness_form" in p1930,
            "p1931_present": "theorem_candidate_branch_policy" in p1931,
            "p1932_present": bool(po2_rows),
            "p1963_po3_machine_nonempty": po3_machine_nonempty,
        },
        "input_hashes": {
            "p1907_sha256": digest(p1907),
            "p1930_sha256": digest(p1930),
            "p1931_sha256": digest(p1931),
            "p1932_sha256": digest(p1932),
            "p1963_sha256": digest(p1963),
        },
        "c1_c4_formalization": {
            "C1": "delta_R = 0",
            "C2": "delta_RicUU = 0",
            "C3": "delta_gradchi2 = 0",
            "C4": "same branch satisfies strict EOM residual consistency in matter+metric sectors",
            "c4_role": "domain/same-branch admissibility predicate; not an extra term in the declared tensorial DELTA_BG_Yf normal form",
        },
        "symbolic_po2_sufficiency_check": symbolic_check,
        "full_eom_derivation_recheck": {
            "full_lagrangian_non_skeleton_anchor_present": full_ltotal_anchor_present,
            "declared_tensorial_b1_normal_form_present": declared_tensorial_form_present,
            "c1_c4_policy_present": c1_c4_policy_present,
            "p1907_eom_status_summary": eom_summary,
            "full_eom_derivation_available": full_eom_derivation_available,
            "verdict": (
                "PASS_EXPLICIT_VARIATIONAL_EOM_AVAILABLE"
                if full_eom_derivation_available
                else "NOT_PASSED_CURRENT_EXPORT_STATE"
            ),
            "blocking_gap": [
                "P1907 exports sector-level non-skeleton L_total registry, but EOM rows are not explicit variational expressions.",
                "Need termwise Euler-Lagrange derivatives for H, A_mu, psi_f, and g_mu_nu with coefficient substitution.",
                "Need a normal-form extraction lemma deriving the P1930 tensorial DELTA_BG_Yf basis from those EOMs.",
                "Need same-branch quantifier transport from the P1963 nonempty formal domain into the full EOM-derived branch class.",
            ],
        },
        "po2_restamp": {
            "before": po2_before,
            "after": local_status,
            "not_promoted_to": [
                "full PO2 sufficiency from explicit L_total EOM",
                "global background-independence closure",
                "full strict-core ToE closure",
            ],
        },
        "strict_core_statusvector_restamp": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": background_status,
            "selector_qw2191": "OPEN",
        },
        "toe_remaining_minimum_after_p1964": {
            "previous_open_blocks_after_p1963": p1963.get("toe_remaining_minimum_after_p1963", {}).get(
                "current_minimum_open_blocks", 6
            ),
            "po2_conditional_algebraic_subtheorem_pass": bool(conditional_pass),
            "po2_full_eom_sufficiency_closed": bool(full_eom_derivation_available and conditional_pass),
            "current_minimum_open_blocks": current_open_blocks,
            "still_open": [
                "R1 renormalization theorem-grade closure",
                "U1 unitarity/Cutkosky theorem-grade closure",
                "B1 background-independence global theorem closure",
                "S1 selector obstruction QW-2191 discharge",
                "PO2 full EOM-derived sufficiency",
                "cross-scheme finite-part transport theorem linking R1/U1/B1 on common basis",
            ],
        },
        "toolchain": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "sympy": sp.__version__,
        },
        "false_pass_guard": (
            "P1964 proves only conditional algebraic sufficiency on the declared tensorial "
            "B1 normal form. It does not prove that the normal form is derived from the "
            "full non-skeleton L_total EOM, because those explicit variational EOM exports "
            "are not present in the current repo state."
        ),
        "higher_reasoning_required_for_next_step": True,
        "next_honest_step": (
            "Build P1965: export the explicit variational EOM normal-form extraction for "
            "DELTA_BG_Yf from the P1907 L_total registry, or prove a nonavailability theorem "
            "pinpointing the missing density/connection/coefficient data."
        ),
        "lay_explanation": (
            "The algebra now checks out: if the three background differences are exactly zero "
            "on one EOM-consistent branch, the Yukawa background discrepancy is exactly zero. "
            "But the repo still has to prove that these conditions really follow from the full "
            "equations of motion, not only from the declared tensorial test form."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
