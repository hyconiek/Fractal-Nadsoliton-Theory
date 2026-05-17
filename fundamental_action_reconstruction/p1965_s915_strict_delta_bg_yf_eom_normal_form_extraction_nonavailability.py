#!/usr/bin/env python3
"""P1965 S915 strict DELTA_BG_Yf EOM-normal-form nonavailability theorem.

P1964 proved the conditional algebraic PO2 subtheorem on the declared
tensorial B1 normal form.  This packet checks whether the current repository
exports enough explicit variational data to derive that normal form from the
full non-skeleton L_total.  It does not.  The result is a narrow current-state
nonavailability theorem, not a no-go theorem about the theory.
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
OUT = GEN / "p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def contains_text(obj: Any, needle: str) -> bool:
    return needle.lower() in json.dumps(obj, sort_keys=True, ensure_ascii=True).lower()


def requirement(
    req_id: str,
    required_object: str,
    source_checks: list[str],
    verdict: str,
    missing_or_insufficient: list[str],
) -> dict[str, Any]:
    return {
        "req_id": req_id,
        "required_object": required_object,
        "source_checks": source_checks,
        "verdict": verdict,
        "missing_or_insufficient": missing_or_insufficient,
    }


def eom_status_summary(p1907: dict[str, Any]) -> dict[str, Any]:
    rows = p1907.get("eom_export_matrix") or []
    statuses = {row.get("field", "UNKNOWN"): row.get("status", "MISSING") for row in rows}
    return {
        "row_count": len(rows),
        "statuses_by_field": statuses,
        "all_open_symbolic": bool(rows) and all(status == "OPEN_SYMBOLIC" for status in statuses.values()),
        "any_explicit_variational_pass": any("PASS" in str(status) for status in statuses.values()),
    }


def underdetermination_witness() -> dict[str, Any]:
    dR, dU, dX = sp.symbols("delta_R delta_RicUU delta_gradchi2")
    C_R, C_U, C_chi = sp.symbols("C_R C_U C_chi")
    omega_missing = sp.symbols("Omega_unexported")
    C_missing = sp.symbols("C_missing")

    declared = C_R * dR + C_U * dU + C_chi * dX
    augmented = declared + C_missing * omega_missing
    constraints = {dR: sp.Integer(0), dU: sp.Integer(0), dX: sp.Integer(0)}

    declared_reduced = sp.simplify(declared.subs(constraints))
    augmented_reduced = sp.simplify(augmented.subs(constraints))

    return {
        "declared_normal_form": sp.sstr(declared),
        "augmented_currently_unexcluded_form": sp.sstr(augmented),
        "c1_c3_substitution": {
            "delta_R": "0",
            "delta_RicUU": "0",
            "delta_gradchi2": "0",
        },
        "declared_after_c1_c3": sp.sstr(declared_reduced),
        "augmented_after_c1_c3": sp.sstr(augmented_reduced),
        "declared_is_zero": bool(declared_reduced == 0),
        "augmented_is_not_identically_zero_without_extra_export": bool(augmented_reduced != 0),
        "meaning": (
            "Current exports prove the declared tensorial algebraic form, but they do "
            "not export a full EOM-derived normal-form theorem excluding an additional "
            "unexported contribution Omega_unexported."
        ),
    }


def main() -> None:
    p1662 = load("p1662_s612_strict_full_lagrangian_explicit_density_summary.json")
    p1708 = load("p1708_s658_strict_nonproxy_covariant_eom_first_explicit_formula_export_checkpoint.json")
    p1848 = load("p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json")
    p1880 = load("p1880_s830_strict_full_lagrangian_term_registry_to_eom_transport_probe.json")
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1930 = load("p1930_s880_strict_b1_invariant_triplet_branch_evaluation_probe.json")
    p1964 = load("p1964_s914_strict_po2_conditional_sufficiency_and_eom_gap_certificate.json")

    density_present = "full_lagrangian_density_explicit" in p1662
    p1708_template = contains_text(p1708, "template-level") or contains_text(p1708, "template")
    p1848_componentwise_gravity = "gravity_componentwise_variation_pack" in p1848
    p1880_eom_registry = "covariant_eom_registry" in p1880
    p1907_eom = eom_status_summary(p1907)
    p1930_tensorial_form = bool(p1930.get("tensorial_b1_witness_form"))
    p1964_conditional_pass = (
        p1964.get("symbolic_po2_sufficiency_check", {}).get("sympy_zero_check") is True
    )

    requirements = [
        requirement(
            "E1",
            "Full non-skeleton L_total density with all relevant sectors",
            [f"P1662 density_present={density_present}", f"P1907 registry_present={'full_lagrangian_term_registry_non_skeleton' in p1907}"],
            "PRESENT_SECTOR_LEVEL",
            [
                "Density is sector-level and formula-level; it is not yet tied to a complete variational derivation for DELTA_BG_Yf.",
            ],
        ),
        requirement(
            "E2",
            "Convention-fixed covariant field calculus for H, A_mu, psi, and g_mu_nu",
            [
                f"P1708 template_level={p1708_template}",
                f"P1880 eom_registry_present={p1880_eom_registry}",
            ],
            "TEMPLATE_ONLY",
            [
                "Need explicit covariant derivative, tetrad/spin connection, metric-density variation, and boundary-term convention.",
            ],
        ),
        requirement(
            "E3",
            "Termwise Euler-Lagrange derivatives from L_total",
            [
                f"P1907 eom_statuses={p1907_eom['statuses_by_field']}",
                f"P1907 all_open_symbolic={p1907_eom['all_open_symbolic']}",
            ],
            "MISSING_EXPLICIT_VARIATIONAL_ROWS",
            [
                "H, A_mu, psi_f, and g_mu_nu rows are OPEN_SYMBOLIC in P1907.",
                "P1693/P1866 provide reduced or 1D proxy EOMs, not full 4D covariant termwise derivatives.",
            ],
        ),
        requirement(
            "E4",
            "Gravity curvature variation pack linked to the background-Yukawa normal form",
            [
                f"P1848 gravity_componentwise_variation_pack={p1848_componentwise_gravity}",
                f"P1930 tensorial_form_present={p1930_tensorial_form}",
            ],
            "PARTIAL_UNLINKED",
            [
                "P1848 exports gravity variation structures, but not a derivation of the P1930 DELTA_BG_Yf basis.",
                "Need projection from metric/Higgs/Yukawa EOMs to {delta_R, delta_RicUU, delta_gradchi2}.",
            ],
        ),
        requirement(
            "E5",
            "Normal-form extraction theorem for DELTA_BG_Yf from the EOM bundle",
            [
                f"P1930 declared_tensorial_form_present={p1930_tensorial_form}",
                f"P1964 conditional_algebra_pass={p1964_conditional_pass}",
            ],
            "DECLARED_NOT_DERIVED",
            [
                "The declared normal form is algebraically certified by P1964.",
                "No current artifact derives that normal form from the full EOM bundle.",
            ],
        ),
        requirement(
            "E6",
            "Same-branch quantifier transport from PO3 formal nonempty domain to EOM-derived branch class",
            [
                f"P1964 current_status={p1964.get('status', 'MISSING')}",
            ],
            "MISSING_TRANSPORT_THEOREM",
            [
                "P1963/P1964 provide formal-domain and conditional algebraic data.",
                "Need proof that the nonempty formal branch lies in the full EOM-derived branch class.",
            ],
        ),
    ]

    failed_or_partial = [
        row["req_id"]
        for row in requirements
        if row["verdict"] not in {"PASS", "PRESENT_COMPLETE"}
    ]

    witness = underdetermination_witness()
    nonavailability_pass = bool(
        density_present
        and p1930_tensorial_form
        and p1964_conditional_pass
        and failed_or_partial
        and witness["augmented_is_not_identically_zero_without_extra_export"]
    )

    out = {
        "packet_id": "P1965",
        "stage_id": "S915",
        "status": (
            "PO2_FULL_EOM_NORMAL_FORM_EXTRACTION_NONAVAILABLE_CURRENT_EXPORT_STATE"
            if nonavailability_pass
            else "PO2_EOM_NORMAL_FORM_EXTRACTION_AUDIT_INCONCLUSIVE"
        ),
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1662_present": density_present,
            "p1708_present": "first_explicit_nonproxy_formula_pack" in p1708,
            "p1848_present": p1848_componentwise_gravity,
            "p1880_present": p1880_eom_registry,
            "p1907_present": "eom_export_matrix" in p1907,
            "p1930_present": p1930_tensorial_form,
            "p1964_present": p1964_conditional_pass,
        },
        "input_hashes": {
            "p1662_sha256": digest(p1662),
            "p1708_sha256": digest(p1708),
            "p1848_sha256": digest(p1848),
            "p1880_sha256": digest(p1880),
            "p1907_sha256": digest(p1907),
            "p1930_sha256": digest(p1930),
            "p1964_sha256": digest(p1964),
        },
        "acceptance_requirements_for_full_po2_eom_derivation": requirements,
        "failed_or_partial_requirement_ids": failed_or_partial,
        "symbolic_underdetermination_witness": witness,
        "formal_nonavailability_theorem": {
            "statement": (
                "On the current repository state, the declared tensorial DELTA_BG_Yf "
                "normal form is not derivable as a full EOM theorem from L_total."
            ),
            "proof_trace": [
                "P1662/P1907 provide a sector-level non-skeleton L_total anchor.",
                "P1708/P1880 provide template or registry EOM forms, not a complete variational derivation.",
                "P1907 marks H, A_mu, psi_f, and g_mu_nu EOM rows OPEN_SYMBOLIC.",
                "P1848 exports gravity variation structures but not the projection theorem to DELTA_BG_Yf.",
                "P1964 proves only conditional algebra on the declared normal form.",
                "The symbolic Omega_unexported witness shows an extra unexported EOM contribution is not excluded by current artifacts.",
            ],
            "not_a_no_go_theorem": True,
        },
        "po2_restamp_after_p1965": {
            "before": p1964.get("po2_restamp", {}).get("after", p1964.get("status", "UNKNOWN")),
            "after": (
                "OPEN_FULL_EOM_NORMAL_FORM_EXTRACTION_NONAVAILABLE_CURRENT_EXPORT_STATE"
                if nonavailability_pass
                else "OPEN_EOM_NORMAL_FORM_AUDIT_INCONCLUSIVE"
            ),
            "conditional_algebra_from_p1964_retained": p1964_conditional_pass,
            "full_po2_sufficiency_closed": False,
        },
        "strict_core_statusvector_restamp": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": "OPEN_PO2_EOM_NORMAL_FORM_EXPORT_REQUIRED",
            "selector_qw2191": "OPEN",
        },
        "toe_remaining_minimum_after_p1965": {
            "previous_open_blocks_after_p1964": p1964.get("toe_remaining_minimum_after_p1964", {}).get(
                "current_minimum_open_blocks", 6
            ),
            "current_minimum_open_blocks": 6,
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
            "P1965 is a current-export nonavailability theorem. It does not refute "
            "PO2 in principle and does not weaken the P1964 conditional algebraic pass."
        ),
        "higher_reasoning_required_for_next_step": True,
        "next_honest_step": (
            "Build P1966: a narrow constructive variational extraction for the Higgs/Yukawa/"
            "curvature subsector sqrt(-g)[(D H)^dagger(D H)-mu_H^2 H^dagger H-"
            "lambda_H(H^dagger H)^2-xi_H H^dagger H R-y_f qbar H f], with frozen "
            "metric-density and integration-by-parts conventions, then test whether "
            "Omega_unexported can be forced to zero."
        ),
        "lay_explanation": (
            "The repo has the recipe and the algebraic lock, but it still lacks the full "
            "mechanical derivation showing that the recipe produces exactly that lock. "
            "There may be an uncomputed leftover term until the variational calculation is "
            "done explicitly."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
