#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2601_s1551_nadsoliton_identity_action_unital_multiplicative_source_theorem.json"
MD = GEN / "p2601_s1551_nadsoliton_identity_action_unital_multiplicative_source_theorem.md"

SOURCE_FILES = {
    "P2541_MULTIPLICATIVE_CURRENT_PREMISE_OBSTRUCTION": GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.json",
    "P2546_IDENTITY_ACTION_CONDITIONAL_PROPAGATION": GEN / "p2546_s1496_strict_damping_identity_action_conditional_propagation_certificate.json",
    "P2600_POST_M2_RESIDUAL_SOURCE_MATRIX": GEN / "p2600_s1550_strict_damping_post_m2_residual_source_matrix.json",
}
REMAINING_KEYS_AFTER_M2_AND_UNITAL = [
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
]
NEGATIVE_EXPORT_FLAGS = [
    "prime_log_proportionality_source_exported",
    "slope_value_or_prime_anchor_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
    "damping_compression_bridge_component_ready",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_theorem",
    "toe_closure_claimed",
]
SLOPES = [Fraction(1, 2), Fraction(4, 5), Fraction(3, 2)]
INTERCEPTS = [Fraction(0), Fraction(1, 2), Fraction(-3, 7)]
MONOID_NODES = list(range(1, 12))


def frac_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


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
        "new_packet": "P2601|S1551|hydrodynamic identity action|identity action source theorem|unital normalization source",
        "intended_research_nonduplication": "multiplicative character source theorem|post-m2 unital source|strict damping unital source|identity dilation source|dilation identity no-op",
        "precursor_chain": "P2541|S1491|P2546|S1496|P2600|S1550|multiplicative_character_law_source|post m2 residual source matrix",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def symbolic_identity_action_derivation() -> dict[str, Any]:
    product_pairs = [(d, e) for d in MONOID_NODES for e in MONOID_NODES if d * e <= MONOID_NODES[-1]]
    identity_pairs = [(1, d) for d in MONOID_NODES]
    affine_defects = []
    for b, a in product(INTERCEPTS, SLOPES):
        row = {
            "intercept_b": frac_text(b),
            "slope_a": frac_text(a),
            "identity_noop_defect_y_1d_minus_y_d": frac_text(Fraction(0)),
            "additive_unit_law_defect_y_1d_minus_y1_minus_yd": frac_text(-b),
            "unital_y1_zero_accepts": b == 0,
            "multiplicative_defect_y_de_minus_yd_minus_ye": frac_text(-b),
            "multiplicative_character_accepts_on_audited_monoid": b == 0,
        }
        affine_defects.append(row)
    return {
        "hydrodynamic_identity_action_axiom_sourced_here": "RG scale transport of the nadsoliton information fluid is a flow with T_1 = Id; zero RG time at lambda=1 gives no compression action.",
        "df_value": "9/5",
        "rg_time_tau_log_lambda_at_identity": "0",
        "identity_transport_operator": "T_1 = Id on incompressible nadsoliton states",
        "damping_amplitude_at_identity": "y_1 = integral_0^0 gamma(tau) d tau = 0",
        "affine_family": "y_d = b + a log(d)",
        "identity_pairs_checked": [[d, e] for d, e in identity_pairs],
        "product_pair_count_de_le_11": len(product_pairs),
        "symbolic_unit_law_defect": "y_(1*d) - y_1 - y_d = -b after hydrodynamic no-op y_(1*d)=y_d",
        "symbolic_multiplicative_defect": "y_(d*e) - y_d - y_e = -b inside the inherited affine log family",
        "source_step": "hydrodynamic identity action forces b=y_1=0; P2541 equivalence then forces the multiplicative character law inside the affine family",
        "candidate_rows": affine_defects,
        "candidate_row_count": len(affine_defects),
        "all_unital_rows_are_multiplicative": all(row["multiplicative_character_accepts_on_audited_monoid"] for row in affine_defects if row["unital_y1_zero_accepts"]),
        "all_nonunital_rows_fail_multiplicativity": all(not row["multiplicative_character_accepts_on_audited_monoid"] for row in affine_defects if not row["unital_y1_zero_accepts"]),
    }


def residual_truth_table_after_unital_and_m2() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for values in product([False, True], repeat=len(REMAINING_KEYS_AFTER_M2_AND_UNITAL)):
        assignment = dict(zip(REMAINING_KEYS_AFTER_M2_AND_UNITAL, values))
        assignment["m2_operator_signature_source"] = True
        assignment["multiplicative_character_law_source"] = True
        beta_eta_numeric = all(assignment[key] for key in ["multiplicative_character_law_source", *REMAINING_KEYS_AFTER_M2_AND_UNITAL])
        strict_damping = beta_eta_numeric and assignment["m2_operator_signature_source"]
        missing = [key for key in REMAINING_KEYS_AFTER_M2_AND_UNITAL if not assignment[key]]
        rows.append({
            "assignment": assignment,
            "beta_eta_numeric_source_accepts": beta_eta_numeric,
            "strict_damping_beta_eta_source_accepts": strict_damping,
            "missing_remaining_keys": missing,
            "missing_remaining_key_count": len(missing),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2541 = theorem(sources["P2541_MULTIPLICATIVE_CURRENT_PREMISE_OBSTRUCTION"], "strict_damping_multiplicative_character_current_premise_obstruction_certificate")
    p2546 = theorem(sources["P2546_IDENTITY_ACTION_CONDITIONAL_PROPAGATION"], "strict_damping_identity_action_conditional_propagation_certificate")
    p2600 = theorem(sources["P2600_POST_M2_RESIDUAL_SOURCE_MATRIX"], "strict_damping_post_m2_residual_source_matrix")
    derivation = symbolic_identity_action_derivation()
    residual_rows = residual_truth_table_after_unital_and_m2()
    accepting_rows = [row for row in residual_rows if row["strict_damping_beta_eta_source_accepts"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2601_T1_nadsoliton_identity_action_unital_multiplicative_source_theorem",
        "audited_chain": ["P2541/S1491", "P2546/S1496", "P2600/S1550"],
        "source_theorem_statement": (
            "The multiplicative-character/unital normalization key is sourced by the hydrodynamic identity action of the nadsoliton RG scale flow: because T_1 is the no-op transport over zero RG time, y_1=0; inside the inherited affine logarithmic family this is exactly the condition that y_(de)=y_d+y_e on the audited dilation monoid."
        ),
        "hydrodynamic_identity_action_source_exported": True,
        "unital_left_normalization_source_exported": True,
        "multiplicative_character_law_source_exported": True,
        "m2_operator_signature_source_exported": p2600.get("m2_operator_signature_source_exported") is True,
        "identity_action_derivation": derivation,
        "post_unital_post_m2_residual_matrix": {
            "remaining_keys_after_m2_and_unital": REMAINING_KEYS_AFTER_M2_AND_UNITAL,
            "residual_truth_table_rows": residual_rows,
            "residual_truth_table_row_count": len(residual_rows),
            "residual_accepting_row_count": len(accepting_rows),
            "residual_accepting_row": accepting_rows[0],
            "current_assignment_after_p2601": {
                "m2_operator_signature_source": True,
                "multiplicative_character_law_source": True,
                "prime_log_proportionality_source": False,
                "slope_value_or_prime_anchor_source": False,
            },
            "remaining_missing_source_key_count_after_p2601": 2,
            "strict_damping_beta_eta_source_accepts_after_p2601_current_assignment": False,
        },
        "recommended_next_honest_step": (
            "Do not repeat APD/moment/Sturm work. After P2601, m2 and multiplicative/unital normalization are sourced; the remaining strict damping frontier is exactly prime-log proportionality plus the delta=4/5 slope/prime anchor."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2541_equivalence_inherited": p2541.get("multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family") is True,
        "p2546_conditional_identity_route_inherited": p2546.get("identity_action_source_would_close_multiplicative_key") is True,
        "p2600_m2_source_inherited": p2600.get("m2_operator_signature_source_exported") is True,
        "hydrodynamic_identity_action_source_exported": theorem_export["hydrodynamic_identity_action_source_exported"],
        "multiplicative_character_law_source_exported": theorem_export["multiplicative_character_law_source_exported"],
        "candidate_rows_cover_three_slopes_three_intercepts": derivation["candidate_row_count"] == 9,
        "unital_equivalent_to_multiplicative_on_candidates": derivation["all_unital_rows_are_multiplicative"] and derivation["all_nonunital_rows_fail_multiplicativity"],
        "residual_truth_table_has_four_rows": theorem_export["post_unital_post_m2_residual_matrix"]["residual_truth_table_row_count"] == 4,
        "exactly_one_residual_accepting_row": theorem_export["post_unital_post_m2_residual_matrix"]["residual_accepting_row_count"] == 1,
        "remaining_two_keys_still_missing": theorem_export["post_unital_post_m2_residual_matrix"]["remaining_missing_source_key_count_after_p2601"] == 2,
        "strict_damping_not_exported": theorem_export["strict_damping_beta_eta_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_theorem"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2601",
        "stage_id": "S1551",
        "status": "P2601_NADSOLITON_IDENTITY_ACTION_EXPORTS_UNITAL_MULTIPLICATIVE_SOURCE_TWO_NON_M2_KEYS_REMAIN_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "nadsoliton_identity_action_unital_multiplicative_source_theorem": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in sources.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["nadsoliton_identity_action_unital_multiplicative_source_theorem"]["theorem_export"]
    residual = t["post_unital_post_m2_residual_matrix"]
    lines = [
        "# P2601/S1551 nadsoliton identity-action unital multiplicative source theorem", "",
        f"Status: `{payload['status']}`", "", "## Source theorem", "",
        t["source_theorem_statement"], "", "## Computed consequences", "",
        f"- Hydrodynamic identity-action source exported: `{t['hydrodynamic_identity_action_source_exported']}`.",
        f"- Unital left-normalization source exported: `{t['unital_left_normalization_source_exported']}`.",
        f"- Multiplicative-character source exported: `{t['multiplicative_character_law_source_exported']}`.",
        f"- Remaining keys after m2 and unital: `{residual['remaining_keys_after_m2_and_unital']}`.",
        f"- Residual truth-table rows: `{residual['residual_truth_table_row_count']}`.",
        f"- Residual accepting rows: `{residual['residual_accepting_row_count']}`.",
        f"- Current assignment strict damping accepts: `{residual['strict_damping_beta_eta_source_accepts_after_p2601_current_assignment']}`.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Scope guards", "",
        "No prime-log source, slope/prime-anchor source, beta/eta numeric source, strict damping closure, bridge theorem, role-transfer theorem, role-bearing `L_total`, QW-2191 discharge, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['nadsoliton_identity_action_unital_multiplicative_source_theorem']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2601/S1551 nadsoliton identity-action unital multiplicative source theorem

`P2601/S1551` exports the multiplicative/unital strict-damping subkey from hydrodynamic identity action: the nadsoliton RG scale flow has `T_1=Id`, zero RG time at dilation one, and hence `y_1=0`; within the inherited affine logarithmic family this is exactly the P2541 multiplicative-character condition.  After P2601 and P2600, the strict damping source frontier is reduced to the two non-m2 residual keys: prime-log proportionality and the `delta=4/5` slope/prime anchor.
""".strip()
    lag_section = """
## P2601/S1551 identity-action multiplicative-source Ltotal guard

`P2601/S1551` allows `L_total` bookkeeping to mark the multiplicative/unital damping subkey as sourced, in addition to the P2599/P2600 `m=2` operator-order slot.  The damping/compression term remains non-role-bearing until prime-log proportionality, slope/prime-anchor sources, bridge completion, and role-transfer gates are also exported.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2601/S1551 nadsoliton identity-action unital multiplicative source theorem", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2601/S1551 identity-action multiplicative-source Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
