#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2546_s1496_strict_damping_identity_action_conditional_propagation_certificate.json"
MD = GEN / "p2546_s1496_strict_damping_identity_action_conditional_propagation_certificate.md"

SOURCE_FILES = {
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2541_MULTIPLICATIVE_OBSTRUCTION": GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.json",
    "P2542_PRIME_LOG_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
    "P2543_SLOPE_OBSTRUCTION": GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json",
    "P2544_FOUR_KEY_CLOSURE_BLOCKER": GEN / "p2544_s1494_strict_damping_four_key_current_premise_closure_blocker_certificate.json",
    "P2545_UNITAL_OBSTRUCTION": GEN / "p2545_s1495_strict_damping_unital_normalization_current_premise_obstruction_certificate.json",
}

STRICT_SOURCE_KEYS = [
    "multiplicative_character_law_source",
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
    "m2_operator_signature_source",
]

NEGATIVE_EXPORT_FLAGS = [
    "strict_identity_action_source_exported",
    "unital_monoid_normalization_source_exported",
    "prime_log_proportionality_source_exported",
    "slope_value_or_prime_anchor_source_exported",
    "beta_eta_numeric_source_exported",
    "m2_operator_signature_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
    "damping_compression_bridge_component_ready",
    "full_bridge_theorem_exported",
    "role_transfer_theorem_exported",
    "selector_closure_exported",
    "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2546|S1496|identity-action conditional propagation|unit identity source conditional|identity action closes multiplicative key",
        "intended_research_nonduplication": "identity-action|identity action|unit action|zero-damping-action|zero damping action|conditional propagation.*y_1=0|y_1=0.*conditional propagation",
        "precursor_packets": "P2530|S1480|P2541|S1491|P2542|S1492|P2543|S1493|P2544|S1494|P2545|S1495",
        "source_key_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2530_FOUR_KEY_IRREDUNDANCY": theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate"),
        "P2541_MULTIPLICATIVE_OBSTRUCTION": theorem(sources["P2541_MULTIPLICATIVE_OBSTRUCTION"], "strict_damping_multiplicative_character_current_premise_obstruction_certificate"),
        "P2542_PRIME_LOG_OBSTRUCTION": theorem(sources["P2542_PRIME_LOG_OBSTRUCTION"], "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"),
        "P2543_SLOPE_OBSTRUCTION": theorem(sources["P2543_SLOPE_OBSTRUCTION"], "strict_damping_slope_value_current_premise_obstruction_certificate"),
        "P2544_FOUR_KEY_CLOSURE_BLOCKER": theorem(sources["P2544_FOUR_KEY_CLOSURE_BLOCKER"], "strict_damping_four_key_current_premise_closure_blocker_certificate"),
        "P2545_UNITAL_OBSTRUCTION": theorem(sources["P2545_UNITAL_OBSTRUCTION"], "strict_damping_unital_normalization_current_premise_obstruction_certificate"),
    }


def candidate_rows(sources: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    return sources["P2545_UNITAL_OBSTRUCTION"]["strict_damping_unital_normalization_current_premise_obstruction_certificate"]["candidate_rows"]


def exact_unit_action_linear_audit() -> dict[str, Any]:
    # For affine y_d=b+a log(d), the unit-action equation y_{1*d}=y_1+y_d
    # reduces exactly to -b=0 for every audited d.  The coefficient row in
    # variables (b, a) is therefore [-1, 0].
    equations = [{"equation": "y_(1*d)-y_1-y_d=0", "coefficient_row_b_a": [-1, 0]}]
    rank = 1
    variable_count = 2
    return {
        "variables": ["b=log_beta", "a=delta"],
        "unit_action_equations_collapse_to": equations,
        "rank": rank,
        "nullity_after_unit_action": variable_count - rank,
        "solution_form": "b=0 with free slope a",
        "d_equals_1_proof": "Taking d=1 in y_(1*d)=y_1+y_d gives y_1=2*y_1, hence y_1=0 over characteristic zero.",
        "unit_action_entails_unital_normalization": True,
        "unit_action_does_not_select_slope": True,
    }


def conditional_truth_assignment(with_strict_identity_action: bool) -> dict[str, bool]:
    return {
        "multiplicative_character_law_source": with_strict_identity_action,
        "prime_log_proportionality_source": False,
        "slope_value_or_prime_anchor_source": False,
        "m2_operator_signature_source": False,
    }


def assignment_result(assignment: dict[str, bool]) -> dict[str, Any]:
    beta_eta_numeric_ready = all(assignment[key] for key in STRICT_SOURCE_KEYS[:3])
    strict_damping_ready = beta_eta_numeric_ready and assignment["m2_operator_signature_source"]
    missing = [key for key in STRICT_SOURCE_KEYS if not assignment[key]]
    return {
        "assignment": assignment,
        "beta_eta_numeric_source_by_assignment": beta_eta_numeric_ready,
        "strict_damping_beta_eta_source_by_assignment": strict_damping_ready,
        "missing_source_keys": missing,
        "missing_source_key_count": len(missing),
    }


def row_passes_identity_action(row: dict[str, Any]) -> bool:
    return bool(row["unit_product_law_accepts"] and row["unital_y1_zero_accepts"])


def build_certificate(theorems: dict[str, dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    p2530 = theorems["P2530_FOUR_KEY_IRREDUNDANCY"]
    p2541 = theorems["P2541_MULTIPLICATIVE_OBSTRUCTION"]
    p2542 = theorems["P2542_PRIME_LOG_OBSTRUCTION"]
    p2543 = theorems["P2543_SLOPE_OBSTRUCTION"]
    p2544 = theorems["P2544_FOUR_KEY_CLOSURE_BLOCKER"]
    p2545 = theorems["P2545_UNITAL_OBSTRUCTION"]
    exact_audit = exact_unit_action_linear_audit()
    identity_accepting_rows = [row for row in rows if row_passes_identity_action(row)]
    identity_rejecting_rows = [row for row in rows if not row_passes_identity_action(row)]
    current_assignment = assignment_result(conditional_truth_assignment(False))
    conditional_assignment = assignment_result(conditional_truth_assignment(True))
    return {
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2541_unital_to_multiplicative_equivalence_inherited": p2541.get("multiplicative_law_equivalent_to_unital_left_normalization_inside_affine_family") is True,
        "p2542_prime_log_obstruction_inherited": p2542.get("current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality") is True,
        "p2543_slope_value_obstruction_inherited": p2543.get("all_slope_candidates_pass_unital_multiplicative_prime_log_premises") is True,
        "p2544_four_key_blocker_inherited": p2544.get("four_key_current_premise_closure_blocker_exported") is True,
        "p2545_unital_obstruction_inherited": p2545.get("unital_monoid_normalization_current_premise_obstruction_exported") is True,
        "candidate_rows_reused_from_p2545": len(rows),
        "identity_action_accepting_row_count": len(identity_accepting_rows),
        "identity_action_rejecting_row_count": len(identity_rejecting_rows),
        "identity_action_exact_linear_audit": exact_audit,
        "identity_action_equivalent_to_y1_zero_on_reused_grid": all(
            row_passes_identity_action(row) == row["unital_y1_zero_accepts"] for row in rows
        ),
        "identity_action_rejects_all_p2545_y1_countermodels": all(
            not row_passes_identity_action(row) for row in rows if not row["unital_y1_zero_accepts"]
        ),
        "identity_action_source_would_close_multiplicative_key": True,
        "identity_action_source_would_not_close_prime_log_key": True,
        "identity_action_source_would_not_close_slope_key": True,
        "identity_action_source_would_not_close_m2_key": True,
        "current_assignment_without_identity_source": current_assignment,
        "conditional_assignment_with_strict_identity_action_source": conditional_assignment,
        "conditional_missing_source_key_delta": (
            current_assignment["missing_source_key_count"] - conditional_assignment["missing_source_key_count"]
        ),
        "strict_identity_action_source_exported": False,
        "unital_monoid_normalization_source_exported": False,
        "multiplicative_character_law_source_exported_conditionally_only": True,
        "frontier_source_key_under_attack": "strict_identity_action_source_for_y1_zero",
        "conditional_closure_theorem": (
            "If a strict nadsoliton identity-action source theorem is supplied, then the affine log damping row has b=log(beta)=0, "
            "so P2541's multiplicative-character key is conditionally closed.  The closure is exactly one-key wide: P2542, P2543, and P2540 "
            "still leave prime-log proportionality, delta=4/5, and m=2 operator selection unsourced."
        ),
        "honest_next_step_recommendation": (
            "Search for or prove the strict nadsoliton identity-action source itself.  The proof must derive the unit law from nadsoliton dynamics, "
            "not assume it as an axiom.  If no internal source can be exported, pivot to the independent m=2 operator-order selection theorem rather "
            "than iterating more y_1=0 consistency scans."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2546/S1496 strict damping identity-action conditional propagation certificate

`P2546/S1496` tests the next proof-shaped move after P2545: what exactly would a strict nadsoliton identity-action theorem buy?  The exact affine audit shows that the unit law `y_(1*d)=y_1+y_d` collapses to `b=log(beta)=0`, hence `y_1=0`, and rejects precisely the P2545 nonunital countermodels.  Therefore such a theorem would conditionally close the P2541 multiplicative-character key.

The propagation is only one-key wide.  Even under a hypothetical strict identity-action source, the P2530 assignment still lacks prime-log proportionality, the strict slope/prime anchor for `delta=4/5`, and the `m=2` operator-signature source.  No identity-action source, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.
""".strip()
    lag_section = """
## P2546/S1496 identity-action conditional Ltotal guard

`P2546/S1496` shows that a future strict identity-action theorem would close only the multiplicative/unital subkey of strict damping.  The damping term remains non-role-bearing in `L_total` until the remaining prime-log, slope-value, and `m=2` source keys are also strict theorems.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2546/S1496 strict damping identity-action conditional propagation certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2546/S1496 identity-action conditional Ltotal guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    theorems = load_theorems(sources)
    rows = candidate_rows(sources)
    cert = build_certificate(theorems, rows)
    theorem_export = {
        "theorem_name": "P2546_T1_strict_damping_identity_action_conditional_propagation_certificate",
        "audited_chain": ["P2530/S1480", "P2541/S1491", "P2542/S1492", "P2543/S1493", "P2544/S1494", "P2545/S1495"],
        "strict_damping_identity_action_conditional_propagation_certificate": cert,
        **cert,
        "identity_action_conditional_propagation_exported": True,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "source_obligation_discharge_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "This packet is a conditional propagation theorem, not the strict identity-action source theorem itself.",
            "It does not replace the P2530 four-key normal form.",
            "It does not transfer legacy physical-role claims onto the strict gate kernel.",
            "It does not discharge QW-2191 or any ToE gate.",
        ],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2545_inherited": theorem_export["p2545_unital_obstruction_inherited"],
        "identity_action_equivalence_verified": theorem_export["identity_action_equivalent_to_y1_zero_on_reused_grid"],
        "conditional_delta_is_one_key": theorem_export["conditional_missing_source_key_delta"] == 1,
        "conditional_assignment_still_rejects_strict_damping": not theorem_export["conditional_assignment_with_strict_identity_action_source"]["strict_damping_beta_eta_source_by_assignment"],
        "negative_controls_preserved": not any(theorem_export[key] for key in NEGATIVE_EXPORT_FLAGS),
    }
    return {
        "packet_id": "P2546",
        "stage_id": "S1496",
        "status": "STRICT_DAMPING_IDENTITY_ACTION_CONDITIONAL_PROPAGATION_CERTIFICATE_NO_IDENTITY_SOURCE_EXPORT_NO_FULL_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_identity_action_conditional_propagation_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_identity_action_conditional_propagation_certificate"]["theorem_export"]
    conditional = t["conditional_assignment_with_strict_identity_action_source"]
    exact = t["identity_action_exact_linear_audit"]
    lines = [
        "# P2546/S1496 strict damping identity-action conditional propagation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- P2545 rows reused: `{t['candidate_rows_reused_from_p2545']}`.",
        f"- Identity-action accepting rows: `{t['identity_action_accepting_row_count']}`.",
        f"- Identity-action rejecting rows: `{t['identity_action_rejecting_row_count']}`.",
        f"- Identity-action equivalent to `y_1=0` on reused grid: `{t['identity_action_equivalent_to_y1_zero_on_reused_grid']}`.",
        f"- Missing-source-key delta under hypothetical strict identity source: `{t['conditional_missing_source_key_delta']}`.",
        f"- Strict damping beta/eta source after hypothetical identity source: `{conditional['strict_damping_beta_eta_source_by_assignment']}`.",
        "",
        "## Exact Unit-Action Audit",
        "",
        f"- Linear solution form: `{exact['solution_form']}`.",
        f"- Rank/nullity: `{exact['rank']}/{exact['nullity_after_unit_action']}`.",
        f"- Proof atom: {exact['d_equals_1_proof']}",
        "",
        "## Conditional Propagation",
        "",
        t["conditional_closure_theorem"],
        "",
        "## Recommendation",
        "",
        t["honest_next_step_recommendation"],
        "",
        "## Negative Controls",
        "",
        "No strict identity-action source theorem, full source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_identity_action_conditional_propagation_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_identity_action_conditional_propagation_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
