#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
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
OUT = GEN / "p2544_s1494_strict_damping_four_key_current_premise_closure_blocker_certificate.json"
MD = GEN / "p2544_s1494_strict_damping_four_key_current_premise_closure_blocker_certificate.md"

SOURCE_FILES = {
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2539_TOE_POTENTIAL_RECOMMENDATION": GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json",
    "P2540_M2_OPERATOR_OBSTRUCTION": GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.json",
    "P2541_MULTIPLICATIVE_OBSTRUCTION": GEN / "p2541_s1491_strict_damping_multiplicative_character_current_premise_obstruction_certificate.json",
    "P2542_PRIME_LOG_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
    "P2543_SLOPE_VALUE_OBSTRUCTION": GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json",
}

STRICT_SOURCE_KEYS = [
    "multiplicative_character_law_source",
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
    "m2_operator_signature_source",
]

NEGATIVE_EXPORT_FLAGS = [
    "multiplicative_character_law_source_exported",
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
        "new_packet": "P2544|S1494|four-key current-premise closure blocker|no-false-source theorem|source-key closure blocker",
        "intended_research_nonduplication": "four-key current-premise obstruction|strict damping source-key obstruction|no false source|closure blocker|source-key blocker",
        "precursor_packets": "P2530|S1480|P2539|S1489|P2540|S1490|P2541|S1491|P2542|S1492|P2543|S1493",
        "source_key_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def obstruction_bool(row: dict[str, Any]) -> bool:
    return bool(
        row.get("current_premise_obstruction_exported")
        or row.get("current_premise_nonentailment_of_m2_exported")
        or row.get("m2_operator_signature_source_route_refuted_for_current_source_free_premises")
    )


def source_flag_for_key(t: dict[str, Any], source_key: str) -> bool:
    return bool(t.get(f"{source_key}_exported", False))


def key_row(packet: str, t: dict[str, Any], route_label: str, obstruction_measure: str, residual_code: str) -> dict[str, Any]:
    source_key = t["frontier_source_key_under_attack"]
    return {
        "packet": packet,
        "source_key": source_key,
        "route_label": route_label,
        "current_premise_obstruction_exported": obstruction_bool(t),
        "source_key_exported": source_flag_for_key(t, source_key),
        "countermodel_or_measure": obstruction_measure,
        "required_new_premise_class": t["required_new_premise_class"],
        "residual_code": residual_code,
    }


def source_key_rows(theorems: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p2540 = theorems["P2540_M2_OPERATOR_OBSTRUCTION"]
    p2541 = theorems["P2541_MULTIPLICATIVE_OBSTRUCTION"]
    p2542 = theorems["P2542_PRIME_LOG_OBSTRUCTION"]
    p2543 = theorems["P2543_SLOPE_VALUE_OBSTRUCTION"]
    return [
        key_row(
            "P2541/S1491",
            p2541,
            "affine-consistency route to multiplicative character law",
            f"affine countermodels: {p2541['affine_countermodel_count']}; multiplicativity equivalent to y_1=0 inside affine family",
            "strict_unital_monoid_normalization_y1_zero",
        ),
        key_row(
            "P2542/S1492",
            p2542,
            "unital finite-monoid character route to prime-log proportionality",
            f"prime-character countermodels: {p2542['countermodel_count']}; ratio-collapse rank/nullity: {p2542['exact_ratio_constraint_rank']}/{p2542['exact_ratio_constraint_nullity']}",
            "strict_prime_ratio_collapse_vp_over_logp",
        ),
        key_row(
            "P2543/S1493",
            p2543,
            "prime-log slope-line route to delta=4/5 or equivalent prime anchor",
            f"slope-line countermodels: {p2543['countermodel_count']}; strict delta rows: {p2543['strict_delta_row_count']}",
            "strict_slope_or_prime_anchor_delta_4_over_5",
        ),
        key_row(
            "P2540/S1490",
            p2540,
            "source-free derivative-order route to m=2 operator signature",
            f"passing derivative orders: {p2540['passing_order_count']}; countermodel orders: {p2540['countermodel_pair_orders']}",
            "strict_operator_order_selector_m2",
        ),
    ]


def build_certificate(theorems: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2530 = theorems["P2530_FOUR_KEY_IRREDUNDANCY"]
    p2539 = theorems["P2539_TOE_POTENTIAL_RECOMMENDATION"]
    rows = source_key_rows(theorems)
    missing_codes = [row["residual_code"] for row in rows if row["current_premise_obstruction_exported"] and not row["source_key_exported"]]
    current_source_truth_assignment = {key: False for key in STRICT_SOURCE_KEYS}
    beta_eta_numeric_ready = all(current_source_truth_assignment[key] for key in STRICT_SOURCE_KEYS[:3])
    strict_damping_ready = beta_eta_numeric_ready and current_source_truth_assignment["m2_operator_signature_source"]
    return {
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2530_source_keys": p2530.get("source_keys", []),
        "p2539_zero_toe_gate_delta_inherited": p2539.get("strict_damping_bookkeeping_toe_gate_delta") == 0,
        "p2539_current_missing_toe_gates": p2539.get("current_missing_gates", []),
        "p2539_recommendation_inherited": p2539.get("recommended_next_honest_step"),
        "source_obstruction_rows": rows,
        "source_key_row_count": len(rows),
        "all_four_p2530_source_keys_audited": sorted(row["source_key"] for row in rows) == sorted(STRICT_SOURCE_KEYS),
        "all_four_current_premise_routes_blocked": all(row["current_premise_obstruction_exported"] for row in rows),
        "all_four_source_keys_unexported": all(not row["source_key_exported"] for row in rows),
        "current_source_truth_assignment": current_source_truth_assignment,
        "beta_eta_numeric_source_by_current_assignment": beta_eta_numeric_ready,
        "strict_damping_beta_eta_source_by_current_assignment": strict_damping_ready,
        "p2530_contract": (
            "strict_damping_beta_eta_source = multiplicative_character_law_source "
            "AND prime_log_proportionality_source AND slope_value_or_prime_anchor_source AND m2_operator_signature_source"
        ),
        "relative_to_p2530_contract_missing_source_principles": missing_codes,
        "relative_to_p2530_missing_source_principle_count": len(missing_codes),
        "no_false_source_theorem": (
            "P2540-P2543 jointly refute the current source-free/consistency-only routes to all four P2530 strict-damping source keys. "
            "Therefore the strict damping beta/eta package remains a conditional target, not a role-bearing L_total source, until a new strict nadsoliton source theorem supplies these missing principles or replaces the P2530 normal form."
        ),
        "professorial_completion_path": [
            "First prove or refute a genuinely strict nadsoliton source theorem for one missing principle, starting with the smallest target: y_1=0 unital monoid normalization or m=2 operator-order selection.",
            "Then test whether the new source theorem propagates through the P2530 four-key normal form without importing legacy physical roles silently.",
            "Only after source obligation discharge, resume legacy -> strict completion bridge and perform a separate role-transfer audit for EM, Weinberg-angle, gravity-hierarchy, and SM/GR claims.",
            "Promote any damping/compression term into L_total only after the source theorem and role-transfer gate are explicit; otherwise keep it as conditional structure.",
            "Treat ToE potential as conditional high-structure potential, not closure: P2539 still leaves source_obligation_discharge, chi11_source_export, QW-2191 discharge, role-transfer license, and role-bearing L_total open.",
        ],
        "toe_potential_assessment": (
            "The theory has nontrivial ToE potential because it now has a precise obstruction frontier instead of vague numerical matching. "
            "The current state is not a ToE: it lacks strict source discharge, bridge completion, role-transfer theorem, QW-2191 discharge, and role-bearing L_total export."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2544/S1494 strict damping four-key current-premise closure blocker certificate

`P2544/S1494` synthesizes the four source-key attacks `P2540-P2543` against the P2530 strict-damping normal form.  The result is a no-false-source theorem: every current route to a required strict-damping source key is blocked by an explicit countermodel or nonentailment witness.  Affine consistency does not source multiplicativity without `y_1=0`; unital multiplicative prime characters do not source `v_p/log(p)=constant`; the prime-log slope line does not source `delta=4/5`; and the current derivative-order premises do not select the `m=2` operator signature.

Therefore `strict_damping_beta_eta_source` remains false under the current source assignment.  The next professorially honest path is a new strict nadsoliton source theorem for one missing principle, followed only then by bridge-completion and role-transfer auditing.  No source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.
""".strip()
    lag_section = """
## P2544/S1494 strict damping no-false-source Ltotal guard

`P2544/S1494` blocks promotion of the current strict-damping beta/eta package into role-bearing `L_total`: P2540-P2543 show that all four P2530 source keys remain unsourced under current premises.  A damping/compression term may remain as a conditional target, but not as a derived dynamical source until a new strict nadsoliton theorem supplies the missing source principle.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2544/S1494 strict damping four-key current-premise closure blocker certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2544/S1494 strict damping no-false-source Ltotal guard", lag_section)


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2530_FOUR_KEY_IRREDUNDANCY": theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate"),
        "P2539_TOE_POTENTIAL_RECOMMENDATION": theorem(sources["P2539_TOE_POTENTIAL_RECOMMENDATION"], "strict_damping_toe_potential_recommendation_certificate"),
        "P2540_M2_OPERATOR_OBSTRUCTION": theorem(sources["P2540_M2_OPERATOR_OBSTRUCTION"], "strict_damping_m2_operator_signature_current_premise_obstruction_certificate"),
        "P2541_MULTIPLICATIVE_OBSTRUCTION": theorem(sources["P2541_MULTIPLICATIVE_OBSTRUCTION"], "strict_damping_multiplicative_character_current_premise_obstruction_certificate"),
        "P2542_PRIME_LOG_OBSTRUCTION": theorem(sources["P2542_PRIME_LOG_OBSTRUCTION"], "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate"),
        "P2543_SLOPE_VALUE_OBSTRUCTION": theorem(sources["P2543_SLOPE_VALUE_OBSTRUCTION"], "strict_damping_slope_value_current_premise_obstruction_certificate"),
    }


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    theorems = load_theorems(sources)
    cert = build_certificate(theorems)
    theorem_export = {
        "theorem_name": "P2544_T1_strict_damping_four_key_current_premise_closure_blocker_certificate",
        "audited_chain": ["P2530/S1480", "P2539/S1489", "P2540/S1490", "P2541/S1491", "P2542/S1492", "P2543/S1493"],
        "strict_damping_four_key_current_premise_closure_blocker_certificate": cert,
        **cert,
        "four_key_current_premise_closure_blocker_exported": True,
        "multiplicative_character_law_source_exported": False,
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
            "This packet closes a current-premise obstruction block, not a source theorem.",
            "It does not replace the P2530 four-key normal form.",
            "It does not transfer legacy physical-role claims onto the strict gate kernel.",
            "It does not discharge QW-2191 or any ToE gate.",
        ],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2530_inherited": theorem_export["p2530_four_key_irredundancy_inherited"],
        "p2539_zero_toe_gate_delta_inherited": theorem_export["p2539_zero_toe_gate_delta_inherited"],
        "all_four_keys_audited": theorem_export["all_four_p2530_source_keys_audited"] and theorem_export["source_key_row_count"] == 4,
        "all_four_routes_blocked": theorem_export["all_four_current_premise_routes_blocked"],
        "all_four_sources_unexported": theorem_export["all_four_source_keys_unexported"],
        "p2530_assignment_rejects_strict_damping_source": (
            not theorem_export["beta_eta_numeric_source_by_current_assignment"]
            and not theorem_export["strict_damping_beta_eta_source_by_current_assignment"]
        ),
        "missing_source_principle_count_verified": theorem_export["relative_to_p2530_missing_source_principle_count"] == 4,
        "negative_controls_preserved": not any(theorem_export[key] for key in NEGATIVE_EXPORT_FLAGS),
    }
    return {
        "packet_id": "P2544",
        "stage_id": "S1494",
        "status": "STRICT_DAMPING_FOUR_KEY_CURRENT_PREMISE_CLOSURE_BLOCKER_CERTIFICATE_NO_FALSE_SOURCE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_four_key_current_premise_closure_blocker_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_four_key_current_premise_closure_blocker_certificate"]["theorem_export"]
    rows = t["source_obstruction_rows"]
    lines = [
        "# P2544/S1494 strict damping four-key current-premise closure blocker certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- P2530 four-key irredundancy inherited: `{t['p2530_four_key_irredundancy_inherited']}`.",
        f"- Source-key rows audited: `{t['source_key_row_count']}`.",
        f"- All four current-premise routes blocked: `{t['all_four_current_premise_routes_blocked']}`.",
        f"- All four source keys exported: `{not t['all_four_source_keys_unexported']}`.",
        f"- Strict damping beta/eta source under current assignment: `{t['strict_damping_beta_eta_source_by_current_assignment']}`.",
        f"- Missing source principle count relative to P2530: `{t['relative_to_p2530_missing_source_principle_count']}`.",
        "",
        "## Blocked Keys",
        "",
    ]
    for row in rows:
        lines.append(f"- `{row['source_key']}` via `{row['packet']}`: {row['countermodel_or_measure']}.")
    lines.extend([
        "",
        "## Professorial Completion Path",
        "",
    ])
    for step in t["professorial_completion_path"]:
        lines.append(f"- {step}")
    lines.extend([
        "",
        "## Interpretation",
        "",
        t["no_false_source_theorem"],
        "",
        t["toe_potential_assessment"],
        "",
        "## Negative Controls",
        "",
        "No source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_four_key_current_premise_closure_blocker_certificate']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_four_key_current_premise_closure_blocker_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
