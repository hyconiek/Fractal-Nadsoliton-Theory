#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation.json"
MD = GEN / "p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation.md"

SOURCE_FILES = {
    "P2611_ROLE_SEMANTICS": GEN / "p2611_s1561_ltotal_role_semantics_acceptance_predicate.json",
    "P2612_GF2_OBSTRUCTION": GEN / "p2612_s1562_p2607_gf2_physical_origin_obstruction.json",
    "P2613_MONOID_UNITAL": GEN / "p2613_s1563_p2601_monoid_action_uniqueness_proof.json",
    "P2614_CONTINUUM_PRIME_LOG": GEN / "p2614_s1564_p2602_continuum_rg_character_prime_log_proof.json",
    "P2615_LINEAR_SLICE_NONBRIDGE": GEN / "p2615_s1565_p2605_linear_slice_nonbridge_obstruction.json",
    "P2608_QUARANTINED_ROLE_TRANSFER": GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.json",
}

ROLE_GATES = [
    "formal_role_semantics_defined",
    "strict_damping_source_revalidated",
    "legacy_to_strict_bridge_revalidated",
    "role_transfer_revalidated",
]

NEGATIVE_EXPORT_FLAGS = [
    "p2608_role_bearing_ltotal_reenabled",
    "role_transfer_theorem_exported",
    "p2607_bridge_completion_revalidated",
    "p2605_full_bridge_revalidated",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
    "apd_source_exported",
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
        "new_packet": "P2616|S1566|P2608.*role.*obstruction|role-bearing.*acceptance.*obstruction",
        "intended_research_nonduplication": "bridge-valid.*role-bearing|non-bridge bookkeeping.*role|Ltotal role semantics.*current assignment|P2608.*quarantine.*predicate",
        "precursor_chain": "P2611|P2612|P2613|P2614|P2615|P2608|role-bearing L_total|GF2|linear-slice non-bridge",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer|QW-2191|ToE closure|APD source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def role_truth_table() -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(ROLE_GATES)):
        assignment = dict(zip(ROLE_GATES, values))
        accepts = all(assignment.values())
        rows.append({
            "assignment": assignment,
            "role_bearing_ltotal_accepted": accepts,
            "missing_gates": [gate for gate, value in assignment.items() if not value],
            "missing_gate_count": sum(1 for value in assignment.values() if not value),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2611 = theorem(payloads["P2611_ROLE_SEMANTICS"], "ltotal_role_semantics_acceptance_predicate")
    p2612 = theorem(payloads["P2612_GF2_OBSTRUCTION"], "p2607_gf2_physical_origin_obstruction")
    p2613 = theorem(payloads["P2613_MONOID_UNITAL"], "p2601_monoid_action_uniqueness_proof")
    p2614 = theorem(payloads["P2614_CONTINUUM_PRIME_LOG"], "p2602_continuum_rg_character_prime_log_proof")
    p2615 = theorem(payloads["P2615_LINEAR_SLICE_NONBRIDGE"], "p2605_linear_slice_nonbridge_obstruction")
    p2608 = theorem(payloads["P2608_QUARANTINED_ROLE_TRANSFER"], "strict_damping_role_transfer_to_ltotal_theorem")
    table = role_truth_table()
    current_assignment = {
        "formal_role_semantics_defined": p2611.get("role_semantics_definition_exported") is True,
        "strict_damping_source_revalidated": p2613.get("p2601_quarantine_lifted_by_p2613") is True and p2614.get("p2602_quarantine_lifted_by_p2614") is True and p2614.get("strict_damping_beta_eta_source_revalidated_under_df_9_5_scope") is True,
        "legacy_to_strict_bridge_revalidated": p2612.get("bridge_completion_revalidated") is True and p2615.get("p2605_full_bridge_revalidated") is True,
        "role_transfer_revalidated": False,
    }
    current_missing = [gate for gate, value in current_assignment.items() if not value]
    proof = {
        "formal_obstruction_statement": "Under the P2611 acceptance predicate, a role-bearing L_total term requires all four gates. After P2613/P2614 the source gate is repaired, but P2612 and P2615 keep the bridge gate false; hence the conjunction is false and P2608 cannot be re-enabled.",
        "proof_steps": [
            "P2611 defines role-bearing L_total acceptance as the conjunction of role semantics, strict source support, bridge validity when imported through legacy, and role-transfer validity.",
            "P2613 and P2614 revalidate only the non-bridge strict damping source subkeys under the retained D_f=9/5 scope.",
            "P2612 obstructs the P2607 GF(2) bridge path, and P2615 obstructs using the P2605 eta=1 linear slice as a nonlinear bridge completion.",
            "Therefore the bridge-valid conjunct remains false even after source repair.",
            "A conjunction with a false conjunct is false, so the P2608 role-bearing L_total export remains rejected.",
        ],
        "only_escape_conditions": [
            "supply a non-quarantined physical-origin bridge theorem for the phase/topological selector data",
            "supply a completion map whose strict-side nonlinear residual additions are physically sourced rather than inferred from the eta=1 slice",
            "then re-run the P2611 role predicate with an explicit variation-to-effect role-transfer theorem",
        ],
    }
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2616_T1_p2608_role_acceptance_obstruction_after_source_revalidation",
        "audited_chain": ["P2608/S1558", "P2611/S1561", "P2612/S1562", "P2613/S1563", "P2614/S1564", "P2615/S1565"],
        "role_acceptance_gates": ROLE_GATES,
        "truth_table": {
            "row_count": len(table),
            "accepting_row_count": sum(1 for row in table if row["role_bearing_ltotal_accepted"]),
            "rows": table,
        },
        "current_assignment_after_p2615": current_assignment,
        "current_missing_gates_after_p2615": current_missing,
        "current_assignment_accepts_role_bearing_ltotal": all(current_assignment.values()),
        "formal_obstruction_proof": proof,
        "p2608_claim_before_revalidation": p2608.get("role_bearing_L_total_exported", p2608.get("role_bearing_ltotal_exported")),
        "p2608_quarantine_retained_by_p2616": not all(current_assignment.values()),
        "strict_damping_bookkeeping_status": "NON_BRIDGE_STRICT_DAMPING_SOURCE_REVALIDATED_BUT_NOT_ROLE_BEARING_LTOTAL",
        "recommended_next_honest_step": "Stop trying to promote P2608 until a physical bridge theorem replaces the P2607 GF(2) obstruction and the P2605 linear-slice obstruction; otherwise keep the damping term as non-role-bearing bookkeeping only.",
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "role_semantics_inherited": current_assignment["formal_role_semantics_defined"],
        "source_gate_repaired": current_assignment["strict_damping_source_revalidated"],
        "bridge_gate_still_false": current_assignment["legacy_to_strict_bridge_revalidated"] is False,
        "role_transfer_gate_still_false": current_assignment["role_transfer_revalidated"] is False,
        "current_assignment_rejects_role_bearing_ltotal": theorem_export["current_assignment_accepts_role_bearing_ltotal"] is False,
        "p2608_quarantine_retained": theorem_export["p2608_quarantine_retained_by_p2616"],
        "exactly_one_truth_table_accepting_row": theorem_export["truth_table"]["accepting_row_count"] == 1,
        "no_p2608_reenable": theorem_export["p2608_role_bearing_ltotal_reenabled"] is False,
        "no_bridge_revalidation": theorem_export["p2607_bridge_completion_revalidated"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2616",
        "stage_id": "S1566",
        "status": "P2616_P2608_ROLE_ACCEPTANCE_REJECTED_AFTER_SOURCE_REVALIDATION_BRIDGE_GATE_STILL_FALSE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2608_role_acceptance_obstruction_after_source_revalidation": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in payloads.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2608_role_acceptance_obstruction_after_source_revalidation"]["theorem_export"]
    proof = t["formal_obstruction_proof"]
    lines = [
        "# P2616/S1566 P2608 role acceptance obstruction after source revalidation", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        proof["formal_obstruction_statement"], "", "## Proof", "",
    ]
    for step in proof["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend([
        "", "## Current assignment", "",
        f"- Formal role semantics defined: `{t['current_assignment_after_p2615']['formal_role_semantics_defined']}`.",
        f"- Strict damping source revalidated: `{t['current_assignment_after_p2615']['strict_damping_source_revalidated']}`.",
        f"- Legacy-to-strict bridge revalidated: `{t['current_assignment_after_p2615']['legacy_to_strict_bridge_revalidated']}`.",
        f"- Role transfer revalidated: `{t['current_assignment_after_p2615']['role_transfer_revalidated']}`.",
        f"- Role-bearing L_total accepted: `{t['current_assignment_accepts_role_bearing_ltotal']}`.",
        f"- Missing gates: `{t['current_missing_gates_after_p2615']}`.", "",
        "## Scope guards", "",
        "P2616 retains the P2608 quarantine. It does not revalidate bridge completion, role transfer, role-bearing L_total, legacy physical-role transfer, QW-2191 discharge, APD sourcehood, or ToE closure.",
        "", "## Fingerprint", "",
        f"`{payload['p2608_role_acceptance_obstruction_after_source_revalidation']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2616/S1566 P2608 role acceptance obstruction after source revalidation

`P2616/S1566` reruns the P2611 role predicate after the P2613/P2614 source repair and the P2615 linear-slice obstruction.  The source gate is now repaired for non-bridge strict damping bookkeeping, but the bridge-valid conjunct remains false because P2612 blocks the GF(2) bridge and P2615 blocks use of the `eta=1` slice as nonlinear completion.  Therefore P2608 role-bearing `L_total` acceptance remains rejected.
""".strip()
    lag_section = """
## P2616/S1566 post-source role-obstruction Ltotal guard

`P2616/S1566` records the exact reason `L_total` role-bearing status is still unavailable: source revalidation is not enough.  The term stays non-role-bearing until a non-quarantined physical bridge and a role-transfer theorem are supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2616/S1566 P2608 role acceptance obstruction after source revalidation", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2616/S1566 post-source role-obstruction Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
