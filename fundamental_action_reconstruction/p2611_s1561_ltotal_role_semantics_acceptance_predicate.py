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
OUT = GEN / "p2611_s1561_ltotal_role_semantics_acceptance_predicate.json"
MD = GEN / "p2611_s1561_ltotal_role_semantics_acceptance_predicate.md"

SOURCE_FILES = {
    "P2610_CRITICAL_REVALIDATION": GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.json",
    "P2608_QUARANTINED_ROLE_TRANSFER": GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.json",
}

ACCEPTANCE_GATES = [
    "formal_role_semantics_defined",
    "strict_damping_source_revalidated_after_p2610",
    "bridge_completion_revalidated_after_p2610",
    "role_transfer_revalidated_after_p2610",
]

ROLE_SEMANTICS_AXIOMS = [
    {
        "axiom_id": "R1_typed_term",
        "statement": "A candidate role-bearing L_total term must have a typed term name, variable domain, and target layer.",
    },
    {
        "axiom_id": "R2_variational_action",
        "statement": "The term must contribute through a variational action slot, not merely through narrative bookkeeping.",
    },
    {
        "axiom_id": "R3_source_predicate",
        "statement": "Every coefficient/operator appearing in the term must have a non-quarantined strict source predicate.",
    },
    {
        "axiom_id": "R4_bridge_predicate",
        "statement": "If the term is imported through the legacy-to-strict bridge, the bridge predicate must be non-quarantined and physically derived.",
    },
    {
        "axiom_id": "R5_explanatory_role_map",
        "statement": "The term's physical role must be an explicit map from term variation to the claimed observable/structural effect.",
    },
    {
        "axiom_id": "R6_guardrail_noninterference",
        "statement": "The role must not silently discharge QW-2191, ToE, APD sourcehood, or unrelated legacy physical roles.",
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "strict_damping_beta_eta_source_reexported",
    "bridge_completion_reexported",
    "role_transfer_theorem_reexported",
    "role_bearing_ltotal_accepted",
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
        "new_packet": "P2611|S1561|role-bearing semantics|role semantics.*L_total|Ltotal role semantics",
        "intended_research_nonduplication": "accepted role-bearing term|candidate bookkeeping target|role-bearing acceptance predicate|P2608.*role semantics|bridge revalidation.*role semantics",
        "p2610_frontier": "P2610|S1560|GF2.*physical-origin|ROLE_SEMANTICS|role-bearing.*quarantine|bridge completion.*quarantine",
        "guardrails": "K_legacy_ont|K_strict_gate|QW-2191|ToE closure|APD source|legacy physical-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def acceptance_truth_table() -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(ACCEPTANCE_GATES)):
        assignment = dict(zip(ACCEPTANCE_GATES, values))
        accepts = all(assignment.values())
        rows.append({
            "assignment": assignment,
            "accepted_role_bearing_ltotal_term": accepts,
            "missing_gates": [gate for gate, value in assignment.items() if not value],
            "missing_gate_count": sum(1 for value in assignment.values() if not value),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2610_payload = load_json(SOURCE_FILES["P2610_CRITICAL_REVALIDATION"])
    p2608_payload = load_json(SOURCE_FILES["P2608_QUARANTINED_ROLE_TRANSFER"])
    p2610 = theorem(p2610_payload, "p2601_p2608_critical_revalidation_audit")
    p2608 = theorem(p2608_payload, "strict_damping_role_transfer_to_ltotal_theorem")
    quarantined = set(p2610.get("quarantined_packet_ids_after_revalidation", []))
    rows = acceptance_truth_table()
    accepting_rows = [row for row in rows if row["accepted_role_bearing_ltotal_term"]]
    current_assignment = {
        "formal_role_semantics_defined": True,
        "strict_damping_source_revalidated_after_p2610": "P2601" not in quarantined and "P2602" not in quarantined,
        "bridge_completion_revalidated_after_p2610": "P2607" not in quarantined,
        "role_transfer_revalidated_after_p2610": "P2608" not in quarantined,
    }
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2611_T1_ltotal_role_semantics_acceptance_predicate",
        "audited_chain": ["P2608/S1558", "P2610/S1560"],
        "role_semantics_definition_exported": True,
        "role_semantics_statement": (
            "A role-bearing L_total term is accepted only when it is a typed variational action term whose coefficients and operators are strictly sourced, whose legacy-to-strict bridge is non-quarantined when applicable, whose physical role is supplied by an explicit variation-to-effect map, and whose use does not silently discharge unrelated selector, APD, legacy-role, or ToE claims."
        ),
        "role_semantics_axioms": ROLE_SEMANTICS_AXIOMS,
        "acceptance_gates": ACCEPTANCE_GATES,
        "acceptance_truth_table": {
            "truth_table_rows": rows,
            "truth_table_row_count": len(rows),
            "accepting_row_count": len(accepting_rows),
            "accepting_row": accepting_rows[0],
            "current_assignment": current_assignment,
            "current_assignment_accepts": all(current_assignment.values()),
        },
        "p2608_claim_before_revalidation": {
            "role_bearing_ltotal_exported_in_p2608": p2608.get("role_bearing_ltotal_exported") is True,
            "term_name": p2608.get("role_bearing_ltotal_term", {}).get("term_name"),
        },
        "p2610_quarantine_inherited": {
            "quarantined_packet_ids_after_revalidation": sorted(quarantined),
            "p2607_bridge_quarantined": "P2607" in quarantined,
            "p2608_role_transfer_quarantined": "P2608" in quarantined,
        },
        "operational_result": (
            "P2611 supplies the formal role-semantics predicate requested by P2610, but it does not re-accept the strict damping/compression term: P2607 bridge completion and P2608 role transfer remain quarantined, and P2601/P2602 source support remains unrevalidated."
        ),
        "recommended_next_honest_step": (
            "Prove the P2607 GF(2) physical-origin and invariance theorem before attempting to re-enable P2608 role-bearing L_total acceptance."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2610_quarantine_inherited": "P2607" in quarantined and "P2608" in quarantined,
        "role_semantics_definition_exported": theorem_export["role_semantics_definition_exported"],
        "truth_table_has_sixteen_rows": theorem_export["acceptance_truth_table"]["truth_table_row_count"] == 16,
        "truth_table_has_one_accepting_row": theorem_export["acceptance_truth_table"]["accepting_row_count"] == 1,
        "current_assignment_rejected": theorem_export["acceptance_truth_table"]["current_assignment_accepts"] is False,
        "no_reexported_source": theorem_export["strict_damping_beta_eta_source_reexported"] is False,
        "no_reexported_bridge": theorem_export["bridge_completion_reexported"] is False,
        "no_reexported_role_transfer": theorem_export["role_transfer_theorem_reexported"] is False,
        "no_accepted_role_bearing_ltotal": theorem_export["role_bearing_ltotal_accepted"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2611",
        "stage_id": "S1561",
        "status": "P2611_ROLE_SEMANTICS_DEFINED_BUT_LTOTAL_ACCEPTANCE_REMAINS_REJECTED_BY_P2610_QUARANTINES_NO_NEW_SOURCE_BRIDGE_ROLE_EXPORT",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "ltotal_role_semantics_acceptance_predicate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2610_CRITICAL_REVALIDATION": sha256_json(p2610_payload),
                "P2608_QUARANTINED_ROLE_TRANSFER": sha256_json(p2608_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["ltotal_role_semantics_acceptance_predicate"]["theorem_export"]
    table = t["acceptance_truth_table"]
    lines = [
        "# P2611/S1561 L_total role-semantics acceptance predicate", "",
        f"Status: `{payload['status']}`", "", "## Role semantics definition", "",
        t["role_semantics_statement"], "", "## Computed acceptance predicate", "",
        f"- Formal role semantics defined: `{t['role_semantics_definition_exported']}`.",
        f"- Truth-table rows: `{table['truth_table_row_count']}`.",
        f"- Accepting rows: `{table['accepting_row_count']}`.",
        f"- Current assignment accepts: `{table['current_assignment_accepts']}`.",
        f"- P2607 bridge quarantined: `{t['p2610_quarantine_inherited']['p2607_bridge_quarantined']}`.",
        f"- P2608 role transfer quarantined: `{t['p2610_quarantine_inherited']['p2608_role_transfer_quarantined']}`.", "",
        "## Operational result", "", t["operational_result"], "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "",
        "## Scope guards", "",
        "P2611 exports no revalidated strict damping source, no bridge completion, no role-transfer theorem, no accepted role-bearing L_total term, no legacy physical-role transfer, no QW-2191 discharge, no APD source, and no ToE closure.", "", "## Fingerprint", "",
        f"`{payload['ltotal_role_semantics_acceptance_predicate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2611/S1561 Ltotal role-semantics acceptance predicate

`P2611/S1561` supplies the formal acceptance semantics missing from P2608: a role-bearing `L_total` term must be typed, variational, strictly sourced, bridge-valid when imported from legacy, equipped with an explicit variation-to-effect role map, and guardrail-safe.  The computed four-gate truth table has one accepting row, but the current P2610-inherited assignment rejects acceptance because P2601/P2602 source support, P2607 bridge completion, and P2608 role transfer remain quarantined.
""".strip()
    lag_section = """
## P2611/S1561 role-semantics Ltotal guard

`P2611/S1561` defines what role-bearing `L_total` acceptance means, but it does not accept the strict damping/compression term.  The term remains a candidate bookkeeping target until the P2607 GF(2) physical-origin/invariance theorem and subsequent P2608 role-transfer revalidation are supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2611/S1561 Ltotal role-semantics acceptance predicate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2611/S1561 role-semantics Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
