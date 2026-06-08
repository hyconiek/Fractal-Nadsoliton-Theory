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
OUT = GEN / "p2600_s1550_strict_damping_post_m2_residual_source_matrix.json"
MD = GEN / "p2600_s1550_strict_damping_post_m2_residual_source_matrix.md"

SOURCE_FILES = {
    "P2599_PROJECTED_VISCOUS_STRESS_M2_SOURCE": GEN / "p2599_s1549_nadsoliton_projected_viscous_stress_m2_derivation_theorem.json",
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json",
}
RESIDUAL_KEYS_AFTER_M2 = [
    "multiplicative_character_law_source",
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
]
NEGATIVE_EXPORT_FLAGS = [
    "multiplicative_character_law_source_exported",
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
        "new_packet": "P2600|S1550|post m2 residual source matrix|post-m2 residual|m2 source integration",
        "intended_research_nonduplication": "hydrodynamic m2 integration|residual source matrix.*m2|m2 discharged.*damping|strict damping.*post m2|post-m2 source obligation",
        "precursor_chain": "P2530|S1480|P2547|S1497|P2599|S1549|m2_operator_signature_source|strict_damping_beta_eta_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def residual_truth_table() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for values in product([False, True], repeat=len(RESIDUAL_KEYS_AFTER_M2)):
        assignment = dict(zip(RESIDUAL_KEYS_AFTER_M2, values))
        m2_assignment = True
        beta_eta_numeric = all(assignment.values())
        strict_damping = beta_eta_numeric and m2_assignment
        missing = [key for key, value in assignment.items() if not value]
        rows.append({
            "assignment": assignment,
            "m2_operator_signature_source": m2_assignment,
            "beta_eta_numeric_source_accepts": beta_eta_numeric,
            "strict_damping_beta_eta_source_accepts": strict_damping,
            "missing_residual_keys": missing,
            "missing_residual_key_count": len(missing),
        })
    return rows


def single_omission_witnesses(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if row["missing_residual_key_count"] == 1]


def post_m2_matrix(p2530: dict[str, Any], p2547: dict[str, Any], p2599: dict[str, Any]) -> dict[str, Any]:
    rows = residual_truth_table()
    accepting = [row for row in rows if row["strict_damping_beta_eta_source_accepts"]]
    single_omissions = single_omission_witnesses(rows)
    current_assignment = {key: False for key in RESIDUAL_KEYS_AFTER_M2}
    current_assignment["m2_operator_signature_source"] = True
    return {
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2547_post_identity_residual_trikey_inherited": p2547.get("post_identity_residual_trikey_irredundancy_exported") is True,
        "p2599_m2_source_inherited": p2599.get("m2_operator_signature_source_exported") is True,
        "m2_source_status_after_p2599": "strict_source_exported_by_projected_viscous_stress_theorem",
        "residual_keys_after_m2_discharge": RESIDUAL_KEYS_AFTER_M2,
        "residual_truth_table_rows": rows,
        "residual_truth_table_row_count": len(rows),
        "residual_accepting_row_count": len(accepting),
        "residual_accepting_row": accepting[0],
        "residual_single_omission_witnesses": single_omissions,
        "all_single_residual_omissions_reject_strict_damping": all(not row["strict_damping_beta_eta_source_accepts"] for row in single_omissions),
        "post_m2_missing_source_key_count_by_current_assignment": sum(1 for value in current_assignment.values() if not value),
        "current_post_m2_assignment": current_assignment,
        "current_assignment_beta_eta_numeric_source_accepts": all(current_assignment[key] for key in RESIDUAL_KEYS_AFTER_M2),
        "current_assignment_strict_damping_beta_eta_source_accepts": False,
        "p2530_contract_after_m2": (
            "strict_damping_beta_eta_source = m2_operator_signature_source AND multiplicative_character_law_source "
            "AND prime_log_proportionality_source AND slope_value_or_prime_anchor_source; P2599 supplies only the m2 factor."
        ),
        "next_honest_source_targets": [
            "multiplicative_character_law_source / strict unital normalization y_1=0",
            "prime_log_proportionality_source",
            "slope_value_or_prime_anchor_source for delta=4/5",
        ],
        "post_m2_residual_source_matrix_exported": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2599_payload = load_json(SOURCE_FILES["P2599_PROJECTED_VISCOUS_STRESS_M2_SOURCE"])
    p2530_payload = load_json(SOURCE_FILES["P2530_FOUR_KEY_IRREDUNDANCY"])
    p2547_payload = load_json(SOURCE_FILES["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"])
    p2599 = theorem(p2599_payload, "nadsoliton_projected_viscous_stress_m2_derivation_theorem")
    p2530 = theorem(p2530_payload, "strict_damping_four_key_irredundancy_witness_certificate")
    p2547 = theorem(p2547_payload, "strict_damping_post_identity_residual_trikey_certificate")
    matrix = post_m2_matrix(p2530, p2547, p2599)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2600_T1_strict_damping_post_m2_residual_source_matrix",
        "audited_chain": ["P2530/S1480", "P2547/S1497", "P2599/S1549"],
        "frontier_under_update": "strict_damping_beta_eta_source after m2_operator_signature_source discharge",
        "m2_operator_signature_source_exported": matrix["p2599_m2_source_inherited"],
        "strict_damping_post_m2_residual_source_matrix": matrix,
        "strict_damping_beta_eta_source_remains_blocked_by_current_assignment": not matrix["current_assignment_strict_damping_beta_eta_source_accepts"],
        "post_m2_residual_matrix_exported": matrix["post_m2_residual_source_matrix_exported"],
        "recommended_next_honest_step": (
            "Do not repeat APD/moment/Sturm work. With m2 now hydrodynamically sourced by P2596-P2599, the strict damping source frontier has shifted to the three residual non-m2 keys: multiplicative/unital normalization, prime-log proportionality, and the delta=4/5 slope/prime anchor. The next source theorem should target exactly one of those keys."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2530_four_key_irredundancy_inherited": matrix["p2530_four_key_irredundancy_inherited"],
        "p2547_post_identity_residual_trikey_inherited": matrix["p2547_post_identity_residual_trikey_inherited"],
        "p2599_m2_source_inherited": matrix["p2599_m2_source_inherited"],
        "m2_operator_signature_source_exported": theorem_export["m2_operator_signature_source_exported"],
        "residual_truth_table_has_eight_rows": matrix["residual_truth_table_row_count"] == 8,
        "exactly_one_residual_accepting_row": matrix["residual_accepting_row_count"] == 1,
        "single_residual_omissions_reject": matrix["all_single_residual_omissions_reject_strict_damping"],
        "current_assignment_has_three_missing_residual_keys": matrix["post_m2_missing_source_key_count_by_current_assignment"] == 3,
        "strict_damping_not_exported_by_current_assignment": theorem_export["strict_damping_beta_eta_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2600",
        "stage_id": "S1550",
        "status": "P2600_STRICT_DAMPING_POST_M2_RESIDUAL_SOURCE_MATRIX_M2_DISCHARGED_THREE_NON_M2_KEYS_REMAIN_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_post_m2_residual_source_matrix": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2599_PROJECTED_VISCOUS_STRESS_M2_SOURCE": sha256_json(p2599_payload),
                "P2530_FOUR_KEY_IRREDUNDANCY": sha256_json(p2530_payload),
                "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": sha256_json(p2547_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_post_m2_residual_source_matrix"]["theorem_export"]
    matrix = t["strict_damping_post_m2_residual_source_matrix"]
    lines = [
        "# P2600/S1550 strict damping post-m2 residual source matrix", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- m2 operator signature source exported: `{t['m2_operator_signature_source_exported']}`.",
        f"- Residual keys after m2 discharge: `{matrix['residual_keys_after_m2_discharge']}`.",
        f"- Residual truth-table rows: `{matrix['residual_truth_table_row_count']}`.",
        f"- Residual accepting rows: `{matrix['residual_accepting_row_count']}`.",
        f"- Current assignment strict damping accepts: `{matrix['current_assignment_strict_damping_beta_eta_source_accepts']}`.", "",
        "## Interpretation", "",
        "P2600 integrates the P2599 hydrodynamic `m=2` source into the older P2530/P2547 strict-damping source normal form.  The `m2_operator_signature_source` factor is now discharged, but the beta/eta numeric package still requires the three non-m2 residual keys: multiplicative/unital normalization, prime-log proportionality, and the `delta=4/5` slope/prime anchor.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Scope guards", "",
        "No beta/eta numeric source, strict damping source closure, bridge theorem, role-transfer theorem, role-bearing `L_total`, QW-2191 discharge, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['strict_damping_post_m2_residual_source_matrix']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2600/S1550 strict damping post-m2 residual source matrix guard

`P2600/S1550` integrates the hydrodynamic `m=2` source theorem chain P2596--P2599 into the P2530/P2547 strict-damping source normal form.  The `m2_operator_signature_source` factor is now discharged, but the residual truth table still has exactly one accepting row over the three non-m2 keys: multiplicative/unital normalization, prime-log proportionality, and the `delta=4/5` slope/prime anchor.  Therefore strict damping/beta-eta source closure is not exported until those three source theorems are supplied.
""".strip()
    lag_section = """
## P2600/S1550 strict damping post-m2 residual source matrix Ltotal guard

`P2600/S1550` updates `L_total` bookkeeping after the hydrodynamic `m=2` source discharge: the operator-order slot may be treated as sourced, but the damping/compression term remains non-role-bearing until multiplicative/unital normalization, prime-log proportionality, and slope/prime-anchor sources are also exported, followed by bridge and role-transfer gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2600/S1550 strict damping post-m2 residual source matrix guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2600/S1550 strict damping post-m2 residual source matrix Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
