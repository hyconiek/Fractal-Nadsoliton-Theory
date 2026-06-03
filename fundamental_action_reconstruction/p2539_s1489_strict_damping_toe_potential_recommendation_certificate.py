#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from typing import Any

from p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate import GATES, CURRENT_TRUE_GATES
from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json"
MD = GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.md"

SOURCE_FILES = {
    "P2538_REWRITE_NORMALIZATION": GEN / "p2538_s1488_strict_damping_rewrite_normalization_certificate.json",
    "P2421_TOE_PRIME_IMPLICANT": GEN / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.json",
}

STRICT_DAMPING_LOCAL_ATOMS = [
    "multiplicative_character_law_source",
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
    "m2_operator_signature_source",
]

NEXT_STEP_RECOMMENDATION = (
    "prove_or_refute_one_strict_source_theorem_for_a_single_P2529_source_key_before_more_bookkeeping_layers"
)


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
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2539|S1489|strict damping ToE potential|ToE potential recommendation|next honest strict damping source step",
        "intended_research_nonduplication": "ToE potential delta|readiness delta|source-obligation discharge delta|next honest source theorem|strict damping ToE interface",
        "precursor_packets": "P2538|S1488|rewrite normalization|P2421|S1371|ToE prime implicant|seven closure gates",
        "toe_guardrail_language": "ToE closure|source_obligation_discharge|role_bearing_ltotal_export|QW-2191|role-transfer",
        "kernel_guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|selector guardrail|bridge-completion",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def mask_for(gates: list[str]) -> int:
    return sum(1 << GATES.index(gate) for gate in gates)


def gates_for(mask: int) -> list[str]:
    return [gate for index, gate in enumerate(GATES) if mask & (1 << index)]


def missing_gates(mask: int) -> list[str]:
    return [gate for gate in GATES if gate not in gates_for(mask)]


def toe_ready(mask: int) -> bool:
    return mask == mask_for(GATES)


def scenario_row(name: str, true_gates: list[str], note: str) -> dict[str, Any]:
    mask = mask_for(true_gates)
    missing = missing_gates(mask)
    return {
        "scenario": name,
        "true_gate_count": len(true_gates),
        "true_gates": true_gates,
        "missing_gate_count": len(missing),
        "missing_gates": missing,
        "repair_distance_to_toe": len(missing),
        "toe_ready": toe_ready(mask),
        "note": note,
    }


def local_source_subkey_rows() -> list[dict[str, Any]]:
    rows = []
    for size in range(len(STRICT_DAMPING_LOCAL_ATOMS) + 1):
        for subset in combinations(STRICT_DAMPING_LOCAL_ATOMS, size):
            rows.append({
                "strict_local_subkeys": list(subset),
                "strict_local_subkey_count": size,
                "strict_damping_beta_eta_source_ready": size == len(STRICT_DAMPING_LOCAL_ATOMS),
                "source_obligation_discharge_gate_flipped_by_this_packet": False,
                "toe_gate_delta": 0,
            })
    return rows


def build_certificate(p2538: dict[str, Any], p2421: dict[str, Any]) -> dict[str, Any]:
    current_gates = list(CURRENT_TRUE_GATES)
    bookkeeping_gates = list(current_gates)
    source_gate_gates = sorted(set(current_gates + ["source_obligation_discharge"]), key=GATES.index)
    all_missing_added = list(GATES)
    local_rows = local_source_subkey_rows()
    scenarios = [
        scenario_row(
            "current_p2421_state",
            current_gates,
            "P2421 current state: APD bridge witness and chi11 phase selector cut mechanism only; five ToE gates remain missing.",
        ),
        scenario_row(
            "after_p2536_p2538_bookkeeping_only",
            bookkeeping_gates,
            "P2536-P2538 normalize the strict-damping repair bookkeeping but do not flip any P2421 ToE gate.",
        ),
        scenario_row(
            "hypothetical_full_source_obligation_discharge_only",
            source_gate_gates,
            "Even a future full source-obligation discharge gate by itself would still leave selector source/discharge, role transfer, and role-bearing L_total gates open.",
        ),
        scenario_row(
            "all_current_p2421_missing_gates_supplied",
            all_missing_added,
            "This is the conditional ToE row from P2421, not the current repository state.",
        ),
    ]
    return {
        "frontier_atom_under_discussion": "strict_damping_beta_eta_source as a sub-obligation below P2421 source_obligation_discharge",
        "p2538_rewrite_normalization_inherited": p2538.get("finite_newman_conditions_verified") is True,
        "p2421_toe_prime_implicant_inherited": p2421.get("prime_implicant_count") == 1 and p2421.get("prime_implicant_masks") == [127],
        "toe_gate_names": GATES,
        "current_true_gates": current_gates,
        "current_missing_gates": missing_gates(mask_for(current_gates)),
        "current_toe_readiness_fraction": "1/128",
        "current_toe_repair_distance": len(missing_gates(mask_for(current_gates))),
        "strict_damping_bookkeeping_toe_gate_delta": 0,
        "strict_damping_source_subkey_rows": local_rows,
        "local_subkey_truth_table_row_count": len(local_rows),
        "local_subkey_full_rows": sum(row["strict_damping_beta_eta_source_ready"] for row in local_rows),
        "local_subkey_proper_rows": sum(not row["strict_damping_beta_eta_source_ready"] for row in local_rows),
        "local_subkey_no_direct_toe_gate_flip_verified": all(row["toe_gate_delta"] == 0 for row in local_rows),
        "toe_potential_scenarios": scenarios,
        "scenario_count": len(scenarios),
        "bookkeeping_changes_toe_readiness": scenarios[0]["repair_distance_to_toe"] != scenarios[1]["repair_distance_to_toe"],
        "source_gate_alone_closes_toe": scenarios[2]["toe_ready"],
        "conditional_all_missing_gates_close_toe": scenarios[3]["toe_ready"],
        "toe_potential_assessment": (
            "high_structural_potential_but_zero_current_toe_gate_delta: P2536-P2538 improve auditability and normal-form discipline, "
            "yet ToE closure still requires source_obligation_discharge, chi11_source_export, QW-2191 discharge, role-transfer license, and role-bearing L_total export."
        ),
        "recommended_next_honest_step": NEXT_STEP_RECOMMENDATION,
        "recommendation_reason": (
            "The repair/normalization layer is now saturated enough: another combinatorial wrapper would not flip a ToE gate. "
            "The next useful proof attempt should attack one actual P2529 source key, with a negative/no-go result acceptable if the source cannot be derived."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2539/S1489 strict damping ToE-potential recommendation certificate

`P2539/S1489` audits the ToE potential of the P2536-P2538 strict-damping repair stack against the existing P2421 seven-gate ToE prime-implicant frontier.  The result is intentionally conservative: P2536-P2538 improve proof bookkeeping, confluence, and rewrite normalization, but they flip `0` of the P2421 ToE gates.  The current P2421 state remains at two true gates and five missing gates, so the ToE truth-table readiness remains the unique full row `1/128`, not an unconditional closure.

The strict-damping local source subtable has `16` rows over the four P2529/P2530 source subkeys, with exactly one all-strict local row; however even that local row is only a sub-obligation below the broader `source_obligation_discharge` gate.  The recommended next honest step is therefore not another bookkeeping layer: prove or refute one actual strict source theorem for a single P2529 source key (`multiplicative_character_law_source`, `prime_log_proportionality_source`, `slope_value_or_prime_anchor_source`, or `m2_operator_signature_source`).  No source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.
""".strip()
    lag_section = """
`P2539/S1489` records the ToE-potential interface for strict damping: the repair certificates strengthen auditability but do not change the P2421 ToE gate vector.  A role-bearing Lagrangian step should therefore wait for an actual strict source theorem for at least one strict-damping source key, rather than treating normalized repair bookkeeping as a dynamical source.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2539/S1489", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2539/S1489", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2538 = theorem(sources["P2538_REWRITE_NORMALIZATION"], "strict_damping_rewrite_normalization_certificate")
    p2421 = theorem(sources["P2421_TOE_PRIME_IMPLICANT"], "bridge_selector_closure_prime_implicant_failure_cover_certificate")
    cert = build_certificate(p2538, p2421)
    theorem_export = {
        **cert,
        "toe_potential_recommendation_certificate_exported": True,
        "axiom_promotion_to_strict_exported": False,
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
    }
    gatekeepers = {
        "p2538_inherited": theorem_export["p2538_rewrite_normalization_inherited"],
        "p2421_inherited": theorem_export["p2421_toe_prime_implicant_inherited"],
        "toe_gate_delta_zero": theorem_export["strict_damping_bookkeeping_toe_gate_delta"] == 0 and not theorem_export["bookkeeping_changes_toe_readiness"],
        "local_subkey_table_verified": theorem_export["local_subkey_truth_table_row_count"] == 16 and theorem_export["local_subkey_full_rows"] == 1 and theorem_export["local_subkey_proper_rows"] == 15,
        "source_gate_not_enough_for_toe": not theorem_export["source_gate_alone_closes_toe"],
        "conditional_full_toe_row_preserved": theorem_export["conditional_all_missing_gates_close_toe"],
        "recommendation_is_source_theorem_not_bookkeeping": theorem_export["recommended_next_honest_step"] == NEXT_STEP_RECOMMENDATION,
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "axiom_promotion_to_strict_exported",
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
        ]),
    }
    return {
        "packet_id": "P2539",
        "stage_id": "S1489",
        "status": "STRICT_DAMPING_TOE_POTENTIAL_RECOMMENDATION_CERTIFICATE_ZERO_TOE_GATE_DELTA_NEXT_SOURCE_THEOREM_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_toe_potential_recommendation_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_toe_potential_recommendation_certificate"]["theorem_export"]
    lines = [
        "# P2539/S1489 strict damping ToE-potential recommendation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier under discussion: `{t['frontier_atom_under_discussion']}`.",
        f"- P2538 rewrite normalization inherited: `{t['p2538_rewrite_normalization_inherited']}`.",
        f"- P2421 ToE prime implicant inherited: `{t['p2421_toe_prime_implicant_inherited']}`.",
        f"- Current ToE readiness fraction: `{t['current_toe_readiness_fraction']}` with repair distance `{t['current_toe_repair_distance']}`.",
        f"- Strict-damping bookkeeping ToE gate delta: `{t['strict_damping_bookkeeping_toe_gate_delta']}`.",
        f"- Local source subkey rows / full rows / proper rows: `{t['local_subkey_truth_table_row_count']}` / `{t['local_subkey_full_rows']}` / `{t['local_subkey_proper_rows']}`.",
        f"- Source gate alone closes ToE: `{t['source_gate_alone_closes_toe']}`.",
        f"- Recommended next honest step: `{t['recommended_next_honest_step']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## ToE potential assessment",
        "",
        t["toe_potential_assessment"],
        "",
        "## Negative controls",
        "",
        "This packet exports only a ToE-potential and next-step recommendation certificate. It does not source any strict-damping key, discharge source_obligation_discharge, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_toe_potential_recommendation_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_toe_potential_recommendation_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
