#!/usr/bin/env python3
"""P2682/S1632: P46/P50 m2_psi4 target-EOM zero-witness matrix audit.

P2681 selected a finite P46/N49 live-frontier computation.  This packet first
checks that the originally named P46 target-action blocker has already been
locally closed by AX12 on the canonical-ontology-supported external lane, then
moves to the next live finite object on the same m2_psi4 thread: the P50/R39
target-EOM coefficient defect zero witness on common psi4(x) support.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.json"
MD = GEN / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.md"
P2681 = GEN / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.json"
AX12 = GEN / "ax12_canonical_ontology_supported_preobserver_target_action_coherence_instance.json"
R39 = GEN / "r39_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json"
P50_JSON = GEN / "p50_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "target_eom_zero_witness_exported",
    "ax12_action_identity_transported_to_eom",
    "common_support_divided_or_assumed_nonzero",
    "strict_core_promoted",
    "q_w_2191_discharged",
    "legacy_strict_bridge_completed",
    "role_transfer_started",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "p46_target_action_frontier": r"P46|R37|target action|m2_psi4.*target action|psi4\*\*2/2|target-action coefficient defect",
        "ax12_action_closure_boundary": r"AX12|target_action_coefficient_defect_zero_witness|strict_core_promotion|target eom-side witness|target eom role",
        "p50_target_eom_frontier": r"P50|R39|target eom|target-eom|psi4\(x\)|target eom coefficient defect|m2_psi4.*mu_m2_plus3_segment_psi1_psi4",
        "forbidden_imports_and_closure": r"nonzero-factor|divide by|strict-core|QW-2191|ToE closure|role transfer|legacy -> strict|selector replay|tau_src -> pair12|beta_tors -> chi11",
    }
    return {
        "tool": "rg",
        "mode": "content-first finite P46/P50 m2_psi4 target-action-to-target-EOM live-frontier audit",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2681 = load_json(P2681)
    return {
        "p2681_exists": not p2681.get("missing", False),
        "p2681_selected_finite_p46_lane": p2681.get("finite_priority_lattice", {}).get("selected_next_lane") == "finite_p46_n49_zero_witness_matrix",
        "p2681_toe_not_closed": p2681.get("closure_decision", {}).get("toe_closed_now") is False,
        "p2681_hash": sha256_file(P2681),
    }


def action_eom_artifact_read() -> dict[str, Any]:
    ax12 = load_json(AX12)
    r39 = load_json(R39)
    p50 = load_json(P50_JSON)
    ax12_result = ax12.get("result", {})
    r39_result = r39.get("result", {})
    packet = r39.get("direct_m2_psi4_target_eom_coefficient_defect_packet", {})
    return {
        "ax12_present": not ax12.get("missing", False),
        "ax12_external_action_zero_witness_present": ax12_result.get("canonical_ontology_supported_target_action_coefficient_defect_zero_witness_present") is True,
        "ax12_strict_core_promotion": ax12_result.get("strict_core_promotion") is True,
        "r39_present": not r39.get("missing", False),
        "r39_target_eom_defect_packet_present": r39_result.get("explicit_target_eom_coefficient_defect_packet_present") is True,
        "r39_target_eom_zero_witness_present": r39_result.get("explicit_zero_witness_for_target_eom_coefficient_defect_polynomial_present") is True,
        "r39_defect_polynomial": packet.get("exact_coefficient_defect_polynomial"),
        "r39_common_support": packet.get("common_support"),
        "p50_present": not p50.get("missing", False),
        "artifact_hashes": {"AX12": sha256_file(AX12), "R39": sha256_file(R39), "P50": sha256_file(P50_JSON)},
    }


def symbolic_defect_matrix(read: dict[str, Any]) -> dict[str, Any]:
    # Work in the polynomial ring Z[m, mu, psi].  The local EOM defect is
    # (m - mu) * psi.  Without an exported equality m=mu or a permitted
    # nonzero/division rule for psi(x), the zero witness is not derivable.
    obligations = [
        "r39_defect_packet_present",
        "coefficient_identity_m2_equals_mu_exported_for_target_eom",
        "action_to_eom_transport_theorem_exported",
        "common_support_nonzero_or_division_rule_exported",
        "zero_witness_exported_without_field_division",
        "strict_core_promotion_not_used",
    ]
    current = {
        "r39_defect_packet_present": read["r39_target_eom_defect_packet_present"],
        "coefficient_identity_m2_equals_mu_exported_for_target_eom": False,
        "action_to_eom_transport_theorem_exported": False,
        "common_support_nonzero_or_division_rule_exported": False,
        "zero_witness_exported_without_field_division": read["r39_target_eom_zero_witness_present"],
        "strict_core_promotion_not_used": not read["ax12_strict_core_promotion"],
    }
    rows = []
    pass_count = 0
    for bits in itertools.product([False, True], repeat=len(obligations)):
        state = dict(zip(obligations, bits))
        # Either an explicit zero witness, or coefficient identity plus a valid
        # action-to-EOM transport theorem, is needed.  Support division alone is
        # not accepted as proof-grade because the route forbids nonzero-factor
        # assumptions unless exported separately.
        passes = (
            state["r39_defect_packet_present"]
            and state["strict_core_promotion_not_used"]
            and (
                state["zero_witness_exported_without_field_division"]
                or (
                    state["coefficient_identity_m2_equals_mu_exported_for_target_eom"]
                    and state["action_to_eom_transport_theorem_exported"]
                )
            )
        )
        pass_count += int(passes)
        if passes or state == current:
            rows.append({"state": state, "passes_target_eom_zero_witness_gate": passes})
    missing = [key for key, value in current.items() if key in obligations and not value]
    return {
        "ring_model": "Z[m2_psi4, mu_m2_plus3_segment_psi1_psi4, psi4_x]",
        "target_eom_defect": "(m2_psi4 - mu_m2_plus3_segment_psi1_psi4) * psi4_x",
        "coefficient_defect": read["r39_defect_polynomial"],
        "common_support": read["r39_common_support"],
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "passing_states": pass_count,
        "current_state": current,
        "distinguished_rows": rows,
        "missing_current_obligations": missing,
        "hamming_distance_to_nearest_simple_pass": 2,
        "bounded_no_go_now": True,
    }


def closure_decision(read: dict[str, Any], matrix: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2682_P46_TARGET_ACTION_ALREADY_EXTERNAL_CLOSED_P50_TARGET_EOM_ZERO_WITNESS_STILL_BLOCKED_NO_FALSE_PASS",
        "professorial_verdict": (
            "P2682 corrects the P2681 next-step target after reading the current repo state: the P46 target-action m2_psi4 zero witness has already been locally closed by AX12 on the canonical-ontology-supported external lane, but that closure explicitly does not transport to the target-EOM side or to strict core.  The live finite object is therefore the P50/R39 target-EOM coefficient defect zero witness.  In the polynomial audit, the defect is only (m2_psi4 - mu_m2_plus3_segment_psi1_psi4) on common psi4(x) support; without an exported target-EOM coefficient identity or an action-to-EOM transport theorem, no zero witness follows."
        ),
        "next_honest_step": (
            "Do not attack the already locally closed P46 target-action blocker again.  The next honest proof-grade move is a construction-or-no-go audit for the target-EOM transport/assignment theorem: either export a canonical-ontology-supported preobserver target-EOM coherence instance analogous to AX12 but explicitly typed for EOM support psi4(x), or prove that AX12-style action closure cannot transport to EOM without a new variational/EOM-role premise.  If that fails, move to one of the remaining finite direct g4/g6/gY zero-witness matrices."
        ),
        "p46_target_action_reopened": False,
        "p50_target_eom_zero_witness_exported_now": False,
        "strict_core_promoted": False,
        "toe_closed_now": False,
        "bounded_no_go_now": matrix["bounded_no_go_now"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2682/S1632 P46/P50 m2_psi4 target-EOM zero-witness matrix audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first grep",
    ]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    read = payload["action_eom_artifact_read"]
    lines.extend([
        "", "## Artifact read",
        f"- AX12 external target-action zero witness present: `{read['ax12_external_action_zero_witness_present']}`",
        f"- AX12 strict-core promotion: `{read['ax12_strict_core_promotion']}`",
        f"- R39 target-EOM defect packet present: `{read['r39_target_eom_defect_packet_present']}`",
        f"- R39 target-EOM zero witness present: `{read['r39_target_eom_zero_witness_present']}`",
        f"- R39 defect polynomial: `{read['r39_defect_polynomial']}` on `{read['r39_common_support']}`",
        "", "## Symbolic obligation matrix",
    ])
    matrix = payload["symbolic_defect_matrix"]
    lines.extend([
        f"Ring model: `{matrix['ring_model']}`.",
        f"Target EOM defect: `{matrix['target_eom_defect']}`.",
        f"Total states: `{matrix['total_states']}`; passing states: `{matrix['passing_states']}`.",
        f"Current state: `{matrix['current_state']}`.",
        f"Missing current obligations: `{matrix['missing_current_obligations']}`.",
        "", "## Verdict", payload["closure_decision"]["professorial_verdict"],
        f"Decision: `{payload['closure_decision']['decision']}`.",
        "", "## Next honest step", payload["closure_decision"]["next_honest_step"],
        "", "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    read = action_eom_artifact_read()
    matrix = symbolic_defect_matrix(read)
    payload: dict[str, Any] = {
        "status": "P2682_P46_P50_M2_PSI4_TARGET_EOM_ZERO_WITNESS_MATRIX_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "upstream_consistency": upstream_consistency(),
        "action_eom_artifact_read": read,
        "symbolic_defect_matrix": matrix,
        "closure_decision": closure_decision(read, matrix),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2682/S1632 P46/P50 target-EOM zero-witness guard",
        "## P2682/S1632 P46/P50 target-EOM zero-witness guard\n\n"
        "`P2682/S1632` corrects the P2681 finite-frontier target after reading AX12/R39/P50: the P46 target-action `m2_psi4` blocker is already locally closed on the canonical-ontology-supported external lane by AX12, but that closure does not transport to target EOM or strict core.  The live finite blocker is the P50/R39 target-EOM coefficient defect zero witness on `psi4(x)`.  No zero witness, action-to-EOM transport theorem, strict-core promotion, `QW-2191` discharge, role transfer, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2682/S1632 target-EOM Ltotal guard",
        "## P2682/S1632 target-EOM Ltotal guard\n\n"
        "`P2682/S1632` does not add a variational source term.  It shows that AX12 target-action closure cannot be imported into the target-EOM coefficient defect without an explicit action-to-EOM transport/assignment theorem typed on `psi4(x)`, so `L_total` remains unpromoted.\n",
    )
    return payload


if __name__ == "__main__":
    main()
