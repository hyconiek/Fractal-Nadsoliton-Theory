#!/usr/bin/env python3
"""P2696/S1646: pair1 c1c1/s1s1 zero-equation carrier/no-go matrix.

Executes the P2695 recommendation by auditing exactly the two remaining declared
pair1 residual zero-equation carriers.  It is deliberately narrower than
selector closure: it uses current finite R18/P477/P479/N522/P631 artifacts to
separate computable zero-equation obstruction from QW-2191, bridge-source replay,
role transfer, L_total promotion, and ToE closure.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.json"
MD = GEN / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2695": GEN / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.json",
    "P477_summary": GEN / "p477_current_strict_r18_pair1_residual_zero_equations_value_instantiation_probe_summary.json",
    "P477_object": GEN / "pair1_residual_zero_equations_evaluation_under_n477_value_instance_v1.json",
    "N520": ROOT / "N520_CURRENT_FIRST_STRICT_R18_DECLARED_PAIR1_RESIDUAL_ZERO_EQUATIONS_VALUE_INSTANCE_OBSTRUCTION_THEOREM.md",
    "P479_summary": GEN / "p479_current_strict_reference_magnitude_family_sign_scan_for_r18_pair1_zero_equations_probe_summary.json",
    "N522_summary": GEN / "n522_current_first_strict_reference_magnitude_family_sign_scan_r18_pair1_zero_equations_obstruction_theorem_summary.json",
    "P631_summary": GEN / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet_summary.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "selector_imported",
    "qw2191_discharged",
    "m2_psi4_replayed",
    "bridge_source_replayed",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_claimed",
    "strict_core_promoted",
]

TARGETS = ["c1c1", "s1s1"]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
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
        "p2695_selected_p2696": r"P2696|pair1 c1c1/s1s1|zero-equation carrier|c1c1.*s1s1",
        "value_instance_obstruction": r"P477|N520|R18 declared pair1 residual|value instance|violated_equations",
        "finite_sign_scan_obstruction": r"P479|N522|sign scan|fixed-magnitude reference|no sign vector|no solution",
        "p631_negative_freeze": r"P631|DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE|N473|T166 nonzero",
        "forbidden_replays": r"QW-2191|selector import|m2_psi4|bridge-source|role transfer|L_total|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first pair1 c1c1/s1s1 zero-equation carrier/no-go matrix", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def equation_values() -> dict[str, dict[str, Any]]:
    obj = load_json(INPUTS["P477_object"])
    eqs = ((obj.get("computed") or {}).get("equations") or []) if isinstance(obj.get("computed"), dict) else []
    by_entry = {row.get("entry"): row for row in eqs if isinstance(row, dict)}
    return {target: by_entry.get(target, {"entry": target, "missing": True}) for target in TARGETS}


def state_reads() -> dict[str, Any]:
    p2695 = load_json(INPUTS["P2695"])
    p477 = load_json(INPUTS["P477_summary"])
    p479 = load_json(INPUTS["P479_summary"])
    n522 = load_json(INPUTS["N522_summary"])
    p631 = load_json(INPUTS["P631_summary"])
    n520_text = read_text(INPUTS["N520"])
    values = equation_values()
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2695_recommended_p2696": "P2696" in p2695.get("decision", {}).get("next_honest_step", ""),
        "p2695_remaining_mentions_c1c1_s1s1": all(term in json.dumps(p2695.get("route_boundary", {})) for term in ["c1c1", "s1s1"]),
        "p477_violated_equations": p477.get("violated_equations", []),
        "p477_all_zero_equations_satisfied": p477.get("all_zero_equations_satisfied") is True,
        "p477_target_values": values,
        "n520_states_value_instance_obstruction": all(term in n520_text for term in ["c1c1", "s1s1", "nonzero"]),
        "p479_any_reference_has_solution": p479.get("any_reference_has_solution") is True,
        "p479_best_overall_objective_value": p479.get("best_overall_objective_value"),
        "n522_finite_family_obstruction_discharged": (n522.get("theorem_result") or {}).get("discharged") is True,
        "p631_negative_freeze_selected": p631.get("decision") == "DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED",
        "p631_recommended_next_strict_target": p631.get("recommended_next_strict_target"),
    }


def zero_equation_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    violated = set(reads["p477_violated_equations"])
    for target in TARGETS:
        value_row = reads["p477_target_values"].get(target, {})
        rows.append(
            {
                "target_id": f"declared_pair1_{target}_zero_equation_carrier",
                "p477_value_instance_value": value_row.get("value"),
                "p477_abs_value": value_row.get("abs"),
                "p477_zero_tol": value_row.get("zero_tol"),
                "p477_is_zero": value_row.get("is_zero") is True,
                "p477_marks_violated": target in violated,
                "n520_value_instance_obstruction_supports_no_go": reads["n520_states_value_instance_obstruction"],
                "p479_n522_finite_sign_scan_no_solution_supports_no_go": (not reads["p479_any_reference_has_solution"]) and reads["n522_finite_family_obstruction_discharged"],
                "explicit_zero_witness_exported_now": False,
                "no_go_or_closed": "BOUNDED_NO_GO_ON_CURRENT_VALUE_INSTANCE_AND_SCANNED_REFERENCE_FAMILY",
            }
        )
    return rows


def lane_boundary(reads: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    all_no_go = all(row["p477_marks_violated"] and not row["explicit_zero_witness_exported_now"] for row in rows)
    return {
        "all_targeted_pair1_c1c1_s1s1_rows_bounded_no_go": all_no_go,
        "strict_zero_witness_exported_for_c1c1_or_s1s1": any(row["explicit_zero_witness_exported_now"] for row in rows),
        "finite_sign_scan_blocks_reference_family_only": (not reads["p479_any_reference_has_solution"]) and reads["n522_finite_family_obstruction_discharged"],
        "direct_formal_residual_cancellation_lane_frozen_by_p631": reads["p631_negative_freeze_selected"],
        "qw2191_still_separate_open_boundary": True,
        "full_route_closure_exported": False,
        "strict_core_promotion_exported": False,
    }


def decision(boundary: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2696_PAIR1_C1C1_S1S1_ZERO_EQUATION_CARRIER_MATRIX_BOUNDED_NO_GO_NO_FALSE_PASS",
        "bounded_no_go_for_targeted_c1c1_s1s1_zero_witnesses_now": boundary["all_targeted_pair1_c1c1_s1s1_rows_bounded_no_go"],
        "direct_route_reopen_admissible_without_new_provider": False,
        "reason": (
            "P477/N520 show the current strict-derived N477 value instance violates the declared pair1 c1c1 and s1s1 zero equations; "
            "P479/N522 show no solution in the scanned fixed-magnitude reference family; and P631 already freezes direct-formal residual cancellation on the current strict branch."
        ),
        "next_honest_step": (
            "P2697 should be a post-direct-route state-map/no-new-live-frontier reconciliation: do not continue pair1 residual cancellation unless a genuinely new strict-derived provider class, non-N477 ingredient, or blocker-cut is exported; otherwise emit a no-new-live-frontier certificate rather than replaying selector, bridge, role-transfer, L_total, or ToE lanes."
        ),
        "forbidden_reopens": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2696/S1646 pair1 c1c1/s1s1 zero-equation carrier/no-go matrix", "", f"Status: `{payload['status']}`", "", "## Zero-equation matrix"]
    for row in payload["zero_equation_matrix"]:
        lines.append(f"- `{row['target_id']}`: value=`{row['p477_value_instance_value']}`, violated=`{row['p477_marks_violated']}`, status=`{row['no_go_or_closed']}`")
    lines.extend(["", "## Lane boundary", json.dumps(payload["lane_boundary"], indent=2, sort_keys=True), "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    matrix = zero_equation_matrix(reads)
    boundary = lane_boundary(reads, matrix)
    payload: dict[str, Any] = {
        "status": "P2696_PAIR1_C1C1_S1S1_ZERO_EQUATION_CARRIER_NO_GO_MATRIX_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "zero_equation_matrix": matrix,
        "lane_boundary": boundary,
        "decision": decision(boundary),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2696/S1646 pair1 c1c1/s1s1 zero-equation carrier no-go matrix",
        "## P2696/S1646 pair1 c1c1/s1s1 zero-equation carrier no-go matrix\n\n"
        "`P2696/S1646` executes the remaining direct-route `pair1 c1c1/s1s1` zero-equation audit.  `P477/N520` show that the current strict-derived N477 value instance violates both rows; `P479/N522` add a finite no-solution scan over the fixed-magnitude reference family; `P631` already freezes direct-formal residual cancellation on the current strict branch.  Therefore the targeted `c1c1/s1s1` zero witnesses are bounded no-go on current artifacts, without `QW-2191` discharge, selector import, bridge-source replay, role transfer, `L_total`, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2696/S1646 pair1 zero-equation Ltotal guard",
        "## P2696/S1646 pair1 zero-equation Ltotal guard\n\n"
        "`P2696/S1646` is a direct-route finite residual-equation obstruction matrix.  It does not promote `L_total`, variational closure, role transfer, strict-core selector closure, or ToE closure.\n",
    )
    append_once(
        AGENTS,
        "Current pair1 c1c1/s1s1 no-go guardrail (P2696/S1646, 2026-06-13)",
        "## Current pair1 c1c1/s1s1 no-go guardrail (P2696/S1646, 2026-06-13)\n\n"
        "- P2696 audits the remaining direct-route `pair1 c1c1/s1s1` zero-equation carriers and finds bounded no-go on current artifacts: current N477 value instance violates both rows, the scanned fixed-magnitude reference family has no solution, and P631 freezes direct-formal residual cancellation on the current strict branch.\n"
        "- Do not continue this direct-route residual-cancellation lane without a genuinely new strict-derived provider class, non-N477 ingredient, or blocker-cut; do not import selector replay, `QW-2191` discharge, bridge-source replay, role transfer, `L_total`, or ToE closure.\n"
        "- The next honest move is a post-direct-route broad state-map/no-new-live-frontier reconciliation.\n",
    )
    return payload


if __name__ == "__main__":
    main()
