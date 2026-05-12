#!/usr/bin/env python3
"""P1223: execute a minimal W1 witness check under strict guardrails."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_or_build_input_payload(input_path: Path, target_id: str, w1: dict) -> dict:
    if input_path.exists():
        existing = _load_json(input_path)
        if isinstance(existing, dict) and existing.get("witness_id") == target_id:
            return existing

    return {
        "packet": "P1223_INPUT",
        "as_of": "2026-05-11",
        "witness_id": target_id,
        "target": w1.get("target"),
        "required_artifact": w1.get("required_artifact"),
        "pass_condition": w1.get("pass_condition"),
        "available_inputs": {
            "strict_selector_source_theorem_exported": False,
            "explicit_symmetry_breaking_axiom_declared": False,
            "scoped_non_strict_tag_if_axiom_used": False,
        },
        "guardrails": {
            "qw_2191_is_real_obstruction": True,
            "strict_closure_claim_allowed": False,
        },
    }




def _classify_discharge(strict_theorem_available: bool, non_strict_axiom_path_available: bool) -> str:
    if strict_theorem_available:
        return "STRICT_PATH_DISCHARGE"
    if non_strict_axiom_path_available:
        return "NON_STRICT_AXIOM_DISCHARGE"
    return "NO_DISCHARGE"

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1222", type=Path, default=GEN / "p1222_next_witness_attack_planner_summary.json")
    parser.add_argument("--p1192", type=Path, default=GEN / "p1192_selector_premise_witness_map_summary.json")
    parser.add_argument("--input", type=Path, default=GEN / "p1223_w1_input_payload.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1223_w1_execution_summary.json")
    args = parser.parse_args()

    p1222 = _load_json(args.p1222)
    p1192 = _load_json(args.p1192)

    selected = p1222.get("selected_open_witness") or {}
    target_id = selected.get("id", "W1_selector_uniqueness_bridge")

    witness_map = {
        w.get("id"): w for w in p1192.get("witness_obligations", []) if isinstance(w, dict) and w.get("id")
    }
    w1 = witness_map.get(target_id, selected)

    input_payload = _load_or_build_input_payload(args.input, target_id, w1)

    args.input.write_text(json.dumps(input_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    available_inputs = input_payload.get("available_inputs", {})
    strict_theorem_available = bool(available_inputs.get("strict_selector_source_theorem_exported", False))
    non_strict_axiom_path_available = (
        bool(available_inputs.get("explicit_symmetry_breaking_axiom_declared", False))
        and bool(available_inputs.get("scoped_non_strict_tag_if_axiom_used", False))
    )

    # W1 allows either strict theorem path OR explicitly tagged non-strict axiom path.
    pass_condition_met = strict_theorem_available or non_strict_axiom_path_available
    discharge_mode = _classify_discharge(strict_theorem_available, non_strict_axiom_path_available)
    witness_status = "DISCHARGED" if pass_condition_met else "OPEN"

    # Strict-core closure may only be discussed when strict theorem path is available.
    strict_core_closure_eligibility = strict_theorem_available

    out = {
        "packet": "P1223",
        "as_of": "2026-05-11",
        "input_artifact": str(args.input),
        "witness_under_test": {
            "id": target_id,
            "target": w1.get("target"),
            "pass_condition": w1.get("pass_condition"),
        },
        "evaluation": {
            "strict_selector_source_theorem_exported": strict_theorem_available,
            "explicit_symmetry_breaking_axiom_declared": bool(available_inputs.get("explicit_symmetry_breaking_axiom_declared", False)),
            "scoped_non_strict_tag_if_axiom_used": bool(available_inputs.get("scoped_non_strict_tag_if_axiom_used", False)),
            "strict_theorem_path_available": strict_theorem_available,
            "non_strict_axiom_path_available": non_strict_axiom_path_available,
            "pass_condition_met": pass_condition_met,
            "discharge_mode": discharge_mode,
            "strict_core_closure_eligibility": strict_core_closure_eligibility,
        },
        "witness_status": witness_status,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Minimal execution checkpoint completed; non-strict axiom discharge never implies strict-core closure.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1223] witness_id={target_id} witness_status={witness_status} wrote {args.out}")


if __name__ == "__main__":
    main()
