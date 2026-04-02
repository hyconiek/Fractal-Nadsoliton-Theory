#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1039 = GENERATED / "p1039_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_probe.json"

OUT_JSON = GENERATED / "p1040_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1040_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_nonexport_audit_probe_summary.json"

ACTUAL_REALIZATION_OBJECT_ID = (
    "nadsoliton_neural_character_support_reference_to_exact_selector_interface_"
    "bridge_target_actual_realization_attempt_exact_bridge_refinement_target_"
    "actual_realization_v1"
)

EXCLUDED_PATHS = {
    "fundamental_action_reconstruction/P1040_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
    "fundamental_action_reconstruction/N873_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
    "fundamental_action_reconstruction/T299_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
    "fundamental_action_reconstruction/N874_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
    "fundamental_action_reconstruction/p1040_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_nonexport_audit_probe.py",
    "fundamental_action_reconstruction/p1041_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_attempt_probe.py",
    "fundamental_action_reconstruction/generated/p1040_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_nonexport_audit_probe.json",
    "fundamental_action_reconstruction/generated/p1040_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_nonexport_audit_probe_summary.json",
    "fundamental_action_reconstruction/generated/p1041_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_attempt_probe.json",
    "fundamental_action_reconstruction/generated/p1041_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_attempt_probe_summary.json",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def search_repo_for_actual_realization() -> list[dict[str, str]]:
    hits: list[dict[str, str]] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        repo_rel = rel(path)
        if repo_rel in EXCLUDED_PATHS:
            continue
        if path.suffix.lower() not in {".md", ".py", ".json", ".tex"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if ACTUAL_REALIZATION_OBJECT_ID in text:
            hits.append({"path": repo_rel, "match": "candidate_specific_t298_actual_realization_marker"})
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1039]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1040",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT_JSON, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1039 = load_json(IN_P1039)
    actual_realization_hits = search_repo_for_actual_realization()

    checks = [
        {
            "id": "p1039_already_exports_exact_t298_target",
            "pass": (
                p1039.get("status")
                == "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_EXPORTED"
                and p1039.get("t298_target_exported_on_current_repo_state") is True
                and p1039.get("t298_target_keeps_strict_selector_interface_open") is True
                and p1039.get("t298_target_keeps_strict_selector_source_open") is True
            ),
            "details": "P1039 already exports the exact T298 bridge-refinement target.",
        },
        {
            "id": "repo_has_no_exact_actual_realization_of_t298_target",
            "pass": len(actual_realization_hits) == 0,
            "details": "Repo scan finds no exact actual-realization marker for the T298 target outside the files being introduced now.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "current_repo_has_exported_actual_realization_of_t298_target": False if all_pass else None,
        "t298_target_still_remains_future_only_not_actual_export": True if all_pass else None,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t298_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1040_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if all_pass
        else "P1040_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1040",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1039_bridge_refinement_target_probe": rel(IN_P1039),
        },
        "theorem_result": theorem_result,
        "actual_realization_hits_outside_current_introduction": actual_realization_hits,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact T298 bridge-refinement target is exported, but still future-only.",
            "No exact actual realization of that target is exported on the current repo state.",
            "The honest next move is therefore one exact actual-realization attempt below that same T298 target.",
        ],
        "recommended_next_probe": {
            "id": "P1041_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_PROBE",
            "goal": "Export one exact actual-realization attempt below the same T298 bridge-refinement target.",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1040",
        "status": status,
        "as_of": AS_OF,
        "current_repo_has_exported_actual_realization_of_t298_target": theorem_result[
            "current_repo_has_exported_actual_realization_of_t298_target"
        ],
        "t298_target_still_remains_future_only_not_actual_export": theorem_result[
            "t298_target_still_remains_future_only_not_actual_export"
        ],
        "recommended_next_probe": artifact["recommended_next_probe"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT_JSON, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
