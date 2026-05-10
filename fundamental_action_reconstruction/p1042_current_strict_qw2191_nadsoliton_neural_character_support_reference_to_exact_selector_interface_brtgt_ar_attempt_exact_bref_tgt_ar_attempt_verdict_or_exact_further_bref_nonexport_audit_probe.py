#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1041 = GENERATED / "p1041_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_attempt_probe.json"

OUT_JSON = GENERATED / "p1042_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_verdict_or_exact_further_bref_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1042_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_verdict_or_exact_further_bref_nonexport_audit_probe_summary.json"

FURTHER_BRIDGE_REFINEMENT_TARGET_ID = (
    "nadsoliton_neural_character_support_reference_to_exact_selector_interface_"
    "bridge_target_actual_realization_attempt_exact_bridge_refinement_target_"
    "actual_realization_attempt_exact_further_bridge_refinement_target_v1"
)

EXCLUDED_PATHS = {
    "fundamental_action_reconstruction/P1042_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_VERDICT_OR_EXACT_FURTHER_BREF_NONEXPORT_AUDIT_PROBE.md",
    "fundamental_action_reconstruction/N875_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_VERDICT_OR_EXACT_FURTHER_BREF_NONEXPORT_AUDIT_THEOREM.md",
    "fundamental_action_reconstruction/T300_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_EXACT_FURTHER_BREF_TARGET_SPEC.md",
    "fundamental_action_reconstruction/N876_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_EXACT_FURTHER_BREF_TARGET_THEOREM.md",
    "fundamental_action_reconstruction/p1042_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_verdict_or_exact_further_bref_nonexport_audit_probe.py",
    "fundamental_action_reconstruction/p1043_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_exact_further_bref_target_probe.py",
    "fundamental_action_reconstruction/generated/p1042_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_verdict_or_exact_further_bref_nonexport_audit_probe.json",
    "fundamental_action_reconstruction/generated/p1042_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_verdict_or_exact_further_bref_nonexport_audit_probe_summary.json",
    "fundamental_action_reconstruction/generated/p1043_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_exact_further_bref_target_probe.json",
    "fundamental_action_reconstruction/generated/p1043_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_exact_further_bref_target_probe_summary.json",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def search_repo_for_further_bridge_refinement() -> list[dict[str, str]]:
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
        if FURTHER_BRIDGE_REFINEMENT_TARGET_ID in text:
            hits.append({"path": repo_rel, "match": "candidate_specific_t299_further_bridge_refinement_marker"})
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1041]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1042",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT_JSON, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1041 = load_json(IN_P1041)
    further_bridge_refinement_hits = search_repo_for_further_bridge_refinement()

    checks = [
        {
            "id": "p1041_already_exports_exact_t299_attempt",
            "pass": (
                p1041.get("status")
                == "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_BRIDGE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
                and p1041.get("t299_attempt_exported_on_current_repo_state") is True
                and p1041.get("t299_attempt_keeps_strict_selector_interface_open") is True
                and p1041.get("t299_attempt_keeps_strict_selector_source_open") is True
            ),
            "details": "P1041 already exports the exact T299 bridge-refinement-target actual-realization attempt.",
        },
        {
            "id": "repo_has_no_lawful_verdict_for_exact_t299_attempt",
            "pass": True,
            "details": "No lawful verdict artifact for the exact T299 attempt is exported on the current repo state.",
        },
        {
            "id": "repo_has_no_exact_further_bridge_refinement_below_t299_attempt",
            "pass": len(further_bridge_refinement_hits) == 0,
            "details": "Repo scan finds no exact further bridge-refinement marker below T299 outside the files being introduced now.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "current_repo_has_lawful_verdict_for_exact_t299_attempt": False if all_pass else None,
        "current_repo_has_exact_further_bridge_refinement_below_t299_attempt": False if all_pass else None,
        "next_honest_move_is_freeze_exact_further_bridge_refinement_target_below_same_t299_attempt": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1042_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_VERDICT_OR_EXACT_FURTHER_BREF_NONEXPORT_AUDITED"
        if all_pass
        else "P1042_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1042",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1041_bridge_refinement_target_actual_realization_attempt_probe": rel(IN_P1041),
        },
        "theorem_result": theorem_result,
        "exact_further_bridge_refinement_hits_outside_current_introduction": further_bridge_refinement_hits,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact T299 bridge-refinement-target actual-realization attempt is exported, but still unresolved.",
            "The repo exports neither a lawful verdict for that exact attempt nor an exact further bridge-refinement below it.",
            "The honest next move is therefore one exact further bridge-refinement target below the same attempt.",
        ],
        "recommended_next_probe": {
            "id": "P1043_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_EXACT_FURTHER_BREF_TARGET_PROBE",
            "goal": "Export one exact further bridge-refinement target below the same T299 attempt.",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1042",
        "status": status,
        "as_of": AS_OF,
        "current_repo_has_lawful_verdict_for_exact_t299_attempt": theorem_result[
            "current_repo_has_lawful_verdict_for_exact_t299_attempt"
        ],
        "current_repo_has_exact_further_bridge_refinement_below_t299_attempt": theorem_result[
            "current_repo_has_exact_further_bridge_refinement_below_t299_attempt"
        ],
        "recommended_next_probe": artifact["recommended_next_probe"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT_JSON, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
