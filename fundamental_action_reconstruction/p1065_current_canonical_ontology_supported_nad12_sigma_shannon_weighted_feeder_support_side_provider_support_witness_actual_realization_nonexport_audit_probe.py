#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F963 = GENERATED / "f963_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_packet_summary.json"
IN_N363 = ROOT / "N363_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md"
IN_F252 = ROOT / "F252_FIRST_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_PACKET.md"

OUT_JSON = GENERATED / "p1065_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_side_provider_support_witness_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1065_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_side_provider_support_witness_actual_realization_nonexport_audit_probe_summary.json"

WITNESS_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_support_side_provider_support_witness_v1"
CURRENT_THEOREM_FILE = (
    "N900_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_"
    "PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1065_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N363_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md",
        "F252_FIRST_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_PACKET.md",
        OUT_JSON.name,
        OUT_SUMMARY.name,
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if WITNESS_NAME in text and "actual_realization_attempt" in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F963, IN_N363, IN_F252]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1065",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f963 = load_json(IN_F963)
    n363_text = load_text(IN_N363)
    f252_text = load_text(IN_F252)
    positive_candidates = scan_positive_actual_realization_candidates()

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    preferred_first_noncyclic_pivot_branch_already_selected = (
        f963.get("status")
        == "F963_EXECUTED_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_EXHAUSTION_AND_FEEDER_PIVOT_PACKET_NO_FALSE_PASS"
        and f963.get("preferred_first_pivot_branch") == WITNESS_NAME
    )

    witness_already_exported_only_at_future_only_strength = (
        WITNESS_NAME in n363_text
        and WITNESS_NAME in f252_text
        and "future-only feeder-support-side provider support witness" in n363_text
        and "future_only_feeder_support_side_provider_support_witness" in f252_text
    )

    current_repo_has_exported_actual_realization_of_witness = len(positive_candidates) > 0

    witness_still_remains_future_only_not_actual_export = (
        preferred_first_noncyclic_pivot_branch_already_selected
        and witness_already_exported_only_at_future_only_strength
        and not current_repo_has_exported_actual_realization_of_witness
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_witness = (
        witness_still_remains_future_only_not_actual_export
    )

    add_check(
        "preferred_first_noncyclic_pivot_branch_already_selected",
        preferred_first_noncyclic_pivot_branch_already_selected,
        True,
        "F963 already selects the feeder-support-side provider support witness as the preferred next pivot branch.",
    )
    add_check(
        "witness_already_exported_only_at_future_only_strength",
        witness_already_exported_only_at_future_only_strength,
        True,
        "N363/F252 already export the feeder-side witness only at future-only strength.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_witness",
        current_repo_has_exported_actual_realization_of_witness,
        False,
        "No stronger actual-realization artifact for this exact feeder-side witness is exported on the current repo state.",
    )
    add_check(
        "witness_still_remains_future_only_not_actual_export",
        witness_still_remains_future_only_not_actual_export,
        True,
        "Therefore the feeder-side witness still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_witness",
        next_honest_move_is_exact_actual_realization_attempt_of_same_witness,
        True,
        "The next honest move is now one exact actual-realization attempt on the same feeder-side witness branch.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and witness_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1065",
        "status": status,
        "as_of": AS_OF,
        "witness_name": WITNESS_NAME,
        "current_repo_has_exported_actual_realization_of_witness": current_repo_has_exported_actual_realization_of_witness,
        "witness_still_remains_future_only_not_actual_export": witness_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_witness": next_honest_move_is_exact_actual_realization_attempt_of_same_witness,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "witness_name": artifact["witness_name"],
        "current_repo_has_exported_actual_realization_of_witness": artifact["current_repo_has_exported_actual_realization_of_witness"],
        "witness_still_remains_future_only_not_actual_export": artifact["witness_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_witness": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_witness"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
