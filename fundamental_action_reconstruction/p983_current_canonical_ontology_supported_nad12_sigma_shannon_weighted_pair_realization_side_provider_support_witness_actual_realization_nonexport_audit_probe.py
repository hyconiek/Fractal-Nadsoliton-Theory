#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_N362 = ROOT / "N362_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md"
IN_F251 = ROOT / "F251_FIRST_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_PACKET.md"

OUT_JSON = GENERATED / "p983_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p983_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_nonexport_probe_summary.json"

WITNESS_NAME = "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1"
CURRENT_THEOREM_FILE = (
    "N816_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_"
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
        "P983_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N362_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md",
        "F251_FIRST_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_PACKET.md",
        "p983_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_nonexport_probe.json",
        "p983_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_nonexport_probe_summary.json",
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

    prerequisites = [IN_F949, IN_N362, IN_F251]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P983",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f949 = load_json(IN_F949)
    n362_text = load_text(IN_N362)
    f251_text = load_text(IN_F251)
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
        f949.get("status")
        == "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
        and f949.get("preferred_first_pivot_branch") == WITNESS_NAME
    )

    witness_already_exported_only_at_future_only_strength = (
        WITNESS_NAME in n362_text
        and WITNESS_NAME in f251_text
        and "future-only pair-realization-side provider support witness" in n362_text
        and "future_only_pair_realization_side_provider_support_witness" in f251_text
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
        "F949 already selects the pair-realization-side provider support witness as the preferred first pivot branch.",
    )
    add_check(
        "witness_already_exported_only_at_future_only_strength",
        witness_already_exported_only_at_future_only_strength,
        True,
        "N362/F251 already export the witness only at future-only strength.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_witness",
        current_repo_has_exported_actual_realization_of_witness,
        False,
        "No stronger actual-realization artifact for this exact witness is exported on the current repo state.",
    )
    add_check(
        "witness_still_remains_future_only_not_actual_export",
        witness_still_remains_future_only_not_actual_export,
        True,
        "Therefore the witness still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_witness",
        next_honest_move_is_exact_actual_realization_attempt_of_same_witness,
        True,
        "The next honest move is now one exact actual-realization attempt on the same witness branch.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and witness_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P983",
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
