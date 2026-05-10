#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P750 = GENERATED / "p750_current_strict_t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_nonexport_audit_probe_summary.json"
IN_F322_MD = ROOT / "F322_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_SUPPORT_REFINEMENT_CANDIDATE_PACKET.md"
IN_N433_MD = ROOT / "N433_CURRENT_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_SUPPORT_REFINEMENT_CANDIDATE_THEOREM.md"
IN_T26_MD = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p751_current_strict_t205_strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p751_current_strict_t205_strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_nonexport_audit_probe_summary.json"

THETA_SUPPORT_OBJECT = "ThetaPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1"
PAIR_INDEXED_POPULATION_ANCHOR_TARGET = "pair_indexed_population_anchor_target_v1"
CURRENT_THEOREM_FILE = (
    "N747_CURRENT_STRICT_T205_STRICT_SOURCE_SHANNON_THETA_SUPPORT_REFINEMENT_TO_MINIMAL_DESIGNATED_PAIR12_THETA_ENTRY_NONEXPORT_BOUNDARY_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_theta_entry_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "f*.py", "n*.py", "t*.py")
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path.name == CURRENT_THEOREM_FILE or path in seen:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if (
                THETA_SUPPORT_OBJECT in text
                and "[pair1, pair2]" in text
                and "theta_1" in text
                and "theta_2" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P750, IN_F322_MD, IN_N433_MD, IN_T26_MD]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P751",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p750 = load_json(IN_P750)
    f322_text = IN_F322_MD.read_text(encoding="utf-8")
    n433_text = IN_N433_MD.read_text(encoding="utf-8")
    t26_text = IN_T26_MD.read_text(encoding="utf-8")
    positive_theta_entry_candidates = scan_positive_theta_entry_candidates()

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

    p750_combined_theta_population_entry_bridge_already_nonexported = (
        bool(p750.get("current_strict_source_shannon_theta_support_refinement_candidate_exported"))
        and bool(p750.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported"))
        and bool(
            p750.get(
                "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry"
            )
        )
        and not bool(
            p750.get("current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry")
        )
        and not bool(p750.get("t204_target_exported_on_current_repo_state"))
    )

    current_minimal_designated_pair12_component2_family_frozen = (
        PAIR_INDEXED_POPULATION_ANCHOR_TARGET in t26_text
        and "strict_internal_orientation_provider_target_v1" in t26_text
        and "[pair1, pair2]" in t26_text
    )

    current_strict_source_shannon_theta_support_refinement_candidate_exported = (
        THETA_SUPPORT_OBJECT in f322_text
        and "actual packaged strict-source Shannon-weighted theta-export support refinement candidate"
        in f322_text
        and THETA_SUPPORT_OBJECT in n433_text
    )

    current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only = (
        current_strict_source_shannon_theta_support_refinement_candidate_exported
        and "pair_indexed_slot_status = present_via_R1" in f322_text
        and "conditional_population_status = present_via_C48_C49" in f322_text
        and "theta_1" in f322_text
        and "theta_2" in f322_text
        and "no actual theta values are exported" in f322_text
        and "no actual population of `u_1`, `u_2` is exported" in f322_text
        and "[pair1, pair2]" not in f322_text
        and "pair1" not in f322_text
        and "pair2" not in f322_text
        and PAIR_INDEXED_POPULATION_ANCHOR_TARGET not in f322_text
        and "tau_src_candidate_v1" not in f322_text
    )

    current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export = (
        current_strict_source_shannon_theta_support_refinement_candidate_exported
        and "still below actual theta export" in f322_text
        and "no actual theta values are exported" in f322_text
    )

    current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry = (
        p750_combined_theta_population_entry_bridge_already_nonexported
        and current_minimal_designated_pair12_component2_family_frozen
        and current_strict_source_shannon_theta_support_refinement_candidate_exported
        and len(positive_theta_entry_candidates) > 0
        and not current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only
        and not current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export
    )

    current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry = (
        p750_combined_theta_population_entry_bridge_already_nonexported
        and current_minimal_designated_pair12_component2_family_frozen
        and current_strict_source_shannon_theta_support_refinement_candidate_exported
        and current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only
        and current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export
        and len(positive_theta_entry_candidates) == 0
        and not current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry
    )

    add_check(
        "p750_combined_theta_population_entry_bridge_already_nonexported",
        {
            "theta_support_refinement_exported": bool(
                p750.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")
            ),
            "pair_population_support_refinement_exported": bool(
                p750.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
            ),
            "combined_theta_population_entry_nonentering": bool(
                p750.get(
                    "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry"
                )
            ),
            "combined_theta_population_entry_exported": bool(
                p750.get(
                    "current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry"
                )
            ),
            "t204_exported": bool(p750.get("t204_target_exported_on_current_repo_state")),
        },
        {
            "theta_support_refinement_exported": True,
            "pair_population_support_refinement_exported": True,
            "combined_theta_population_entry_nonentering": True,
            "combined_theta_population_entry_exported": False,
            "t204_exported": False,
        },
        "P750 already freezes that the combined minimal designated [pair1,pair2] theta/population entry bridge is still absent on current exports.",
    )
    add_check(
        "current_minimal_designated_pair12_component2_family_frozen",
        {
            "pair_indexed_population_anchor_target_present": PAIR_INDEXED_POPULATION_ANCHOR_TARGET in t26_text,
            "strict_internal_orientation_provider_target_present": (
                "strict_internal_orientation_provider_target_v1" in t26_text
            ),
            "minimal_designated_pair_family_pair12_present": "[pair1, pair2]" in t26_text,
        },
        {
            "pair_indexed_population_anchor_target_present": True,
            "strict_internal_orientation_provider_target_present": True,
            "minimal_designated_pair_family_pair12_present": True,
        },
        "T26 already freezes component 2 on at least the minimal designated pair family [pair1, pair2].",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_candidate_exported",
        {
            "theta_support_object_present_in_f322": THETA_SUPPORT_OBJECT in f322_text,
            "theta_support_phrase_present": (
                "actual packaged strict-source Shannon-weighted theta-export support refinement candidate"
                in f322_text
            ),
            "theta_support_object_present_in_n433": THETA_SUPPORT_OBJECT in n433_text,
        },
        {
            "theta_support_object_present_in_f322": True,
            "theta_support_phrase_present": True,
            "theta_support_object_present_in_n433": True,
        },
        "F322/N433 already export one real strict-source Shannon theta-support refinement candidate.",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only",
        {
            "pair_indexed_slot_status_present": "pair_indexed_slot_status = present_via_R1" in f322_text,
            "conditional_population_status_present": "conditional_population_status = present_via_C48_C49"
            in f322_text,
            "theta_slot_values_mentioned": ("theta_1" in f322_text and "theta_2" in f322_text),
            "actual_theta_values_exported": "no actual theta values are exported" not in f322_text,
            "actual_u1_u2_population_exported": "no actual population of `u_1`, `u_2` is exported"
            not in f322_text,
            "minimal_designated_pair_family_pair12_mentioned": "[pair1, pair2]" in f322_text,
            "pair1_or_pair2_mentioned": ("pair1" in f322_text or "pair2" in f322_text),
            "pair_indexed_population_anchor_target_mentioned": PAIR_INDEXED_POPULATION_ANCHOR_TARGET in f322_text,
            "tau_src_candidate_v1_mentioned": "tau_src_candidate_v1" in f322_text,
        },
        {
            "pair_indexed_slot_status_present": True,
            "conditional_population_status_present": True,
            "theta_slot_values_mentioned": True,
            "actual_theta_values_exported": False,
            "actual_u1_u2_population_exported": False,
            "minimal_designated_pair_family_pair12_mentioned": False,
            "pair1_or_pair2_mentioned": False,
            "pair_indexed_population_anchor_target_mentioned": False,
            "tau_src_candidate_v1_mentioned": False,
        },
        "The strongest strict-source Shannon theta-support refinement lane still carries only generic pair-indexed candidate theta-slot values, not one actual minimal designated [pair1,pair2] theta-entry packet.",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export",
        current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export,
        True,
        "That strongest strict-source Shannon theta-support refinement lane still remains explicitly below actual theta export on current exports.",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry",
        current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry,
        True,
        "So the strongest strict-source Shannon theta-support refinement route still remains nonentering for one actual minimal designated [pair1,pair2] theta-entry packet on current exports.",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry",
        current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry,
        False,
        "No current export turns the strongest strict-source Shannon theta-support refinement route into one actual minimal designated [pair1,pair2] theta-entry packet.",
    )

    status = (
        "PASS_STRICT_SOURCE_SHANNON_THETA_SUPPORT_REFINEMENT_TO_MINIMAL_DESIGNATED_PAIR12_THETA_ENTRY_NONEXPORT_AUDITED"
        if not blocking
        and current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry
        and not current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry
        else "FAIL_STRICT_SOURCE_SHANNON_THETA_SUPPORT_REFINEMENT_TO_MINIMAL_DESIGNATED_PAIR12_THETA_ENTRY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P751",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t205_target_name": "StrictSourceShannonThetaSupportRefinementToMinimalDesignatedPair12ThetaEntry_strict_v1",
            "t205_target_exported_on_current_repo_state": current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry,
            "p750_combined_theta_population_entry_bridge_already_nonexported": p750_combined_theta_population_entry_bridge_already_nonexported,
            "current_minimal_designated_pair12_component2_family_frozen": current_minimal_designated_pair12_component2_family_frozen,
            "current_strict_source_shannon_theta_support_refinement_candidate_exported": current_strict_source_shannon_theta_support_refinement_candidate_exported,
            "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only": current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only,
            "current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export": current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export,
            "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry": current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry,
            "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry": current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry,
            "positive_theta_entry_candidates": positive_theta_entry_candidates,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P751",
        "status": status,
        "as_of": AS_OF,
        "t205_target_name": artifact["theorem_result"]["t205_target_name"],
        "t205_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t205_target_exported_on_current_repo_state"
        ],
        "p750_combined_theta_population_entry_bridge_already_nonexported": artifact[
            "theorem_result"
        ]["p750_combined_theta_population_entry_bridge_already_nonexported"],
        "current_minimal_designated_pair12_component2_family_frozen": artifact[
            "theorem_result"
        ]["current_minimal_designated_pair12_component2_family_frozen"],
        "current_strict_source_shannon_theta_support_refinement_candidate_exported": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_candidate_exported"],
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only"],
        "current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export"],
        "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry"],
        "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
