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

IN_P749 = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"
IN_F322_MD = ROOT / "F322_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_SUPPORT_REFINEMENT_CANDIDATE_PACKET.md"
IN_F323_MD = ROOT / "F323_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_SUPPORT_REFINEMENT_CANDIDATE_PACKET.md"
IN_N433_MD = ROOT / "N433_CURRENT_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_SUPPORT_REFINEMENT_CANDIDATE_THEOREM.md"
IN_N434_MD = ROOT / "N434_CURRENT_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_SUPPORT_REFINEMENT_CANDIDATE_THEOREM.md"
IN_T26_MD = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p750_current_strict_t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p750_current_strict_t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_nonexport_audit_probe_summary.json"

THETA_SUPPORT_OBJECT = "ThetaPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1"
PAIR_POPULATION_SUPPORT_OBJECT = (
    "BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1"
)
PAIR_INDEXED_POPULATION_ANCHOR_TARGET = "pair_indexed_population_anchor_target_v1"
CURRENT_THEOREM_FILE = (
    "N746_CURRENT_STRICT_T204_STRICT_SOURCE_SHANNON_SUPPORT_REFINEMENT_TO_MINIMAL_DESIGNATED_PAIR12_THETA_POPULATION_ENTRY_NONEXPORT_BOUNDARY_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_entry_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "f*.py", "n*.py", "t*.py")
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path.name == CURRENT_THEOREM_FILE or path in seen:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            mentions_support = (
                THETA_SUPPORT_OBJECT in text or PAIR_POPULATION_SUPPORT_OBJECT in text
            )
            mentions_pair12 = "[pair1, pair2]" in text or ("pair1" in text and "pair2" in text)
            mentions_entry = (
                "theta_1" in text
                or "theta_2" in text
                or "u_1" in text
                or "u_2" in text
                or "populated_instance" in text
            )
            if mentions_support and mentions_pair12 and mentions_entry:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P749, IN_F322_MD, IN_F323_MD, IN_N433_MD, IN_N434_MD, IN_T26_MD]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P750",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p749 = load_json(IN_P749)
    f322_text = IN_F322_MD.read_text(encoding="utf-8")
    f323_text = IN_F323_MD.read_text(encoding="utf-8")
    n433_text = IN_N433_MD.read_text(encoding="utf-8")
    n434_text = IN_N434_MD.read_text(encoding="utf-8")
    t26_text = IN_T26_MD.read_text(encoding="utf-8")
    positive_entry_candidates = scan_positive_entry_candidates()

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

    p749_pair_indexed_population_anchor_entry_bridge_already_nonexported = (
        bool(p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported"))
        and bool(
            p749.get(
                "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor"
            )
        )
        and not bool(
            p749.get(
                "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry"
            )
        )
        and not bool(p749.get("t203_target_exported_on_current_repo_state"))
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

    current_strict_source_shannon_pair_population_support_refinement_candidate_exported = (
        PAIR_POPULATION_SUPPORT_OBJECT in f323_text
        and "actual packaged strict-source Shannon-weighted pair-population support refinement candidate"
        in f323_text
        and PAIR_POPULATION_SUPPORT_OBJECT in n434_text
    )

    current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only = (
        current_strict_source_shannon_theta_support_refinement_candidate_exported
        and "pair_indexed_slot_status = present_via_R1" in f322_text
        and "conditional_population_status = present_via_C48_C49" in f322_text
        and "theta_1" in f322_text
        and "theta_2" in f322_text
        and "no actual theta values are exported" in f322_text
        and "no actual population of `u_1`, `u_2` is exported" in f322_text
        and PAIR_INDEXED_POPULATION_ANCHOR_TARGET not in f322_text
        and "[pair1, pair2]" not in f322_text
        and "pair1" not in f322_text
        and "pair2" not in f322_text
        and "tau_src_candidate_v1" not in f322_text
    )

    current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only = (
        current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and "pair_indexed_slot_status = present_via_R1" in f323_text
        and "conditional_population_status = present_via_C48_C49" in f323_text
        and "u_1: cos(theta_1^{cand,sh,sup,strict})c_1 + sin(theta_1^{cand,sh,sup,strict})s_1"
        in f323_text
        and "u_2: cos(theta_2^{cand,sh,sup,strict})c_2 + sin(theta_2^{cand,sh,sup,strict})s_2"
        in f323_text
        and "no actual populated basis-pair instance is exported" in f323_text
        and PAIR_INDEXED_POPULATION_ANCHOR_TARGET not in f323_text
        and "[pair1, pair2]" not in f323_text
        and "pair1" not in f323_text
        and "pair2" not in f323_text
        and "tau_src_candidate_v1" not in f323_text
    )

    current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry = (
        p749_pair_indexed_population_anchor_entry_bridge_already_nonexported
        and current_minimal_designated_pair12_component2_family_frozen
        and current_strict_source_shannon_theta_support_refinement_candidate_exported
        and current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and len(positive_entry_candidates) > 0
        and not current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only
        and not current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only
    )

    current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry = (
        p749_pair_indexed_population_anchor_entry_bridge_already_nonexported
        and current_minimal_designated_pair12_component2_family_frozen
        and current_strict_source_shannon_theta_support_refinement_candidate_exported
        and current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only
        and current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only
        and len(positive_entry_candidates) == 0
        and not current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry
    )

    add_check(
        "p749_pair_indexed_population_anchor_entry_bridge_already_nonexported",
        {
            "strongest_support_refinement_exported": bool(
                p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
            ),
            "anchor_entry_nonentering": bool(
                p749.get(
                    "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor"
                )
            ),
            "anchor_entry_exported": bool(
                p749.get(
                    "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry"
                )
            ),
            "t203_exported": bool(p749.get("t203_target_exported_on_current_repo_state")),
        },
        {
            "strongest_support_refinement_exported": True,
            "anchor_entry_nonentering": True,
            "anchor_entry_exported": False,
            "t203_exported": False,
        },
        "P749 already freezes that even the strongest current strict-source Shannon support-refinement route still does not enter T26 component 2.",
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
        "F322/N433 already export one real strict-source Shannon theta-support refinement candidate on this route.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported",
        {
            "pair_population_support_object_present_in_f323": PAIR_POPULATION_SUPPORT_OBJECT in f323_text,
            "pair_population_support_phrase_present": (
                "actual packaged strict-source Shannon-weighted pair-population support refinement candidate"
                in f323_text
            ),
            "pair_population_support_object_present_in_n434": PAIR_POPULATION_SUPPORT_OBJECT in n434_text,
        },
        {
            "pair_population_support_object_present_in_f323": True,
            "pair_population_support_phrase_present": True,
            "pair_population_support_object_present_in_n434": True,
        },
        "F323/N434 already export one real strict-source Shannon pair-population support-refinement candidate on this route.",
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
            "pair12_designated_family_mentioned": "[pair1, pair2]" in f322_text,
            "pair1_or_pair2_mentioned": ("pair1" in f322_text or "pair2" in f322_text),
            "anchor_target_mentioned": PAIR_INDEXED_POPULATION_ANCHOR_TARGET in f322_text,
            "tau_src_candidate_v1_mentioned": "tau_src_candidate_v1" in f322_text,
        },
        {
            "pair_indexed_slot_status_present": True,
            "conditional_population_status_present": True,
            "theta_slot_values_mentioned": True,
            "actual_theta_values_exported": False,
            "actual_u1_u2_population_exported": False,
            "pair12_designated_family_mentioned": False,
            "pair1_or_pair2_mentioned": False,
            "anchor_target_mentioned": False,
            "tau_src_candidate_v1_mentioned": False,
        },
        "The strongest strict-source Shannon theta-support refinement lane still carries only generic pair-indexed candidate theta-slot values, not one actual [pair1,pair2] entry packet.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only",
        {
            "pair_indexed_slot_status_present": "pair_indexed_slot_status = present_via_R1" in f323_text,
            "conditional_population_status_present": "conditional_population_status = present_via_C48_C49"
            in f323_text,
            "u1_u2_population_syntax_present": (
                "u_1: cos(theta_1^{cand,sh,sup,strict})c_1 + sin(theta_1^{cand,sh,sup,strict})s_1"
                in f323_text
                and "u_2: cos(theta_2^{cand,sh,sup,strict})c_2 + sin(theta_2^{cand,sh,sup,strict})s_2"
                in f323_text
            ),
            "actual_populated_basis_pair_exported": "no actual populated basis-pair instance is exported"
            not in f323_text,
            "pair12_designated_family_mentioned": "[pair1, pair2]" in f323_text,
            "pair1_or_pair2_mentioned": ("pair1" in f323_text or "pair2" in f323_text),
            "anchor_target_mentioned": PAIR_INDEXED_POPULATION_ANCHOR_TARGET in f323_text,
            "tau_src_candidate_v1_mentioned": "tau_src_candidate_v1" in f323_text,
        },
        {
            "pair_indexed_slot_status_present": True,
            "conditional_population_status_present": True,
            "u1_u2_population_syntax_present": True,
            "actual_populated_basis_pair_exported": False,
            "pair12_designated_family_mentioned": False,
            "pair1_or_pair2_mentioned": False,
            "anchor_target_mentioned": False,
            "tau_src_candidate_v1_mentioned": False,
        },
        "The strongest strict-source Shannon pair-population support-refinement lane still carries only generic pair-indexed populated-instance syntax, not one actual [pair1,pair2] entry packet.",
    )
    add_check(
        "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry",
        current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry,
        True,
        "So the strongest current strict-source Shannon support-refinement route still remains nonentering for one minimal designated [pair1,pair2] theta/population entry packet on current exports.",
    )
    add_check(
        "current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry",
        current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry,
        False,
        "No current export turns the strongest strict-source Shannon support-refinement route into one actual minimal designated [pair1,pair2] theta/population entry packet.",
    )

    status = (
        "PASS_STRICT_SOURCE_SHANNON_SUPPORT_REFINEMENT_TO_MINIMAL_DESIGNATED_PAIR12_THETA_POPULATION_ENTRY_NONEXPORT_AUDITED"
        if not blocking
        and current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry
        and not current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry
        else "FAIL_STRICT_SOURCE_SHANNON_SUPPORT_REFINEMENT_TO_MINIMAL_DESIGNATED_PAIR12_THETA_POPULATION_ENTRY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P750",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t204_target_name": "StrictSourceShannonSupportRefinementToMinimalDesignatedPair12ThetaPopulationEntry_strict_v1",
            "t204_target_exported_on_current_repo_state": current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry,
            "p749_pair_indexed_population_anchor_entry_bridge_already_nonexported": p749_pair_indexed_population_anchor_entry_bridge_already_nonexported,
            "current_minimal_designated_pair12_component2_family_frozen": current_minimal_designated_pair12_component2_family_frozen,
            "current_strict_source_shannon_theta_support_refinement_candidate_exported": current_strict_source_shannon_theta_support_refinement_candidate_exported,
            "current_strict_source_shannon_pair_population_support_refinement_candidate_exported": current_strict_source_shannon_pair_population_support_refinement_candidate_exported,
            "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only": current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only,
            "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only": current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only,
            "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry": current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry,
            "current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry": current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry,
            "positive_entry_candidates": positive_entry_candidates,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P750",
        "status": status,
        "as_of": AS_OF,
        "t204_target_name": artifact["theorem_result"]["t204_target_name"],
        "t204_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t204_target_exported_on_current_repo_state"
        ],
        "p749_pair_indexed_population_anchor_entry_bridge_already_nonexported": artifact[
            "theorem_result"
        ]["p749_pair_indexed_population_anchor_entry_bridge_already_nonexported"],
        "current_minimal_designated_pair12_component2_family_frozen": artifact[
            "theorem_result"
        ]["current_minimal_designated_pair12_component2_family_frozen"],
        "current_strict_source_shannon_theta_support_refinement_candidate_exported": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_candidate_exported"],
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_pair_population_support_refinement_candidate_exported"],
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only"],
        "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only"],
        "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry"],
        "current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
