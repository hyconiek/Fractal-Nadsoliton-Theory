#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p72_legacy_weinberg_replaced_successor_subbranch_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p72_legacy_weinberg_replaced_successor_subbranch_probe_summary.json"
)

STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
    "t1_nonanchor_observables_input_qw2085_2086.json",
]

LEGACY_MARKERS = [
    "legacy weinberg",
    "legacy weinberg-angle role",
    "alpha_geo/12",
    "sin^2(theta_W)=alpha_geo/12",
]
OBJECT_MARKERS = ["sin2_theta_w_mz"]
METHOD_MARKERS = ["qw2098_sin2_from_nonanchor_ew_pole_chain"]
REPLACEMENT_MARKERS = [
    "replaced",
    "replacement",
    "replaced by",
    "superseded",
    "supersession",
    "strict successor semantics",
    "successor semantics",
]


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def object_replaced_verdict_present(text: str) -> bool:
    return (
        any(marker in text for marker in LEGACY_MARKERS)
        and any(marker in text for marker in OBJECT_MARKERS)
        and any(marker in text for marker in REPLACEMENT_MARKERS)
    )


def method_replaced_verdict_present(text: str) -> bool:
    return (
        any(marker in text for marker in LEGACY_MARKERS)
        and any(marker in text for marker in METHOD_MARKERS)
        and any(marker in text for marker in REPLACEMENT_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f10 = load_json(
        "fundamental_action_reconstruction/generated/f10_legacy_weinberg_replaced_branch_refinement_packet_summary.json"
    )

    per_source_object = {}
    per_source_method = {}
    any_object_verdict = False
    any_method_verdict = False
    for source in STRICT_SIDE_SOURCES:
        text = load_text(source)
        object_present = object_replaced_verdict_present(text)
        method_present = method_replaced_verdict_present(text)
        per_source_object[source] = object_present
        per_source_method[source] = method_present
        any_object_verdict = any_object_verdict or object_present
        any_method_verdict = any_method_verdict or method_present

    checks_spec = [
        {
            "id": "f10_object_candidate_present",
            "actual": f10["candidate_state"]["strict_object_candidate_present"],
            "expected": True,
            "meaning": "F10 already confirms that the strict side exports sin2_theta_w_mz",
        },
        {
            "id": "f10_method_candidate_present",
            "actual": f10["candidate_state"]["strict_method_candidate_present"],
            "expected": True,
            "meaning": "F10 already confirms that the strict side exports qw2098_sin2_from_nonanchor_ew_pole_chain",
        },
        {
            "id": "object_replaced_verdict_present",
            "actual": any_object_verdict,
            "expected": False,
            "meaning": "the current repo exports no explicit object-successor replaced verdict for the legacy Weinberg role",
        },
        {
            "id": "method_replaced_verdict_present",
            "actual": any_method_verdict,
            "expected": False,
            "meaning": "the current repo exports no explicit method-successor-semantics replaced verdict for the legacy Weinberg role",
        },
    ]

    checks = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "P72",
        "lane": "legacy_weinberg_replaced_successor_subbranch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_either_the_object_successor_or_method_successor_replaced_subbranch_for_the_legacy_weinberg_role",
        "status": "CURRENT_REPO_EXPORTS_STRICT_WEINBERG_SUCCESSOR_CANDIDATES_BUT_NEITHER_OBJECT_NOR_METHOD_REPLACED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P72",
        "reason": "the strict side exports both sin2_theta_w_mz and qw2098_sin2_from_nonanchor_ew_pole_chain as real successor-candidate structures, but no current source upgrades either one into an explicit replaced verdict for the legacy Weinberg role",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "per_source_object_replaced_verdict_presence": per_source_object,
        "per_source_method_replaced_verdict_presence": per_source_method,
        "object_replaced_verdict_present": any_object_verdict,
        "method_replaced_verdict_present": any_method_verdict,
        "remaining_missing_objects": [
            "explicit_object_successor_verdict_identifying_sin2_theta_w_mz_as_the_strict_side_successor_object_replacing_the_legacy_weinberg_angle_role",
            "explicit_method_successor_semantics_verdict_identifying_qw2098_sin2_from_nonanchor_ew_pole_chain_as_the_strict_side_successor_semantics_replacing_the_legacy_weinberg_angle_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P72",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_object_replaced_verdict_presence": per_source_object,
        "per_source_method_replaced_verdict_presence": per_source_method,
        "object_replaced_verdict_present": any_object_verdict,
        "method_replaced_verdict_present": any_method_verdict,
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
