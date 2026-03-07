#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p73_strict_side_object_successor_subbranch_probe_for_legacy_weinberg_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p73_strict_side_object_successor_subbranch_probe_for_legacy_weinberg_role_summary.json"
)

STRICT_SIDE_OBJECT_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
]

LEGACY_MARKERS = [
    "legacy weinberg",
    "legacy weinberg-angle role",
    "alpha_geo/12",
    "sin^2(theta_W)=alpha_geo/12",
]
OBJECT_MARKERS = ["sin2_theta_w_mz"]
TEXTUAL_SUCCESSOR_MARKERS = [
    "successor object",
    "strict-side successor object",
    "replacing the legacy weinberg-angle role",
    "replacing the legacy weinberg role",
    "replaced by",
    "replacement semantics",
]
LINEAGE_MARKERS = [
    "strict-derived",
    "target parameter",
    "qw2098_sin2_from_nonanchor_ew_pole_chain",
    "matches_qw2098_method",
    "matches_qw2098_status",
]
UPGRADE_MARKERS = [
    "upgraded to replacement",
    "replacement semantics",
    "successor object",
    "replaced by",
    "legacy-to-strict",
]


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def textual_object_successor_verdict_present(text: str) -> bool:
    return (
        any(marker in text for marker in LEGACY_MARKERS)
        and any(marker in text for marker in OBJECT_MARKERS)
        and any(marker in text for marker in TEXTUAL_SUCCESSOR_MARKERS)
    )


def object_lineage_upgrade_verdict_present(text: str) -> bool:
    return (
        any(marker in text for marker in LEGACY_MARKERS)
        and any(marker in text for marker in OBJECT_MARKERS)
        and any(marker in text for marker in LINEAGE_MARKERS)
        and any(marker in text for marker in UPGRADE_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f11 = load_json(
        "fundamental_action_reconstruction/generated/f11_legacy_weinberg_object_successor_refinement_packet_summary.json"
    )

    per_source_textual = {}
    per_source_lineage = {}
    any_textual = False
    any_lineage = False
    for source in STRICT_SIDE_OBJECT_SOURCES:
        text = load_text(source)
        textual = textual_object_successor_verdict_present(text)
        lineage = object_lineage_upgrade_verdict_present(text)
        per_source_textual[source] = textual
        per_source_lineage[source] = lineage
        any_textual = any_textual or textual
        any_lineage = any_lineage or lineage

    checks_spec = [
        {
            "id": "f11_object_chain_present",
            "actual": f11["candidate_state"]["object_chain_present"],
            "expected": True,
            "meaning": "F11 already confirms that the repo exports the sin2_theta_w_mz object chain",
        },
        {
            "id": "f11_object_chain_consistency_present",
            "actual": f11["candidate_state"]["object_chain_consistency_present"],
            "expected": True,
            "meaning": "F11 already confirms that QW-2094 supplies object-chain consistency markers",
        },
        {
            "id": "textual_object_successor_verdict_present",
            "actual": any_textual,
            "expected": False,
            "meaning": "the current repo exports no explicit textual object-successor verdict for the legacy Weinberg role",
        },
        {
            "id": "object_lineage_upgrade_verdict_present",
            "actual": any_lineage,
            "expected": False,
            "meaning": "the current repo exports no explicit object-lineage-upgrade verdict for the legacy Weinberg role",
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
        "stage": "P73",
        "lane": "strict_side_object_successor_subbranch_probe_for_legacy_weinberg_role_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_either_the_textual_object_successor_or_object_lineage_upgrade_subbranch_for_the_legacy_weinberg_role",
        "status": "CURRENT_REPO_EXPORTS_THE_STRICT_SIN2_THETA_W_MZ_OBJECT_CHAIN_BUT_NEITHER_TEXTUAL_OBJECT_SUCCESSOR_NOR_OBJECT_LINEAGE_UPGRADE_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P73",
        "reason": "the repo exports the sin2_theta_w_mz object chain and object-side consistency markers, but no current source upgrades that chain into replacement semantics for the legacy Weinberg role",
        "strict_side_sources_checked": STRICT_SIDE_OBJECT_SOURCES,
        "per_source_textual_object_successor_verdict_presence": per_source_textual,
        "per_source_object_lineage_upgrade_verdict_presence": per_source_lineage,
        "textual_object_successor_verdict_present": any_textual,
        "object_lineage_upgrade_verdict_present": any_lineage,
        "remaining_missing_objects": [
            "explicit_textual_object_successor_verdict_identifying_sin2_theta_w_mz_as_the_strict_side_successor_object_replacing_the_legacy_weinberg_angle_role",
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_sin2_theta_w_mz_candidate_chain_into_replacement_semantics_for_the_legacy_weinberg_angle_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P73",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_textual_object_successor_verdict_presence": per_source_textual,
        "per_source_object_lineage_upgrade_verdict_presence": per_source_lineage,
        "textual_object_successor_verdict_present": any_textual,
        "object_lineage_upgrade_verdict_present": any_lineage,
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
