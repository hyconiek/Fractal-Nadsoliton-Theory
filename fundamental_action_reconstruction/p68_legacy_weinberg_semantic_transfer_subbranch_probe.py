#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p68_legacy_weinberg_semantic_transfer_subbranch_probe.json"
)
OUT_SUMMARY = (
    GENERATED / "p68_legacy_weinberg_semantic_transfer_subbranch_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def textual_successor_verdict_present(text: str) -> bool:
    return (
        "sin2_theta_w_mz" in text
        and any(
            marker in text
            for marker in [
                "retained successor",
                "retained strict-side successor",
                "same legacy Weinberg role",
                "legacy Weinberg-angle role",
                "successor semantics",
            ]
        )
    )


def lineage_upgrade_verdict_present(text: str) -> bool:
    return (
        "alpha_geo" in text
        and ("sin2_theta_w_mz" in text or "sin2_theta_w_eff" in text)
        and any(
            marker in text
            for marker in [
                "semantic transfer",
                "role-equivalence",
                "role equivalence",
                "retained role",
                "upgraded to retained",
                "legacy-to-strict",
            ]
        )
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f9 = load_json(
        "fundamental_action_reconstruction/generated/f9_legacy_weinberg_semantic_transfer_refinement_packet_summary.json"
    )
    source_map = {
        "RELEASE_4_9_TEXTBOOK_EN_PL.md": load_text("RELEASE_4_9_TEXTBOOK_EN_PL.md"),
        "report_qw2069_full_sm_gr_derivation_package.json": (
            REPO
            / "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
        ).read_text(encoding="utf-8"),
        "RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md": load_text(
            "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md"
        ),
        "report_qw2098_ew_secondary_nonanchor_closure_gate.json": (
            REPO
            / "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json"
        ).read_text(encoding="utf-8"),
        "t1_nonanchor_observables_input_qw2085_2086.json": (
            REPO / "t1_nonanchor_observables_input_qw2085_2086.json"
        ).read_text(encoding="utf-8"),
    }

    per_source_textual_successor = {
        source: textual_successor_verdict_present(text)
        for source, text in source_map.items()
    }
    per_source_lineage_upgrade = {
        source: lineage_upgrade_verdict_present(text)
        for source, text in source_map.items()
    }
    textual_successor_present = any(per_source_textual_successor.values())
    lineage_upgrade_present = any(per_source_lineage_upgrade.values())

    checks_spec = [
        {
            "id": "f9_candidate_present",
            "actual": f9["semantic_transfer_state"]["strict_side_candidate_object_present"],
            "expected": True,
            "meaning": "F9 already confirms the strict-side candidate object exists",
        },
        {
            "id": "f9_lineage_touchpoint_present",
            "actual": f9["semantic_transfer_state"]["qw2093_alpha_geo_lineage_touchpoint_present"],
            "expected": True,
            "meaning": "F9 already confirms the QW-2093 output-side alpha_geo lineage touchpoint exists",
        },
        {
            "id": "textual_successor_verdict_present",
            "actual": textual_successor_present,
            "expected": False,
            "meaning": "no explicit textual retained-successor verdict is currently exported",
        },
        {
            "id": "lineage_upgrade_verdict_present",
            "actual": lineage_upgrade_present,
            "expected": False,
            "meaning": "no explicit lineage-upgrade verdict is currently exported",
        },
    ]

    checks: list[dict[str, Any]] = []
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
        "stage": "P68",
        "lane": "legacy_weinberg_semantic_transfer_subbranch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_either_the_textual_successor_or_the_lineage_upgrade_subbranch_of_retained_weinberg_semantic_transfer",
        "status": "CURRENT_REPO_EXPORTS_NEITHER_TEXTUAL_SUCCESSOR_NOR_LINEAGE_UPGRADE_TRANSFER_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P68",
        "reason": "the strict-side candidate object and the QW-2093 alpha_geo lineage touchpoint are both present, but neither one is upgraded into an explicit retained semantic-transfer verdict",
        "per_source_textual_successor_verdict_presence": per_source_textual_successor,
        "per_source_lineage_upgrade_verdict_presence": per_source_lineage_upgrade,
        "textual_successor_verdict_present": textual_successor_present,
        "lineage_upgrade_verdict_present": lineage_upgrade_present,
        "remaining_missing_objects": [
            "explicit_textual_retained_successor_verdict_identifying_sin2_theta_w_mz_as_the_retained_strict_side_successor_of_the_legacy_weinberg_angle_role",
            "explicit_lineage_upgrade_verdict_elevating_the_qw2093_alpha_geo_touchpoint_into_retained_strict_side_weinberg_role_transfer",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P68",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "textual_successor_verdict_present": textual_successor_present,
        "lineage_upgrade_verdict_present": lineage_upgrade_present,
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
