#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f8_legacy_weinberg_role_equivalence_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED / "f8_legacy_weinberg_role_equivalence_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def any_entry_with_id(obj: Any, target_id: str) -> bool:
    if isinstance(obj, dict):
        if obj.get("id") == target_id:
            return True
        return any(any_entry_with_id(value, target_id) for value in obj.values())
    if isinstance(obj, list):
        return any(any_entry_with_id(value, target_id) for value in obj)
    return False


def role_transfer_verdict_present(text: str) -> bool:
    legacy_markers = [
        "alpha_geo/12",
        "sin^2(theta_W)=alpha_geo/12",
        "legacy weinberg",
        "legacy weinberg-angle role",
    ]
    strict_markers = ["sin2_theta_w_mz", "sin2_theta_w_eff"]
    transfer_markers = [
        "role-equivalence",
        "role equivalence",
        "retained successor",
        "same role",
        "semantic transfer",
        "successor semantics",
        "retained role",
    ]
    return (
        any(marker in text for marker in legacy_markers)
        and any(marker in text for marker in strict_markers)
        and any(marker in text for marker in transfer_markers)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p66 = load_json(
        "fundamental_action_reconstruction/generated/p66_strict_side_literal_retention_probe_for_legacy_weinberg_formula_summary.json"
    )
    release_text = load_text("RELEASE_4_9_TEXTBOOK_EN_PL.md")
    qw2068 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json"
    )
    qw2069 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
    )
    qw2098_md = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md"
    )
    qw2098_json = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json"
    )

    release_promotes_candidate = "sin2_theta_w_mz" in release_text
    registry_has_candidate = any_entry_with_id(qw2068, "sin2_theta_w_mz")
    package_has_candidate = any_entry_with_id(qw2069, "sin2_theta_w_mz")
    gate_has_candidate = any_entry_with_id(qw2098_json, "sin2_theta_w_mz") or (
        "sin2_theta_w_mz" in qw2098_md
    )
    candidate_object_present = (
        release_promotes_candidate
        and registry_has_candidate
        and package_has_candidate
        and gate_has_candidate
    )

    semantic_transfer_sources = {
        "RELEASE_4_9_TEXTBOOK_EN_PL.md": release_text,
        "report_qw2069_full_sm_gr_derivation_package.json": json.dumps(
            qw2069, ensure_ascii=True
        ),
        "RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md": qw2098_md,
        "report_qw2098_ew_secondary_nonanchor_closure_gate.json": json.dumps(
            qw2098_json, ensure_ascii=True
        ),
    }
    per_source_semantic_transfer = {
        source: role_transfer_verdict_present(text)
        for source, text in semantic_transfer_sources.items()
    }
    explicit_semantic_transfer_verdict_present = any(
        per_source_semantic_transfer.values()
    )

    checks_spec = [
        {
            "id": "p66_literal_retention_path_closed",
            "actual": p66["literal_retention_present"],
            "expected": False,
            "meaning": "P66 already closed literal retention negatively on the current repo state",
        },
        {
            "id": "strict_side_candidate_object_present",
            "actual": candidate_object_present,
            "expected": True,
            "meaning": "the strict side exports a real Weinberg candidate object via sin2_theta_w_mz",
        },
        {
            "id": "explicit_semantic_transfer_verdict_present",
            "actual": explicit_semantic_transfer_verdict_present,
            "expected": False,
            "meaning": "the current repo still exports no explicit legacy-to-strict semantic-transfer verdict for that candidate object",
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
        "stage": "F8",
        "lane": "legacy_weinberg_role_equivalence_refinement_current_repo_state_only",
        "goal": "refine_the_remaining_retained_side_weinberg_role_equivalence_frontier_into_candidate_object_presence_versus_explicit_legacy_to_strict_semantic_transfer_verdict",
        "status": "F8_EXECUTED_LEGACY_WEINBERG_ROLE_EQUIVALENCE_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "the strict side already exports sin2_theta_w_mz as a real candidate object, but that is still weaker than an explicit verdict that this object retains the legacy Weinberg-angle role",
        "role_equivalence_state": {
            "strict_side_candidate_object_id": "sin2_theta_w_mz",
            "strict_side_candidate_object_present": candidate_object_present,
            "explicit_semantic_transfer_verdict_present": explicit_semantic_transfer_verdict_present,
        },
        "per_source_semantic_transfer_verdict_presence": per_source_semantic_transfer,
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_semantic_transfer_verdict_identifying_sin2_theta_w_mz_as_the_retained_strict_side_successor_of_the_legacy_weinberg_angle_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F8",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "role_equivalence_state": artifact["role_equivalence_state"],
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
