#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p82_strict_side_role_equivalence_probe_for_legacy_fine_structure_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p82_strict_side_role_equivalence_probe_for_legacy_fine_structure_role_summary.json"
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
        "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "legacy fine-structure",
        "legacy fine-structure role",
    ]
    strict_markers = ["alpha_em_inv_mz"]
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

    f15 = load_json(
        "fundamental_action_reconstruction/generated/f15_legacy_fine_structure_role_equivalence_refinement_packet_summary.json"
    )
    qw2068 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json"
    )
    qw2069 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
    )
    qw2094 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json"
    )
    qw2098_md = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md"
    )
    qw2098_json = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json"
    )

    candidate_chain = {
        "q2068_registry_has_alpha_em_inv_mz": any_entry_with_id(
            qw2068, "alpha_em_inv_mz"
        ),
        "q2069_package_has_alpha_em_inv_mz": any_entry_with_id(
            qw2069, "alpha_em_inv_mz"
        ),
        "q2098_gate_has_alpha_em_inv_mz": any_entry_with_id(
            qw2098_json, "alpha_em_inv_mz"
        )
        or ("alpha_em_inv_mz" in qw2098_md),
        "q2094_status_matches_qw2098": bool(
            qw2094["checks"].get("QW-2069_alpha_em_inv_mz_matches_qw2098_status")
        ),
        "q2094_method_matches_qw2098": bool(
            qw2094["checks"].get("QW-2069_alpha_em_inv_mz_matches_qw2098_method")
        ),
    }
    candidate_object_present = all(candidate_chain.values())

    semantic_transfer_sources = {
        "report_qw2069_full_sm_gr_derivation_package.json": json.dumps(
            qw2069, ensure_ascii=True
        ),
        "report_qw2094_strict_rigor_defect_sweep.json": json.dumps(
            qw2094, ensure_ascii=True
        ),
        "RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md": qw2098_md,
        "report_qw2098_ew_secondary_nonanchor_closure_gate.json": json.dumps(
            qw2098_json, ensure_ascii=True
        ),
    }
    per_source_role_equivalence_verdict = {
        source: role_transfer_verdict_present(text)
        for source, text in semantic_transfer_sources.items()
    }
    explicit_role_equivalence_verdict_present = any(
        per_source_role_equivalence_verdict.values()
    )

    checks_spec = [
        {
            "id": "f15_candidate_present",
            "actual": f15["role_equivalence_state"]["strict_side_candidate_object_present"],
            "expected": True,
            "meaning": "F15 already identified a real strict-side candidate object",
        },
        {
            "id": "strict_side_candidate_object_present",
            "actual": candidate_object_present,
            "expected": True,
            "meaning": "the current strict-side sources export a real fine-structure candidate object via alpha_em_inv_mz",
        },
        {
            "id": "explicit_role_equivalence_verdict_present",
            "actual": explicit_role_equivalence_verdict_present,
            "expected": False,
            "meaning": "the current strict-side sources still export no explicit legacy-to-strict role-equivalence verdict",
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
        "stage": "P82",
        "lane": "strict_side_role_equivalence_probe_for_legacy_fine_structure_role_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_an_explicit_legacy_to_strict_role_equivalence_verdict_for_the_legacy_fine_structure_role",
        "status": "CURRENT_REPO_EXPORTS_STRICT_SIDE_ALPHA_EM_INV_MZ_CANDIDATE_BUT_NO_EXPLICIT_LEGACY_FINE_STRUCTURE_ROLE_EQUIVALENCE_VERDICT_AFTER_P82",
        "reason": "the strict side exports alpha_em_inv_mz as a real derived candidate object, but still no explicit semantic-transfer verdict identifies that object as the retained role-equivalent successor of the legacy fine-structure role",
        "strict_side_candidate_object": {
            "id": "alpha_em_inv_mz",
            "present": candidate_object_present,
            "candidate_chain": candidate_chain,
        },
        "per_source_role_equivalence_verdict_presence": per_source_role_equivalence_verdict,
        "explicit_role_equivalence_verdict_present": explicit_role_equivalence_verdict_present,
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_semantic_transfer_verdict_identifying_alpha_em_inv_mz_as_the_retained_strict_side_successor_of_the_legacy_fine_structure_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P82",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "strict_side_candidate_object": artifact["strict_side_candidate_object"],
        "explicit_role_equivalence_verdict_present": explicit_role_equivalence_verdict_present,
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
