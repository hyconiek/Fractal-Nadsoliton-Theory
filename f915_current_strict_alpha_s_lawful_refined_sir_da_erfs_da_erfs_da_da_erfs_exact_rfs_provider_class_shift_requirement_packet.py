#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P915 = GENERATED / "p915_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_rfs_class_support_reuse_as_same_lane_exact_rfs_source_audit_probe.json"
IN_F914 = GENERATED / "f914_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_continuation_boundary_packet.json"

OUT = GENERATED / "f915_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_exact_rfs_provider_class_shift_requirement_packet.json"
OUT_SUMMARY = GENERATED / "f915_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_exact_rfs_provider_class_shift_requirement_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P915, IN_F914]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F915",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p915 = load_json(IN_P915)
    f914 = load_json(IN_F914)

    if (
        p915.get("status")
        == "P915_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_REQUIRED_FORM_STATEMENT_CLASS_SUPPORT_REUSE_AS_SAME_LANE_EXACT_REQUIRED_FORM_STATEMENT_SOURCE_NEGATIVE_PROVIDER_SHIFT_REQUIRED"
        and (p915.get("theorem_result") or {}).get(
            "current_repo_exports_no_genuinely_new_same_lane_exact_required_form_statement_source_for_alpha_s"
        )
        is True
        and (p915.get("theorem_result") or {}).get(
            "next_honest_move_requires_required_form_statement_provider_class_shift"
        )
        is True
        and f914.get("status")
        == "F914_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    ):
        status = "F915_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
    else:
        status = "F915_REQUIRES_REVIEW"

    rejected_candidate = p915.get("rejected_same_lane_source_candidate_class") or {}
    continuation_boundary = f914.get("exported_object") or {}

    artifact = {
        "stage": "F915",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedSirDaErfsDaErfsDaDaErfsExactRfsProviderClassShiftRequirement_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p915_same_lane_source_audit_probe": rel(IN_P915),
            "f914_continuation_boundary_packet": rel(IN_F914),
        },
        "exported_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_provider_class_shift_requirement_v1",
            "goal": "Export the current-repo-state requirement that deeper lawful refined alpha_s continuation now proceed by provider-class shift after same-lane neighboring statement/form reuse fails on the current exact-required-form-statement domain-admission lane.",
            "continuation_boundary_ref": continuation_boundary.get("object_id"),
            "rejected_same_lane_source_candidate_id": rejected_candidate.get("candidate_id"),
            "remaining_admitted_move_class": "shift_to_a_different_required_form_statement_provider_class_lane",
            "candidate_shift_lane_hint": "different_required_form_statement_provider_lane_not_yet_current_lawful_refined_deeper_domain_admission_exact_required_form_statement_audited",
            "forbidden_move_clause": "no_same_lane_verbal_promotion_of_neighboring_statement_or_form_scaffolding_into_current_lawful_refined_deeper_alpha_s_domain_admission_exact_required_form_statement_source",
            "scope": "strict_alpha_s_current_lawful_refined_deeper_domain_admission_exact_required_form_statement_provider_class_shift_requirement_only",
        },
        "current_honest_reading": [
            "The repo now exports an explicit provider-class-shift requirement for the current deeper lawful refined domain-admission exact required-form statement lane.",
            "Same-lane reuse of neighboring statement/form scaffolding has been audited and rejected as a current exact required-form statement source.",
            "The remaining honest continuation class is a provider shift, not local rhetorical promotion of neighboring supports.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether a different provider-class lane can serve as an admissible shift candidate for the current deeper lawful refined alpha_s exact required-form-statement domain-admission lane without silent domain identification.",
        "hard_limits": [
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim that neighboring statement/form scaffolding already supplies alpha_s semantics.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F915",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "remaining_admitted_move_class": artifact["exported_object"]["remaining_admitted_move_class"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
