#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def main() -> None:
    files = {
        "B4": "B4_MINIMAL_SIGMA_INT_CANDIDATE.md",
        "B6": "B6_SIGMA_TO_SELECTOR_FACTORIZED_BRIDGE.md",
        "B8": "B8_SELECTOR_TRACK_ANTI_OVERCLAIM_AUDIT.md",
        "C36": "C36_AXIOM_BRANCH_TO_STRICT_TRACK_BRIDGE_AUDIT.md",
        "C37": "C37_RESIDUAL_ORIENTATION_DATUM_INTERNALIZATION_CANDIDATE_AUDIT.md",
        "C38": "C38_SIGMA_INT_RESIDUAL_DATUM_THEOREM_SPEC_AUDIT.md",
        "C39": "C39_SIGMA_INT_ACCEPTANCE_SKELETON_AUDIT.md",
    }

    texts = {key: read_text(path) for key, path in files.items()}

    field_list = {
        "candidate_object": (
            "sigma_int_candidate" in texts["B4"]
            and "sigma_int_candidate" in texts["C37"]
        ),
        "target_slot_or_target_datum": (
            "orientation_sign_convention" in texts["B6"]
            and "residual orientation datum" in texts["C37"]
        ),
        "current_support_lane": (
            "overlay lane" in texts["C37"].lower()
            and "strict-core bridge absent" in texts["C36"].lower()
        ),
        "strict_absence_claim": (
            "theorem-spec" in texts["C38"].lower()
            and "acceptance skeleton" in texts["C39"].lower()
        ),
        "forbidden_overclaim_set": (
            "anti-overclaim" in texts["B8"].lower()
            and "candidate-fit" in texts["C38"].lower()
        ),
    }

    summary = {
        "step": "C40",
        "status": "C40_EXECUTED_MINIMAL_FIELD_LIST_AUDIT_NO_FALSE_PASS",
        "goal": "Check whether strict core already contains the minimal semantic field list needed for a future acceptance skeleton identifying sigma_int_candidate with the residual orientation datum.",
        "inputs": files,
        "field_list": field_list,
        "frontier_after_c40": {
            "C40_B1": "no_explicit_assembled_acceptance_artifact_built_from_the_already_present_minimal_field_list_for_identifying_sigma_int_candidate_with_the_residual_orientation_datum",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that field-list presence equals acceptance skeleton",
            "no claim that theorem-spec or export-spec exists",
            "no claim that QW-2191 is discharged",
        ],
    }

    out = ROOT / "generated" / "c40_minimal_field_list_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
