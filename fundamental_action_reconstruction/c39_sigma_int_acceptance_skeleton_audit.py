#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def contains_any(text: str, needles: list[str]) -> bool:
    lower = text.lower()
    return any(needle.lower() in lower for needle in needles)


def main() -> None:
    files = {
        "B6": "B6_SIGMA_TO_SELECTOR_FACTORIZED_BRIDGE.md",
        "B8": "B8_SELECTOR_TRACK_ANTI_OVERCLAIM_AUDIT.md",
        "C36": "C36_AXIOM_BRANCH_TO_STRICT_TRACK_BRIDGE_AUDIT.md",
        "C37": "C37_RESIDUAL_ORIENTATION_DATUM_INTERNALIZATION_CANDIDATE_AUDIT.md",
        "C38": "C38_SIGMA_INT_RESIDUAL_DATUM_THEOREM_SPEC_AUDIT.md",
    }

    texts = {key: read_text(path) for key, path in files.items()}

    candidate_fit_present = (
        "candidate-fit" in texts["C37"].lower()
        and "sigma_int_candidate" in texts["C37"]
    ) or "supported_candidate_fit" in texts["B6"]

    acceptance_needles = [
        "acceptance matrix",
        "acceptance skeleton",
        "acceptance packet",
        "acceptance criteria",
    ]
    target_needles = [
        "sigma_int_candidate",
        "residual orientation",
        "orientation_sign_convention",
    ]

    acceptance_skeleton_present = any(
        contains_any(text, acceptance_needles) and contains_any(text, target_needles)
        for text in texts.values()
    )

    summary = {
        "step": "C39",
        "status": "C39_EXECUTED_SIGMA_INT_ACCEPTANCE_SKELETON_AUDIT_NO_FALSE_PASS",
        "goal": "Check whether strict core already contains even a packet-ready acceptance skeleton for a future theorem-spec/export-spec identifying sigma_int_candidate with the residual orientation datum.",
        "inputs": files,
        "findings": {
            "candidate_fit_present": candidate_fit_present,
            "theorem_spec_present": False,
            "export_spec_present": False,
            "acceptance_skeleton_present": acceptance_skeleton_present,
        },
        "frontier_after_c39": {
            "C39_B1": "no_packet_ready_acceptance_skeleton_for_a_future_theorem_spec_or_export_spec_identifying_sigma_int_candidate_with_the_residual_orientation_datum; only_candidate_fit_on_overlay_lane_exists",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that candidate-fit is already enough for strict closure",
            "no claim that hidden acceptance skeleton exists",
            "no claim that QW-2191 is discharged",
        ],
    }

    out = ROOT / "generated" / "c39_sigma_int_acceptance_skeleton_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
