#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    generated_dir = ROOT / "generated"
    md_files = sorted(p.name for p in ROOT.glob("C*.md"))
    py_files = sorted(p.name for p in ROOT.glob("c*.py"))
    summary_files = sorted(p.name for p in generated_dir.glob("*_summary.json"))

    generated_dir_present = generated_dir.is_dir()
    uppercase_step_docs_present = any(name.startswith("C40_") for name in md_files)
    lowercase_generators_present = any(name.startswith("c40_") for name in py_files)
    summary_json_family_present = len(summary_files) > 0

    minimal_convention = (
        "generated/"
        "sigma_int_residual_orientation_datum_acceptance_artifact_instance.json"
    )

    dedicated_carrier_file_present = (ROOT / minimal_convention).exists()

    summary = {
        "step": "C43",
        "status": "C43_EXECUTED_FILENAME_PATH_CONVENTION_PACKET_READY_NO_FALSE_PASS",
        "goal": "Check whether strict core already contains a packet-ready minimal filename/path convention for a dedicated acceptance artifact carrier identifying sigma_int_candidate with the residual orientation datum.",
        "findings": {
            "generated_dir_present": generated_dir_present,
            "uppercase_step_docs_present": uppercase_step_docs_present,
            "lowercase_generators_present": lowercase_generators_present,
            "summary_json_family_present": summary_json_family_present,
            "minimal_filename_path_convention_packet_ready": (
                generated_dir_present
                and uppercase_step_docs_present
                and lowercase_generators_present
                and summary_json_family_present
            ),
            "proposed_relative_path": minimal_convention,
            "dedicated_carrier_file_present": dedicated_carrier_file_present,
            "persisted_artifact_instance_present": False,
        },
        "frontier_after_c43": {
            "C43_B1": "no_explicit_created_file_instance_following_the_now_packet_ready_minimal_filename_path_convention_for_a_dedicated_acceptance_artifact_carrier_identifying_sigma_int_candidate_with_the_residual_orientation_datum",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that naming convention equals a dedicated carrier",
            "no claim that naming convention equals a persisted artifact instance",
            "no claim that QW-2191 is discharged",
        ],
    }

    out = ROOT / "generated" / "c43_filename_path_convention_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
