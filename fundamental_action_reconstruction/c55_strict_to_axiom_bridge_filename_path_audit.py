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
    uppercase_step_docs_present = any(name.startswith("C54_") for name in md_files)
    lowercase_generators_present = any(name.startswith("c54_") for name in py_files)
    summary_json_family_present = len(summary_files) > 0

    minimal_convention = (
        "generated/"
        "strict_to_axiom_sigma_int_residual_orientation_datum_bridge_artifact_instance.json"
    )

    dedicated_bridge_carrier_file_present = (ROOT / minimal_convention).exists()

    summary = {
        "step": "C55",
        "status": "C55_EXECUTED_STRICT_TO_AXIOM_BRIDGE_FILENAME_PATH_PACKET_READY_NO_FALSE_PASS",
        "goal": "Check whether strict core already contains a packet-ready minimal filename/path convention for a dedicated strict-to-axiom bridge carrier reducing C50_B1.",
        "findings": {
            "generated_dir_present": generated_dir_present,
            "uppercase_step_docs_present": generated_dir_present and uppercase_step_docs_present,
            "lowercase_generators_present": generated_dir_present and lowercase_generators_present,
            "summary_json_family_present": summary_json_family_present,
            "minimal_filename_path_convention_packet_ready": (
                generated_dir_present
                and uppercase_step_docs_present
                and lowercase_generators_present
                and summary_json_family_present
            ),
            "proposed_relative_path": minimal_convention,
            "dedicated_bridge_carrier_file_present": dedicated_bridge_carrier_file_present,
            "persisted_bridge_artifact_instance_present": False,
        },
        "frontier_after_C55": {
            "C55_B1": "no_explicit_created_file_instance_following_the_now_packet_ready_minimal_filename_path_convention_for_a_dedicated_strict_to_axiom_bridge_carrier_reducing_C50_B1",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that naming convention equals a dedicated bridge carrier",
            "no claim that naming convention equals a persisted bridge artifact instance",
            "no claim that QW-2191 is discharged",
        ],
        "next_step": "C56",
    }

    out = ROOT / "generated" / "c55_strict_to_axiom_bridge_filename_path_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
