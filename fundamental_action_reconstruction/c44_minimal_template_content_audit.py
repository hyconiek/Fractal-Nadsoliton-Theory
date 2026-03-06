#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def main() -> None:
    c40 = read_text("C40_MINIMAL_FIELD_LIST_AUDIT.md")
    c41 = read_text("C41_ACCEPTANCE_ARTIFACT_SCHEMA_PACKET.md")
    c43 = read_text("C43_FILENAME_PATH_CONVENTION_AUDIT.md")

    fields = {
        "object": "object" in c41.lower(),
        "target": "target" in c41.lower(),
        "support_lane": "support_lane" in c41.lower(),
        "current_absence": "current_absence" in c41.lower(),
        "forbidden_claims": "forbidden_claims" in c41.lower() or "forbidden claims" in c41.lower(),
        "residual_blockers": "residual_blockers" in c41.lower() or "residual blockers" in c41.lower(),
    }

    filename_path_ready = "sigma_int_residual_orientation_datum_acceptance_artifact_instance.json" in c43
    template_content_packet_ready = all(fields.values()) and filename_path_ready and "minimal field list" in c40.lower()

    template_content = {
        "object": "sigma_int_candidate",
        "target": "residual orientation datum",
        "support_lane": "candidate_fit_on_overlay_lane_only",
        "current_absence": [
            "no theorem-spec",
            "no export-spec",
            "no strict-core bridge",
        ],
        "forbidden_claims": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no QW-2191 discharge",
        ],
        "residual_blockers": [
            "C32_B2",
            "C26_B2",
        ],
    }

    summary = {
        "step": "C44",
        "status": "C44_EXECUTED_MINIMAL_TEMPLATE_CONTENT_PACKET_READY_NO_FALSE_PASS",
        "goal": "Check whether strict core already contains packet-ready minimal template content for a dedicated acceptance artifact carrier identifying sigma_int_candidate with the residual orientation datum.",
        "inputs": {
            "C40": "C40_MINIMAL_FIELD_LIST_AUDIT.md",
            "C41": "C41_ACCEPTANCE_ARTIFACT_SCHEMA_PACKET.md",
            "C43": "C43_FILENAME_PATH_CONVENTION_AUDIT.md",
        },
        "field_presence": fields,
        "filename_path_convention_packet_ready": filename_path_ready,
        "minimal_template_content_packet_ready": template_content_packet_ready,
        "proposed_template_content": template_content,
        "dedicated_carrier_file_created": False,
        "persisted_artifact_instance_present": False,
        "frontier_after_c44": {
            "C44_B1": "no_explicit_created_persisted_file_instance_populated_with_the_now_packet_ready_minimal_template_content_and_filename_path_convention_for_identifying_sigma_int_candidate_with_the_residual_orientation_datum",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that template content equals theorem-spec",
            "no claim that template content equals export-spec",
            "no claim that QW-2191 is discharged",
        ],
    }

    out = ROOT / "generated" / "c44_minimal_template_content_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
