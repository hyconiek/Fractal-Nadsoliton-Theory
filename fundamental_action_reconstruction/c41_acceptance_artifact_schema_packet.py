#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def main() -> None:
    files = {
        "B6": "B6_SIGMA_TO_SELECTOR_FACTORIZED_BRIDGE.md",
        "B8": "B8_SELECTOR_TRACK_ANTI_OVERCLAIM_AUDIT.md",
        "C36": "C36_AXIOM_BRANCH_TO_STRICT_TRACK_BRIDGE_AUDIT.md",
        "C37": "C37_RESIDUAL_ORIENTATION_DATUM_INTERNALIZATION_CANDIDATE_AUDIT.md",
        "C38": "C38_SIGMA_INT_RESIDUAL_DATUM_THEOREM_SPEC_AUDIT.md",
        "C39": "C39_SIGMA_INT_ACCEPTANCE_SKELETON_AUDIT.md",
        "C40": "C40_MINIMAL_FIELD_LIST_AUDIT.md",
    }

    texts = {k: read_text(v) for k, v in files.items()}

    field_presence = {
        "object": "sigma_int_candidate" in texts["C40"],
        "target": "residual orientation datum" in texts["C40"],
        "support_lane": "candidate_fit_on_overlay_lane_only" in texts["C40"],
        "current_absence": "no theorem-spec" in texts["C40"] and "no export-spec" in texts["C40"],
        "forbidden_claims": "no theorem-level PASS" in texts["C40"] and "no QW-2191 discharge" in texts["C40"],
    }

    schema_packet_ready = all(field_presence.values())

    summary = {
        "step": "C41",
        "status": "C41_EXECUTED_ACCEPTANCE_ARTIFACT_SCHEMA_PACKET_NO_FALSE_PASS",
        "goal": "Check whether the already present minimal field list can now be assembled into a packet-ready acceptance artifact schema for identifying sigma_int_candidate with the residual orientation datum.",
        "inputs": files,
        "field_presence": field_presence,
        "schema_packet_ready": schema_packet_ready,
        "frontier_after_c41": {
            "C41_B1": "no_explicit_persisted_acceptance_artifact_instance_for_identifying_sigma_int_candidate_with_the_residual_orientation_datum_even_though_a_minimal_schema_is_now_packet_ready",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that acceptance artifact schema is a theorem-spec",
            "no claim that acceptance artifact schema is an export-spec",
            "no claim that QW-2191 is discharged",
        ],
    }

    out = ROOT / "generated" / "c41_acceptance_artifact_schema_packet_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
