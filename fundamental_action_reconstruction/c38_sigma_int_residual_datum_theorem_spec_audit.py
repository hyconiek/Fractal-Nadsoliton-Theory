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
        "B7": "B7_FACTORIZED_SELECTOR_MODE_SCAFFOLD_COMPATIBILITY_AUDIT.md",
        "B8": "B8_SELECTOR_TRACK_ANTI_OVERCLAIM_AUDIT.md",
        "C36": "C36_AXIOM_BRANCH_TO_STRICT_TRACK_BRIDGE_AUDIT.md",
        "C37": "C37_RESIDUAL_ORIENTATION_DATUM_INTERNALIZATION_CANDIDATE_AUDIT.md",
    }

    texts = {key: read_text(path) for key, path in files.items()}

    candidate_fit_present = (
        "sigma_int_candidate  ~  residual orientation_sign_convention slot" in texts["C37"]
        or "fits residual `Z2` orientation slot" in texts["B6"]
    )

    overlay_lane_present = contains_any(
        texts["B7"] + "\n" + texts["C36"],
        ["overlay compatibility", "control-route overlay"],
    )

    theorem_spec_needles = [
        "theorem-spec gate",
        "theorem spec gate",
        "minimal lemma dag",
        "acceptance matrix",
        "assumption map",
        "target theorem",
    ]
    export_spec_needles = [
        "export-spec",
        "export spec",
        "export theorem",
        "export obligations",
        "attachment spec",
    ]
    target_needles = [
        "sigma_int_candidate",
        "residual orientation",
        "orientation_sign_convention",
    ]

    theorem_spec_present = any(
        contains_any(text, theorem_spec_needles) and contains_any(text, target_needles)
        for text in texts.values()
    )
    export_spec_present = any(
        contains_any(text, export_spec_needles) and contains_any(text, target_needles)
        for text in texts.values()
    )

    summary = {
        "step": "C38",
        "status": "C38_EXECUTED_SIGMA_INT_RESIDUAL_DATUM_THEOREM_SPEC_AUDIT_NO_FALSE_PASS",
        "goal": "Check whether strict core already contains a packet-ready theorem-spec or export-spec for identifying sigma_int_candidate with the residual orientation datum, beyond candidate-fit on the overlay lane.",
        "inputs": files,
        "findings": {
            "candidate_fit_present": candidate_fit_present,
            "overlay_lane_present": overlay_lane_present,
            "strict_core_theorem_spec_present": theorem_spec_present,
            "strict_core_export_spec_present": export_spec_present,
        },
        "frontier_after_c38": {
            "C38_B1": "no_packet_ready_strict_core_theorem_spec_or_export_spec_for_identifying_sigma_int_candidate_with_the_residual_orientation_datum; only_candidate_fit_on_overlay_lane_exists",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that candidate-fit is already a theorem-spec",
            "no claim that overlay compatibility is already a strict-core bridge",
            "no claim that QW-2191 is discharged",
        ],
    }

    out = ROOT / "generated" / "c38_sigma_int_residual_datum_theorem_spec_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
