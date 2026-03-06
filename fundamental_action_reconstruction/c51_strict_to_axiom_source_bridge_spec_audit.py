from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "c51_strict_to_axiom_source_bridge_spec_audit_summary.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    c35 = load(ROOT / "generated" / "c35_actual_phase_source_branch_audit_summary.json")
    c36 = load(ROOT / "generated" / "c36_axiom_branch_to_strict_track_bridge_audit_summary.json")
    c50 = load(ROOT / "generated" / "c50_actual_phase_source_skeleton_audit_summary.json")

    summary = {
        "step": "C51",
        "status": "C51_EXECUTED_STRICT_TO_AXIOM_SOURCE_BRIDGE_SPEC_AUDIT_NO_FALSE_PASS",
        "goal": "Reduce C50_B1 by checking whether the repo already contains a packet-ready strict-to-axiom source bridge specification for actual theta_1, theta_2, or whether only a fallback branch citation to QW-2192/QW-2193 is available.",
        "sources": {
            "C35": "actual phase source branch exists only on axiom-augmented lane",
            "C36": "bridge to selector track exists only as control-route overlay",
            "C50": c50["frontier_after_C50"]["C50_B1"],
            "QW_2192_QW_2193": "axiom-augmented fallback lane for theta_i*=0 mod 2pi",
            "A10": "anti-overclaim boundary"
        },
        "findings": {
            "strict_core_actual_phase_source_skeleton": c50["findings"]["strict_core_minimal_source_skeleton"],
            "axiom_augmented_fallback_lane": c50["findings"]["axiom_augmented_source_branch"],
            "strict_to_axiom_source_bridge_spec": "not_shown",
            "selector_track_overlay_bridge": c36["result"]["bridge_from_axiom_branch_to_selector_track"]
        },
        "frontier_after_C51": {
            "C51_B1": "no_packet_ready_strict_to_axiom_source_bridge_spec_for_reducing_C50_B1; only_fallback_branch_citation_to_QW_2192_QW_2193_is_available",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12"
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_qw_2192_qw_2193_become_strict_core_source",
            "no_claim_that_c50_b1_is_discharged",
            "no_claim_that_qw_2191_is_discharged"
        ],
        "next_step": "C52"
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
