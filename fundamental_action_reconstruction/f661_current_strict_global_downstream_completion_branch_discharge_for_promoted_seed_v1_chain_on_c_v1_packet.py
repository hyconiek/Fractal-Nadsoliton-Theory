#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
GLOBAL_STATE_PROJECTIVE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"

F658_OBJECT = GENERATED / "selector_bridge_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
F659_OBJECT = GENERATED / "selector_reduction_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
F660_OBJECT = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"

OUT_OBJECT = (
    GENERATED / "selector_downstream_completion_branch_global_c_v1_seed_v1_promoted_strict_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f661_current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    atlas = load_json(GLOBAL_ATLAS)
    transition = load_json(GLOBAL_TRANSITION)
    state = load_json(GLOBAL_STATE_PROJECTIVE)
    b_global = load_json(F658_OBJECT)
    r_global = load_json(F659_OBJECT)
    o_global = load_json(F660_OBJECT)

    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]

    obj: dict[str, Any] = {
        "object": "SelectorDownstreamCompletionBranch_global_C_v1_seed_v1_promoted_strict_v1",
        "stage": "F661",
        "status": "actual_exported_global_downstream_completion_branch_bundle_for_seed_v1_promoted_selector_operator_chain_on_C_v1__no_false_pass",
        "as_of": "2026-03-17",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit global downstream-completion branch discharge bundle for the promoted seed-v1 selector operator chain "
            "on the declared strict configuration space C_v1, by packaging references to the already exported global atlas/transition/state "
            "objects (F469/F470) and the already exported promoted operator-chain objects (F658/F659/F660), while keeping the N512 "
            "no-operator-level-groupoid-promotion boundary and residual sign-gauge discipline explicit."
        ),
        "domain": atlas.get("domain"),
        "inputs": {
            "global_selector_atlas": str(GLOBAL_ATLAS.relative_to(REPO)),
            "global_selector_transition": str(GLOBAL_TRANSITION.relative_to(REPO)),
            "global_projective_selector_state": str(GLOBAL_STATE_PROJECTIVE.relative_to(REPO)),
            "global_selector_bridge_operator": str(F658_OBJECT.relative_to(REPO)),
            "global_selector_reduction_operator": str(F659_OBJECT.relative_to(REPO)),
            "global_selector_output_operator_and_channels": str(F660_OBJECT.relative_to(REPO)),
            "boundary": "N512 (no operator-level groupoid promotion; gluing remains projector/section-level; residual sign gauge explicit where axis-only transport is used)",
        },
        "charts": expected_charts,
        "assembled_chain": {
            "selector_atlas": {"object": atlas.get("object"), "ref": str(GLOBAL_ATLAS.relative_to(REPO))},
            "selector_transition": {
                "object": transition.get("object"),
                "ref": str(GLOBAL_TRANSITION.relative_to(REPO)),
            },
            "selector_state_projective": {
                "object": state.get("object"),
                "ref": str(GLOBAL_STATE_PROJECTIVE.relative_to(REPO)),
            },
            "selector_bridge_operator": {
                "object": b_global.get("object"),
                "ref": str(F658_OBJECT.relative_to(REPO)),
            },
            "selector_reduction_operator": {
                "object": r_global.get("object"),
                "ref": str(F659_OBJECT.relative_to(REPO)),
            },
            "selector_output_operator": {
                "object": o_global.get("object"),
                "ref": str(F660_OBJECT.relative_to(REPO)),
            },
            "chartwise_output_channels": {
                "definition": "Y_sel(pair_m) := O_sel ∘ R_sel(pair_m) : pair_m -> Q_out_v1",
                "stored_in": str(F660_OBJECT.relative_to(REPO)),
            },
        },
        "hard_limits": [
            "no_admissible_S_sel_int",
            "no_strict_core_selector_closure",
            "no_global_QW2191_discharge",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_physical_sign_orientation_claim",
            "no_emergent_observer_construction",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary: dict[str, Any] = {
        "stage": "F661",
        "lane": "current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_only",
        "goal": "export_one_global_downstream_completion_branch_bundle_for_promoted_seed_v1_selector_operator_chain_on_C_v1",
        "status": "F661_EXECUTED_CURRENT_STRICT_GLOBAL_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_PROMOTED_SEED_V1_CHAIN_ON_C_V1_PACKET_NO_FALSE_PASS",
        "exported_object": obj["object"],
        "exported_object_file": str(OUT_OBJECT.relative_to(REPO)),
        "referenced_objects": {
            "selector_atlas_global": str(GLOBAL_ATLAS.relative_to(REPO)),
            "selector_transition_global": str(GLOBAL_TRANSITION.relative_to(REPO)),
            "selector_state_projective_global": str(GLOBAL_STATE_PROJECTIVE.relative_to(REPO)),
            "selector_bridge_operator_global": str(F658_OBJECT.relative_to(REPO)),
            "selector_reduction_operator_global": str(F659_OBJECT.relative_to(REPO)),
            "selector_output_operator_global": str(F660_OBJECT.relative_to(REPO)),
        },
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "global_QW2191_discharge": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_OBJECT)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

