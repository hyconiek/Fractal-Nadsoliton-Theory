#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P480 = GENERATED / "p480_current_strict_p16_route_negative_freeze_decision_packet_summary.json"
IN_P631 = GENERATED / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet_summary.json"
IN_P632 = GENERATED / "p632_current_strict_directed_continuation_decision_packet_summary.json"

IN_P119 = GENERATED / "p119_first_source_seed_construction_target_probe_summary.json"
IN_N676 = GENERATED / "n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json"

OUT_JSON = GENERATED / "p633_current_strict_source_seed_route_selection_decision_packet.json"
OUT_SUMMARY = GENERATED / "p633_current_strict_source_seed_route_selection_decision_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = (IN_P480, IN_P631, IN_P632, IN_P119)
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P633",
            "date": AS_OF,
            "goal": "declare_professorial_next_move_after_freezing_P16_and_direct_formal_lanes__shift_to_source_seed_frontier__no_false_pass",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_SOURCE_SEED_ROUTE_DECISION",
            "missing": missing,
            "recommended_next_strict_target": {
                "target": "P480" if not IN_P480.exists() else ("P631" if not IN_P631.exists() else ("P632" if not IN_P632.exists() else "P119")),
                "note": "This decision requires explicit freeze/continuation decision packets (P480/P631/P632) and the source-seed target probe (P119).",
            },
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": "P633",
                    "status": artifact["status"],
                    "recommended_next_strict_target": artifact["recommended_next_strict_target"]["target"],
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p480 = load_json(IN_P480)
    p631 = load_json(IN_P631)
    p632 = load_json(IN_P632)
    p119 = load_json(IN_P119)
    n676 = load_json(IN_N676) if IN_N676.exists() else {}

    p480_selected = str(p480.get("decision") or "") == "P16_ROUTE_NEGATIVE_FREEZE_SELECTED"
    p631_selected = str(p631.get("decision") or "") == "DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED"
    p632_directed = (str(p632.get("decision") or "") == "DIRECTED_CONTINUATION_SELECTED") or (
        str(p632.get("selected_continuation") or "") == "directed"
    )

    # If the prereq decisions are not selected, this packet is a no-op and should not override dashboards.
    if not (p480_selected and p631_selected and p632_directed):
        artifact = {
            "stage": "P633",
            "date": AS_OF,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "goal": "declare_professorial_next_move_after_freezing_P16_and_direct_formal_lanes__shift_to_source_seed_frontier__no_false_pass",
            "status": "NOT_APPLICABLE_PREREQUISITES_NOT_SELECTED_FOR_SOURCE_SEED_ROUTE_DECISION",
            "prerequisites": {
                "p480_selected": p480_selected,
                "p631_selected": p631_selected,
                "p632_directed_selected": p632_directed,
            },
            "recommended_next_strict_target": {
                "target": str(p632.get("recommended_next_strict_target") or "P11"),
                "note": "Source-seed route selection applies only after P480+P631 freezes and directed continuation selection; otherwise defer to the currently exported continuation decision.",
            },
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": "P633",
                    "status": artifact["status"],
                    "recommended_next_strict_target": artifact["recommended_next_strict_target"]["target"],
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    decision = "STRICT_CORE_SOURCE_SEED_ROUTE_SELECTED"
    seed_exported = bool(
        n676.get("theorem_result", {}).get("admissible_S_sel_int_source_object_in_F34_sense")
    )
    recommended_next = "T172" if seed_exported else "P119"
    meaning = (
        "Source-seed lane is now materially exported: an admissible strict-core source object for S_sel_int exists (F34 sense) and downstream "
        "operators are available in the declared scope. Next strict bottleneck shifts to global strict selector closure and QW-2191 discipline (T172), "
        "without implying any kernel-alone discharge."
        if seed_exported
        else (
            "Shift the next strict move to the genuinely-new strict-core source-seed construction frontier for S_sel_int. "
            "P119 packages the first explicit construction target; later branches (E_orient / B_sel->R_sel->O_sel) remain open."
        )
    )
    artifact = {
        "stage": "P633",
        "date": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "declare_professorial_next_move_after_freezing_P16_and_direct_formal_lanes__shift_to_source_seed_frontier__no_false_pass",
        "status": "PASS_DECISION_DECLARED_STRICT_CORE_SOURCE_SEED_ROUTE_SELECTED",
        "decision": decision,
        "decision_basis": {
            "p16_lane_frozen_negative": True,
            "p16_freeze_packet": str(IN_P480.relative_to(REPO)),
            "direct_formal_lane_frozen_negative_on_t166_nonzero_branch": True,
            "direct_formal_freeze_packet": str(IN_P631.relative_to(REPO)),
            "directed_continuation_selected": True,
            "directed_continuation_packet": str(IN_P632.relative_to(REPO)),
            "target_probe_present": True,
            "target_probe": str(IN_P119.relative_to(REPO)),
            "strict_boundary_note": "This is a professorial routing decision only. It does not assert existence of S_sel_int nor any strict-core selector closure.",
        },
        "next_move": {
            "recommended_next_strict_target": recommended_next,
            "meaning": meaning,
            "target_state_excerpt": p119.get("target_state"),
        },
        "hard_limits": [
            "No strict-core selector closure / admissible S_sel_int claim.",
            "No global discharge of QW-2191 claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(
        json.dumps(
            {
                "stage": "P633",
                "status": artifact["status"],
                "decision": decision,
                "recommended_next_strict_target": recommended_next,
                "no_false_pass": True,
            },
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="ascii",
    )

    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
