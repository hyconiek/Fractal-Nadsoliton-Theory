#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_OBJECT = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f660_current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_packet_summary.json"
)

GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
F657_SUMMARY = (
    GENERATED
    / "f657_first_exported_s_sel_int_strict_core_source_object_selector_output_operator_packet_summary.json"
)
F659_OBJECT = GENERATED / "selector_reduction_operator_global_c_v1_seed_v1_promoted_strict_v1.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_2x2_matrix(v: Any) -> bool:
    return (
        isinstance(v, list)
        and len(v) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in v)
    )


def mat_mul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [sum(float(a[i][k]) * float(b[k][j]) for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def main() -> None:
    atlas = load_json(GLOBAL_ATLAS)
    f657 = load_json(F657_SUMMARY)
    r_global = load_json(F659_OBJECT)

    o_sel = (f657.get("selector_output_operator", {}) or {}).get("matrix")
    if not is_2x2_matrix(o_sel):
        raise ValueError("missing 2x2 selector output operator matrix in F657 summary")
    o_sel = [[float(o_sel[i][j]) for j in range(2)] for i in range(2)]

    charts = list((r_global.get("charts", {}) or {}).keys())
    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    if sorted(charts) != sorted(expected_charts):
        raise ValueError(f"unexpected charts in global R_sel object: {charts}")

    chart_payload: dict[str, Any] = {}
    for chart in expected_charts:
        entry = (r_global.get("charts", {}) or {}).get(chart) or {}
        r = entry.get("R_sel_matrix_in_c_s_to_q")
        if not is_2x2_matrix(r):
            raise ValueError(f"missing R_sel_matrix_in_c_s_to_q for {chart} in global R_sel object")
        r = [[float(r[i][j]) for j in range(2)] for i in range(2)]
        y = mat_mul(o_sel, r)
        chart_payload[chart] = {
            "chart_id": chart,
            "domain_basis": [f"c{chart[-1]}", f"s{chart[-1]}"],
            "selector_sector_basis": ["q_+", "q_-"],
            "output_basis": ["o_+", "o_-"],
            "R_sel_ref": str(F659_OBJECT.relative_to(REPO)),
            "O_sel_ref": str(F657_SUMMARY.relative_to(REPO)),
            "Y_sel_matrix_in_c_s_to_o": y,
        }

    obj: dict[str, Any] = {
        "object": "SelectorOutputOperator_global_C_v1_seed_v1_promoted_strict_v1",
        "stage": "F660",
        "status": "actual_exported_global_selector_output_operator_object_and_induced_chartwise_output_channels__no_false_pass",
        "as_of": "2026-03-17",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Promote the seed-v1 local selector output operator O_sel : Q_sel_v1 -> Q_out_v1 to a global typed object "
            "and package the induced chartwise output channels Y_sel(pair_m) := O_sel ∘ R_sel(pair_m) for m=1..5, "
            "without implying emergent observer construction, strict-core selector closure, or global QW-2191 discharge."
        ),
        "domain": atlas.get("domain"),
        "inputs": {
            "seed_v1_local_O_sel_summary": str(F657_SUMMARY.relative_to(REPO)),
            "global_selector_atlas": str(GLOBAL_ATLAS.relative_to(REPO)),
            "global_selector_reduction_operator": str(F659_OBJECT.relative_to(REPO)),
            "boundary": "N512 (no operator-level groupoid promotion; overlap gluing remains section-level and sign-gauge-explicit)",
        },
        "output_operator": {
            "object": "O_sel_s_sel_int_source_object_v1",
            "domain_basis": ["q_+", "q_-"],
            "codomain_basis": ["o_+", "o_-"],
            "matrix_q_to_o": o_sel,
        },
        "charts": chart_payload,
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
        "stage": "F660",
        "lane": "current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_only",
        "goal": "export_one_global_selector_output_operator_object_and_induced_chartwise_output_channels_on_C_v1",
        "status": "F660_EXECUTED_CURRENT_STRICT_GLOBAL_SELECTOR_OUTPUT_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_PACKET_NO_FALSE_PASS",
        "exported_object": obj["object"],
        "exported_object_file": str(OUT_OBJECT.relative_to(REPO)),
        "chart_count": len(expected_charts),
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJECT)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

