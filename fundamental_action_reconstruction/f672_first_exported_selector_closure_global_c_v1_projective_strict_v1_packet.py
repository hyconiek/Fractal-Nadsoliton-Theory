#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_OBJECT = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f672_first_exported_selector_closure_global_c_v1_projective_strict_v1_packet_summary.json"
)

GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
GLOBAL_PROJECTIVE_STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
GLOBAL_OUTPUT_OPERATOR = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"


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


def mat_t(a: list[list[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*a)]


def mat_max_abs_diff(a: list[list[float]], b: list[list[float]]) -> float:
    return max(abs(float(a[i][j]) - float(b[i][j])) for i in range(2) for j in range(2))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    atlas = load_json(GLOBAL_ATLAS)
    transition = load_json(GLOBAL_TRANSITION)
    state = load_json(GLOBAL_PROJECTIVE_STATE)
    y_global = load_json(GLOBAL_OUTPUT_OPERATOR)

    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    charts = list((y_global.get("charts", {}) or {}).keys())
    if sorted(charts) != sorted(expected_charts):
        raise ValueError(f"unexpected charts in global output operator object: {charts}")

    chart_payload: dict[str, Any] = {}
    b_ref: list[list[float]] | None = None
    max_diff_to_ref = 0.0
    tol = 1e-12

    for chart in expected_charts:
        y_entry = (y_global.get("charts", {}) or {}).get(chart) or {}
        y = y_entry.get("Y_sel_matrix_in_c_s_to_o")
        if not is_2x2_matrix(y):
            raise ValueError(f"missing 2x2 Y_sel_matrix_in_c_s_to_o for {chart} in global output operator object")
        y = [[float(y[i][j]) for j in range(2)] for i in range(2)]

        m = int(chart[-1])
        a_path = GENERATED / f"a_{m}_pair{m}_orientation_projector_operator_strict_core_v1.json"
        a_obj = load_json(a_path)
        a = (a_obj.get("data", {}) or {}).get(f"A_{m}_pair{m}_matrix_in_c{m}_s{m}")
        if not is_2x2_matrix(a):
            raise ValueError(f"missing 2x2 A_{m}(pair{m}) matrix in {a_path.name}")
        a = [[float(a[i][j]) for j in range(2)] for i in range(2)]

        b_out = mat_mul(y, mat_mul(a, mat_t(y)))
        if b_ref is None:
            b_ref = b_out
        diff = mat_max_abs_diff(b_out, b_ref)
        max_diff_to_ref = max(max_diff_to_ref, diff)

        chart_payload[chart] = {
            "chart_id": chart,
            "domain_basis": y_entry.get("domain_basis"),
            "output_basis": y_entry.get("output_basis"),
            "A_m_ref": str(a_path.relative_to(REPO)),
            "Y_sel_ref": str(GLOBAL_OUTPUT_OPERATOR.relative_to(REPO)),
            "B_out_matrix_in_o_plus_o_minus": b_out,
            "max_abs_diff_to_reference_chart_output_projector": diff,
        }

    if b_ref is None:
        raise RuntimeError("internal error: missing reference projector")

    obj: dict[str, Any] = {
        "object": "SelectorClosure_global_C_v1_projective_strict_v1",
        "stage": "F672",
        "status": "actual_exported_global_projective_selector_closure_object_on_C_v1__no_false_pass",
        "as_of": "2026-03-17",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit global selector closure object on C_v1 in projective (ray-level) scope by packaging "
            "the exported global atlas/transition and projective selector state objects with the globally promoted "
            "seed-v1 selector output operator. The closure observable is the chartwise output projector "
            "B_out(pair_m) := Y_sel(pair_m) · A_m(pair_m) · Y_sel(pair_m)^T in the (o_+,o_-) basis, with a "
            "projector-level well-definedness certificate (chart-independence within tolerance). This export does not "
            "claim strict-core selector closure, global QW-2191 discharge, or ToE closure."
        ),
        "domain": atlas.get("domain"),
        "closure_scope": {
            "level": "projective_ray_state",
            "closure_type": "global_output_projector_section_on_output_space",
            "sign_gauge": "projector/span semantics (u and -u identified at state level)",
        },
        "inputs": {
            "global_selector_atlas": str(GLOBAL_ATLAS.relative_to(REPO)),
            "global_selector_transition": str(GLOBAL_TRANSITION.relative_to(REPO)),
            "global_projective_selector_state": str(GLOBAL_PROJECTIVE_STATE.relative_to(REPO)),
            "global_selector_output_operator": str(GLOBAL_OUTPUT_OPERATOR.relative_to(REPO)),
            "boundary": "N512 (no operator-level transition groupoid promotion; closure is projector/section-level only)",
        },
        "global_atlas_object": atlas.get("object"),
        "global_transition_object": transition.get("object"),
        "global_projective_state_object": state.get("object"),
        "global_output_operator_object": y_global.get("object"),
        "output_observable": {
            "basis": ["o_+", "o_-"],
            "output_projector_matrix_in_o_plus_o_minus": b_ref,
            "meaning": "Global chart-independent output projector obtained from the promoted output channels and the exported projector-level state section.",
        },
        "charts": chart_payload,
        "well_definedness_certificate": {
            "criterion": "max_abs_diff_across_chartwise_output_projectors <= tolerance",
            "tolerance": tol,
            "max_abs_diff_to_reference_chart_output_projector": max_diff_to_ref,
            "certificate_pass": bool(max_diff_to_ref <= tol),
            "note": "Projector-level equality is sufficient for projective closure scope; this does not promote to operator-level groupoid identities (N512).",
        },
        "hard_limits": [
            "projective_closure_only (no directed/sign-sensitive physical claim)",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_strict_core_selector_closure",
            "no_global_QW2191_discharge",
            "no_emergent_observer_construction",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary: dict[str, Any] = {
        "stage": "F672",
        "lane": "first_exported_selector_closure_global_c_v1_projective_strict_v1_packet_only",
        "goal": "export_one_global_projective_selector_closure_object_and_chart_independence_certificate_on_C_v1",
        "status": "F672_EXECUTED_FIRST_EXPORTED_SELECTOR_CLOSURE_GLOBAL_C_V1_PROJECTIVE_STRICT_V1_PACKET_NO_FALSE_PASS",
        "exported_object": obj["object"],
        "exported_object_file": str(OUT_OBJECT.relative_to(REPO)),
        "chart_count": len(expected_charts),
        "certificate_pass": obj["well_definedness_certificate"]["certificate_pass"],
        "max_abs_diff_to_reference_chart_output_projector": max_diff_to_ref,
        "tolerance": tol,
        "strict_core_selector_closure": False,
        "QW2191_discharge": False,
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJECT)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

