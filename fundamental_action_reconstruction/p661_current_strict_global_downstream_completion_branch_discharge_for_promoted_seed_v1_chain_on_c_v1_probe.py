#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

F661_OBJECT = (
    GENERATED / "selector_downstream_completion_branch_global_c_v1_seed_v1_promoted_strict_v1.json"
)
F661_SUMMARY = (
    GENERATED
    / "f661_current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_packet_summary.json"
)

F658_OBJECT = GENERATED / "selector_bridge_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
F659_OBJECT = GENERATED / "selector_reduction_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
F660_OBJECT = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"

OUT = (
    GENERATED
    / "p661_current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_2x2_matrix(v: Any) -> bool:
    return (
        isinstance(v, list)
        and len(v) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in v)
    )


def mat_t(a: list[list[float]]) -> list[list[float]]:
    return [[float(a[j][i]) for j in range(len(a))] for i in range(len(a[0]))]


def mat_mul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [sum(float(a[i][k]) * float(b[k][j]) for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[float(a[i][j]) + float(b[i][j]) for j in range(len(a[0]))] for i in range(len(a))]


def mat_sub(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[float(a[i][j]) - float(b[i][j]) for j in range(len(a[0]))] for i in range(len(a))]


def max_abs(a: list[list[float]]) -> float:
    return max(abs(float(v)) for row in a for v in row)


def extract_so2_transport(op: dict[str, Any]) -> list[list[float]]:
    outputs = op.get("outputs") or {}
    if not isinstance(outputs, dict):
        raise ValueError("operator outputs missing or not a dict")
    for key, value in outputs.items():
        if isinstance(key, str) and key.startswith("G") and key.endswith("_so2") and is_2x2_matrix(value):
            return [[float(value[i][j]) for j in range(2)] for i in range(2)]
    alpha_candidates: list[float] = []
    for key, value in outputs.items():
        if not isinstance(key, str):
            continue
        if "alpha" in key and ("mod_pi" in key or "mod_2pi" in key):
            try:
                alpha_candidates.append(float(value))
            except Exception:
                pass
    if alpha_candidates:
        alpha = alpha_candidates[0]
        c = float(math.cos(alpha))
        s = float(math.sin(alpha))
        return [[c, -s], [s, c]]
    raise ValueError("no 2x2 SO(2) transport (G*_so2 or alpha*_mod_*) found in operator outputs")


def parse_edge_key(edge_key: str) -> tuple[str, str]:
    parts = edge_key.split("_to_")
    if len(parts) != 2:
        raise ValueError(f"unexpected transition key: {edge_key}")
    return parts[0], parts[1]


def main() -> None:
    tol = 1.0e-12
    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]

    f661 = load_json(F661_OBJECT)
    f661_summary = load_json(F661_SUMMARY)
    b_global = load_json(F658_OBJECT)
    r_global = load_json(F659_OBJECT)
    o_global = load_json(F660_OBJECT)
    transition = load_json(GLOBAL_TRANSITION)

    charts_list = f661.get("charts")
    charts_ok = isinstance(charts_list, list) and sorted(charts_list) == sorted(expected_charts)

    # Load per-chart matrices.
    b_by_pair: dict[str, list[list[float]]] = {}
    r_by_pair: dict[str, list[list[float]]] = {}
    y_by_pair: dict[str, list[list[float]]] = {}

    for chart in expected_charts:
        b = ((b_global.get("charts", {}) or {}).get(chart) or {}).get("B_sel_matrix_in_c_s")
        r = ((r_global.get("charts", {}) or {}).get(chart) or {}).get("R_sel_matrix_in_c_s_to_q")
        y = ((o_global.get("charts", {}) or {}).get(chart) or {}).get("Y_sel_matrix_in_c_s_to_o")
        if is_2x2_matrix(b):
            b_by_pair[chart] = [[float(b[i][j]) for j in range(2)] for i in range(2)]
        if is_2x2_matrix(r):
            r_by_pair[chart] = [[float(r[i][j]) for j in range(2)] for i in range(2)]
        if is_2x2_matrix(y):
            y_by_pair[chart] = [[float(y[i][j]) for j in range(2)] for i in range(2)]

    o_sel = ((o_global.get("output_operator", {}) or {}).get("matrix_q_to_o")) or None
    if not is_2x2_matrix(o_sel):
        o_sel_m: list[list[float]] | None = None
    else:
        o_sel_m = [[float(o_sel[i][j]) for j in range(2)] for i in range(2)]

    # End-to-end checks per chart:
    # - B = R^T D R
    # - Y = O R
    d = [[1.0, 0.0], [0.0, -1.0]]
    per_chart_checks = []
    per_chart_fail = []
    for chart in expected_charts:
        if chart not in b_by_pair or chart not in r_by_pair or chart not in y_by_pair or o_sel_m is None:
            per_chart_fail.append(chart)
            continue
        b = b_by_pair[chart]
        r = r_by_pair[chart]
        y = y_by_pair[chart]

        b_from_r = mat_mul(mat_mul(mat_t(r), d), r)
        b_diff = max_abs(mat_sub(b_from_r, b))

        y_expected = mat_mul(o_sel_m, r)
        y_diff = max_abs(mat_sub(y_expected, y))

        ok = b_diff <= 5 * tol and y_diff <= 5 * tol
        per_chart_checks.append(
            {
                "chart": chart,
                "alignment_max_abs_diff_B_vs_RtDR": b_diff,
                "alignment_max_abs_diff_Y_vs_O_R": y_diff,
                "pass": ok,
            }
        )
        if not ok:
            per_chart_fail.append(chart)

    # Overlap transport consistency across all exported transition edges:
    # - B_dst ?= G B_src G^T
    # - R_dst ?= ± R_src G^T
    # - Y_dst ?= ± Y_src G^T
    ops = transition.get("transition_operators") or {}
    overlap_checks = []
    overlap_fail = []
    for edge_key, info in ops.items():
        if not isinstance(edge_key, str) or not isinstance(info, dict):
            continue
        src, dst = parse_edge_key(edge_key)
        op_ref = info.get("operator_ref")
        if not isinstance(op_ref, str):
            overlap_fail.append(edge_key)
            continue
        if src not in b_by_pair or dst not in b_by_pair or src not in r_by_pair or dst not in r_by_pair:
            overlap_fail.append(edge_key)
            continue
        if src not in y_by_pair or dst not in y_by_pair:
            overlap_fail.append(edge_key)
            continue

        op_json = load_json(REPO / op_ref)
        g = extract_so2_transport(op_json)
        gt = mat_t(g)

        b_lhs = b_by_pair[dst]
        b_rhs = mat_mul(mat_mul(g, b_by_pair[src]), gt)
        b_diff = max_abs(mat_sub(b_lhs, b_rhs))

        r_lhs = r_by_pair[dst]
        r_rhs = mat_mul(r_by_pair[src], gt)
        r_diff_pos = max_abs(mat_sub(r_lhs, r_rhs))
        r_diff_neg = max_abs(mat_add(r_lhs, r_rhs))
        r_best = min(r_diff_pos, r_diff_neg)

        y_lhs = y_by_pair[dst]
        y_rhs = mat_mul(y_by_pair[src], gt)
        y_diff_pos = max_abs(mat_sub(y_lhs, y_rhs))
        y_diff_neg = max_abs(mat_add(y_lhs, y_rhs))
        y_best = min(y_diff_pos, y_diff_neg)

        ok = b_diff <= 5 * tol and r_best <= 5 * tol and y_best <= 5 * tol
        overlap_checks.append(
            {
                "edge": edge_key,
                "operator_ref": op_ref,
                "B_max_abs_diff": b_diff,
                "R_max_abs_diff_best_of_plus_minus": r_best,
                "Y_max_abs_diff_best_of_plus_minus": y_best,
                "pass": ok,
            }
        )
        if not ok:
            overlap_fail.append(edge_key)

    no_false_pass_ok = (
        bool(f661_summary.get("no_false_pass"))
        and not bool(f661_summary.get("strict_core_selector_closure"))
        and not bool(f661_summary.get("full_closure_pass"))
    )

    checks = [
        {
            "id": "charts_present_pair1_to_pair5",
            "actual": charts_list,
            "expected": expected_charts,
            "pass": charts_ok,
        },
        {
            "id": "per_chart_chain_alignment_B_and_Y",
            "actual": per_chart_checks,
            "expected": "all pass",
            "pass": len(per_chart_fail) == 0,
        },
        {
            "id": "overlap_transport_consistency_across_transition_edges",
            "actual": overlap_checks,
            "expected": "all pass",
            "pass": len(overlap_fail) == 0,
        },
        {
            "id": "no_false_pass_flags",
            "actual": no_false_pass_ok,
            "expected": True,
            "pass": no_false_pass_ok,
        },
    ]

    blocking = [c["id"] for c in checks if not c["pass"]]

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_GLOBAL_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_PROMOTED_SEED_V1_CHAIN_ON_C_V1_AFTER_P661"
    )

    if blocking:
        summary = {
            "stage": "P661",
            "lane": "current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_only",
            "status": "P661_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_GLOBAL_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_STATE",
            "checks": checks,
            "blocking_mismatches": blocking,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P661",
            "lane": "current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_only",
            "status": expected_status,
            "exported_object": f661.get("object"),
            "exported_object_ref": str(F661_OBJECT.relative_to(REPO)),
            "checks": checks,
            "projector_section_level_only": True,
            "residual_sign_gauge_explicit": True,
            "strict_core_selector_closure": False,
            "QW2191_discharge": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

