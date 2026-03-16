#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

F660_OBJECT = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
F660_SUMMARY = (
    GENERATED
    / "f660_current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_packet_summary.json"
)
F657_SUMMARY = (
    GENERATED
    / "f657_first_exported_s_sel_int_strict_core_source_object_selector_output_operator_packet_summary.json"
)
F659_OBJECT = GENERATED / "selector_reduction_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"

OUT = (
    GENERATED
    / "p660_current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_probe_summary.json"
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
    f660 = load_json(F660_OBJECT)
    f660_summary = load_json(F660_SUMMARY)
    f657 = load_json(F657_SUMMARY)
    r_global = load_json(F659_OBJECT)
    transition = load_json(GLOBAL_TRANSITION)

    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    charts = list((f660.get("charts", {}) or {}).keys())
    charts_ok = sorted(charts) == sorted(expected_charts)

    # Seed alignment: O_sel matrix matches F657.
    o_seed = (f657.get("selector_output_operator", {}) or {}).get("matrix")
    o_global = (f660.get("output_operator", {}) or {}).get("matrix_q_to_o")
    seed_o_ok = False
    seed_o_diff = None
    if is_2x2_matrix(o_seed) and is_2x2_matrix(o_global):
        o_seed_m = [[float(o_seed[i][j]) for j in range(2)] for i in range(2)]
        o_global_m = [[float(o_global[i][j]) for j in range(2)] for i in range(2)]
        seed_o_diff = max_abs(mat_sub(o_global_m, o_seed_m))
        seed_o_ok = seed_o_diff <= 5 * tol
    else:
        o_seed_m = None
        o_global_m = None

    # Load global R_sel per chart.
    r_by_pair: dict[str, list[list[float]]] = {}
    for chart in expected_charts:
        entry = (r_global.get("charts", {}) or {}).get(chart) or {}
        r = entry.get("R_sel_matrix_in_c_s_to_q")
        if not is_2x2_matrix(r):
            continue
        r_by_pair[chart] = [[float(r[i][j]) for j in range(2)] for i in range(2)]

    # Per chart: check Y_sel = O_sel ∘ R_sel.
    chart_checks = []
    chart_fail = []
    for chart in expected_charts:
        entry = (f660.get("charts", {}) or {}).get(chart) or {}
        y = entry.get("Y_sel_matrix_in_c_s_to_o")
        if not is_2x2_matrix(y) or chart not in r_by_pair or o_global_m is None:
            chart_fail.append(chart)
            continue
        y_m = [[float(y[i][j]) for j in range(2)] for i in range(2)]
        y_expected = mat_mul(o_global_m, r_by_pair[chart])
        diff = max_abs(mat_sub(y_m, y_expected))
        chart_checks.append({"chart": chart, "max_abs_diff": diff, "pass": diff <= 5 * tol})
        if diff > 5 * tol:
            chart_fail.append(chart)

    # Overlap consistency across all transition edges, up to residual sign gauge:
    # Y_dst ?= ± Y_src * G^T
    ops = transition.get("transition_operators") or {}
    overlap_checks = []
    overlap_fail = []
    y_by_pair: dict[str, list[list[float]]] = {}
    for chart in expected_charts:
        entry = (f660.get("charts", {}) or {}).get(chart) or {}
        y = entry.get("Y_sel_matrix_in_c_s_to_o")
        if is_2x2_matrix(y):
            y_by_pair[chart] = [[float(y[i][j]) for j in range(2)] for i in range(2)]

    for edge_key, info in ops.items():
        if not isinstance(edge_key, str) or not isinstance(info, dict):
            continue
        src, dst = parse_edge_key(edge_key)
        op_ref = info.get("operator_ref")
        if not isinstance(op_ref, str):
            continue
        if src not in y_by_pair or dst not in y_by_pair:
            overlap_fail.append(edge_key)
            continue
        op_json = load_json(REPO / op_ref)
        g = extract_so2_transport(op_json)
        lhs = y_by_pair[dst]
        rhs = mat_mul(y_by_pair[src], mat_t(g))
        diff_pos = max_abs(mat_sub(lhs, rhs))
        diff_neg = max_abs(mat_add(lhs, rhs))
        best = min(diff_pos, diff_neg)
        overlap_checks.append(
            {
                "edge": edge_key,
                "operator_ref": op_ref,
                "max_abs_diff_best_of_plus_minus": best,
                "pass": best <= 5 * tol,
            }
        )
        if best > 5 * tol:
            overlap_fail.append(edge_key)

    checks = [
        {
            "id": "charts_present_pair1_to_pair5",
            "actual": charts,
            "expected": expected_charts,
            "pass": charts_ok,
        },
        {
            "id": "seed_alignment_O_sel",
            "actual": {"max_abs_diff": seed_o_diff, "pass": seed_o_ok},
            "expected": "<= tolerance",
            "pass": seed_o_ok,
        },
        {
            "id": "per_chart_channel_construction",
            "actual": chart_checks,
            "expected": "all pass",
            "pass": len(chart_fail) == 0,
        },
        {
            "id": "overlap_transport_consistency_up_to_residual_sign_gauge",
            "actual": overlap_checks,
            "expected": "all pass",
            "pass": len(overlap_fail) == 0,
        },
        {
            "id": "no_false_pass_flags",
            "actual": bool(f660_summary.get("no_false_pass")) and not bool(f660_summary.get("strict_core_selector_closure")) and not bool(f660_summary.get("emergent_observer_constructed")),
            "expected": True,
            "pass": bool(f660_summary.get("no_false_pass")) and not bool(f660_summary.get("strict_core_selector_closure")) and not bool(f660_summary.get("emergent_observer_constructed")),
        },
    ]

    blocking = []
    for item in checks:
        if not item["pass"]:
            blocking.append(item["id"])

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_GLOBAL_SELECTOR_OUTPUT_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_AFTER_P660"
    )

    if blocking:
        summary = {
            "stage": "P660",
            "lane": "current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_only",
            "status": "P660_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_GLOBAL_SELECTOR_OUTPUT_OPERATOR_PROMOTION_STATE",
            "checks": checks,
            "blocking_mismatches": blocking,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P660",
            "lane": "current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_only",
            "status": expected_status,
            "exported_object": f660.get("object"),
            "exported_object_ref": str(F660_OBJECT.relative_to(REPO)),
            "checks": checks,
            "projector_section_level_only": True,
            "residual_sign_gauge_explicit": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

