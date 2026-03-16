#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

F659_OBJECT = GENERATED / "selector_reduction_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
F659_SUMMARY = (
    GENERATED
    / "f659_current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_packet_summary.json"
)
GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
GLOBAL_STATE_PROJECTIVE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
F656_SUMMARY = (
    GENERATED
    / "f656_first_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_packet_summary.json"
)
F658_OBJECT = GENERATED / "selector_bridge_operator_global_c_v1_seed_v1_promoted_strict_v1.json"

OUT = (
    GENERATED
    / "p659_current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_probe_summary.json"
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


def identity2() -> list[list[float]]:
    return [[1.0, 0.0], [0.0, 1.0]]


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


def extract_projector_matrix_2x2(a_op: dict[str, Any]) -> list[list[float]]:
    data = a_op.get("data") or {}
    if not isinstance(data, dict):
        raise ValueError("A_* operator missing data dict")
    mats = []
    for key, value in data.items():
        if isinstance(key, str) and "matrix_in_c" in key and is_2x2_matrix(value):
            mats.append([[float(value[i][j]) for j in range(2)] for i in range(2)])
    if len(mats) != 1:
        raise ValueError(f"expected exactly one 2x2 projector matrix candidate, got {len(mats)}")
    return mats[0]


def parse_edge_key(edge_key: str) -> tuple[str, str]:
    parts = edge_key.split("_to_")
    if len(parts) != 2:
        raise ValueError(f"unexpected transition key: {edge_key}")
    return parts[0], parts[1]


def main() -> None:
    tol = 1.0e-12
    f659 = load_json(F659_OBJECT)
    f659_summary = load_json(F659_SUMMARY)
    transition = load_json(GLOBAL_TRANSITION)
    state = load_json(GLOBAL_STATE_PROJECTIVE)
    f656 = load_json(F656_SUMMARY)
    b_global = load_json(F658_OBJECT)

    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    charts = list((f659.get("charts", {}) or {}).keys())
    charts_ok = sorted(charts) == sorted(expected_charts)

    # Load chartwise R matrices.
    r_by_pair: dict[str, list[list[float]]] = {}
    for chart in expected_charts:
        entry = (f659.get("charts", {}) or {}).get(chart) or {}
        r = entry.get("R_sel_matrix_in_c_s_to_q")
        if not is_2x2_matrix(r):
            continue
        r_by_pair[chart] = [[float(r[i][j]) for j in range(2)] for i in range(2)]

    # Load chartwise B matrices from the already exported global bridge operator object.
    b_by_pair: dict[str, list[list[float]]] = {}
    for chart in expected_charts:
        entry = (b_global.get("charts", {}) or {}).get(chart) or {}
        b = entry.get("B_sel_matrix_in_c_s")
        if not is_2x2_matrix(b):
            continue
        b_by_pair[chart] = [[float(b[i][j]) for j in range(2)] for i in range(2)]

    i2 = identity2()
    d = [[1.0, 0.0], [0.0, -1.0]]
    p_q_plus = [[1.0, 0.0], [0.0, 0.0]]

    # Per-chart checks: orthogonality and alignment vs B_sel and A_m.
    chart_checks = []
    chart_fail = []
    for chart in expected_charts:
        r = r_by_pair.get(chart)
        b = b_by_pair.get(chart)
        if r is None or b is None:
            chart_fail.append(chart)
            continue
        orth = max_abs(mat_sub(mat_mul(r, mat_t(r)), i2))
        b_from_r = mat_mul(mat_mul(mat_t(r), d), r)
        b_diff = max_abs(mat_sub(b_from_r, b))

        a_ref = state.get("charts", {}).get(chart, {}).get("local_operator", {}).get("ref")
        if not isinstance(a_ref, str):
            chart_fail.append(chart)
            continue
        a_op = load_json(REPO / a_ref)
        a_mat = extract_projector_matrix_2x2(a_op)
        p_plus = mat_mul(mat_mul(mat_t(r), p_q_plus), r)
        a_diff = max_abs(mat_sub(p_plus, a_mat))

        ok = orth <= 5 * tol and b_diff <= 5 * tol and a_diff <= 5 * tol
        chart_checks.append(
            {
                "chart": chart,
                "R_orthogonality_inf_norm": orth,
                "alignment_max_abs_diff_B_vs_RtDR": b_diff,
                "alignment_max_abs_diff_P_plus_vs_A": a_diff,
                "pass": ok,
            }
        )
        if not ok:
            chart_fail.append(chart)

    # Overlap transport consistency across all exported transition edges:
    # R_dst ?= ± R_src * G^T (residual sign gauge explicit on axis-only edges).
    ops = transition.get("transition_operators") or {}
    overlap_checks = []
    overlap_fail = []
    for edge_key, info in ops.items():
        if not isinstance(edge_key, str) or not isinstance(info, dict):
            continue
        src, dst = parse_edge_key(edge_key)
        op_ref = info.get("operator_ref")
        if not isinstance(op_ref, str):
            continue
        op_json = load_json(REPO / op_ref)
        g = extract_so2_transport(op_json)
        if src not in r_by_pair or dst not in r_by_pair:
            overlap_fail.append(edge_key)
            continue
        lhs = r_by_pair[dst]
        rhs = mat_mul(r_by_pair[src], mat_t(g))
        diff_pos = max_abs(mat_sub(lhs, rhs))
        diff_neg = max_abs(mat_add(lhs, rhs))
        best = min(diff_pos, diff_neg)
        overlap_checks.append(
            {
                "edge": edge_key,
                "operator_ref": op_ref,
                "max_abs_diff_best_of_plus_minus": best,
                "max_abs_diff_pos": diff_pos,
                "max_abs_diff_neg": diff_neg,
                "pass": best <= 5 * tol,
            }
        )
        if best > 5 * tol:
            overlap_fail.append(edge_key)

    # Seed alignment on pair1.
    r1_seed = f656.get("selector_reduction_operator", {}).get("matrix_in_c1_s1")
    seed_alignment_ok = False
    seed_alignment_diff = None
    if is_2x2_matrix(r1_seed) and "pair1" in r_by_pair:
        r_seed = [[float(r1_seed[i][j]) for j in range(2)] for i in range(2)]
        seed_alignment_diff = max_abs(mat_sub(r_by_pair["pair1"], r_seed))
        seed_alignment_ok = seed_alignment_diff <= 5 * tol

    checks = [
        {
            "id": "charts_present_pair1_to_pair5",
            "actual": charts,
            "expected": expected_charts,
            "pass": charts_ok,
        },
        {
            "id": "per_chart_orthogonality_and_alignment",
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
            "id": "seed_pair1_alignment",
            "actual": {
                "max_abs_diff": seed_alignment_diff,
                "pass": seed_alignment_ok,
            },
            "expected": "<= tolerance",
            "pass": seed_alignment_ok,
        },
        {
            "id": "no_false_pass_flags",
            "actual": bool(f659_summary.get("no_false_pass")) and not bool(f659_summary.get("strict_core_selector_closure")),
            "expected": True,
            "pass": bool(f659_summary.get("no_false_pass")) and not bool(f659_summary.get("strict_core_selector_closure")),
        },
    ]

    blocking = []
    for item in checks:
        if not item["pass"]:
            blocking.append(item["id"])

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_GLOBAL_SELECTOR_REDUCTION_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_AFTER_P659"
    )

    if blocking:
        summary = {
            "stage": "P659",
            "lane": "current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_only",
            "status": "P659_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_GLOBAL_SELECTOR_REDUCTION_OPERATOR_PROMOTION_STATE",
            "checks": checks,
            "blocking_mismatches": blocking,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P659",
            "lane": "current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_only",
            "status": expected_status,
            "exported_object": f659.get("object"),
            "exported_object_ref": str(F659_OBJECT.relative_to(REPO)),
            "checks": checks,
            "projector_section_level_only": True,
            "residual_sign_gauge_explicit": True,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

