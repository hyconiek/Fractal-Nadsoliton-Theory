#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_OBJECT = GENERATED / "selector_reduction_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f659_current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_packet_summary.json"
)

GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
GLOBAL_STATE_PROJECTIVE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
F656_SUMMARY = (
    GENERATED
    / "f656_first_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_packet_summary.json"
)
F658_OBJECT = GENERATED / "selector_bridge_operator_global_c_v1_seed_v1_promoted_strict_v1.json"


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


def mat_scale(alpha: float, a: list[list[float]]) -> list[list[float]]:
    return [[float(alpha) * float(a[i][j]) for j in range(len(a[0]))] for i in range(len(a))]


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


def extract_projector_matrix_2x2(a_op: dict[str, Any]) -> tuple[str, list[list[float]]]:
    data = a_op.get("data") or {}
    if not isinstance(data, dict):
        raise ValueError("A_* operator missing data dict")
    candidates: list[tuple[str, list[list[float]]]] = []
    for key, value in data.items():
        if isinstance(key, str) and "matrix_in_c" in key and is_2x2_matrix(value):
            mat = [[float(value[i][j]) for j in range(2)] for i in range(2)]
            candidates.append((key, mat))
    if len(candidates) != 1:
        raise ValueError(f"expected exactly one 2x2 projector matrix candidate, got {len(candidates)}")
    return candidates[0]


def parse_edge_key(edge_key: str) -> tuple[str, str]:
    parts = edge_key.split("_to_")
    if len(parts) != 2:
        raise ValueError(f"unexpected transition key: {edge_key}")
    return parts[0], parts[1]


def iter_edges(transition: dict[str, Any]) -> Iterable[tuple[str, str, list[list[float]], str]]:
    ops = transition.get("transition_operators") or {}
    if not isinstance(ops, dict):
        raise ValueError("transition_operators missing or not a dict")
    for edge_key, info in ops.items():
        if not isinstance(edge_key, str) or not isinstance(info, dict):
            continue
        src, dst = parse_edge_key(edge_key)
        op_ref = info.get("operator_ref")
        if not isinstance(op_ref, str):
            raise ValueError(f"missing operator_ref for edge {edge_key}")
        op_json = load_json(REPO / op_ref)
        g = extract_so2_transport(op_json)
        yield src, dst, g, op_ref


def main() -> None:
    f656 = load_json(F656_SUMMARY)
    r1 = f656.get("selector_reduction_operator", {}).get("matrix_in_c1_s1")
    if not is_2x2_matrix(r1):
        raise ValueError("missing seed R_sel matrix_in_c1_s1 in F656 summary")
    r1 = [[float(r1[i][j]) for j in range(2)] for i in range(2)]

    transition = load_json(GLOBAL_TRANSITION)
    state = load_json(GLOBAL_STATE_PROJECTIVE)
    atlas = load_json(GLOBAL_ATLAS)
    b_global = load_json(F658_OBJECT)

    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    charts = list((state.get("charts", {}) or {}).keys())
    if sorted(charts) != sorted(expected_charts):
        raise ValueError(f"unexpected global state charts: {charts}")

    # Build a bidirectional graph of SO(2) transport matrices between the pair planes.
    g_edge: dict[tuple[str, str], list[list[float]]] = {}
    edge_ref: dict[tuple[str, str], str] = {}
    for src, dst, g, op_ref in iter_edges(transition):
        g_edge[(src, dst)] = g
        edge_ref[(src, dst)] = op_ref
        g_edge[(dst, src)] = mat_t(g)
        edge_ref[(dst, src)] = op_ref

    # Promote R_sel from pair1 to all charts by transport:
    # R_dst := R_src * G_{src->dst}^T, so that R_dst(x_dst)=R_src(x_src) for x_dst=G x_src.
    r_by_pair: dict[str, list[list[float]]] = {"pair1": r1}
    q: deque[str] = deque(["pair1"])
    while q:
        src = q.popleft()
        for (i, j), g in list(g_edge.items()):
            if i != src:
                continue
            if j in r_by_pair:
                continue
            r_by_pair[j] = mat_mul(r_by_pair[src], mat_t(g))
            q.append(j)

    missing = [c for c in expected_charts if c not in r_by_pair]
    if missing:
        raise ValueError(f"could not transport R_sel to charts: {missing}")

    # Consistency witnesses vs already exported global B_sel and projective state operators A_m.
    i2 = identity2()
    d = [[1.0, 0.0], [0.0, -1.0]]
    p_q_plus = [[1.0, 0.0], [0.0, 0.0]]

    chart_payload: dict[str, Any] = {}
    r_orthogonality_max_abs_diff: dict[str, float] = {}
    b_alignment_max_abs_diff: dict[str, float] = {}
    plus_projector_alignment_max_abs_diff: dict[str, float] = {}

    for chart in expected_charts:
        r = r_by_pair[chart]

        # Orthogonality check (exporter-side witness; probe repeats).
        rrt = mat_mul(r, mat_t(r))
        r_orthogonality_max_abs_diff[chart] = max_abs(mat_sub(rrt, i2))

        b_from_r = mat_mul(mat_mul(mat_t(r), d), r)
        b_ref = (
            (b_global.get("charts", {}) or {}).get(chart, {}) or {}
        ).get("B_sel_matrix_in_c_s")
        if not is_2x2_matrix(b_ref):
            raise ValueError(f"missing B_sel_matrix_in_c_s for {chart} in F658 object")
        b_ref = [[float(b_ref[i][j]) for j in range(2)] for i in range(2)]
        b_alignment_max_abs_diff[chart] = max_abs(mat_sub(b_from_r, b_ref))

        p_plus_from_r = mat_mul(mat_mul(mat_t(r), p_q_plus), r)

        a_ref = state["charts"][chart]["local_operator"]["ref"]
        a_op = load_json(REPO / a_ref)
        a_key, a_mat = extract_projector_matrix_2x2(a_op)
        plus_projector_alignment_max_abs_diff[chart] = max_abs(mat_sub(p_plus_from_r, a_mat))

        chart_payload[chart] = {
            "chart_id": chart,
            "domain_basis": [f"c{chart[-1]}", f"s{chart[-1]}"],
            "codomain_basis": ["q_+", "q_-"],
            "R_sel_matrix_in_c_s_to_q": r,
            "implied_e_parallel_coords_in_c_s": [float(r[0][0]), float(r[0][1])],
            "implied_e_transverse_coords_in_c_s": [float(r[1][0]), float(r[1][1])],
            "aligned_global_B_sel_ref": str(F658_OBJECT.relative_to(REPO)),
            "alignment_max_abs_diff_B_vs_RtDR": b_alignment_max_abs_diff[chart],
            "aligned_projective_state_local_operator_ref": a_ref,
            "aligned_projective_state_local_operator_matrix_key": a_key,
            "alignment_max_abs_diff_P_plus_vs_A": plus_projector_alignment_max_abs_diff[chart],
        }

    obj: dict[str, Any] = {
        "object": "SelectorReductionOperator_global_C_v1_seed_v1_promoted_strict_v1",
        "stage": "F659",
        "status": "actual_exported_global_selector_reduction_operator_object__chartwise__projector_section_level__sign_gauge_explicit__no_false_pass",
        "as_of": "2026-03-17",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Promote the seed-v1 local selector reduction operator R_sel on pair1 to a global C_v1-typed chartwise "
            "selector reduction operator family on {pair1..pair5}, glued via the exported global selector transition data "
            "(up to residual Z2 sign gauge on axis-only edges), and aligned with the already exported global selector bridge operator "
            "and projective selector state operators."
        ),
        "domain": atlas.get("domain"),
        "inputs": {
            "seed_v1_local_R_sel_summary": str(F656_SUMMARY.relative_to(REPO)),
            "global_selector_atlas": str(GLOBAL_ATLAS.relative_to(REPO)),
            "global_selector_transition": str(GLOBAL_TRANSITION.relative_to(REPO)),
            "global_projective_selector_state": str(GLOBAL_STATE_PROJECTIVE.relative_to(REPO)),
            "global_selector_bridge_operator": str(F658_OBJECT.relative_to(REPO)),
            "boundary": "N512 (no operator-level groupoid promotion; projector/section-level gluing only)",
        },
        "construction": {
            "seed_pair1_operator": "R_sel_s_sel_int_source_object_v1 (pair1 chart; from F656)",
            "promotion_rule": "R_dst := R_src * G_{src->dst}^T on pair planes (SO(2) transport from the exported transition object)",
            "residual_sign_note": (
                "Some transport edges are axis-only (alpha mod pi), so the underlying plane transport is defined only up to a global sign. "
                "Therefore overlap consistency is asserted only up to residual Z2 sign gauge (projector/section sense), "
                "and no directed physical sign datum is claimed."
            ),
        },
        "charts": chart_payload,
        "exporter_witnesses": {
            "R_orthogonality_max_abs_diff_by_chart": r_orthogonality_max_abs_diff,
            "B_alignment_max_abs_diff_by_chart": b_alignment_max_abs_diff,
            "projective_state_alignment_max_abs_diff_P_plus_vs_A_by_chart": plus_projector_alignment_max_abs_diff,
        },
        "hard_limits": [
            "no_admissible_S_sel_int",
            "no_strict_core_selector_closure",
            "no_global_QW2191_discharge",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_physical_sign_orientation_claim",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary: dict[str, Any] = {
        "stage": "F659",
        "lane": "current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_only",
        "goal": "export_one_global_selector_reduction_operator_object_promoted_from_seed_v1_local_R_sel_on_pair1_to_C_v1",
        "status": "F659_EXECUTED_CURRENT_STRICT_GLOBAL_SELECTOR_REDUCTION_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_PACKET_NO_FALSE_PASS",
        "exported_object": obj["object"],
        "exported_object_file": str(OUT_OBJECT.relative_to(REPO)),
        "chart_count": len(expected_charts),
        "R_orthogonality_max_abs_diff_by_chart": r_orthogonality_max_abs_diff,
        "B_alignment_max_abs_diff_by_chart": b_alignment_max_abs_diff,
        "projective_state_alignment_max_abs_diff_P_plus_vs_A_by_chart": plus_projector_alignment_max_abs_diff,
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

