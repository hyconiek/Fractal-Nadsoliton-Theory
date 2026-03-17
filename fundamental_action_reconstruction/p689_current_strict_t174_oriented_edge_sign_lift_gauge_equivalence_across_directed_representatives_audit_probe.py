#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Tuple

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F688_LIFT = (
    GENERATED
    / "selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_w_break_rooted_directed_state_strict_convention_v1.json"
)
IN_F469_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
IN_DIRECTED_STATE_PREMISE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"

OUT = (
    GENERATED
    / "p689_current_strict_t174_oriented_edge_sign_lift_gauge_equivalence_across_directed_representatives_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p689_current_strict_t174_oriented_edge_sign_lift_gauge_equivalence_across_directed_representatives_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def matvec12(m: list[list[float]], v: list[float]) -> list[float]:
    return [sum(float(m[i][j]) * float(v[j]) for j in range(12)) for i in range(12)]


def max_abs_diff(a: list[float], b: list[float]) -> float:
    return max(abs(float(x) - float(y)) for x, y in zip(a, b))


def extract_12x12_matrix_from_outputs(obj: dict[str, Any], *, path: Path) -> Tuple[str, list[list[float]]]:
    outs = obj.get("outputs")
    if not isinstance(outs, dict):
        raise ValueError(f"{path.relative_to(REPO)}: missing dict 'outputs'")

    for key, value in outs.items():
        if isinstance(value, list) and len(value) == 12 and all(is_numeric_list_len(row, 12) for row in value):
            return str(key), [[float(x) for x in row] for row in value]
    raise ValueError(f"{path.relative_to(REPO)}: outputs contains no 12x12 matrix")


def parse_edge_id(edge_id: str) -> tuple[str, str]:
    parts = edge_id.split("_to_")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(f"Invalid edge id: {edge_id!r} (expected 'pairi_to_pairj')")
    return parts[0], parts[1]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F688_LIFT, IN_F469_TRANSITION, IN_DIRECTED_STATE_PREMISE]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P689",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    tol_match = 1e-9

    lift_a = load_json(IN_F688_LIFT)
    trans = load_json(IN_F469_TRANSITION)
    state_b = load_json(IN_DIRECTED_STATE_PREMISE)

    # Load directed vectors (premise-based representative B).
    u_outs_b = ((state_b.get("outputs") or {}).get("u_vectors_directed") or {})
    if not isinstance(u_outs_b, dict):
        artifact = {
            "stage": "P689",
            "status": "INVALID_DIRECTED_STATE_SHAPE",
            "as_of": AS_OF,
            "error": "Premise-based directed state must contain outputs.u_vectors_directed dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    u_by_pair_b: dict[str, list[float]] = {}
    for m in range(1, 6):
        key = f"u_{m}"
        vec = u_outs_b.get(key)
        if not is_numeric_list_len(vec, 12):
            artifact = {
                "stage": "P689",
                "status": "INVALID_DIRECTED_STATE_VECTOR_SHAPE",
                "as_of": AS_OF,
                "error": f"Premise-based directed state outputs.u_vectors_directed.{key} must be length-12 numeric list",
                "no_false_pass": True,
            }
            OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT)
            return
        u_by_pair_b[f"pair{m}"] = [float(v) for v in vec]

    # Load edge sign-lift pattern A from F688.
    edges_a = ((lift_a.get("edge_sign_lift") or {}).get("transition_edges") or {})
    if not isinstance(edges_a, dict):
        artifact = {
            "stage": "P689",
            "status": "INVALID_LIFT_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "F688 lift object must contain edge_sign_lift.transition_edges dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    s_a: dict[str, int] = {}
    s_a_errors: list[dict[str, Any]] = []
    for edge_id, edge_spec in edges_a.items():
        if not isinstance(edge_spec, dict):
            s_a_errors.append({"edge": str(edge_id), "error": "edge_spec is not a dict"})
            continue
        s = edge_spec.get("oriented_sign_on_operator")
        if s not in (-1, 1):
            s_a_errors.append({"edge": str(edge_id), "error": "missing or invalid oriented_sign_on_operator (expected ±1)"})
            continue
        s_a[str(edge_id)] = int(s)

    # Compute edge sign-lift pattern B induced by the premise-based directed representative under base axis-only transitions.
    trans_ops = trans.get("transition_operators") or {}
    if not isinstance(trans_ops, dict) or not trans_ops:
        artifact = {
            "stage": "P689",
            "status": "INVALID_TRANSITION_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "Transition object must contain a non-empty dict transition_operators keyed by edge_id",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    # Restrict to edges present in pattern A, to keep comparison strict.
    edge_ids = sorted(set(s_a.keys()))

    missing_in_trans = [eid for eid in edge_ids if eid not in trans_ops]
    if missing_in_trans:
        artifact = {
            "stage": "P689",
            "status": "NOT_COMPUTABLE_EDGE_SET_MISMATCH",
            "as_of": AS_OF,
            "error": "Some F688 lift edges are missing in the base transition object",
            "missing_edges_in_transition_object": missing_in_trans,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    operator_cache: dict[str, tuple[str, list[list[float]]]] = {}

    def load_operator_matrix(operator_ref: str) -> tuple[str, list[list[float]]]:
        if operator_ref in operator_cache:
            return operator_cache[operator_ref]
        op_path = (REPO / operator_ref).resolve()
        if not op_path.exists():
            raise FileNotFoundError(f"Missing operator file: {operator_ref}")
        op_obj = load_json(op_path)
        matrix_key, O = extract_12x12_matrix_from_outputs(op_obj, path=op_path)
        operator_cache[operator_ref] = (matrix_key, O)
        return matrix_key, O

    s_b: dict[str, int] = {}
    edge_audit_b: dict[str, Any] = {}
    edge_errors: list[dict[str, Any]] = []

    for edge_id in edge_ids:
        try:
            trans_edge = trans_ops.get(edge_id)
            if not isinstance(trans_edge, dict):
                raise ValueError("transition_operators edge entry is not a dict")
            operator_ref = trans_edge.get("operator_ref")
            if not isinstance(operator_ref, str) or not operator_ref:
                raise ValueError("Missing operator_ref for base transition edge")

            src, dst = parse_edge_id(edge_id)
            if src not in u_by_pair_b or dst not in u_by_pair_b:
                raise ValueError(f"Missing directed vector for src={src!r} or dst={dst!r}")

            matrix_key, O = load_operator_matrix(operator_ref)
            v = matvec12(O, u_by_pair_b[src])
            u_dst = u_by_pair_b[dst]

            diff_plus = float(max_abs_diff(v, u_dst))
            diff_minus = float(max_abs_diff(v, [-x for x in u_dst]))

            # Choose sign that eliminates the sign flip: (s * O) u_src ≈ u_dst.
            s = 1 if diff_plus <= diff_minus else -1
            best = diff_plus if s == 1 else diff_minus

            s_b[edge_id] = int(s)
            edge_audit_b[edge_id] = {
                "edge_id": edge_id,
                "src": src,
                "dst": dst,
                "base_operator_ref": operator_ref,
                "base_operator_matrix_key": matrix_key,
                "s_b_oriented_sign_on_operator": int(s),
                "max_abs_diff_to_uj": best,
                "max_abs_diff_to_uj_without_sign": diff_plus,
                "max_abs_diff_to_minus_uj_without_sign": diff_minus,
                "edge_ok": bool(best <= tol_match),
            }
        except Exception as exc:
            edge_errors.append({"edge": edge_id, "error": str(exc)})

    bad_edges_b = [eid for eid, e in edge_audit_b.items() if not bool(e.get("edge_ok"))]

    # Gauge check: does there exist a chart-level 0-cochain t_i ∈ {±1} such that
    # s_b(e) = t(src) * s_a(e) * t(dst)?
    r_edges: dict[str, int] = {}
    for edge_id in edge_ids:
        if edge_id in s_b:
            r_edges[edge_id] = int(s_b[edge_id] * s_a[edge_id])

    nodes: set[str] = set()
    for edge_id in edge_ids:
        src, dst = parse_edge_id(edge_id)
        nodes.add(src)
        nodes.add(dst)

    # BFS solve r_ij = t_i * t_j.
    t_by_pair: dict[str, int] = {}
    gauge_mismatches: list[dict[str, Any]] = []

    if nodes:
        root = "pair1" if "pair1" in nodes else sorted(nodes)[0]
        t_by_pair[root] = 1
        q: deque[str] = deque([root])

        adj: dict[str, list[tuple[str, str]]] = {n: [] for n in nodes}
        for edge_id in edge_ids:
            src, dst = parse_edge_id(edge_id)
            adj[src].append((dst, edge_id))
            adj[dst].append((src, edge_id))

        while q:
            src = q.popleft()
            for dst, edge_id in adj[src]:
                if edge_id not in r_edges:
                    continue
                implied = int(r_edges[edge_id] * t_by_pair[src])
                if dst not in t_by_pair:
                    t_by_pair[dst] = implied
                    q.append(dst)
                else:
                    if t_by_pair[dst] != implied:
                        gauge_mismatches.append(
                            {
                                "edge": edge_id,
                                "src": src,
                                "dst": dst,
                                "r_edge": int(r_edges[edge_id]),
                                "t_src": int(t_by_pair[src]),
                                "t_dst_existing": int(t_by_pair[dst]),
                                "t_dst_implied": int(implied),
                            }
                        )

    disconnected_nodes = sorted([n for n in nodes if n not in t_by_pair])

    # Final verification on all edges.
    gauge_ok = (len(gauge_mismatches) == 0) and (len(disconnected_nodes) == 0)
    if gauge_ok:
        for edge_id in edge_ids:
            src, dst = parse_edge_id(edge_id)
            if int(t_by_pair[src] * t_by_pair[dst]) != int(r_edges[edge_id]):
                gauge_ok = False
                gauge_mismatches.append(
                    {
                        "edge": edge_id,
                        "src": src,
                        "dst": dst,
                        "r_edge": int(r_edges[edge_id]),
                        "t_src": int(t_by_pair.get(src, 0)),
                        "t_dst": int(t_by_pair.get(dst, 0)),
                        "error": "postcheck_failed_r_edge_not_equal_t_src_times_t_dst",
                    }
                )

    edges_where_patterns_differ = sorted([eid for eid in edge_ids if eid in s_b and s_b[eid] != s_a[eid]])

    all_edges_ok_b = (len(edge_ids) > 0) and (len(edge_errors) == 0) and (len(bad_edges_b) == 0)

    if s_a_errors:
        status = "NOT_COMPUTABLE_INVALID_F688_EDGE_SIGN_LIFT_PATTERN"
    elif edge_errors:
        status = "NOT_COMPUTABLE_EDGE_AUDIT_ERRORS_PRESENT"
    elif not all_edges_ok_b:
        status = "FAIL_PREMISE_BASED_DIRECTED_REPRESENTATIVE_NOT_EDGE_COHERENT_UNDER_ANY_SIGN_LIFT"
    elif gauge_ok:
        status = "PASS_ORIENTED_EDGE_SIGN_LIFT_PATTERNS_GAUGE_EQUIVALENT_ACROSS_DIRECTED_REPRESENTATIVES"
    else:
        status = "FAIL_ORIENTED_EDGE_SIGN_LIFT_PATTERNS_NOT_GAUGE_EQUIVALENT_ACROSS_DIRECTED_REPRESENTATIVES"

    artifact = {
        "stage": "P689",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_that_the_T174_oriented_edge_sign_lift_pattern_depends_on_directed_representative_only_by_chart_level_Z2_gauge__no_false_pass",
        "inputs": {
            "f688_oriented_edge_sign_lift_ref": str(IN_F688_LIFT.relative_to(REPO)),
            "base_transition_object_ref": str(IN_F469_TRANSITION.relative_to(REPO)),
            "premise_based_directed_state_ref": str(IN_DIRECTED_STATE_PREMISE.relative_to(REPO)),
        },
        "status": status,
        "tolerances": {"tol_match_max_abs_diff": tol_match},
        "pattern_A_from_F688": {
            "edge_count": len(edge_ids),
            "s_a_by_edge": s_a,
            "errors": s_a_errors,
        },
        "pattern_B_induced_by_premise_based_directed_state": {
            "edge_count": len(edge_ids),
            "s_b_by_edge": s_b,
            "edge_errors": edge_errors,
            "bad_edges": bad_edges_b,
            "edges": edge_audit_b,
            "all_edges_ok": all_edges_ok_b,
        },
        "gauge_equivalence_check": {
            "meaning": "Check if s_b(e) = t(src)*s_a(e)*t(dst) for some chart-level t ∈ {±1}^charts.",
            "gauge_equivalent": gauge_ok,
            "t_by_pair": {k: int(v) for k, v in sorted(t_by_pair.items())},
            "disconnected_nodes": disconnected_nodes,
            "mismatches": gauge_mismatches,
            "edges_where_patterns_differ": edges_where_patterns_differ,
        },
        "counts_as_strict_physical_orientation_datum": False,
        "hard_limits": [
            "Gauge/convention audit only; does not claim a strict-core physical sign datum.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim Aut(Z_12)-invariant sign canonicity (N462).",
            "Does not promote section-level/vector-level coherence into operator-level groupoid identities (N512 boundary).",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P689",
        "status": status,
        "edge_count": len(edge_ids),
        "premise_based_state_edge_coherent_under_induced_sign_lift": all_edges_ok_b,
        "patterns_gauge_equivalent": gauge_ok,
        "t_by_pair": {k: int(v) for k, v in sorted(t_by_pair.items())},
        "edges_where_patterns_differ": edges_where_patterns_differ,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

