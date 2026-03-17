#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Tuple

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
IN_STATE = (
    GENERATED
    / "selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)

OUT = (
    GENERATED
    / "p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def dot(a: list[float], b: list[float]) -> float:
    return sum(float(x) * float(y) for x, y in zip(a, b))


def matvec12(m: list[list[float]], v: list[float]) -> list[float]:
    return [sum(float(m[i][j]) * float(v[j]) for j in range(12)) for i in range(12)]


def l2_norm(v: list[float]) -> float:
    return math.sqrt(sum(float(x) * float(x) for x in v))


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


def parse_edge_id(edge_id: str) -> Tuple[str, str]:
    parts = edge_id.split("_to_")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(f"Invalid edge id: {edge_id!r} (expected 'pairi_to_pairj')")
    return parts[0], parts[1]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_TRANSITION, IN_STATE]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P686",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    transition = load_json(IN_TRANSITION)
    state = load_json(IN_STATE)

    edges: dict[str, Any] = {}
    edge_errors: list[dict[str, Any]] = []

    tol_match = 1e-9
    tol_line_dot = 1e-9

    # Load directed vectors u_m from the exported directed state representative.
    u_outs = ((state.get("outputs") or {}).get("u_vectors_directed") or {})
    if not isinstance(u_outs, dict):
        artifact = {
            "stage": "P686",
            "status": "INVALID_DIRECTED_STATE_SHAPE",
            "as_of": AS_OF,
            "error": "F684 directed state must contain outputs.u_vectors_directed dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    u_by_pair: dict[str, list[float]] = {}
    for m in range(1, 6):
        key = f"u_{m}"
        vec = u_outs.get(key)
        if not is_numeric_list_len(vec, 12):
            artifact = {
                "stage": "P686",
                "status": "INVALID_DIRECTED_STATE_VECTOR_SHAPE",
                "as_of": AS_OF,
                "error": f"F684 directed state outputs.u_vectors_directed.{key} must be length-12 numeric list",
                "no_false_pass": True,
            }
            OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT)
            return
        u_by_pair[f"pair{m}"] = [float(v) for v in vec]

    # Iterate all exported transition edges.
    trans_ops = (transition.get("transition_operators") or {})
    if not isinstance(trans_ops, dict):
        artifact = {
            "stage": "P686",
            "status": "INVALID_TRANSITION_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "Global transition object must contain transition_operators dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    for edge_id, edge_spec in trans_ops.items():
        if not isinstance(edge_spec, dict):
            edge_errors.append({"edge": edge_id, "error": "edge_spec is not a dict"})
            continue

        try:
            src, dst = parse_edge_id(str(edge_id))
            if src not in u_by_pair or dst not in u_by_pair:
                raise ValueError(f"Missing directed vector for src={src!r} or dst={dst!r} in F684 outputs")

            operator_ref = edge_spec.get("operator_ref")
            if not isinstance(operator_ref, str) or not operator_ref:
                raise ValueError("Missing operator_ref")

            op_path = (REPO / operator_ref).resolve()
            if not op_path.exists():
                raise FileNotFoundError(f"Missing operator file: {operator_ref}")

            op_obj = load_json(op_path)
            matrix_key, O = extract_12x12_matrix_from_outputs(op_obj, path=op_path)

            v = matvec12(O, u_by_pair[src])
            u_dst = u_by_pair[dst]

            dot_v_u_dst = float(dot(v, u_dst))
            dot_abs = abs(dot_v_u_dst)
            line_ok = bool(dot_abs >= 1.0 - tol_line_dot)

            diff_same = float(max_abs_diff(v, u_dst))
            diff_flip = float(max_abs_diff(v, [-x for x in u_dst]))

            best_sign = None
            if diff_same <= tol_match:
                best_sign = 1
            elif diff_flip <= tol_match:
                best_sign = -1

            edge_ok = best_sign is not None

            edges[str(edge_id)] = {
                "edge_id": str(edge_id),
                "src": src,
                "dst": dst,
                "operator_ref": operator_ref,
                "operator_matrix_key": matrix_key,
                "dot_Oij_ui_with_uj": dot_v_u_dst,
                "abs_dot": dot_abs,
                "line_compatible_by_dot": line_ok,
                "max_abs_diff_to_uj": diff_same,
                "max_abs_diff_to_minus_uj": diff_flip,
                "best_match_sign": best_sign,
                "edge_compatible_up_to_sign": edge_ok,
                "sign_compatible_directed": best_sign == 1,
                "sign_flip_present": best_sign == -1,
                "l2_norm_Oij_ui": float(l2_norm(v)),
                "l2_norm_uj": float(l2_norm(u_dst)),
            }
        except Exception as exc:
            edge_errors.append({"edge": str(edge_id), "error": str(exc)})

    edge_count = len(edges)
    bad_edges = [eid for eid, e in edges.items() if not bool(e.get("edge_compatible_up_to_sign"))]
    sign_flip_edges = [eid for eid, e in edges.items() if bool(e.get("sign_flip_present"))]
    sign_mismatch_edges = [eid for eid, e in edges.items() if not bool(e.get("sign_compatible_directed"))]

    all_edges_up_to_sign = (edge_count > 0) and (len(bad_edges) == 0) and (len(edge_errors) == 0)
    all_edges_directed_sign = all_edges_up_to_sign and (len(sign_mismatch_edges) == 0)

    if edge_errors:
        status = "NOT_COMPUTABLE_EDGE_AUDIT_ERRORS_PRESENT"
    elif not edges:
        status = "NOT_COMPUTABLE_NO_EDGES_FOUND_IN_TRANSITION_OBJECT"
    elif not all_edges_up_to_sign:
        status = "FAIL_DIRECTED_STATE_NOT_COMPATIBLE_WITH_GLOBAL_TRANSITION_OPERATORS_UP_TO_SIGN"
    elif all_edges_directed_sign:
        status = "PASS_DIRECTED_STATE_EDGE_COMPATIBLE_WITH_GLOBAL_TRANSITION_OPERATORS_WITHOUT_SIGN_FLIPS"
    else:
        status = "PASS_DIRECTED_STATE_EDGE_COMPATIBLE_UP_TO_SIGN_BUT_SIGN_FLIPS_PRESENT"

    artifact = {
        "stage": "P686",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_full_edgewise_compatibility_of_the_exported_w_break_rooted_directed_state_representative_with_the_exported_global_selector_transition_object_on_C_v1__no_false_pass",
        "inputs": {
            "directed_state_ref": str(IN_STATE.relative_to(REPO)),
            "transition_object_ref": str(IN_TRANSITION.relative_to(REPO)),
        },
        "status": status,
        "tolerances": {
            "tol_match_max_abs_diff": tol_match,
            "tol_line_abs_dot": tol_line_dot,
        },
        "edge_audit": {
            "edge_count": edge_count,
            "edge_errors": edge_errors,
            "edges": edges,
        },
        "result": {
            "all_edges_compatible_up_to_sign": all_edges_up_to_sign,
            "all_edges_compatible_with_directed_sign": all_edges_directed_sign,
            "bad_edges": bad_edges,
            "sign_flip_edges": sign_flip_edges,
            "sign_flip_count": len(sign_flip_edges),
        },
        "counts_as_strict_physical_orientation_datum": False,
        "hard_limits": [
            "Convention/section audit only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
            "Does not promote projector-level cocycle into an operator-level groupoid (N512 boundary).",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P686",
        "status": status,
        "edge_count": edge_count,
        "all_edges_compatible_up_to_sign": all_edges_up_to_sign,
        "all_edges_compatible_with_directed_sign": all_edges_directed_sign,
        "sign_flip_count": len(sign_flip_edges),
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

