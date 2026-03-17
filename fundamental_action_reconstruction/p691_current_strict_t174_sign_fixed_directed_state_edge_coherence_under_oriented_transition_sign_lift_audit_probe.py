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

IN_LIFT = (
    GENERATED
    / "selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_sign_fixed_directed_state_strict_convention_v1.json"
)
IN_STATE = (
    GENERATED
    / "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)

OUT = (
    GENERATED
    / "p691_current_strict_t174_sign_fixed_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p691_current_strict_t174_sign_fixed_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe_summary.json"
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

    prereq = [IN_LIFT, IN_STATE]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P691",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    lift = load_json(IN_LIFT)
    state = load_json(IN_STATE)

    tol_match = 1e-9

    u_outs = ((state.get("outputs") or {}).get("u_vectors_directed_sign_fixed") or {})
    if not isinstance(u_outs, dict):
        artifact = {
            "stage": "P691",
            "status": "INVALID_SIGN_FIXED_STATE_SHAPE",
            "as_of": AS_OF,
            "error": "Sign-fixed state must contain outputs.u_vectors_directed_sign_fixed dict",
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
                "stage": "P691",
                "status": "INVALID_SIGN_FIXED_STATE_VECTOR_SHAPE",
                "as_of": AS_OF,
                "error": f"Sign-fixed state outputs.u_vectors_directed_sign_fixed.{key} must be length-12 numeric list",
                "no_false_pass": True,
            }
            OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT)
            return
        u_by_pair[f"pair{m}"] = [float(v) for v in vec]

    edges_obj = ((lift.get("edge_sign_lift") or {}).get("transition_edges") or {})
    if not isinstance(edges_obj, dict):
        artifact = {
            "stage": "P691",
            "status": "INVALID_LIFT_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "Lift object must contain edge_sign_lift.transition_edges dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    edge_results: dict[str, Any] = {}
    edge_errors: list[dict[str, Any]] = []

    for edge_id, edge_spec in edges_obj.items():
        if not isinstance(edge_spec, dict):
            edge_errors.append({"edge": str(edge_id), "error": "edge_spec is not a dict"})
            continue

        try:
            src, dst = parse_edge_id(str(edge_id))
            if src not in u_by_pair or dst not in u_by_pair:
                raise ValueError(f"Missing directed vector for src={src!r} or dst={dst!r}")

            operator_ref = edge_spec.get("base_operator_ref")
            if not isinstance(operator_ref, str) or not operator_ref:
                raise ValueError("Missing base_operator_ref")

            sign_on_operator = edge_spec.get("oriented_sign_on_operator")
            if sign_on_operator not in (-1, 1):
                raise ValueError("Missing or invalid oriented_sign_on_operator (expected ±1)")

            op_path = (REPO / operator_ref).resolve()
            if not op_path.exists():
                raise FileNotFoundError(f"Missing operator file: {operator_ref}")

            op_obj = load_json(op_path)
            matrix_key, O = extract_12x12_matrix_from_outputs(op_obj, path=op_path)

            v = matvec12(O, u_by_pair[src])
            v_oriented = [float(sign_on_operator) * x for x in v]
            u_dst = u_by_pair[dst]

            diff = float(max_abs_diff(v_oriented, u_dst))
            ok = bool(diff <= tol_match)

            edge_results[str(edge_id)] = {
                "edge_id": str(edge_id),
                "src": src,
                "dst": dst,
                "base_operator_ref": operator_ref,
                "base_operator_matrix_key": matrix_key,
                "oriented_sign_on_operator": int(sign_on_operator),
                "max_abs_diff_to_uj": diff,
                "edge_ok_no_sign_flip": ok,
            }
        except Exception as exc:
            edge_errors.append({"edge": str(edge_id), "error": str(exc)})

    edge_count = len(edge_results)
    bad_edges = [eid for eid, e in edge_results.items() if not bool(e.get("edge_ok_no_sign_flip"))]

    all_edges_ok = (edge_count > 0) and (len(bad_edges) == 0) and (len(edge_errors) == 0)

    if edge_errors:
        status = "NOT_COMPUTABLE_EDGE_AUDIT_ERRORS_PRESENT"
    elif edge_count == 0:
        status = "NOT_COMPUTABLE_NO_EDGES_FOUND_IN_LIFT_OBJECT"
    elif all_edges_ok:
        status = "PASS_DIRECTED_STATE_EDGE_COHERENT_UNDER_ORIENTED_TRANSITION_SIGN_LIFT_WITHOUT_SIGN_FLIPS"
    else:
        status = "FAIL_DIRECTED_STATE_NOT_EDGE_COHERENT_UNDER_ORIENTED_TRANSITION_SIGN_LIFT"

    artifact = {
        "stage": "P691",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_full_edgewise_coherence_of_the_exported_sign_fixed_directed_state_under_the_exported_oriented_transition_edge_sign_lift__no_false_pass",
        "inputs": {
            "sign_fixed_directed_state_ref": str(IN_STATE.relative_to(REPO)),
            "oriented_transition_edge_sign_lift_ref": str(IN_LIFT.relative_to(REPO)),
        },
        "status": status,
        "tolerances": {"tol_match_max_abs_diff": tol_match},
        "edge_audit": {
            "edge_count": edge_count,
            "edge_errors": edge_errors,
            "edges": edge_results,
        },
        "result": {
            "all_edges_ok_without_sign_flips": all_edges_ok,
            "bad_edges": bad_edges,
            "bad_edge_count": len(bad_edges),
        },
        "counts_as_strict_physical_orientation_datum": False,
        "hard_limits": [
            "Convention/section audit only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
            "Does not promote section-level/vector-level coherence into operator-level groupoid identities (N512 boundary).",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P691",
        "status": status,
        "edge_count": edge_count,
        "all_edges_ok_without_sign_flips": all_edges_ok,
        "bad_edge_count": len(bad_edges),
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

