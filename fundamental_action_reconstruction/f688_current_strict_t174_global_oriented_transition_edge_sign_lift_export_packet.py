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

IN_CV1 = GENERATED / "c_v1_void_configuration_space_in_local_b_tilde_1_sector_v1.json"
IN_TRANSITION_AXIS_ONLY = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
IN_DIRECTED_STATE = (
    GENERATED
    / "selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)

OUT = (
    GENERATED
    / "selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_w_break_rooted_directed_state_strict_convention_v1.json"
)
OUT_SUMMARY = GENERATED / "f688_current_strict_t174_global_oriented_transition_edge_sign_lift_export_packet_summary.json"


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


def find_alpha_in_edge_spec(edge_spec: dict[str, Any]) -> tuple[str | None, float | None]:
    # Prefer alpha_mod_2pi if present; otherwise accept alpha_mod_pi.
    best_key = None
    best_val: float | None = None
    for key in sorted(edge_spec.keys()):
        if not key.startswith("alpha"):
            continue
        val = edge_spec.get(key)
        if not isinstance(val, (int, float)) or not math.isfinite(float(val)):
            continue
        if key.endswith("_mod_2pi"):
            return key, float(val)
        if key.endswith("_mod_pi") and best_key is None:
            best_key, best_val = key, float(val)
    return best_key, best_val


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_CV1, IN_TRANSITION_AXIS_ONLY, IN_DIRECTED_STATE]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "F688",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "expected": [
                "F469 exported global axis-only transition object on C_v1",
                "F684 exported w_break-rooted directed state representative section on C_v1 (strict_convention scope)",
            ],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    cv1 = load_json(IN_CV1)
    transition_axis_only = load_json(IN_TRANSITION_AXIS_ONLY)
    directed_state = load_json(IN_DIRECTED_STATE)

    # Load directed vectors u_m from the exported directed state representative.
    u_outs = ((directed_state.get("outputs") or {}).get("u_vectors_directed") or {})
    if not isinstance(u_outs, dict):
        raise SystemExit("Invalid directed state: expected outputs.u_vectors_directed dict")

    u_by_pair: dict[str, list[float]] = {}
    for m in range(1, 6):
        key = f"u_{m}"
        vec = u_outs.get(key)
        if not is_numeric_list_len(vec, 12):
            raise SystemExit(f"Invalid directed state: outputs.u_vectors_directed.{key} must be length-12 numeric list")
        u_by_pair[f"pair{m}"] = [float(v) for v in vec]

    trans_ops = (transition_axis_only.get("transition_operators") or {})
    if not isinstance(trans_ops, dict):
        raise SystemExit("Invalid transition object: expected transition_operators dict")

    tol_match = 1e-9
    edges: dict[str, Any] = {}
    edge_errors: list[dict[str, Any]] = []

    for edge_id, edge_spec in trans_ops.items():
        if not isinstance(edge_spec, dict):
            edge_errors.append({"edge": str(edge_id), "error": "edge_spec is not a dict"})
            continue

        try:
            src, dst = parse_edge_id(str(edge_id))
            if src not in u_by_pair or dst not in u_by_pair:
                raise ValueError(f"Missing directed vectors for src={src!r} or dst={dst!r}")

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

            diff_same = float(max_abs_diff(v, u_dst))
            diff_flip = float(max_abs_diff(v, [-x for x in u_dst]))

            sign_on_operator = None
            if diff_same <= tol_match:
                sign_on_operator = 1
            elif diff_flip <= tol_match:
                sign_on_operator = -1

            if sign_on_operator is None:
                raise ValueError(
                    f"Base operator does not map u_src to ±u_dst within tol={tol_match}; "
                    f"diff_to_uj={diff_same}, diff_to_minus_uj={diff_flip}"
                )

            alpha_key, alpha_val = find_alpha_in_edge_spec(edge_spec)
            oriented_alpha_mod_2pi = None
            if alpha_key is not None and alpha_val is not None:
                oriented_alpha_mod_2pi = (float(alpha_val) + (0.0 if sign_on_operator == 1 else math.pi)) % (
                    2.0 * math.pi
                )

            edges[str(edge_id)] = {
                "edge_id": str(edge_id),
                "src": src,
                "dst": dst,
                "base_operator_ref": operator_ref,
                "base_operator": edge_spec.get("operator"),
                "base_operator_matrix_key": matrix_key,
                "base_alpha_source_key": alpha_key,
                "base_alpha_value": alpha_val,
                "oriented_sign_on_operator": int(sign_on_operator),
                "oriented_alpha_mod_2pi": oriented_alpha_mod_2pi,
                "checks": {
                    "tol_match_max_abs_diff": tol_match,
                    "max_abs_diff_to_uj": diff_same,
                    "max_abs_diff_to_minus_uj": diff_flip,
                },
            }
        except Exception as exc:
            edge_errors.append({"edge": str(edge_id), "error": str(exc)})

    edge_count = len(edges)
    sign_lift_edges = [eid for eid, e in edges.items() if int(e.get("oriented_sign_on_operator") or 0) == -1]

    if edge_errors:
        status = "NOT_COMPUTABLE_EDGE_SIGN_LIFT_ERRORS_PRESENT"
    elif edge_count == 0:
        status = "NOT_COMPUTABLE_NO_EDGES_FOUND_IN_BASE_TRANSITION_OBJECT"
    else:
        status = "PASS_EXPORTED_GLOBAL_ORIENTED_TRANSITION_EDGE_SIGN_LIFT_OBJECT"

    artifact = {
        "object": "SelectorTransition_global_C_v1_oriented_mod_2pi_edge_sign_lift_from_w_break_rooted_directed_state_strict_convention_v1",
        "stage": "F688",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export an explicit edgewise oriented (α mod 2π) sign-lift of the exported global strict transition object on C_v1 by "
            "choosing per-edge signs s_ij ∈ {±1} such that (s_ij * O_ij) transports the already exported w_break-rooted directed "
            "state representative u_i to u_j without sign flips on every edge. This is a strict_convention layer: it does not claim any "
            "strict-core physical sign datum, does not claim Aut(Z_12)-invariant canonicity, does not promote to operator-level groupoid "
            "identities (N512), and does not imply kernel-alone/global QW-2191 discharge."
        ),
        "domain": {
            "configuration_space_object_ref": str(IN_CV1.relative_to(REPO)),
            "configuration_space_object": str(cv1.get("object") or "C_v1"),
            "charts_cover": "U_pair1 = ... = U_pair5 = C_v1 (declared global cover; convention layer)",
        },
        "depends_on": {
            "base_axis_only_transition_object_ref": str(IN_TRANSITION_AXIS_ONLY.relative_to(REPO)),
            "base_axis_only_transition_object": str(transition_axis_only.get("object") or ""),
            "directed_state_ref": str(IN_DIRECTED_STATE.relative_to(REPO)),
            "directed_state_object": str(directed_state.get("object") or ""),
        },
        "charts": ["pair1", "pair2", "pair3", "pair4", "pair5"],
        "edge_sign_lift": {
            "meaning": "For each overlap edge, export s_ij ∈ {±1} defining an oriented lift O_ij^(oriented) := s_ij * O_ij^(axis-only).",
            "transition_edges": edges,
            "edge_count": edge_count,
            "edge_errors": edge_errors,
            "sign_lift_edges": sign_lift_edges,
            "sign_lift_count": len(sign_lift_edges),
        },
        "cocycle_discipline": {
            "level": "vector_section_level_only",
            "supports": [
                "P686/N686: base global transport is compatible only up to sign on directed representatives",
                "P687/N687: no per-chart sign relift eliminates edge sign flips under fixed axis-only reps",
                "N512: no operator-level groupoid promotion",
            ],
            "explicit_non_promotion": [
                "Do not infer operator-level identities O_jk O_ij = O_ik on the full carrier from edgewise sign-lift data.",
                "Do not infer a strict-core physical sign datum from this convention-layer oriented lift.",
            ],
        },
        "hard_limits": [
            "Convention layer only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim Aut(Z_12)-invariant sign canonicity (N462 discipline).",
            "Does not promote section-level/vector-level consistency into operator-level groupoid identities (N512 boundary).",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F688",
        "status": status,
        "exported": {
            "oriented_transition_edge_sign_lift_ref": str(OUT.relative_to(REPO)),
        },
        "edge_count": edge_count,
        "sign_lift_count": len(sign_lift_edges),
        "sign_lift_edges": sign_lift_edges,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

