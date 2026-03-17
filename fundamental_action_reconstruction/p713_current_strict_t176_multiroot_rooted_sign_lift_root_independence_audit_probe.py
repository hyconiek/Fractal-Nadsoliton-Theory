#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"
IN_Y_SEL = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
IN_F684 = GENERATED / "f684_first_exported_selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_packet_summary.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p713_current_strict_t176_multiroot_rooted_sign_lift_root_independence_audit_probe.json"
OUT_SUMMARY = GENERATED / "p713_current_strict_t176_multiroot_rooted_sign_lift_root_independence_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def is_numeric_list_len(value: Any, length: int) -> bool:
    return (
        isinstance(value, list)
        and len(value) == length
        and all(isinstance(entry, (int, float)) and math.isfinite(float(entry)) for entry in value)
    )


def sign(value: float) -> int:
    return 1 if value > 0.0 else -1


def dot(lhs: list[float], rhs: list[float]) -> float:
    return sum(float(x) * float(y) for x, y in zip(lhs, rhs))


def matvec2(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [
        float(matrix[0][0]) * float(vector[0]) + float(matrix[0][1]) * float(vector[1]),
        float(matrix[1][0]) * float(vector[0]) + float(matrix[1][1]) * float(vector[1]),
    ]


def matvec12(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(float(matrix[row_idx][col_idx]) * float(vector[col_idx]) for col_idx in range(12)) for row_idx in range(12)]


def l2_norm(vector: list[float]) -> float:
    return math.sqrt(sum(float(entry) * float(entry) for entry in vector))


def max_abs_pair_diff(lhs: list[float], rhs: list[float]) -> float:
    return max(abs(float(lhs[0]) - float(rhs[0])), abs(float(lhs[1]) - float(rhs[1])))


def load_transport_matrix(path: Path) -> list[list[float]]:
    obj = load_json(path)
    outputs = obj.get("outputs") or {}
    for key in ("O12", "O13", "O14", "O15", "O23", "O24", "O25", "O34", "O35", "O45", "O"):
        matrix = outputs.get(key)
        if isinstance(matrix, list) and len(matrix) == 12 and all(is_numeric_list_len(row, 12) for row in matrix):
            return [[float(entry) for entry in row] for row in matrix]
    raise ValueError(f"{path} does not expose a 12x12 transport matrix under expected output keys")


def pair_name(pair_index: int) -> str:
    return f"pair{pair_index}"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F647, IN_ATLAS, IN_Y_SEL, IN_F684] + [IN_A[pair_index] for pair_index in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P713",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f647 = load_json(IN_F647)
    atlas = load_json(IN_ATLAS)
    y_sel = load_json(IN_Y_SEL)
    f684_summary = load_json(IN_F684)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break = payload.get("w_break_by_x")
    if not is_numeric_list_len(w_break, 12):
        artifact = {
            "stage": "P713",
            "status": "INVALID_W_BREAK_INPUT_SHAPE",
            "as_of": AS_OF,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return
    w_break_by_x = [float(entry) for entry in w_break]

    vectors: dict[int, list[float]] = {}
    coords: dict[int, list[float]] = {}
    for pair_index in sorted(IN_A):
        operator_packet = load_json(IN_A[pair_index])
        data = operator_packet.get("data") or {}
        vector_key = "u_1" if pair_index == 1 else f"u_{pair_index}"
        coords_key = "u_1_coords_in_c1_s1" if pair_index == 1 else f"u_{pair_index}_coords_in_c{pair_index}_s{pair_index}"
        vector_value = data.get(vector_key)
        coords_value = data.get(coords_key)
        if not is_numeric_list_len(vector_value, 12):
            raise SystemExit(f"Invalid {IN_A[pair_index].relative_to(REPO)}: data.{vector_key} must be length-12 numeric list")
        if not is_numeric_list_len(coords_value, 2):
            raise SystemExit(f"Invalid {IN_A[pair_index].relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")
        vectors[pair_index] = [float(entry) for entry in vector_value]
        coords[pair_index] = [float(entry) for entry in coords_value]

    transition_refs: dict[tuple[int, int], Path] = {}
    transitions = atlas.get("transitions") or {}
    for source_index in range(1, 6):
        for target_index in range(source_index + 1, 6):
            transition_key = f"{pair_name(source_index)}_to_{pair_name(target_index)}"
            transition_obj = transitions.get(transition_key) or {}
            operator_ref = transition_obj.get("operator_ref")
            if not isinstance(operator_ref, str):
                raise SystemExit(f"Missing operator_ref for atlas transition {transition_key}")
            operator_path = REPO / operator_ref
            if not operator_path.exists():
                raise SystemExit(f"Missing operator file referenced by atlas: {operator_ref}")
            transition_refs[(source_index, target_index)] = operator_path

    transport_matrices = {
        pair_key: load_transport_matrix(path)
        for pair_key, path in transition_refs.items()
    }

    tolerance_dot_nonzero = 1e-12
    tolerance_alignment = 1e-9
    tolerance_output = 1e-9

    rooted_results: dict[str, Any] = {}
    supported_roots: list[str] = []
    unsupported_roots: list[str] = []

    for root_index in range(1, 6):
        root_pair_name = pair_name(root_index)
        root_dot = dot(w_break_by_x, vectors[root_index])
        root_supported = abs(root_dot) > tolerance_dot_nonzero

        if not root_supported:
            unsupported_roots.append(root_pair_name)
            rooted_results[root_pair_name] = {
                "root_pair": root_pair_name,
                "dot_w_break_u_root": float(root_dot),
                "root_supported": False,
            }
            continue

        root_sign = sign(root_dot)
        signed_root_vector = [float(root_sign) * value for value in vectors[root_index]]

        sign_vector_by_pair: dict[str, int] = {root_pair_name: int(root_sign)}
        coords_star_by_pair: dict[str, list[float]] = {
            root_pair_name: [float(root_sign) * coords[root_index][0], float(root_sign) * coords[root_index][1]]
        }
        pair_audits: dict[str, Any] = {}
        alignment_ok = True

        for target_index in range(1, 6):
            if target_index == root_index:
                pair_audits[root_pair_name] = {
                    "pair": root_pair_name,
                    "root_chart": True,
                    "dot_w_break_u_root": float(root_dot),
                    "root_sign": int(root_sign),
                }
                continue

            pair_key = (min(root_index, target_index), max(root_index, target_index))
            operator_path = transition_refs[pair_key]
            transport_matrix = transport_matrices[pair_key]
            transported_vector = matvec12(transport_matrix, signed_root_vector)
            target_dot = dot(transported_vector, vectors[target_index])
            if abs(target_dot) <= tolerance_dot_nonzero:
                alignment_ok = False
                pair_audits[pair_name(target_index)] = {
                    "pair": pair_name(target_index),
                    "transport_ref": str(operator_path.relative_to(REPO)),
                    "dot_transport_root_with_u_target": float(target_dot),
                    "supported": False,
                }
                continue

            target_sign = sign(target_dot)
            sign_vector_by_pair[pair_name(target_index)] = int(target_sign)
            coords_star_by_pair[pair_name(target_index)] = [
                float(target_sign) * coords[target_index][0],
                float(target_sign) * coords[target_index][1],
            ]

            alignment_residual = l2_norm(
                [
                    transported_vector[site_index] - float(target_sign) * vectors[target_index][site_index]
                    for site_index in range(12)
                ]
            )
            if alignment_residual > tolerance_alignment:
                alignment_ok = False

            pair_audits[pair_name(target_index)] = {
                "pair": pair_name(target_index),
                "transport_ref": str(operator_path.relative_to(REPO)),
                "dot_transport_root_with_u_target": float(target_dot),
                "propagated_sign": int(target_sign),
                "alignment_l2_residual": float(alignment_residual),
                "supported": True,
            }

        output_vectors_by_pair: dict[str, list[float]] = {}
        output_signs_by_pair: dict[str, str | None] = {}
        for pair_index in range(1, 6):
            current_pair_name = pair_name(pair_index)
            y_chart = ((y_sel.get("charts") or {}).get(current_pair_name) or {})
            y_matrix = y_chart.get("Y_sel_matrix_in_c_s_to_o")
            if not (isinstance(y_matrix, list) and len(y_matrix) == 2 and all(is_numeric_list_len(row, 2) for row in y_matrix)):
                raise SystemExit(
                    f"Invalid {IN_Y_SEL.relative_to(REPO)}: charts.{current_pair_name}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
                )
            if current_pair_name not in coords_star_by_pair:
                alignment_ok = False
                continue

            matrix_2x2 = [
                [float(y_matrix[0][0]), float(y_matrix[0][1])],
                [float(y_matrix[1][0]), float(y_matrix[1][1])],
            ]
            output_vector = matvec2(matrix_2x2, coords_star_by_pair[current_pair_name])
            output_vectors_by_pair[current_pair_name] = output_vector
            output_signs_by_pair[current_pair_name] = (
                None if abs(float(output_vector[0])) <= tolerance_output else ("+" if float(output_vector[0]) > 0 else "-")
            )

        root_result_supported = alignment_ok and len(sign_vector_by_pair) == 5 and len(output_vectors_by_pair) == 5
        if root_result_supported:
            supported_roots.append(root_pair_name)
        else:
            unsupported_roots.append(root_pair_name)

        rooted_results[root_pair_name] = {
            "root_pair": root_pair_name,
            "dot_w_break_u_root": float(root_dot),
            "root_sign": int(root_sign),
            "root_supported": root_result_supported,
            "transport_alignment_ok": alignment_ok,
            "sign_vector_by_pair": sign_vector_by_pair,
            "output_vectors_by_pair": output_vectors_by_pair,
            "output_o_plus_sign_by_pair": output_signs_by_pair,
            "pair_audits": pair_audits,
        }

    all_roots_supported = len(supported_roots) == 5

    supported_roots_sign_vector_agree = False
    supported_roots_output_vectors_agree = False
    root_independent_sign_vector = False
    root_independent_output_vectors = False
    max_output_difference_across_roots = None
    reference_root = supported_roots[0] if supported_roots else None

    if reference_root is not None:
        reference_sign_vector = rooted_results[reference_root]["sign_vector_by_pair"]
        reference_output_vectors = rooted_results[reference_root]["output_vectors_by_pair"]

        supported_roots_sign_vector_agree = all(
            rooted_results[root_pair]["sign_vector_by_pair"] == reference_sign_vector for root_pair in supported_roots
        )

        max_output_difference_across_roots = 0.0
        for root_pair in supported_roots:
            current_output_vectors = rooted_results[root_pair]["output_vectors_by_pair"]
            for current_pair_name in sorted(reference_output_vectors):
                current_diff = max_abs_pair_diff(current_output_vectors[current_pair_name], reference_output_vectors[current_pair_name])
                max_output_difference_across_roots = max(float(max_output_difference_across_roots), float(current_diff))
        supported_roots_output_vectors_agree = bool(float(max_output_difference_across_roots) <= tolerance_output)

        root_independent_sign_vector = bool(all_roots_supported and supported_roots_sign_vector_agree)
        root_independent_output_vectors = bool(all_roots_supported and supported_roots_output_vectors_agree)

    same_as_f684_pair1_root = None
    if isinstance(f684_summary, dict):
        f684_signs = f684_summary.get("signs_by_pair")
        if isinstance(f684_signs, dict) and "pair1" in rooted_results and rooted_results["pair1"].get("root_supported"):
            same_as_f684_pair1_root = rooted_results["pair1"]["sign_vector_by_pair"] == f684_signs

    if all_roots_supported and root_independent_sign_vector and root_independent_output_vectors:
        status = "PASS_MULTIROOT_ROOTED_SIGN_LIFT_ROOT_INDEPENDENT_CONVENTION_LAYER_ONLY"
    elif len(supported_roots) >= 2 and supported_roots_sign_vector_agree and supported_roots_output_vectors_agree:
        status = "PARTIAL_SUPPORTED_ROOT_CORRIDOR_ONLY_WITH_AGREEING_SECTION"
    elif not all_roots_supported:
        status = "FAIL_NOT_ALL_ROOTS_SUPPORT_ROOTED_SIGN_LIFT"
    elif not root_independent_sign_vector:
        status = "FAIL_MULTIROOT_SIGN_VECTOR_DEPENDENCE_DETECTED"
    elif not root_independent_output_vectors:
        status = "FAIL_MULTIROOT_OUTPUT_VECTOR_DEPENDENCE_DETECTED"

    artifact = {
        "stage": "P713",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the currently surviving rooted directed-section route is genuinely pair1-root dependent, "
            "or whether every chart root recovers the same sign vector and output-space section once the local root sign is fixed by w_break "
            "and propagated through the exported oriented transport family"
        ),
        "inputs": {
            "F647": str(IN_F647.relative_to(REPO)),
            "F467_atlas": str(IN_ATLAS.relative_to(REPO)),
            "Y_sel": str(IN_Y_SEL.relative_to(REPO)),
            "F684_summary": str(IN_F684.relative_to(REPO)),
            "A_m_refs": {str(pair_index): str(IN_A[pair_index].relative_to(REPO)) for pair_index in sorted(IN_A)},
        },
        "tolerances": {
            "dot_nonzero_abs_tol": tolerance_dot_nonzero,
            "alignment_l2_tol": tolerance_alignment,
            "output_max_abs_diff_tol": tolerance_output,
        },
        "rooted_results": rooted_results,
        "result": {
            "all_roots_supported": all_roots_supported,
            "supported_roots": supported_roots,
            "unsupported_roots": unsupported_roots,
            "reference_root": reference_root,
            "supported_roots_sign_vector_agree": supported_roots_sign_vector_agree,
            "supported_roots_output_vectors_agree": supported_roots_output_vectors_agree,
            "root_independent_sign_vector": root_independent_sign_vector,
            "root_independent_output_vectors": root_independent_output_vectors,
            "max_output_difference_across_roots": max_output_difference_across_roots,
            "pair1_root_matches_existing_F684_sign_vector": same_as_f684_pair1_root,
            "counts_as_strict_physical_orientation_datum": False,
            "implies_t176_discharge": False,
            "surviving_candidate_class": (
                "root-independent convention-layer directed section from strict w_break anchors plus exported oriented transport"
                if root_independent_sign_vector and root_independent_output_vectors
                else (
                    "supported-root corridor convention-layer section from strict w_break anchors plus exported oriented transport"
                    if len(supported_roots) >= 2 and supported_roots_sign_vector_agree and supported_roots_output_vectors_agree
                    else None
                )
            ),
            "remaining_gap_after_positive_result": (
                "strict-core/global upgrade of the transport/provider layer beyond convention/lane-scoped status"
                if root_independent_sign_vector and root_independent_output_vectors
                else (
                    "current w_break payload anchors only a supported-root corridor, not all chart roots; no full global provider anchor exists"
                    if len(supported_roots) >= 2 and supported_roots_sign_vector_agree and supported_roots_output_vectors_agree
                    else "root dependence remains unresolved"
                )
            ),
        },
        "hard_limits": [
            "No strict-core directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No T176 discharge claim.",
            "No promotion of convention-layer oriented transport into strict-core/global provider.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P713",
        "status": status,
        "as_of": AS_OF,
        "all_roots_supported": all_roots_supported,
        "supported_roots": supported_roots,
        "unsupported_roots": unsupported_roots,
        "supported_roots_sign_vector_agree": supported_roots_sign_vector_agree,
        "supported_roots_output_vectors_agree": supported_roots_output_vectors_agree,
        "root_independent_sign_vector": root_independent_sign_vector,
        "root_independent_output_vectors": root_independent_output_vectors,
        "pair1_root_matches_existing_F684_sign_vector": same_as_f684_pair1_root,
        "counts_as_strict_physical_orientation_datum": False,
        "implies_t176_discharge": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
