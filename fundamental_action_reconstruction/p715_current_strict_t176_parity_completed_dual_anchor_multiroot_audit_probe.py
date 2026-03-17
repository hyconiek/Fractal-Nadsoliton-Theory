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
IN_P714 = GENERATED / "p714_current_strict_t176_w_break_parity_root_support_profile_audit_probe_summary.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json"
OUT_SUMMARY = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe_summary.json"


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


def negate_sign_vector(sign_vector: dict[str, int]) -> dict[str, int]:
    return {pair: int(-value) for pair, value in sign_vector.items()}


def sign_vector_relation(reference: dict[str, int], current: dict[str, int]) -> str:
    if current == reference:
        return "same"
    if current == negate_sign_vector(reference):
        return "global_negation"
    return "other"


def output_vector_relation(reference: dict[str, list[float]], current: dict[str, list[float]], tolerance: float) -> tuple[str, float, float]:
    same_diff = 0.0
    negated_diff = 0.0
    for pair in sorted(reference):
        same_diff = max(same_diff, max_abs_pair_diff(current[pair], reference[pair]))
        negated_diff = max(
            negated_diff,
            max_abs_pair_diff(current[pair], [-float(reference[pair][0]), -float(reference[pair][1])]),
        )
    if same_diff <= tolerance:
        return "same", float(same_diff), float(negated_diff)
    if negated_diff <= tolerance:
        return "global_negation", float(same_diff), float(negated_diff)
    return "other", float(same_diff), float(negated_diff)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F647, IN_ATLAS, IN_Y_SEL, IN_P714] + [IN_A[pair_index] for pair_index in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P715",
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
    p714 = load_json(IN_P714)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break = payload.get("w_break_by_x")
    w_ref = payload.get("w_ref_unnormalized_by_x")
    if not is_numeric_list_len(w_break, 12):
        artifact = {
            "stage": "P715",
            "status": "INVALID_W_BREAK_INPUT_SHAPE",
            "as_of": AS_OF,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return
    if not is_numeric_list_len(w_ref, 12):
        artifact = {
            "stage": "P715",
            "status": "INVALID_W_REF_INPUT_SHAPE",
            "as_of": AS_OF,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    w_break_by_x = [float(entry) for entry in w_break]
    w_ref_by_x = [float(entry) for entry in w_ref]

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
    anchor_source_by_root: dict[str, str | None] = {}

    for root_index in range(1, 6):
        root_pair_name = pair_name(root_index)
        root_dot_break = dot(w_break_by_x, vectors[root_index])
        root_dot_ref = dot(w_ref_by_x, vectors[root_index])

        if abs(root_dot_break) > tolerance_dot_nonzero:
            root_supported = True
            root_sign = sign(root_dot_break)
            root_anchor_source = "w_break"
            root_anchor_scalar = float(root_dot_break)
        elif abs(root_dot_ref) > tolerance_dot_nonzero:
            root_supported = True
            root_sign = sign(root_dot_ref)
            root_anchor_source = "w_ref_unnormalized"
            root_anchor_scalar = float(root_dot_ref)
        else:
            root_supported = False
            root_sign = None
            root_anchor_source = None
            root_anchor_scalar = None

        anchor_source_by_root[root_pair_name] = root_anchor_source

        if not root_supported:
            unsupported_roots.append(root_pair_name)
            rooted_results[root_pair_name] = {
                "root_pair": root_pair_name,
                "dot_w_break_u_root": float(root_dot_break),
                "dot_w_ref_u_root": float(root_dot_ref),
                "root_anchor_source": root_anchor_source,
                "root_anchor_scalar": root_anchor_scalar,
                "root_supported": False,
            }
            continue

        signed_root_vector = [float(root_sign) * value for value in vectors[root_index]]

        sign_vector_by_pair: dict[str, int] = {root_pair_name: int(root_sign)}
        coords_star_by_pair: dict[str, list[float]] = {
            root_pair_name: [float(root_sign) * coords[root_index][0], float(root_sign) * coords[root_index][1]]
        }
        pair_audits: dict[str, Any] = {}
        alignment_ok = True

        for target_index in range(1, 6):
            target_pair_name = pair_name(target_index)
            if target_index == root_index:
                pair_audits[root_pair_name] = {
                    "pair": root_pair_name,
                    "root_chart": True,
                    "dot_w_break_u_root": float(root_dot_break),
                    "dot_w_ref_u_root": float(root_dot_ref),
                    "root_anchor_source": root_anchor_source,
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
                pair_audits[target_pair_name] = {
                    "pair": target_pair_name,
                    "transport_ref": str(operator_path.relative_to(REPO)),
                    "dot_transport_root_with_u_target": float(target_dot),
                    "supported": False,
                }
                continue

            target_sign = sign(target_dot)
            sign_vector_by_pair[target_pair_name] = int(target_sign)
            coords_star_by_pair[target_pair_name] = [
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

            pair_audits[target_pair_name] = {
                "pair": target_pair_name,
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
            "dot_w_break_u_root": float(root_dot_break),
            "dot_w_ref_u_root": float(root_dot_ref),
            "root_anchor_source": root_anchor_source,
            "root_anchor_scalar": root_anchor_scalar,
            "root_sign": int(root_sign),
            "root_supported": root_result_supported,
            "transport_alignment_ok": alignment_ok,
            "sign_vector_by_pair": sign_vector_by_pair,
            "output_vectors_by_pair": output_vectors_by_pair,
            "output_o_plus_sign_by_pair": output_signs_by_pair,
            "pair_audits": pair_audits,
        }

    all_roots_supported = len(supported_roots) == 5

    exact_sign_vector_agree = False
    exact_output_vectors_agree = False
    projective_sign_orbit_agree = False
    projective_output_orbit_agree = False
    exact_root_independent_sign_vector = False
    exact_root_independent_output_vectors = False
    projective_root_independent_sign_orbit = False
    projective_root_independent_output_orbit = False
    reference_root = supported_roots[0] if supported_roots else None
    sign_vector_orbit_relation_by_root: dict[str, str] = {}
    output_vector_orbit_relation_by_root: dict[str, str] = {}
    same_orbit_roots: list[str] = []
    negated_orbit_roots: list[str] = []
    other_orbit_roots: list[str] = []
    max_output_difference_same_orbit = None
    max_output_difference_negated_orbit = None

    if reference_root is not None:
        reference_sign_vector = rooted_results[reference_root]["sign_vector_by_pair"]
        reference_output_vectors = rooted_results[reference_root]["output_vectors_by_pair"]

        exact_sign_vector_agree = all(
            rooted_results[root_pair]["sign_vector_by_pair"] == reference_sign_vector for root_pair in supported_roots
        )

        exact_output_vectors_agree = True
        max_output_difference_same_orbit = 0.0
        max_output_difference_negated_orbit = 0.0
        for root_pair in supported_roots:
            current_output_vectors = rooted_results[root_pair]["output_vectors_by_pair"]
            relation, same_diff, neg_diff = output_vector_relation(reference_output_vectors, current_output_vectors, tolerance_output)
            output_vector_orbit_relation_by_root[root_pair] = relation
            max_output_difference_same_orbit = max(float(max_output_difference_same_orbit), float(same_diff))
            max_output_difference_negated_orbit = max(float(max_output_difference_negated_orbit), float(neg_diff))
            if relation not in {"same", "global_negation"}:
                exact_output_vectors_agree = False

        if exact_output_vectors_agree:
            exact_output_vectors_agree = all(
                output_vector_orbit_relation_by_root[root_pair] == "same" for root_pair in supported_roots
            )

        for root_pair in supported_roots:
            relation = sign_vector_relation(reference_sign_vector, rooted_results[root_pair]["sign_vector_by_pair"])
            sign_vector_orbit_relation_by_root[root_pair] = relation
            if relation == "same":
                same_orbit_roots.append(root_pair)
            elif relation == "global_negation":
                negated_orbit_roots.append(root_pair)
            else:
                other_orbit_roots.append(root_pair)

        projective_sign_orbit_agree = len(other_orbit_roots) == 0
        projective_output_orbit_agree = all(
            relation in {"same", "global_negation"} for relation in output_vector_orbit_relation_by_root.values()
        )

        exact_root_independent_sign_vector = bool(all_roots_supported and exact_sign_vector_agree)
        exact_root_independent_output_vectors = bool(all_roots_supported and exact_output_vectors_agree)
        projective_root_independent_sign_orbit = bool(all_roots_supported and projective_sign_orbit_agree)
        projective_root_independent_output_orbit = bool(all_roots_supported and projective_output_orbit_agree)

    if all_roots_supported and exact_root_independent_sign_vector and exact_root_independent_output_vectors:
        status = "PASS_DUAL_ANCHOR_ALL_ROOT_SUPPORT_WITH_EXACT_SECTION_CONVENTION_LAYER_ONLY"
    elif all_roots_supported and projective_root_independent_sign_orbit and projective_root_independent_output_orbit:
        status = "PARTIAL_DUAL_ANCHOR_ALL_ROOT_SUPPORT_PROJECTIVE_ORBIT_ONLY"
    elif not all_roots_supported:
        status = "FAIL_DUAL_ANCHOR_NOT_ALL_ROOTS_SUPPORTED"
    elif not projective_root_independent_sign_orbit:
        status = "FAIL_DUAL_ANCHOR_PROJECTIVE_SIGN_ORBIT_DEPENDENCE_DETECTED"
    else:
        status = "FAIL_DUAL_ANCHOR_PROJECTIVE_OUTPUT_ORBIT_DEPENDENCE_DETECTED"

    dual_anchor_adds_even_component = bool(p714.get("current_w_break_explains_supported_root_corridor"))

    artifact = {
        "stage": "P715",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the minimal parity-completed strict-core dual-anchor rule "
            "(w_break first, w_ref_unnormalized fallback) removes the root-support gap and, if so, "
            "whether the recovered rooted section is exact or only projective up to one global Z2 sign"
        ),
        "inputs": {
            "F647": str(IN_F647.relative_to(REPO)),
            "F467_atlas": str(IN_ATLAS.relative_to(REPO)),
            "Y_sel": str(IN_Y_SEL.relative_to(REPO)),
            "P714_summary": str(IN_P714.relative_to(REPO)),
            "A_m_refs": {str(pair_index): str(IN_A[pair_index].relative_to(REPO)) for pair_index in sorted(IN_A)},
        },
        "root_anchor_rule": {
            "rule_name": "w_break_first_then_w_ref_unnormalized_fallback",
            "definition": (
                "sign(dot(w_break,u_r)) if nonzero; else sign(dot(w_ref_unnormalized,u_r)) if nonzero; else undefined"
            ),
            "uses_only_exported_strict_core_payloads": True,
            "introduces_hidden_chart_slot": False,
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
            "anchor_source_by_root": anchor_source_by_root,
            "dual_anchor_adds_even_component_beyond_P714_odd_corridor": dual_anchor_adds_even_component,
            "exact_sign_vector_agree_across_roots": exact_sign_vector_agree,
            "exact_output_vectors_agree_across_roots": exact_output_vectors_agree,
            "exact_root_independent_sign_vector": exact_root_independent_sign_vector,
            "exact_root_independent_output_vectors": exact_root_independent_output_vectors,
            "projective_sign_orbit_agree_across_roots": projective_sign_orbit_agree,
            "projective_output_orbit_agree_across_roots": projective_output_orbit_agree,
            "projective_root_independent_sign_orbit": projective_root_independent_sign_orbit,
            "projective_root_independent_output_orbit": projective_root_independent_output_orbit,
            "sign_vector_orbit_relation_by_root": sign_vector_orbit_relation_by_root,
            "output_vector_orbit_relation_by_root": output_vector_orbit_relation_by_root,
            "same_orbit_roots_relative_to_reference": same_orbit_roots,
            "negated_orbit_roots_relative_to_reference": negated_orbit_roots,
            "other_orbit_roots_relative_to_reference": other_orbit_roots,
            "max_output_difference_same_orbit": max_output_difference_same_orbit,
            "max_output_difference_negated_orbit": max_output_difference_negated_orbit,
            "all_roots_supported_but_exact_section_split_by_global_z2": bool(
                all_roots_supported
                and projective_root_independent_sign_orbit
                and projective_root_independent_output_orbit
                and not exact_root_independent_sign_vector
                and not exact_root_independent_output_vectors
            ),
            "counts_as_strict_physical_orientation_datum": False,
            "implies_t176_discharge": False,
            "surviving_candidate_class": (
                "all-root supported convention-layer section class defined only up to one global sign orbit"
                if all_roots_supported and projective_root_independent_sign_orbit and projective_root_independent_output_orbit
                else (
                    "exact all-root convention-layer section from the parity-completed dual-anchor rule"
                    if all_roots_supported and exact_root_independent_sign_vector and exact_root_independent_output_vectors
                    else None
                )
            ),
            "remaining_gap_after_positive_result": (
                "residual global Z2 orbit split remains: the exact directed section is not root-independent even though the projective orbit is"
                if all_roots_supported and projective_root_independent_sign_orbit and projective_root_independent_output_orbit and not exact_root_independent_sign_vector
                else (
                    "strict-core/global upgrade of the transport/provider layer beyond convention/lane-scoped status"
                    if all_roots_supported and exact_root_independent_sign_vector and exact_root_independent_output_vectors
                    else "dual-anchor route still fails before an all-root projective orbit is recovered"
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
        "stage": "P715",
        "status": status,
        "as_of": AS_OF,
        "all_roots_supported": all_roots_supported,
        "supported_roots": supported_roots,
        "unsupported_roots": unsupported_roots,
        "anchor_source_by_root": anchor_source_by_root,
        "exact_root_independent_sign_vector": exact_root_independent_sign_vector,
        "exact_root_independent_output_vectors": exact_root_independent_output_vectors,
        "projective_root_independent_sign_orbit": projective_root_independent_sign_orbit,
        "projective_root_independent_output_orbit": projective_root_independent_output_orbit,
        "same_orbit_roots_relative_to_reference": same_orbit_roots,
        "negated_orbit_roots_relative_to_reference": negated_orbit_roots,
        "counts_as_strict_physical_orientation_datum": False,
        "implies_t176_discharge": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
