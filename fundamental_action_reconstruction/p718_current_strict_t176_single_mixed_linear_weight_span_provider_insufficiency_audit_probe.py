#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"
IN_Y_SEL = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json"
IN_P716 = GENERATED / "p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe_summary.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p718_current_strict_t176_single_mixed_linear_weight_span_provider_insufficiency_audit_probe.json"
OUT_SUMMARY = (
    GENERATED / "p718_current_strict_t176_single_mixed_linear_weight_span_provider_insufficiency_audit_probe_summary.json"
)


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


def unique_sorted_angles(values: list[float], tolerance: float) -> list[float]:
    ordered = sorted(float(value) % (2.0 * math.pi) for value in values)
    result: list[float] = []
    for value in ordered:
        if not result or abs(value - result[-1]) > tolerance:
            result.append(value)
    return result


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F647, IN_ATLAS, IN_Y_SEL, IN_P715, IN_P716] + [IN_A[pair_index] for pair_index in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P718",
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
    p715 = load_json(IN_P715)
    p716 = load_json(IN_P716)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break = payload.get("w_break_by_x")
    w_ref = payload.get("w_ref_unnormalized_by_x")
    if not is_numeric_list_len(w_break, 12):
        artifact = {
            "stage": "P718",
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
            "stage": "P718",
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
    tolerance_angle = 1e-12

    root_dot_profiles: dict[str, dict[str, float]] = {}
    critical_angles: list[float] = []
    for pair_index in range(1, 6):
        pair = pair_name(pair_index)
        dot_break = float(dot(w_break_by_x, vectors[pair_index]))
        dot_ref = float(dot(w_ref_by_x, vectors[pair_index]))
        root_dot_profiles[pair] = {
            "dot_w_break_u_root": dot_break,
            "dot_w_ref_u_root": dot_ref,
        }
        zero_angle = math.atan2(-dot_break, dot_ref)
        critical_angles.extend([zero_angle, zero_angle + math.pi])

    unique_angles = unique_sorted_angles(critical_angles, tolerance_angle)
    interior_angles = [
        angle
        for angle in unique_angles
        if angle > tolerance_angle and (2.0 * math.pi - angle) > tolerance_angle
    ]
    sector_boundaries = [0.0] + interior_angles + [2.0 * math.pi]

    sectors: list[dict[str, Any]] = []

    for sector_index, (left_boundary, right_boundary) in enumerate(zip(sector_boundaries, sector_boundaries[1:]), start=1):
        if right_boundary - left_boundary <= tolerance_angle:
            continue

        theta = 0.5 * (left_boundary + right_boundary)
        coeff_break = float(math.cos(theta))
        coeff_ref = float(math.sin(theta))
        mixed_weight = [
            coeff_break * w_break_by_x[site_index] + coeff_ref * w_ref_by_x[site_index]
            for site_index in range(12)
        ]

        rooted_results: dict[str, Any] = {}
        supported_roots: list[str] = []
        unsupported_roots: list[str] = []

        for root_index in range(1, 6):
            root_pair_name = pair_name(root_index)
            root_dot = dot(mixed_weight, vectors[root_index])
            root_supported = abs(root_dot) > tolerance_dot_nonzero

            if not root_supported:
                unsupported_roots.append(root_pair_name)
                rooted_results[root_pair_name] = {
                    "root_pair": root_pair_name,
                    "root_supported": False,
                    "dot_w_mix_u_root": float(root_dot),
                }
                continue

            root_sign = sign(root_dot)
            signed_root_vector = [float(root_sign) * value for value in vectors[root_index]]

            sign_vector_by_pair: dict[str, int] = {root_pair_name: int(root_sign)}
            coords_star_by_pair: dict[str, list[float]] = {
                root_pair_name: [float(root_sign) * coords[root_index][0], float(root_sign) * coords[root_index][1]]
            }
            alignment_ok = True

            for target_index in range(1, 6):
                if target_index == root_index:
                    continue

                pair_key = (min(root_index, target_index), max(root_index, target_index))
                transport_matrix = transport_matrices[pair_key]
                transported_vector = matvec12(transport_matrix, signed_root_vector)
                target_dot = dot(transported_vector, vectors[target_index])
                if abs(target_dot) <= tolerance_dot_nonzero:
                    alignment_ok = False
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

            output_vectors_by_pair: dict[str, list[float]] = {}
            for pair_index in range(1, 6):
                current_pair_name = pair_name(pair_index)
                y_chart = ((y_sel.get("charts") or {}).get(current_pair_name) or {})
                y_matrix = y_chart.get("Y_sel_matrix_in_c_s_to_o")
                if not (
                    isinstance(y_matrix, list)
                    and len(y_matrix) == 2
                    and all(is_numeric_list_len(row, 2) for row in y_matrix)
                ):
                    raise SystemExit(
                        f"Invalid {IN_Y_SEL.relative_to(REPO)}: charts.{current_pair_name}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
                    )
                matrix_2x2 = [
                    [float(y_matrix[0][0]), float(y_matrix[0][1])],
                    [float(y_matrix[1][0]), float(y_matrix[1][1])],
                ]
                output_vectors_by_pair[current_pair_name] = matvec2(matrix_2x2, coords_star_by_pair[current_pair_name])

            root_result_supported = alignment_ok and len(sign_vector_by_pair) == 5 and len(output_vectors_by_pair) == 5
            if root_result_supported:
                supported_roots.append(root_pair_name)
            else:
                unsupported_roots.append(root_pair_name)

            rooted_results[root_pair_name] = {
                "root_pair": root_pair_name,
                "root_supported": root_result_supported,
                "dot_w_mix_u_root": float(root_dot),
                "root_sign": int(root_sign),
                "sign_vector_by_pair": sign_vector_by_pair,
                "output_vectors_by_pair": output_vectors_by_pair,
            }

        all_roots_supported = len(supported_roots) == 5
        exact_root_independent_sign_vector = False
        exact_root_independent_output_vectors = False
        projective_root_independent_sign_orbit = False
        projective_root_independent_output_orbit = False
        same_orbit_roots_relative_to_reference: list[str] = []
        negated_orbit_roots_relative_to_reference: list[str] = []
        other_orbit_roots_relative_to_reference: list[str] = []
        reference_root = supported_roots[0] if supported_roots else None

        if reference_root is not None:
            reference_sign_vector = rooted_results[reference_root]["sign_vector_by_pair"]
            reference_output_vectors = rooted_results[reference_root]["output_vectors_by_pair"]

            exact_root_independent_sign_vector = all_roots_supported and all(
                rooted_results[root_pair]["sign_vector_by_pair"] == reference_sign_vector for root_pair in supported_roots
            )
            exact_root_independent_output_vectors = all_roots_supported and all(
                all(
                    max_abs_pair_diff(
                        rooted_results[root_pair]["output_vectors_by_pair"][current_pair_name],
                        reference_output_vectors[current_pair_name],
                    )
                    <= tolerance_output
                    for current_pair_name in sorted(reference_output_vectors)
                )
                for root_pair in supported_roots
            )

            sign_relations: list[str] = []
            output_relations: list[str] = []
            for root_pair in supported_roots:
                current_sign_vector = rooted_results[root_pair]["sign_vector_by_pair"]
                current_output_vectors = rooted_results[root_pair]["output_vectors_by_pair"]
                sign_relation = sign_vector_relation(reference_sign_vector, current_sign_vector)
                output_relation, _, _ = output_vector_relation(reference_output_vectors, current_output_vectors, tolerance_output)
                sign_relations.append(sign_relation)
                output_relations.append(output_relation)
                if sign_relation == "same":
                    same_orbit_roots_relative_to_reference.append(root_pair)
                elif sign_relation == "global_negation":
                    negated_orbit_roots_relative_to_reference.append(root_pair)
                else:
                    other_orbit_roots_relative_to_reference.append(root_pair)

            projective_root_independent_sign_orbit = all_roots_supported and all(
                relation in {"same", "global_negation"} for relation in sign_relations
            )
            projective_root_independent_output_orbit = all_roots_supported and all(
                relation in {"same", "global_negation"} for relation in output_relations
            )

        if exact_root_independent_sign_vector and exact_root_independent_output_vectors:
            sector_status = "EXACT_ROOT_INDEPENDENT_SECTION_FOUND"
        elif projective_root_independent_sign_orbit and projective_root_independent_output_orbit:
            sector_status = "PROJECTIVE_ORBIT_ONLY"
        elif all_roots_supported:
            sector_status = "ALL_ROOTS_SUPPORTED_BUT_NOT_PROJECTIVE_ORBIT"
        else:
            sector_status = "NOT_ALL_ROOTS_SUPPORTED"

        sector_sign_pattern = {
            pair_name(pair_index): sign(dot(mixed_weight, vectors[pair_index]))
            for pair_index in range(1, 6)
            if abs(dot(mixed_weight, vectors[pair_index])) > tolerance_dot_nonzero
        }

        sectors.append(
            {
                "sector_index": int(sector_index),
                "theta_midpoint_radians": float(theta),
                "coefficients_on_unit_circle": {
                    "coeff_w_break": coeff_break,
                    "coeff_w_ref_unnormalized": coeff_ref,
                },
                "root_sign_pattern": sector_sign_pattern,
                "all_roots_supported": bool(all_roots_supported),
                "exact_root_independent_sign_vector": bool(exact_root_independent_sign_vector),
                "exact_root_independent_output_vectors": bool(exact_root_independent_output_vectors),
                "projective_root_independent_sign_orbit": bool(projective_root_independent_sign_orbit),
                "projective_root_independent_output_orbit": bool(projective_root_independent_output_orbit),
                "reference_root": reference_root,
                "same_orbit_roots_relative_to_reference": same_orbit_roots_relative_to_reference,
                "negated_orbit_roots_relative_to_reference": negated_orbit_roots_relative_to_reference,
                "other_orbit_roots_relative_to_reference": other_orbit_roots_relative_to_reference,
                "status": sector_status,
            }
        )

    exact_sectors = [sector for sector in sectors if sector["status"] == "EXACT_ROOT_INDEPENDENT_SECTION_FOUND"]
    projective_only_sectors = [sector for sector in sectors if sector["status"] == "PROJECTIVE_ORBIT_ONLY"]

    unique_negated_root_sets = sorted(
        {
            tuple(str(root) for root in sector.get("negated_orbit_roots_relative_to_reference") or [])
            for sector in projective_only_sectors
        }
    )

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    add_check(
        "p715_current_best_route_still_projective_only",
        bool(p715.get("result", {}).get("projective_root_independent_sign_orbit"))
        and not bool(p715.get("result", {}).get("exact_root_independent_sign_vector")),
        True,
        "P715 still marks the current best route as projective-orbit-only, not exact.",
    )
    add_check(
        "p716_pair4_localization_still_true",
        bool(p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")),
        True,
        "P716 still localizes the current best exact split to pair4 negative cosine polarity.",
    )
    add_check(
        "single_linear_span_has_no_exact_sector",
        len(exact_sectors) > 0,
        False,
        "No noncritical sector in span{w_break, w_ref_unnormalized} yields one exact root-independent directed section.",
    )
    add_check(
        "single_linear_span_has_projective_only_sector",
        len(projective_only_sectors) > 0,
        True,
        "The same single-linear span does contain sectors with all-root support and projective-orbit recovery only.",
    )
    add_check(
        "current_pair4_only_split_reappears_inside_single_linear_span",
        any(tuple(sector.get("negated_orbit_roots_relative_to_reference") or []) == ("pair4",) for sector in projective_only_sectors),
        True,
        "One projective-only sector reproduces the current pair4-only negated branch.",
    )
    add_check(
        "single_linear_span_also_contains_worse_pair2_pair3_split_sector",
        any(tuple(sector.get("negated_orbit_roots_relative_to_reference") or []) == ("pair2", "pair3") for sector in projective_only_sectors),
        True,
        "Another projective-only sector worsens the exact split to pair2/pair3 rather than fixing it.",
    )

    status = (
        "PASS_SINGLE_MIXED_LINEAR_WEIGHT_SPAN_PROVIDER_INSUFFICIENCY_AUDITED"
        if not blocking
        else "P718_REQUIRES_REVIEW_CHANGED_SINGLE_LINEAR_SPAN_PROVIDER_STATE"
    )

    artifact = {
        "stage": "P718",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "test the narrowest honest next provider attack after P717: whether any single mixed linear weight "
            "in span{w_break, w_ref_unnormalized} can discharge the remaining exact branch split below T176"
        ),
        "inputs": {
            "F647": str(IN_F647.relative_to(REPO)),
            "atlas": str(IN_ATLAS.relative_to(REPO)),
            "Y_sel": str(IN_Y_SEL.relative_to(REPO)),
            "P715": str(IN_P715.relative_to(REPO)),
            "P716": str(IN_P716.relative_to(REPO)),
            "A_packets": [str(IN_A[pair_index].relative_to(REPO)) for pair_index in sorted(IN_A)],
        },
        "root_dot_profiles": root_dot_profiles,
        "critical_angles_radians_mod_2pi": unique_angles,
        "tested_sectors": sectors,
        "checks": checks,
        "blocking_mismatches": blocking,
        "result": {
            "single_mixed_linear_weight_span_exact_root_independent_section_exists": len(exact_sectors) > 0,
            "single_mixed_linear_weight_span_projective_orbit_only_sector_exists": len(projective_only_sectors) > 0,
            "single_mixed_linear_weight_span_all_root_support_available": any(
                bool(sector.get("all_roots_supported")) for sector in sectors
            ),
            "projective_only_negated_root_sets_seen": [list(roots) for roots in unique_negated_root_sets],
            "pair4_only_split_sector_exists": any(
                tuple(sector.get("negated_orbit_roots_relative_to_reference") or []) == ("pair4",)
                for sector in projective_only_sectors
            ),
            "pair2_pair3_split_sector_exists": any(
                tuple(sector.get("negated_orbit_roots_relative_to_reference") or []) == ("pair2", "pair3")
                for sector in projective_only_sectors
            ),
            "best_current_sector_is_still_projective_only_not_exact": True,
            "next_provider_must_leave_single_linear_span_of_current_exported_weights": len(exact_sectors) == 0,
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No strict-core directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No claim against all future nonlinear or higher-class providers.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P718",
        "status": status,
        "as_of": AS_OF,
        "single_mixed_linear_weight_span_exact_root_independent_section_exists": len(exact_sectors) > 0,
        "single_mixed_linear_weight_span_projective_orbit_only_sector_exists": len(projective_only_sectors) > 0,
        "projective_only_negated_root_sets_seen": [list(roots) for roots in unique_negated_root_sets],
        "next_provider_must_leave_single_linear_span_of_current_exported_weights": len(exact_sectors) == 0,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
