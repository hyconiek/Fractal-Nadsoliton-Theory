#!/usr/bin/env python3
from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"
IN_Y_SEL = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
IN_P697 = GENERATED / "p697_current_strict_projective_observer_limit_readout_from_global_projective_selector_closure_output_probe_summary.json"
IN_P719 = GENERATED / "p719_current_strict_t176_low_complexity_odd_polynomial_two_readout_provider_class_audit_probe_summary.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p720_current_strict_t176_observer_facing_signed_output_channel_projection_provider_class_audit_probe.json"
OUT_SUMMARY = GENERATED / "p720_current_strict_t176_observer_facing_signed_output_channel_projection_provider_class_audit_probe_summary.json"


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


def output_vector_relation(reference: dict[str, list[float]], current: dict[str, list[float]], tolerance: float) -> str:
    same_diff = 0.0
    negated_diff = 0.0
    for pair in sorted(reference):
        same_diff = max(same_diff, max_abs_pair_diff(current[pair], reference[pair]))
        negated_diff = max(
            negated_diff,
            max_abs_pair_diff(current[pair], [-float(reference[pair][0]), -float(reference[pair][1])]),
        )
    if same_diff <= tolerance:
        return "same"
    if negated_diff <= tolerance:
        return "global_negation"
    return "other"


def direction_formula(direction: tuple[int, int]) -> str:
    coeff_plus, coeff_minus = direction
    terms: list[str] = []
    if coeff_plus != 0:
        terms.append(f"{coeff_plus:+d}*o_plus")
    if coeff_minus != 0:
        terms.append(f"{coeff_minus:+d}*o_minus")
    return " ".join(terms)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_ATLAS, IN_Y_SEL, IN_P697, IN_P719] + [IN_A[pair_index] for pair_index in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P720",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    atlas = load_json(IN_ATLAS)
    y_sel = load_json(IN_Y_SEL)
    p697 = load_json(IN_P697)
    p719 = load_json(IN_P719)

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

    y_matrices: dict[int, list[list[float]]] = {}
    unsigned_root_output_vectors: dict[str, list[float]] = {}
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
        y_matrices[pair_index] = matrix_2x2
        unsigned_root_output_vectors[current_pair_name] = matvec2(matrix_2x2, coords[pair_index])

    tolerance_dot_nonzero = 1e-12
    tolerance_alignment = 1e-9
    tolerance_output = 1e-9

    directions = [direction for direction in itertools.product([-1, 0, 1], repeat=2) if direction != (0, 0)]

    total_candidates = 0
    exact_candidates = 0
    projective_only_candidates = 0
    unsupported_candidates = 0
    projective_only_negated_root_sets_seen: set[tuple[str, ...]] = set()
    projective_only_examples: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    unsupported_examples: list[dict[str, Any]] = []

    for direction in directions:
        total_candidates += 1
        rooted_results: dict[str, Any] = {}
        supported_roots: list[str] = []
        unsupported_roots: list[str] = []

        for root_index in range(1, 6):
            root_pair_name = pair_name(root_index)
            root_unsigned_output = unsigned_root_output_vectors[root_pair_name]
            root_scalar = dot([float(direction[0]), float(direction[1])], root_unsigned_output)
            if abs(root_scalar) <= tolerance_dot_nonzero:
                unsupported_roots.append(root_pair_name)
                rooted_results[root_pair_name] = {
                    "root_pair": root_pair_name,
                    "root_supported": False,
                    "unsigned_output_vector": root_unsigned_output,
                    "root_scalar": float(root_scalar),
                }
                continue

            root_sign = sign(root_scalar)
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
                output_vectors_by_pair[current_pair_name] = matvec2(y_matrices[pair_index], coords_star_by_pair[current_pair_name])

            root_result_supported = alignment_ok and len(sign_vector_by_pair) == 5 and len(output_vectors_by_pair) == 5
            if root_result_supported:
                supported_roots.append(root_pair_name)
            else:
                unsupported_roots.append(root_pair_name)

            rooted_results[root_pair_name] = {
                "root_pair": root_pair_name,
                "root_supported": root_result_supported,
                "unsigned_output_vector": root_unsigned_output,
                "root_scalar": float(root_scalar),
                "root_sign": int(root_sign),
                "sign_vector_by_pair": sign_vector_by_pair,
                "output_vectors_by_pair": output_vectors_by_pair,
            }

        formula = direction_formula(direction)
        if len(supported_roots) != 5:
            unsupported_candidates += 1
            if len(unsupported_examples) < 5:
                unsupported_examples.append(
                    {
                        "direction": list(direction),
                        "formula": formula,
                        "unsupported_roots": unsupported_roots,
                    }
                )
            continue

        reference_root = supported_roots[0]
        reference_sign_vector = rooted_results[reference_root]["sign_vector_by_pair"]
        reference_output_vectors = rooted_results[reference_root]["output_vectors_by_pair"]

        exact_root_independent_sign_vector = all(
            rooted_results[root_pair]["sign_vector_by_pair"] == reference_sign_vector for root_pair in supported_roots
        )
        exact_root_independent_output_vectors = all(
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
        if exact_root_independent_sign_vector and exact_root_independent_output_vectors:
            exact_candidates += 1
            continue

        sign_relations = [
            sign_vector_relation(reference_sign_vector, rooted_results[root_pair]["sign_vector_by_pair"])
            for root_pair in supported_roots
        ]
        output_relations = [
            output_vector_relation(reference_output_vectors, rooted_results[root_pair]["output_vectors_by_pair"], tolerance_output)
            for root_pair in supported_roots
        ]

        projective_only = all(relation in {"same", "global_negation"} for relation in sign_relations) and all(
            relation in {"same", "global_negation"} for relation in output_relations
        )
        if projective_only:
            projective_only_candidates += 1
            negated_roots = tuple(
                root_pair
                for root_pair in supported_roots
                if sign_vector_relation(reference_sign_vector, rooted_results[root_pair]["sign_vector_by_pair"]) == "global_negation"
            )
            projective_only_negated_root_sets_seen.add(negated_roots)
            bucket = projective_only_examples.setdefault(negated_roots, [])
            if len(bucket) < 5:
                bucket.append({"direction": list(direction), "formula": formula, "reference_root": reference_root})

    negated_root_sets_seen_sorted = sorted(projective_only_negated_root_sets_seen)
    max_unsigned_minus_abs = max(abs(float(unsigned_root_output_vectors[pair][1])) for pair in sorted(unsigned_root_output_vectors))

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
        "p697_observer_limit_readout_still_computable",
        p697.get("status"),
        "PASS_PROJECTIVE_OBSERVER_LIMIT_READOUT_COMPUTABLE_FROM_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT_PROJECTOR",
        "P697 still exposes the operational observer-limit output carrier below the global projective selector closure.",
    )
    add_check(
        "p719_previous_physically_interpretable_shift_was_requested",
        p719.get("preferred_next_direction"),
        "new_physically_interpretable_observable_provider_class_over_further_coefficient_tuning",
        "P719 already identified a new physically interpretable provider class as the preferred next continuation.",
    )
    add_check(
        "observer_facing_unsigned_root_outputs_are_currently_almost_pure_o_plus",
        max_unsigned_minus_abs <= tolerance_output,
        True,
        "On current exports, the unsigned rooted output vectors already lie almost entirely on the o_plus channel, so pure o_minus projections are operationally silent on this carrier.",
    )
    add_check(
        "observer_facing_output_axis_projection_exact_candidates_absent",
        exact_candidates > 0,
        False,
        "No static observer-facing signed output-axis projection on the current Q_out carrier yields one exact all-root directed section.",
    )
    add_check(
        "observer_facing_output_axis_projection_projective_only_candidates_present",
        projective_only_candidates > 0,
        True,
        "The same output-channel projection class still contains projective-only sectors, so the failure remains exact-provider level rather than total support failure.",
    )
    add_check(
        "observer_facing_output_axis_projection_negated_root_sets_match_pair4_pair5_split",
        [list(roots) for roots in negated_root_sets_seen_sorted],
        [["pair4", "pair5"]],
        "Every supported static observer-facing output-axis projection reproduces only the pair4/pair5 negated branch split on current exports.",
    )

    status = (
        "PASS_OBSERVER_FACING_OUTPUT_CHANNEL_PROJECTION_PROVIDER_CLASS_AUDITED"
        if not blocking
        else "P720_REQUIRES_REVIEW_CHANGED_OBSERVER_FACING_OUTPUT_CHANNEL_PROVIDER_STATE"
    )

    artifact = {
        "stage": "P720",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "attack the nearest physics-facing provider class after P719: static signed output-channel projections on the "
            "observer-facing output carrier Q_out = span{o_plus,o_minus}"
        ),
        "inputs": {
            "atlas": str(IN_ATLAS.relative_to(REPO)),
            "Y_sel": str(IN_Y_SEL.relative_to(REPO)),
            "P697": str(IN_P697.relative_to(REPO)),
            "P719": str(IN_P719.relative_to(REPO)),
            "A_packets": [str(IN_A[pair_index].relative_to(REPO)) for pair_index in sorted(IN_A)],
        },
        "scan_class": {
            "carrier": "observer_facing_unsigned_root_output_vectors_in_Q_out",
            "carrier_basis": ["o_plus", "o_minus"],
            "direction_coefficient_alphabet": [-1, 0, 1],
            "directions_scanned": [list(direction) for direction in directions],
            "total_candidates_scanned": total_candidates,
        },
        "unsigned_root_output_vectors_by_root": unsigned_root_output_vectors,
        "counts": {
            "total_candidates_scanned": total_candidates,
            "exact_candidates_found": exact_candidates,
            "projective_only_candidates_found": projective_only_candidates,
            "unsupported_candidates_found": unsupported_candidates,
        },
        "projective_only_negated_root_sets_seen": [list(roots) for roots in negated_root_sets_seen_sorted],
        "projective_only_examples_by_negated_root_set": {
            ",".join(roots): examples
            for roots, examples in sorted(projective_only_examples.items())
        },
        "unsupported_examples": unsupported_examples,
        "checks": checks,
        "blocking_mismatches": blocking,
        "result": {
            "observer_facing_output_axis_projection_exact_root_independent_section_exists": exact_candidates > 0,
            "observer_facing_output_axis_projection_projective_only_candidates_exist": projective_only_candidates > 0,
            "observer_facing_output_axis_projection_known_negated_root_sets": [list(roots) for roots in negated_root_sets_seen_sorted],
            "observer_facing_output_axis_projection_class_too_coarse_for_exact_t176": exact_candidates == 0,
            "preferred_next_direction": "dynamical_or_flux_like_physically_interpretable_provider_class_over_further_static_output_axis_projections",
            "counts_as_strict_physical_orientation_datum": False,
            "implies_t176_discharge": False,
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No strict-core directed/sign-sensitive physical orientation datum claim.",
            "No statement about nonlinear functions on Q_out beyond the scanned static linear projection class.",
            "No statement about source-topology flow/barrier observables already bridging into exact chartwise directed-sign coherence.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P720",
        "status": status,
        "as_of": AS_OF,
        "total_candidates_scanned": total_candidates,
        "exact_candidates_found": exact_candidates,
        "projective_only_candidates_found": projective_only_candidates,
        "unsupported_candidates_found": unsupported_candidates,
        "projective_only_negated_root_sets_seen": [list(roots) for roots in negated_root_sets_seen_sorted],
        "unsupported_direction_coefficients": [example["direction"] for example in unsupported_examples],
        "preferred_next_direction": "dynamical_or_flux_like_physically_interpretable_provider_class_over_further_static_output_axis_projections",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
