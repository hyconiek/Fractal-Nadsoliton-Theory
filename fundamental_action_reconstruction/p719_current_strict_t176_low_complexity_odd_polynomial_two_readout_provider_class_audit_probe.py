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

IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"
IN_Y_SEL = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe_summary.json"
IN_P716 = GENERATED / "p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe_summary.json"
IN_P718 = GENERATED / "p718_current_strict_t176_single_mixed_linear_weight_span_provider_insufficiency_audit_probe_summary.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p719_current_strict_t176_low_complexity_odd_polynomial_two_readout_provider_class_audit_probe.json"
OUT_SUMMARY = GENERATED / "p719_current_strict_t176_low_complexity_odd_polynomial_two_readout_provider_class_audit_probe_summary.json"


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


def candidate_formula(coeffs: tuple[int, ...], monomials: list[tuple[str, tuple[int, int]]]) -> str:
    pieces: list[str] = []
    for coeff, (name, _) in zip(coeffs, monomials):
        if coeff == 0:
            continue
        pieces.append(f"{coeff:+d}*{name}")
    return " ".join(pieces)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F647, IN_ATLAS, IN_Y_SEL, IN_P715, IN_P716, IN_P718] + [IN_A[pair_index] for pair_index in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P719",
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
    p718 = load_json(IN_P718)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break = payload.get("w_break_by_x")
    w_ref = payload.get("w_ref_unnormalized_by_x")
    if not is_numeric_list_len(w_break, 12):
        artifact = {
            "stage": "P719",
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
            "stage": "P719",
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
    readouts_by_root: dict[str, dict[str, float]] = {}
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
        readouts_by_root[pair_name(pair_index)] = {
            "dot_w_break_u_root": float(dot(w_break_by_x, vectors[pair_index])),
            "dot_w_ref_u_root": float(dot(w_ref_by_x, vectors[pair_index])),
        }

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

    odd_monomials: list[tuple[str, tuple[int, int]]] = [
        ("b", (1, 0)),
        ("r", (0, 1)),
        ("b3", (3, 0)),
        ("b2r", (2, 1)),
        ("br2", (1, 2)),
        ("r3", (0, 3)),
    ]

    total_candidates = 0
    projective_only_candidates = 0
    exact_candidates = 0
    not_all_roots_supported_candidates = 0
    projective_only_negated_root_sets_seen: set[tuple[str, ...]] = set()
    projective_only_examples: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    unsupported_examples: list[dict[str, Any]] = []

    for coeffs in itertools.product([-1, 0, 1], repeat=len(odd_monomials)):
        if all(coeff == 0 for coeff in coeffs):
            continue

        total_candidates += 1

        def scalar_provider(dot_break_root: float, dot_ref_root: float, coeffs_local: tuple[int, ...] = coeffs) -> float:
            total = 0.0
            for coeff, (_, (power_break, power_ref)) in zip(coeffs_local, odd_monomials):
                if coeff == 0:
                    continue
                total += coeff * (dot_break_root**power_break) * (dot_ref_root**power_ref)
            return float(total)

        rooted_results: dict[str, Any] = {}
        supported_roots: list[str] = []
        unsupported_roots: list[str] = []

        for root_index in range(1, 6):
            root_pair_name = pair_name(root_index)
            root_dot_break = readouts_by_root[root_pair_name]["dot_w_break_u_root"]
            root_dot_ref = readouts_by_root[root_pair_name]["dot_w_ref_u_root"]
            root_scalar = scalar_provider(root_dot_break, root_dot_ref)
            if abs(root_scalar) <= tolerance_dot_nonzero:
                unsupported_roots.append(root_pair_name)
                rooted_results[root_pair_name] = {
                    "root_pair": root_pair_name,
                    "root_supported": False,
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
                "root_scalar": float(root_scalar),
                "root_sign": int(root_sign),
                "sign_vector_by_pair": sign_vector_by_pair,
                "output_vectors_by_pair": output_vectors_by_pair,
            }

        formula = candidate_formula(coeffs, odd_monomials)
        if len(supported_roots) != 5:
            not_all_roots_supported_candidates += 1
            if len(unsupported_examples) < 5:
                unsupported_examples.append(
                    {
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
                bucket.append({"formula": formula, "reference_root": reference_root})

    negated_root_sets_seen_sorted = sorted(projective_only_negated_root_sets_seen)

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
        bool(p715.get("projective_root_independent_sign_orbit")) and not bool(p715.get("exact_root_independent_sign_vector")),
        True,
        "P715 still marks the current best all-root route as projective-only, not exact.",
    )
    add_check(
        "p716_pair4_localization_still_true",
        bool(p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")),
        True,
        "P716 still localizes the current best exact split to pair4 negative cosine polarity.",
    )
    add_check(
        "p718_single_linear_span_still_insufficient",
        bool(p718.get("single_mixed_linear_weight_span_exact_root_independent_section_exists")),
        False,
        "P718 still rules out the full single linear span of current exported weights.",
    )
    add_check(
        "low_complexity_odd_polynomial_two_readout_exact_candidates_absent",
        exact_candidates > 0,
        False,
        "No low-complexity untuned odd polynomial candidate on the current two-readout carrier yields one exact all-root directed section.",
    )
    add_check(
        "low_complexity_odd_polynomial_two_readout_projective_only_candidates_present",
        projective_only_candidates > 0,
        True,
        "The same nonlinear class still contains many projective-only candidates, so the failure remains exact-provider level rather than total support failure.",
    )
    add_check(
        "low_complexity_odd_polynomial_negated_root_sets_match_current_frontier_patterns",
        [list(roots) for roots in negated_root_sets_seen_sorted],
        [["pair2", "pair3"], ["pair4"]],
        "This minimal nonlinear class reproduces only the already-known exact split patterns: pair4-only or pair2/pair3-only negated branches.",
    )

    status = (
        "PASS_LOW_COMPLEXITY_ODD_POLYNOMIAL_TWO_READOUT_PROVIDER_CLASS_AUDITED"
        if not blocking
        else "P719_REQUIRES_REVIEW_CHANGED_LOW_COMPLEXITY_ODD_POLYNOMIAL_PROVIDER_STATE"
    )

    artifact = {
        "stage": "P719",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "attack the next honest nonlinear provider class after P718: untuned odd polynomials on the current two-readout "
            "carrier (<w_break,u>, <w_ref,u>) with total degree <= 3 and coefficient alphabet {-1,0,1}"
        ),
        "inputs": {
            "F647": str(IN_F647.relative_to(REPO)),
            "atlas": str(IN_ATLAS.relative_to(REPO)),
            "Y_sel": str(IN_Y_SEL.relative_to(REPO)),
            "P715": str(IN_P715.relative_to(REPO)),
            "P716": str(IN_P716.relative_to(REPO)),
            "P718": str(IN_P718.relative_to(REPO)),
            "A_packets": [str(IN_A[pair_index].relative_to(REPO)) for pair_index in sorted(IN_A)],
        },
        "scan_class": {
            "carrier": "two_readout_root_local_scalars_(dot_w_break_u_root,dot_w_ref_u_root)",
            "odd_monomial_basis_total_degree_le_3": [name for name, _ in odd_monomials],
            "coefficient_alphabet": [-1, 0, 1],
            "total_candidates_scanned": total_candidates,
        },
        "root_readout_profiles": readouts_by_root,
        "counts": {
            "total_candidates_scanned": total_candidates,
            "exact_candidates_found": exact_candidates,
            "projective_only_candidates_found": projective_only_candidates,
            "not_all_roots_supported_candidates": not_all_roots_supported_candidates,
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
            "low_complexity_odd_polynomial_two_readout_exact_root_independent_section_exists": exact_candidates > 0,
            "low_complexity_odd_polynomial_two_readout_projective_only_candidates_exist": projective_only_candidates > 0,
            "known_negated_root_sets_exhaust_current_scan": [list(roots) for roots in negated_root_sets_seen_sorted],
            "next_provider_must_either_raise_complexity_or_leave_current_two_readout_carrier": exact_candidates == 0,
            "preferred_next_direction": "new_physically_interpretable_observable_provider_class_over_further_coefficient_tuning",
            "counts_as_strict_physical_orientation_datum": False,
            "implies_t176_discharge": False,
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No strict-core directed/sign-sensitive physical orientation datum claim.",
            "No statement about tuned higher-degree coefficients outside the scanned class.",
            "No statement about provider classes using new observables beyond the current two-readout carrier.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P719",
        "status": status,
        "as_of": AS_OF,
        "total_candidates_scanned": total_candidates,
        "exact_candidates_found": exact_candidates,
        "projective_only_candidates_found": projective_only_candidates,
        "not_all_roots_supported_candidates": not_all_roots_supported_candidates,
        "projective_only_negated_root_sets_seen": [list(roots) for roots in negated_root_sets_seen_sorted],
        "next_provider_must_either_raise_complexity_or_leave_current_two_readout_carrier": exact_candidates == 0,
        "preferred_next_direction": "new_physically_interpretable_observable_provider_class_over_further_coefficient_tuning",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
