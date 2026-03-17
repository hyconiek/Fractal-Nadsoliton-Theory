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

IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json"
IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe.json"
OUT_SUMMARY = GENERATED / "p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe_summary.json"


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


def dot(lhs: list[float], rhs: list[float]) -> float:
    return sum(float(x) * float(y) for x, y in zip(lhs, rhs))


def classify_axis(coords: list[float], tolerance: float) -> str:
    c_coord, s_coord = float(coords[0]), float(coords[1])
    if abs(c_coord) <= tolerance and abs(abs(s_coord) - 1.0) <= tolerance:
        return "sine_axis"
    if abs(s_coord) <= tolerance and abs(abs(c_coord) - 1.0) <= tolerance:
        return "cosine_axis"
    return "mixed_axis"


def polarity_label(coords: list[float], tolerance: float) -> str:
    c_coord, s_coord = float(coords[0]), float(coords[1])
    if abs(s_coord) <= tolerance and abs(abs(c_coord) - 1.0) <= tolerance:
        return "positive_cosine" if c_coord > 0.0 else "negative_cosine"
    if abs(c_coord) <= tolerance and abs(abs(s_coord) - 1.0) <= tolerance:
        return "positive_sine" if s_coord > 0.0 else "negative_sine"
    return "mixed"


def negate_sign_vector(sign_vector: dict[str, int]) -> dict[str, int]:
    return {pair: int(-value) for pair, value in sign_vector.items()}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P715, IN_F647] + [IN_A[pair_index] for pair_index in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P716",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p715 = load_json(IN_P715)
    f647 = load_json(IN_F647)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_ref = payload.get("w_ref_unnormalized_by_x")
    w_break = payload.get("w_break_by_x")
    if not is_numeric_list_len(w_ref, 12):
        artifact = {
            "stage": "P716",
            "status": "INVALID_W_REF_INPUT_SHAPE",
            "as_of": AS_OF,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return
    if not is_numeric_list_len(w_break, 12):
        artifact = {
            "stage": "P716",
            "status": "INVALID_W_BREAK_INPUT_SHAPE",
            "as_of": AS_OF,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    w_ref_by_x = [float(entry) for entry in w_ref]
    w_break_by_x = [float(entry) for entry in w_break]

    tolerance_axis = 1e-9

    per_pair: dict[str, Any] = {}
    positive_cosine_pairs: list[str] = []
    negative_cosine_pairs: list[str] = []

    for pair_index in sorted(IN_A):
        pair_packet = load_json(IN_A[pair_index])
        data = pair_packet.get("data") or {}
        pair = f"pair{pair_index}"
        vector_key = "u_1" if pair_index == 1 else f"u_{pair_index}"
        coords_key = "u_1_coords_in_c1_s1" if pair_index == 1 else f"u_{pair_index}_coords_in_c{pair_index}_s{pair_index}"
        vector_value = data.get(vector_key)
        coords_value = data.get(coords_key)
        if not is_numeric_list_len(vector_value, 12):
            raise SystemExit(f"Invalid {IN_A[pair_index].relative_to(REPO)}: data.{vector_key} must be length-12 numeric list")
        if not is_numeric_list_len(coords_value, 2):
            raise SystemExit(f"Invalid {IN_A[pair_index].relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")

        vector = [float(entry) for entry in vector_value]
        coords = [float(entry) for entry in coords_value]
        axis_type = classify_axis(coords, tolerance_axis)
        polarity_type = polarity_label(coords, tolerance_axis)
        dot_ref = dot(w_ref_by_x, vector)
        dot_break = dot(w_break_by_x, vector)

        if polarity_type == "positive_cosine":
            positive_cosine_pairs.append(pair)
        if polarity_type == "negative_cosine":
            negative_cosine_pairs.append(pair)

        per_pair[pair] = {
            "pair": pair,
            "coords": coords,
            "axis_type": axis_type,
            "polarity_type": polarity_type,
            "dot_w_ref_u_m": float(dot_ref),
            "dot_w_break_u_m": float(dot_break),
        }

    result = p715.get("result") or {}
    rooted_results = p715.get("rooted_results") or {}
    reference_root = result.get("reference_root")
    reference_sign_vector = (rooted_results.get(reference_root) or {}).get("sign_vector_by_pair") if isinstance(reference_root, str) else None
    pair4_sign_vector = (rooted_results.get("pair4") or {}).get("sign_vector_by_pair")

    pair4_is_global_negation = False
    if isinstance(reference_sign_vector, dict) and isinstance(pair4_sign_vector, dict):
        pair4_is_global_negation = pair4_sign_vector == negate_sign_vector(reference_sign_vector)

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
        "p715_all_roots_supported",
        bool(p715.get("result", {}).get("all_roots_supported")),
        True,
        "P715 already supports all roots under the dual-anchor rule.",
    )
    add_check(
        "p715_projective_orbit_survives",
        bool(p715.get("result", {}).get("projective_root_independent_sign_orbit"))
        and bool(p715.get("result", {}).get("projective_root_independent_output_orbit")),
        True,
        "P715 already recovers one common projective rooted section orbit.",
    )
    add_check(
        "p715_exact_section_still_absent",
        bool(p715.get("result", {}).get("exact_root_independent_sign_vector"))
        or bool(p715.get("result", {}).get("exact_root_independent_output_vectors")),
        False,
        "P715 still does not recover one exact directed section across all roots.",
    )
    add_check(
        "pair4_is_unique_negated_orbit_root",
        result.get("negated_orbit_roots_relative_to_reference"),
        ["pair4"],
        "The residual exact split is localized to one negated-orbit root: pair4.",
    )
    add_check(
        "pair4_is_unique_negative_cosine_axis_chart",
        negative_cosine_pairs,
        ["pair4"],
        "Among current cosine-axis charts, pair4 is the unique negative-polarity representative.",
    )
    add_check(
        "positive_cosine_axis_charts_are_pair2_pair3",
        positive_cosine_pairs,
        ["pair2", "pair3"],
        "The remaining cosine-axis charts pair2 and pair3 are the positive-polarity representatives.",
    )
    add_check(
        "pair4_uses_even_anchor_fallback",
        (rooted_results.get("pair4") or {}).get("root_anchor_source"),
        "w_ref_unnormalized",
        "The pair4 root is anchored by the even fallback weight.",
    )
    add_check(
        "pair2_pair3_use_even_anchor_fallback",
        [
            (rooted_results.get("pair2") or {}).get("root_anchor_source"),
            (rooted_results.get("pair3") or {}).get("root_anchor_source"),
        ],
        ["w_ref_unnormalized", "w_ref_unnormalized"],
        "The other cosine-axis roots pair2/pair3 also use the same even fallback weight.",
    )
    add_check(
        "pair4_even_anchor_scalar_is_negative",
        bool(float((rooted_results.get("pair4") or {}).get("dot_w_ref_u_root", 0.0)) < 0.0),
        True,
        "The even-anchor scalar on pair4 is negative.",
    )
    add_check(
        "pair2_pair3_even_anchor_scalars_are_positive",
        [
            bool(float((rooted_results.get("pair2") or {}).get("dot_w_ref_u_root", 0.0)) > 0.0),
            bool(float((rooted_results.get("pair3") or {}).get("dot_w_ref_u_root", 0.0)) > 0.0),
        ],
        [True, True],
        "The even-anchor scalars on pair2/pair3 are positive.",
    )
    add_check(
        "pair4_sign_vector_is_global_negation_of_reference",
        pair4_is_global_negation,
        True,
        "The pair4 rooted sign vector is exactly the global negation of the reference sign vector.",
    )

    status = (
        "PASS_PAIR4_NEGATIVE_COSINE_POLARITY_EXPLAINS_CURRENT_GLOBAL_Z2_ORBIT_SPLIT"
        if not blocking
        else "P716_REQUIRES_REVIEW_CHANGED_ORBIT_SPLIT_PROFILE"
    )

    artifact = {
        "stage": "P716",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the residual exact global Z2 split left by the dual-anchor route is already fully explained "
            "by the unique negative cosine-axis role of pair4 under the current even-anchor fallback"
        ),
        "inputs": {
            "P715": str(IN_P715.relative_to(REPO)),
            "F647": str(IN_F647.relative_to(REPO)),
            "A_m_refs": {str(pair_index): str(IN_A[pair_index].relative_to(REPO)) for pair_index in sorted(IN_A)},
        },
        "per_pair_profile": per_pair,
        "checks": checks,
        "blocking_mismatches": blocking,
        "result": {
            "pair4_is_unique_negative_cosine_axis_chart": negative_cosine_pairs == ["pair4"],
            "positive_cosine_axis_charts": positive_cosine_pairs,
            "negative_cosine_axis_charts": negative_cosine_pairs,
            "pair4_negated_orbit_relative_to_reference": result.get("negated_orbit_roots_relative_to_reference") == ["pair4"],
            "pair4_even_anchor_scalar_negative": bool(float((rooted_results.get("pair4") or {}).get("dot_w_ref_u_root", 0.0)) < 0.0),
            "pair2_pair3_even_anchor_scalars_positive": [
                bool(float((rooted_results.get("pair2") or {}).get("dot_w_ref_u_root", 0.0)) > 0.0),
                bool(float((rooted_results.get("pair3") or {}).get("dot_w_ref_u_root", 0.0)) > 0.0),
            ],
            "pair4_sign_vector_is_global_negation_of_reference": pair4_is_global_negation,
            "current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity": len(blocking) == 0,
            "counts_as_strict_physical_orientation_datum": False,
            "implies_t176_discharge": False,
            "remaining_gap_after_positive_result": (
                "an exact global branch-fixing rule beyond the current pair4 negative-cosine polarity split, with no hidden selector slots"
            ),
        },
        "hard_limits": [
            "No strict-core directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No T176 discharge claim.",
            "No claim that the residual split is already gauge-irrelevant in strict core.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P716",
        "status": status,
        "as_of": AS_OF,
        "pair4_is_unique_negative_cosine_axis_chart": negative_cosine_pairs == ["pair4"],
        "pair4_negated_orbit_relative_to_reference": result.get("negated_orbit_roots_relative_to_reference") == ["pair4"],
        "current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity": len(blocking) == 0,
        "counts_as_strict_physical_orientation_datum": False,
        "implies_t176_discharge": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
