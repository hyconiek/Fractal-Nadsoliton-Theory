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
IN_P713 = GENERATED / "p713_current_strict_t176_multiroot_rooted_sign_lift_root_independence_audit_probe_summary.json"
IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p714_current_strict_t176_w_break_parity_root_support_profile_audit_probe.json"
OUT_SUMMARY = GENERATED / "p714_current_strict_t176_w_break_parity_root_support_profile_audit_probe_summary.json"


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


def reflection_index(x: int, n: int) -> int:
    return (-x) % n


def parity_residuals(vector: list[float]) -> tuple[float, float]:
    n = len(vector)
    even_residual = 0.0
    odd_residual = 0.0
    for idx in range(n):
        reflected = reflection_index(idx, n)
        even_residual = max(even_residual, abs(float(vector[idx]) - float(vector[reflected])))
        odd_residual = max(odd_residual, abs(float(vector[idx]) + float(vector[reflected])))
    return float(even_residual), float(odd_residual)


def classify_parity(vector: list[float], tol: float) -> tuple[str, float, float]:
    even_residual, odd_residual = parity_residuals(vector)
    if even_residual <= tol and odd_residual > tol:
        return "even", even_residual, odd_residual
    if odd_residual <= tol and even_residual > tol:
        return "odd", even_residual, odd_residual
    if even_residual <= tol and odd_residual <= tol:
        return "degenerate", even_residual, odd_residual
    return "mixed", even_residual, odd_residual


def classify_axis(coords: list[float], tol: float) -> str:
    c_coord, s_coord = float(coords[0]), float(coords[1])
    if abs(c_coord) <= tol and abs(abs(s_coord) - 1.0) <= tol:
        return "sine_axis"
    if abs(s_coord) <= tol and abs(abs(c_coord) - 1.0) <= tol:
        return "cosine_axis"
    return "mixed_axis"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F647, IN_P713] + [IN_A[m] for m in sorted(IN_A)]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P714",
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
    p713 = load_json(IN_P713)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break = payload.get("w_break_by_x")
    if not is_numeric_list_len(w_break, 12):
        artifact = {
            "stage": "P714",
            "status": "INVALID_W_BREAK_INPUT_SHAPE",
            "as_of": AS_OF,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    tolerance_parity = 1e-12
    tolerance_axis = 1e-9
    tolerance_dot_nonzero = 1e-12

    w_break_vector = [float(entry) for entry in w_break]
    w_break_parity, w_break_even_residual, w_break_odd_residual = classify_parity(w_break_vector, tolerance_parity)

    per_pair: dict[str, Any] = {}
    odd_pairs: list[str] = []
    even_pairs: list[str] = []
    nonzero_anchor_pairs: list[str] = []
    zero_anchor_pairs: list[str] = []

    for pair_index in sorted(IN_A):
        pair_packet = load_json(IN_A[pair_index])
        data = pair_packet.get("data") or {}
        pair_name = f"pair{pair_index}"
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
        parity_label, even_residual, odd_residual = classify_parity(vector, tolerance_parity)
        axis_label = classify_axis(coords, tolerance_axis)
        anchor_value = dot(w_break_vector, vector)
        anchor_nonzero = abs(anchor_value) > tolerance_dot_nonzero

        if parity_label == "odd":
            odd_pairs.append(pair_name)
        if parity_label == "even":
            even_pairs.append(pair_name)
        if anchor_nonzero:
            nonzero_anchor_pairs.append(pair_name)
        else:
            zero_anchor_pairs.append(pair_name)

        per_pair[pair_name] = {
            "pair": pair_name,
            "coords": coords,
            "axis_type": axis_label,
            "parity_type": parity_label,
            "even_reflection_residual": float(even_residual),
            "odd_reflection_residual": float(odd_residual),
            "dot_w_break_with_u_m": float(anchor_value),
            "anchor_nonzero": bool(anchor_nonzero),
        }

    expected_supported_pairs = ["pair1", "pair5"]
    reported_supported_pairs = p713.get("supported_roots") if isinstance(p713.get("supported_roots"), list) else None

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
        "w_break_is_odd_under_reflection",
        w_break_parity,
        "odd",
        "The exported strict-core payload w_break is odd under reflection on Z_12.",
    )
    add_check(
        "odd_pairs_match_endpoint_sine_axis_corridor",
        odd_pairs,
        ["pair1", "pair5"],
        "The current exported odd-axis representatives are exactly pair1 and pair5.",
    )
    add_check(
        "even_pairs_match_interior_cosine_axis_corridor",
        even_pairs,
        ["pair2", "pair3", "pair4"],
        "The current exported even-axis representatives are exactly pair2, pair3, pair4.",
    )
    add_check(
        "nonzero_w_break_anchors_match_odd_pairs",
        nonzero_anchor_pairs,
        ["pair1", "pair5"],
        "Nonzero linear anchors from w_break occur only on the odd-axis pair corridor.",
    )
    add_check(
        "zero_w_break_anchors_match_even_pairs",
        zero_anchor_pairs,
        ["pair2", "pair3", "pair4"],
        "Linear anchors from w_break vanish on the even-axis interior corridor.",
    )
    add_check(
        "p713_supported_roots_match_parity_profile",
        reported_supported_pairs,
        expected_supported_pairs,
        "The supported roots reported by P713 match the parity-typed root-support profile.",
    )

    status = (
        "PASS_W_BREAK_PARITY_ROOT_SUPPORT_PROFILE_AUDITED"
        if not blocking
        else "P714_REQUIRES_REVIEW_CHANGED_W_BREAK_PARITY_SUPPORT_PROFILE"
    )

    artifact = {
        "stage": "P714",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the current supported-root corridor of the exported w_break payload is explained exactly by a parity profile "
            "on Z_12: odd w_break, odd-axis anchor support on pair1/pair5, and forced vanishing on even-axis charts pair2/pair3/pair4"
        ),
        "inputs": {
            "F647": str(IN_F647.relative_to(REPO)),
            "P713": str(IN_P713.relative_to(REPO)),
            "A_m_refs": {str(m): str(IN_A[m].relative_to(REPO)) for m in sorted(IN_A)},
        },
        "tolerances": {
            "parity_residual_tol": tolerance_parity,
            "axis_coord_tol": tolerance_axis,
            "anchor_nonzero_abs_tol": tolerance_dot_nonzero,
        },
        "w_break_profile": {
            "parity_type": w_break_parity,
            "even_reflection_residual": float(w_break_even_residual),
            "odd_reflection_residual": float(w_break_odd_residual),
        },
        "per_pair": per_pair,
        "checks": checks,
        "blocking_mismatches": blocking,
        "result": {
            "odd_pairs": odd_pairs,
            "even_pairs": even_pairs,
            "nonzero_anchor_pairs": nonzero_anchor_pairs,
            "zero_anchor_pairs": zero_anchor_pairs,
            "p713_supported_roots": reported_supported_pairs,
            "current_w_break_anchor_profile": "odd_weight_supports_odd_axis_pairs_only" if not blocking else None,
            "current_w_break_explains_supported_root_corridor": not blocking,
            "current_single_w_break_all_chart_anchor_exists": False,
            "next_honest_requirement": (
                "future strict-core provider must add an even-or-mixed parity anchor component, or use a genuinely different non-linear/non-scalar observable class"
            ),
        },
        "hard_limits": [
            "No impossibility-in-principle claim about future all-chart strict-core anchors.",
            "No T176 discharge claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No directed/sign-sensitive physical orientation datum claim in strict core.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P714",
        "status": status,
        "as_of": AS_OF,
        "w_break_parity_type": w_break_parity,
        "odd_pairs": odd_pairs,
        "even_pairs": even_pairs,
        "nonzero_anchor_pairs": nonzero_anchor_pairs,
        "zero_anchor_pairs": zero_anchor_pairs,
        "p713_supported_roots": reported_supported_pairs,
        "current_w_break_explains_supported_root_corridor": not blocking,
        "current_single_w_break_all_chart_anchor_exists": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
