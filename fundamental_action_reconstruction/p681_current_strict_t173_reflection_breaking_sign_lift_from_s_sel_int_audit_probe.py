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
IN_Y_SEL = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT = GENERATED / "p681_current_strict_t173_reflection_breaking_sign_lift_from_s_sel_int_audit_probe.json"
OUT_SUMMARY = GENERATED / "p681_current_strict_t173_reflection_breaking_sign_lift_from_s_sel_int_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def sign(x: float) -> int:
    return 1 if x > 0 else -1


def matvec2(m: list[list[float]], v: list[float]) -> list[float]:
    return [
        float(m[0][0]) * float(v[0]) + float(m[0][1]) * float(v[1]),
        float(m[1][0]) * float(v[0]) + float(m[1][1]) * float(v[1]),
    ]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F647, IN_Y_SEL] + [IN_A[m] for m in sorted(IN_A)]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P681",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f647 = load_json(IN_F647)
    y_sel = load_json(IN_Y_SEL)

    # Extract weight candidates from the exported strict-core source-object witness provider payload.
    try:
        payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
        w_ref = payload.get("w_ref_unnormalized_by_x")
        w_break = payload.get("w_break_by_x")
        if not is_numeric_list_len(w_ref, 12):
            raise ValueError(
                "F647.constructed_source_object.exported_payload.w_ref_unnormalized_by_x must be length-12 numeric list"
            )
        if not is_numeric_list_len(w_break, 12):
            raise ValueError("F647.constructed_source_object.exported_payload.w_break_by_x must be length-12 numeric list")
        w_ref_by_x = [float(v) for v in w_ref]
        w_break_by_x = [float(v) for v in w_break]
    except Exception as exc:
        artifact = {
            "stage": "P681",
            "status": "INVALID_W_BREAK_INPUT_SHAPE",
            "as_of": AS_OF,
            "error": str(exc),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n = 12

    def reflect(x: int) -> int:
        return (-x) % n

    tol_reflection_even = 1e-12
    tol_dot_nonzero = 1e-12
    tol_output = 1e-9

    w_ref_reflection_even_max_abs = max(
        abs(w_ref_by_x[x] - w_ref_by_x[reflect(x)]) for x in range(n)
    )
    w_ref_reflection_even = bool(w_ref_reflection_even_max_abs <= tol_reflection_even)

    w_break_reflection_even_max_abs = max(
        abs(w_break_by_x[x] - w_break_by_x[reflect(x)]) for x in range(n)
    )
    w_break_reflection_even = bool(w_break_reflection_even_max_abs <= tol_reflection_even)

    per_pair: dict[str, Any] = {}
    sign_lift_supported_pairs: list[str] = []
    sign_lift_missing_pairs: list[str] = []

    output_vectors_after_sign_lift: dict[str, list[float]] = {}

    for m in sorted(IN_A):
        pair = f"pair{m}"
        a = load_json(IN_A[m])
        data = a.get("data") or {}

        u_key = f"u_{m}"
        coords_key = f"u_{m}_coords_in_c{m}_s{m}"
        if m == 1:
            u_key = "u_1"
            coords_key = "u_1_coords_in_c1_s1"

        u = data.get(u_key)
        coords = data.get(coords_key)
        if not is_numeric_list_len(u, 12):
            raise SystemExit(f"Invalid {IN_A[m].relative_to(REPO)}: data.{u_key} must be length-12 numeric list")
        if not is_numeric_list_len(coords, 2):
            raise SystemExit(f"Invalid {IN_A[m].relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")

        u_by_x = [float(v) for v in u]
        u_coords = [float(v) for v in coords]

        dot_ref = sum(w_ref_by_x[x] * u_by_x[x] for x in range(n))
        dot_break = sum(w_break_by_x[x] * u_by_x[x] for x in range(n))
        dot_ref_nonzero = bool(abs(dot_ref) > tol_dot_nonzero)
        dot_break_nonzero = bool(abs(dot_break) > tol_dot_nonzero)
        if dot_ref_nonzero:
            s_m = sign(dot_ref)  # choose sign so dot(w_ref, s_m*u) > 0
            sign_source = "w_ref"
        elif dot_break_nonzero:
            s_m = sign(dot_break)  # choose sign so dot(w_break, s_m*u) > 0
            sign_source = "w_break"
        else:
            s_m = None
            sign_source = None

        if s_m is None:
            sign_lift_missing_pairs.append(pair)
        else:
            sign_lift_supported_pairs.append(pair)

        # Optional sanity: output vector in (o_+,o_-) after any proposed sign lift.
        y_chart = ((y_sel.get("charts") or {}).get(pair) or {})
        Y = y_chart.get("Y_sel_matrix_in_c_s_to_o")
        if not (isinstance(Y, list) and len(Y) == 2 and all(is_numeric_list_len(row, 2) for row in Y)):
            raise SystemExit(
                f"Invalid {IN_Y_SEL.relative_to(REPO)}: charts.{pair}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
            )
        Y2 = [[float(Y[0][0]), float(Y[0][1])], [float(Y[1][0]), float(Y[1][1])]]

        oriented_coords = u_coords if s_m is None else [float(s_m) * u_coords[0], float(s_m) * u_coords[1]]
        y_vec = matvec2(Y2, oriented_coords)
        output_vectors_after_sign_lift[pair] = y_vec

        per_pair[pair] = {
            "pair": pair,
            "u_ref": str(IN_A[m].relative_to(REPO)),
            "dot_w_ref_u_m": float(dot_ref),
            "dot_w_break_u_m": float(dot_break),
            "dot_w_ref_nonzero": dot_ref_nonzero,
            "dot_w_break_nonzero": dot_break_nonzero,
            "proposed_sign_lift_s_m": s_m,
            "proposed_sign_lift_source": sign_source,
            "output_vector_o_plus_o_minus_after_sign_lift": y_vec,
            "output_o_plus_component_sign_after_sign_lift": (
                None
                if abs(float(y_vec[0])) <= tol_output
                else ("+" if float(y_vec[0]) > 0 else "-")
            ),
        }

    all_pairs_supported = len(sign_lift_missing_pairs) == 0
    any_pair_supported = len(sign_lift_supported_pairs) > 0 or all_pairs_supported

    # Sanity: if we can lift sign everywhere, check directed output sign consistency across charts.
    directed_output_sign_consistent = None
    if all_pairs_supported:
        signs = [
            per_pair[p]["output_o_plus_component_sign_after_sign_lift"]
            for p in sorted(per_pair)
        ]
        directed_output_sign_consistent = len(set(signs)) == 1

    if not any_pair_supported and not all_pairs_supported:
        status = "FAIL_NO_DETERMINISTIC_SIGN_LIFT_FROM_EXPORTED_S_SEL_INT_WEIGHT_PAYLOAD"
    elif not all_pairs_supported:
        status = "PARTIAL_SIGN_LIFT_SUPPORTED_ONLY_ON_SUBSET_OF_CHARTS"
    elif bool(directed_output_sign_consistent):
        status = "PASS_DETERMINISTIC_SIGN_LIFT_SUPPORTED_ON_ALL_CHARTS_AND_OUTPUT_SIGN_CONSISTENT"
    else:
        status = "FAIL_OUTPUT_SIGN_MISMATCH_UNDER_DETERMINISTIC_SIGN_LIFT_CANDIDATE"

    artifact = {
        "stage": "P681",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the exported strict-core internal selector-source seed witness provider already contains "
            "a reflection-breaking weight ingredient that can deterministically lift residual Z2 sign across pair1..pair5 "
            "via a scalar dot(w_break,u_m), and whether that supports a global directed output-sign consistency check"
        ),
        "inputs": {
            "f647_ref": str(IN_F647.relative_to(REPO)),
            "y_sel_ref": str(IN_Y_SEL.relative_to(REPO)),
            "a_m_refs": {str(m): str(IN_A[m].relative_to(REPO)) for m in sorted(IN_A)},
            "tolerances": {
                "reflection_even_max_abs_tol": tol_reflection_even,
                "dot_nonzero_abs_tol": tol_dot_nonzero,
                "output_component_zero_tol": tol_output,
            },
        },
        "w_break": {
            "w_ref": {
                "reflection_even_max_abs": float(w_ref_reflection_even_max_abs),
                "reflection_even": w_ref_reflection_even,
                "note": "even reference weights can distinguish sign on even (cos-like) representatives but cannot distinguish sign on odd (sin-like) representatives.",
            },
            "w_break": {
                "reflection_even_max_abs": float(w_break_reflection_even_max_abs),
                "reflection_even": w_break_reflection_even,
                "note": "reflection_even=false indicates reflection-breaking weights (directionful); this is necessary for sign distinction on odd (sin-like) axes but is not sufficient for a globally consistent directed closure outcome.",
            },
        },
        "per_pair": per_pair,
        "result": {
            "status": status,
            "sign_lift_supported_pairs": sign_lift_supported_pairs,
            "sign_lift_missing_pairs": sign_lift_missing_pairs,
            "all_pairs_supported": all_pairs_supported,
            "directed_output_sign_consistent_if_all_pairs_supported": directed_output_sign_consistent,
        },
        "hard_limits": [
            "Probe-only: does not discharge a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
            "Does not upgrade premise-based directed objects (T164/T171) into strict core.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P681",
        "status": status,
        "w_ref_reflection_even": w_ref_reflection_even,
        "w_break_reflection_even": w_break_reflection_even,
        "sign_lift_supported_pairs": sign_lift_supported_pairs,
        "sign_lift_missing_pairs": sign_lift_missing_pairs,
        "all_pairs_supported": all_pairs_supported,
        "directed_output_sign_consistent_if_all_pairs_supported": directed_output_sign_consistent,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
