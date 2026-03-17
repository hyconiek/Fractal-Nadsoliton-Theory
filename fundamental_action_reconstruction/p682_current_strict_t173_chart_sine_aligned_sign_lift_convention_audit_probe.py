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
IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_Y_SEL = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

OUT = GENERATED / "p682_current_strict_t173_chart_sine_aligned_sign_lift_convention_audit_probe.json"
OUT_SUMMARY = GENERATED / "p682_current_strict_t173_chart_sine_aligned_sign_lift_convention_audit_probe_summary.json"


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

    prereq = [IN_F647, IN_THETA_PAIR, IN_Y_SEL] + [IN_A[m] for m in sorted(IN_A)]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P682",
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
    theta_pair = load_json(IN_THETA_PAIR)
    y_sel = load_json(IN_Y_SEL)

    try:
        payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
        w_ref = payload.get("w_ref_unnormalized_by_x")
        if not is_numeric_list_len(w_ref, 12):
            raise ValueError(
                "F647.constructed_source_object.exported_payload.w_ref_unnormalized_by_x must be length-12 numeric list"
            )
        w_ref_by_x = [float(v) for v in w_ref]
    except Exception as exc:
        artifact = {
            "stage": "P682",
            "status": "INVALID_W_REF_INPUT_SHAPE",
            "as_of": AS_OF,
            "error": str(exc),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    try:
        sigma_int = (theta_pair.get("inputs") or {}).get("sigma_int") or {}
        sigma_int_value = int(sigma_int.get("value"))
        if sigma_int_value not in (-1, 1):
            raise ValueError("sigma_int must be ±1")
    except Exception as exc:
        artifact = {
            "stage": "P682",
            "status": "INVALID_SIGMA_INT_INPUT",
            "as_of": AS_OF,
            "error": str(exc),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n = 12
    tol_dot_nonzero = 1e-12
    tol_output = 1e-9

    per_pair: dict[str, Any] = {}
    sign_lift_supported_pairs: list[str] = []
    sign_lift_missing_pairs: list[str] = []

    for m in sorted(IN_A):
        pair = f"pair{m}"
        a = load_json(IN_A[m])
        data = a.get("data") or {}

        u_key = "u_1" if m == 1 else f"u_{m}"
        coords_key = "u_1_coords_in_c1_s1" if m == 1 else f"u_{m}_coords_in_c{m}_s{m}"

        u = data.get(u_key)
        coords = data.get(coords_key)
        if not is_numeric_list_len(u, 12):
            raise SystemExit(f"Invalid {IN_A[m].relative_to(REPO)}: data.{u_key} must be length-12 numeric list")
        if not is_numeric_list_len(coords, 2):
            raise SystemExit(f"Invalid {IN_A[m].relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")

        u_by_x = [float(v) for v in u]
        u_coords = [float(v) for v in coords]

        dot_ref = sum(w_ref_by_x[x] * u_by_x[x] for x in range(n))
        dot_ref_nonzero = bool(abs(dot_ref) > tol_dot_nonzero)

        # Convention-layer reflection-breaking weight: w_break_m(x) = sigma_int * w_ref(x) * s_m(x).
        # NOTE: s_m depends on a chart embedding convention for Z_12 -> U(1); it is not Aut(Z_12)-invariant.
        s_m_by_x = [math.sin(2 * math.pi * m * x / n) for x in range(n)]
        w_break_m_by_x = [float(sigma_int_value) * w_ref_by_x[x] * s_m_by_x[x] for x in range(n)]
        dot_break = sum(w_break_m_by_x[x] * u_by_x[x] for x in range(n))
        dot_break_nonzero = bool(abs(dot_break) > tol_dot_nonzero)

        if dot_ref_nonzero:
            s_lift = sign(dot_ref)
            sign_source = "w_ref"
        elif dot_break_nonzero:
            s_lift = sign(dot_break)
            sign_source = "w_break_m(chart_sine_aligned)"
        else:
            s_lift = None
            sign_source = None

        if s_lift is None:
            sign_lift_missing_pairs.append(pair)
        else:
            sign_lift_supported_pairs.append(pair)

        y_chart = ((y_sel.get("charts") or {}).get(pair) or {})
        Y = y_chart.get("Y_sel_matrix_in_c_s_to_o")
        if not (isinstance(Y, list) and len(Y) == 2 and all(is_numeric_list_len(row, 2) for row in Y)):
            raise SystemExit(
                f"Invalid {IN_Y_SEL.relative_to(REPO)}: charts.{pair}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
            )
        Y2 = [[float(Y[0][0]), float(Y[0][1])], [float(Y[1][0]), float(Y[1][1])]]

        oriented_coords = u_coords if s_lift is None else [float(s_lift) * u_coords[0], float(s_lift) * u_coords[1]]
        y_vec = matvec2(Y2, oriented_coords)

        per_pair[pair] = {
            "pair": pair,
            "u_ref": str(IN_A[m].relative_to(REPO)),
            "dot_w_ref_u_m": float(dot_ref),
            "dot_w_break_m_u_m": float(dot_break),
            "dot_w_ref_nonzero": dot_ref_nonzero,
            "dot_w_break_m_nonzero": dot_break_nonzero,
            "proposed_sign_lift_s_m": s_lift,
            "proposed_sign_lift_source": sign_source,
            "output_vector_o_plus_o_minus_after_sign_lift": y_vec,
            "output_o_plus_component_sign_after_sign_lift": (
                None if abs(float(y_vec[0])) <= tol_output else ("+" if float(y_vec[0]) > 0 else "-")
            ),
        }

    all_pairs_supported = len(sign_lift_missing_pairs) == 0

    directed_output_sign_consistent = None
    if all_pairs_supported:
        signs = [per_pair[p]["output_o_plus_component_sign_after_sign_lift"] for p in sorted(per_pair)]
        directed_output_sign_consistent = len(set(signs)) == 1

    if not all_pairs_supported:
        status = "PARTIAL_CONVENTION_LEVEL_SIGN_LIFT_SUPPORTED_ONLY_ON_SUBSET_OF_CHARTS"
    elif bool(directed_output_sign_consistent):
        status = "PASS_CONVENTION_LEVEL_SIGN_LIFT_SUPPORTED_ON_ALL_CHARTS_AND_OUTPUT_SIGN_CONSISTENT"
    else:
        status = "FAIL_OUTPUT_SIGN_MISMATCH_UNDER_CHART_SINE_ALIGNED_SIGN_LIFT_CONVENTION"

    artifact = {
        "stage": "P682",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit a chart-sine-aligned (non-Aut(Z12)-invariant) sign-lift convention "
            "w_break_m(x)=sigma_int*w_ref(x)*s_m(x) and check whether it yields chartwise output-sign consistency "
            "for the exported global selector output channels"
        ),
        "inputs": {
            "f647_ref": str(IN_F647.relative_to(REPO)),
            "theta_pair_ref": str(IN_THETA_PAIR.relative_to(REPO)),
            "sigma_int_value": sigma_int_value,
            "y_sel_ref": str(IN_Y_SEL.relative_to(REPO)),
            "a_m_refs": {str(m): str(IN_A[m].relative_to(REPO)) for m in sorted(IN_A)},
            "tolerances": {
                "dot_nonzero_abs_tol": tol_dot_nonzero,
                "output_component_zero_tol": tol_output,
            },
        },
        "per_pair": per_pair,
        "result": {
            "status": status,
            "sign_lift_supported_pairs": sign_lift_supported_pairs,
            "sign_lift_missing_pairs": sign_lift_missing_pairs,
            "all_pairs_supported": all_pairs_supported,
            "directed_output_sign_consistent_if_all_pairs_supported": directed_output_sign_consistent,
            "counts_as_strict_physical_orientation_datum": False,
            "note": (
                "PASS here means only: a deterministic chart-convention sign-lift exists under the declared s_m embedding. "
                "It does not upgrade directed physical orientation into strict core."
            ),
        },
        "hard_limits": [
            "Convention-level only: depends on chart sine embedding s_m(x)=sin(2π m x/12), not Aut(Z_12)-invariant.",
            "Does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
            "Does not upgrade premise-based directed objects (T164/T171) into strict core.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P682",
        "status": status,
        "sigma_int_value": sigma_int_value,
        "sign_lift_supported_pairs": sign_lift_supported_pairs,
        "sign_lift_missing_pairs": sign_lift_missing_pairs,
        "all_pairs_supported": all_pairs_supported,
        "directed_output_sign_consistent_if_all_pairs_supported": directed_output_sign_consistent,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

