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

IN_SIGN_FIXED_STATE_F690 = (
    GENERATED
    / "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)
IN_OUTPUT_F660 = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"

OUT_OBJ = (
    GENERATED
    / "selector_closure_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f692_current_strict_t175_sign_fixed_directed_closure_export_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def dot(a: list[float], b: list[float]) -> float:
    return float(sum(float(x) * float(y) for x, y in zip(a, b, strict=True)))


def matvec2(m: list[list[float]], v: list[float]) -> list[float]:
    return [
        float(m[0][0]) * float(v[0]) + float(m[0][1]) * float(v[1]),
        float(m[1][0]) * float(v[0]) + float(m[1][1]) * float(v[1]),
    ]


def max_abs_diff2(a: list[float], b: list[float]) -> float:
    return max(abs(float(a[0] - b[0])), abs(float(a[1] - b[1])))


def fourier_c_s_basis_vectors_on_z12(k: int) -> tuple[list[float], list[float]]:
    if k <= 0 or k >= 6:
        raise ValueError("This export expects k in {1,2,3,4,5} for the Z_12 Fourier-degenerate pairs.")
    n = 12
    norm = math.sqrt(2.0 / n)
    c = [norm * math.cos(2.0 * math.pi * k * x / n) for x in range(n)]
    s = [norm * math.sin(2.0 * math.pi * k * x / n) for x in range(n)]
    return c, s


def is_2x2_numeric_matrix(obj: Any) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in obj)
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for row in obj for v in row)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_SIGN_FIXED_STATE_F690, IN_OUTPUT_F660]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "object": "SelectorClosure_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1",
            "stage": "F692",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_OBJ.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJ)
        return

    state = load_json(IN_SIGN_FIXED_STATE_F690)
    outop = load_json(IN_OUTPUT_F660)

    u_outs = ((state.get("outputs") or {}).get("u_vectors_directed_sign_fixed") or {})
    if not isinstance(u_outs, dict):
        raise SystemExit("Invalid F690: outputs.u_vectors_directed_sign_fixed must be an object")

    charts_out = outop.get("charts") or {}
    if not isinstance(charts_out, dict):
        raise SystemExit("Invalid F660: charts must be an object")

    tol_glue = 1e-9
    tol_sign = 1e-12

    charts: dict[str, Any] = {}
    v_out_by_pair: dict[str, list[float]] = {}
    v_out_raw_by_pair: dict[str, list[float]] = {}
    output_sign_lift_by_pair: dict[str, int] = {}
    coords_by_pair: dict[str, list[float]] = {}
    coords_norm_by_pair: dict[str, float] = {}
    span_residual_inf_by_pair: dict[str, float] = {}

    for k in range(1, 6):
        pair = f"pair{k}"
        u_key = f"u_{k}"
        u = u_outs.get(u_key)
        if not is_numeric_list_len(u, 12):
            raise SystemExit(f"Invalid F690: outputs.u_vectors_directed_sign_fixed.{u_key} must be length-12 numeric list")
        u12 = [float(x) for x in u]

        chart = charts_out.get(pair) or {}
        Y = chart.get("Y_sel_matrix_in_c_s_to_o")
        if not is_2x2_numeric_matrix(Y):
            raise SystemExit(f"Invalid F660: charts.{pair}.Y_sel_matrix_in_c_s_to_o must be 2x2 numeric list")
        Y2 = [[float(Y[0][0]), float(Y[0][1])], [float(Y[1][0]), float(Y[1][1])]]

        c_vec, s_vec = fourier_c_s_basis_vectors_on_z12(k)
        c_coord = dot(c_vec, u12)
        s_coord = dot(s_vec, u12)
        coords = [float(c_coord), float(s_coord)]

        # Sanity: u should lie in span{c_k,s_k}. Track residual to avoid silent basis mismatch.
        recon = [coords[0] * c_vec[i] + coords[1] * s_vec[i] for i in range(12)]
        span_residual = max(abs(u12[i] - recon[i]) for i in range(12))

        v_out_raw = matvec2(Y2, coords)
        o_plus = float(v_out_raw[0])
        if abs(o_plus) <= tol_sign:
            raise SystemExit(
                f"Directed closure output sign is not computable: pair={pair} has |o_+| <= {tol_sign}."
            )
        s_out = 1 if o_plus > 0.0 else -1
        v_out = [float(s_out) * float(v_out_raw[0]), float(s_out) * float(v_out_raw[1])]

        coords_by_pair[pair] = coords
        coords_norm_by_pair[pair] = float(math.sqrt(coords[0] * coords[0] + coords[1] * coords[1]))
        span_residual_inf_by_pair[pair] = float(span_residual)
        output_sign_lift_by_pair[pair] = int(s_out)
        v_out_raw_by_pair[pair] = v_out_raw
        v_out_by_pair[pair] = v_out

        charts[pair] = {
            "chart_id": pair,
            "k": k,
            "u_key": u_key,
            "coords_in_c_s": coords,
            "coords_norm": coords_norm_by_pair[pair],
            "span_residual_inf_norm": span_residual_inf_by_pair[pair],
            "Y_sel_matrix_in_c_s_to_o": Y2,
            "v_out_raw_in_o_plus_o_minus": v_out_raw,
            "raw_o_plus_sign": int(1 if v_out_raw[0] > 0.0 else -1),
            "output_sign_lift_applied": int(s_out),
            "v_out_in_o_plus_o_minus": v_out,
        }

    ref_vec = v_out_by_pair["pair1"]
    max_abs_diff_to_pair1 = max(max_abs_diff2(v_out_by_pair[p], ref_vec) for p in sorted(v_out_by_pair))
    glued = bool(max_abs_diff_to_pair1 <= tol_glue)

    raw_ref_vec = v_out_raw_by_pair["pair1"]
    raw_max_abs_diff_to_pair1 = max(
        max_abs_diff2(v_out_raw_by_pair[p], raw_ref_vec) for p in sorted(v_out_raw_by_pair)
    )
    raw_glued = bool(raw_max_abs_diff_to_pair1 <= tol_glue)

    status = "PASS_EXPORTED" if glued else "FAIL_DIRECTED_CLOSURE_OUTPUT_NOT_GLUED_WITHIN_TOLERANCE"

    artifact = {
        "object": "SelectorClosure_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1",
        "stage": "F692",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit global directed closure object on C_v1 in the declared strict_convention scope by composing the "
            "exported sign-fixed directed representative (F690) with the promoted global output channels Y_sel(pair_m) (F660), "
            "yielding chartwise output vectors that glue in the fixed (o_+,o_-) basis. "
            "This is a convention-layer construction only and does not claim any strict physical sign/orientation datum."
        ),
        "domain": {
            "configuration_space_object_ref": "fundamental_action_reconstruction/generated/c_v1_void_configuration_space_in_local_b_tilde_1_sector_v1.json",
            "configuration_space_object": "C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1",
            "charts_cover": "U_pair1 = ... = U_pair5 = C_v1 (declared global cover)",
        },
        "depends_on": {
            "sign_fixed_directed_state_ref": str(IN_SIGN_FIXED_STATE_F690.relative_to(REPO)),
            "output_operator_ref": str(IN_OUTPUT_F660.relative_to(REPO)),
            "note": "Sign fixing depends on exported strict-core payload weights (F647) via F690; no Aut(Z_12)-invariant sign claim.",
        },
        "output_sign_lift": {
            "meaning": "Explicit per-chart sign-lift required to define a directed global closure outcome under fixed output bases.",
            "rule": "s_out(pair_m) := sign(o_+(pair_m)) with tolerance tol_sign; v_out := s_out · v_out_raw makes o_+ > 0 by construction.",
            "tolerance": tol_sign,
            "signs_by_pair": output_sign_lift_by_pair,
            "note": "This is an explicit section choice; it is not claimed to be Aut(Z_12)-invariant and does not imply any kernel-alone/global QW-2191 discharge.",
        },
        "output_observable": {
            "basis": ["o_+", "o_-"],
            "reference_chart": "pair1",
            "output_vector_in_o_plus_o_minus": ref_vec,
            "gluing": {"max_abs_diff_to_pair1": float(max_abs_diff_to_pair1), "tol": tol_glue, "glued": glued},
        },
        "raw_output_observable": {
            "basis": ["o_+", "o_-"],
            "reference_chart": "pair1",
            "output_vector_raw_in_o_plus_o_minus": raw_ref_vec,
            "gluing": {
                "max_abs_diff_to_pair1": float(raw_max_abs_diff_to_pair1),
                "tol": tol_glue,
                "glued": raw_glued,
            },
        },
        "charts": charts,
        "sanity": {
            "coords_by_pair": coords_by_pair,
            "coords_norm_by_pair": coords_norm_by_pair,
            "span_residual_inf_norm_by_pair": span_residual_inf_by_pair,
        },
        "hard_limits": [
            "Convention layer only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim Aut(Z_12)-invariant sign canonicity (N462 discipline).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F692",
        "status": status,
        "object": artifact["object"],
        "raw_glued": raw_glued,
        "glued": glued,
        "max_abs_diff_to_pair1": float(max_abs_diff_to_pair1),
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT_OBJ.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJ)


if __name__ == "__main__":
    main()
