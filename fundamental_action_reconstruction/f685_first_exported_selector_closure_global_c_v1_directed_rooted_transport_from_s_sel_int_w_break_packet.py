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

IN_STATE_F684 = (
    GENERATED / "selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)
IN_OUTPUT_F660 = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"

OUT_OBJ = (
    GENERATED
    / "selector_closure_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f685_first_exported_selector_closure_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def matvec2(m: list[list[float]], v: list[float]) -> list[float]:
    return [
        float(m[0][0]) * float(v[0]) + float(m[0][1]) * float(v[1]),
        float(m[1][0]) * float(v[0]) + float(m[1][1]) * float(v[1]),
    ]


def max_abs_diff2(a: list[float], b: list[float]) -> float:
    return max(abs(float(a[0] - b[0])), abs(float(a[1] - b[1])))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_STATE_F684, IN_OUTPUT_F660]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "object": "SelectorClosure_global_C_v1_directed_rooted_transport_from_S_sel_int_w_break_strict_convention_v1",
            "stage": "F685",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_OBJ.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJ)
        return

    state = load_json(IN_STATE_F684)
    out = load_json(IN_OUTPUT_F660)

    coords_by_pair = (state.get("outputs") or {}).get("coords_in_c_s_after_sign_lift") or {}
    if not isinstance(coords_by_pair, dict):
        raise SystemExit("Invalid F684: outputs.coords_in_c_s_after_sign_lift must be an object")

    charts_out = out.get("charts") or {}
    if not isinstance(charts_out, dict):
        raise SystemExit("Invalid F660: charts must be an object")

    tol_glue = 1e-9

    charts: dict[str, Any] = {}
    v_out_by_pair: dict[str, list[float]] = {}
    for m in range(1, 6):
        pair = f"pair{m}"
        coords = coords_by_pair.get(pair)
        if not is_numeric_list_len(coords, 2):
            raise SystemExit(f"Invalid F684: outputs.coords_in_c_s_after_sign_lift.{pair} must be length-2 numeric list")
        y_chart = charts_out.get(pair) or {}
        Y = y_chart.get("Y_sel_matrix_in_c_s_to_o")
        if not (isinstance(Y, list) and len(Y) == 2 and all(is_numeric_list_len(row, 2) for row in Y)):
            raise SystemExit(
                f"Invalid F660: charts.{pair}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
            )
        Y2 = [[float(Y[0][0]), float(Y[0][1])], [float(Y[1][0]), float(Y[1][1])]]
        v_out = matvec2(Y2, [float(coords[0]), float(coords[1])])
        v_out_by_pair[pair] = v_out
        charts[pair] = {
            "chart_id": pair,
            "coords_in_c_s_after_sign_lift": [float(coords[0]), float(coords[1])],
            "Y_sel_matrix_in_c_s_to_o": Y2,
            "v_out_in_o_plus_o_minus": v_out,
        }

    ref_vec = v_out_by_pair["pair1"]
    max_abs_diff_to_pair1 = max(max_abs_diff2(v_out_by_pair[p], ref_vec) for p in sorted(v_out_by_pair))
    glued = bool(max_abs_diff_to_pair1 <= tol_glue)

    status = "PASS_EXPORTED" if glued else "FAIL_DIRECTED_CLOSURE_OUTPUT_NOT_GLUED_WITHIN_TOLERANCE"

    artifact = {
        "object": "SelectorClosure_global_C_v1_directed_rooted_transport_from_S_sel_int_w_break_strict_convention_v1",
        "stage": "F685",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit global directed closure object on C_v1 in a tracked convention scope by composing the exported "
            "w_break-rooted directed representative section (F684) with the promoted global output channels Y_sel(pair_m) (F660), "
            "yielding chartwise output vectors that glue in the fixed (o_+,o_-) basis. "
            "This does not claim any strict physical orientation datum."
        ),
        "domain": {
            "configuration_space_object_ref": "fundamental_action_reconstruction/generated/c_v1_void_configuration_space_in_local_b_tilde_1_sector_v1.json",
            "configuration_space_object": "C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1",
            "charts_cover": "U_pair1 = ... = U_pair5 = C_v1 (declared global cover)",
        },
        "depends_on": {
            "directed_state_ref": str(IN_STATE_F684.relative_to(REPO)),
            "output_operator_ref": str(IN_OUTPUT_F660.relative_to(REPO)),
            "note": "Sign is fixed only in the declared strict_convention scope of F684; no physical sign datum claim.",
        },
        "output_observable": {
            "basis": ["o_+", "o_-"],
            "reference_chart": "pair1",
            "output_vector_in_o_plus_o_minus": ref_vec,
            "gluing": {"max_abs_diff_to_pair1": float(max_abs_diff_to_pair1), "tol": tol_glue, "glued": glued},
        },
        "charts": charts,
        "hard_limits": [
            "Convention/section choice only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim Aut(Z_12)-invariant sign canonicity.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F685",
        "status": status,
        "object": artifact["object"],
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

