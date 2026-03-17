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
IN_OUTPUT_SIGN_LIFT_F683 = (
    GENERATED
    / "selector_output_sign_lift_global_c_v1_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)

OUT_JSON = GENERATED / "p685_current_strict_t173_w_break_rooted_directed_closure_audit_probe.json"
OUT_SUMMARY = GENERATED / "p685_current_strict_t173_w_break_rooted_directed_closure_audit_probe_summary.json"


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

    prereq = [IN_STATE_F684, IN_OUTPUT_F660, IN_OUTPUT_SIGN_LIFT_F683]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P685",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_JSON)
        return

    state = load_json(IN_STATE_F684)
    out = load_json(IN_OUTPUT_F660)
    ref = load_json(IN_OUTPUT_SIGN_LIFT_F683)

    coords_by_pair = (state.get("outputs") or {}).get("coords_in_c_s_after_sign_lift") or {}
    if not isinstance(coords_by_pair, dict):
        raise SystemExit("Invalid F684: outputs.coords_in_c_s_after_sign_lift must be an object")

    charts_out = out.get("charts") or {}
    if not isinstance(charts_out, dict):
        raise SystemExit("Invalid F660: charts must be an object")

    tol_glue = 1e-9
    tol_cross_check = 1e-9

    v_out_by_pair: dict[str, list[float]] = {}
    cross_check_max_abs_diff = 0.0

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

        # Cross-check with F683 chartwise output if present.
        ref_chart = ((ref.get("charts") or {}).get(pair) or {})
        v_ref = ref_chart.get("v_out_in_o_plus_o_minus")
        if is_numeric_list_len(v_ref, 2):
            cross_check_max_abs_diff = max(cross_check_max_abs_diff, max_abs_diff2(v_out, [float(v_ref[0]), float(v_ref[1])]))

    ref_vec = v_out_by_pair["pair1"]
    max_abs_diff_to_pair1 = max(max_abs_diff2(v_out_by_pair[p], ref_vec) for p in sorted(v_out_by_pair))
    glued = bool(max_abs_diff_to_pair1 <= tol_glue)
    cross_check_ok = bool(cross_check_max_abs_diff <= tol_cross_check)

    status = "PASS_W_BREAK_ROOTED_DIRECTED_CLOSURE_OUTPUT_GLUED"
    if not glued:
        status = "FAIL_DIRECTED_CLOSURE_OUTPUT_NOT_GLUED_WITHIN_TOLERANCE"
    elif not cross_check_ok:
        status = "FAIL_F684_F660_OUTPUT_MISMATCHES_F683_REFERENCE"

    artifact = {
        "stage": "P685",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_directed_closure_output_from_F684_directed_state_and_F660_output_channels_on_C_v1",
        "inputs": {
            "directed_state_ref": str(IN_STATE_F684.relative_to(REPO)),
            "output_operator_ref": str(IN_OUTPUT_F660.relative_to(REPO)),
            "reference_output_sign_lift_ref": str(IN_OUTPUT_SIGN_LIFT_F683.relative_to(REPO)),
        },
        "output_observable": {
            "basis": ["o_+", "o_-"],
            "reference_chart": "pair1",
            "output_vector_in_o_plus_o_minus": ref_vec,
            "gluing": {"max_abs_diff_to_pair1": float(max_abs_diff_to_pair1), "tol": tol_glue, "glued": glued},
        },
        "cross_checks": {
            "matches_F683_chartwise_outputs_if_present": cross_check_ok,
            "max_abs_diff_to_F683": float(cross_check_max_abs_diff),
            "tol": tol_cross_check,
        },
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
        "stage": "P685",
        "status": status,
        "glued": glued,
        "max_abs_diff_to_pair1": float(max_abs_diff_to_pair1),
        "matches_F683_if_present": cross_check_ok,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

