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

IN_O_1M = {
    2: GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json",
    3: GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json",
    4: GENERATED / "o14_pair1_pair4_selector_chart_transport_operator_axis_only_alpha14_mod_pi_strict_core_v1.json",
    5: GENERATED / "o15_pair1_pair5_selector_chart_transport_operator_axis_only_alpha15_mod_pi_strict_core_v1.json",
}

OUT_OBJ = (
    GENERATED
    / "selector_output_sign_lift_global_c_v1_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f683_first_exported_selector_output_sign_lift_global_c_v1_rooted_transport_from_s_sel_int_w_break_packet_summary.json"
)


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


def dot(a: list[float], b: list[float]) -> float:
    return sum(float(x) * float(y) for x, y in zip(a, b))


def matvec2(m: list[list[float]], v: list[float]) -> list[float]:
    return [
        float(m[0][0]) * float(v[0]) + float(m[0][1]) * float(v[1]),
        float(m[1][0]) * float(v[0]) + float(m[1][1]) * float(v[1]),
    ]


def matvec12(m: list[list[float]], v: list[float]) -> list[float]:
    return [sum(float(m[i][j]) * float(v[j]) for j in range(12)) for i in range(12)]


def load_transport_matrix(path: Path) -> list[list[float]]:
    obj = load_json(path)
    outs = obj.get("outputs") or {}
    for key in ("O12", "O13", "O14", "O15", "O"):
        M = outs.get(key)
        if isinstance(M, list) and len(M) == 12 and all(is_numeric_list_len(row, 12) for row in M):
            return [[float(x) for x in row] for row in M]
    raise ValueError(f"{path} outputs does not contain a 12x12 transport matrix under keys O12/O13/O14/O15/O")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F647, IN_Y_SEL] + [IN_A[m] for m in sorted(IN_A)] + [IN_O_1M[m] for m in sorted(IN_O_1M)]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "F683",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_OBJ.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJ)
        return

    f647 = load_json(IN_F647)
    y_sel = load_json(IN_Y_SEL)

    # Extract w_break from the exported seed payload.
    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break = payload.get("w_break_by_x")
    if not is_numeric_list_len(w_break, 12):
        raise SystemExit("Invalid F647: constructed_source_object.exported_payload.w_break_by_x must be length-12 numeric list")
    w_break_by_x = [float(v) for v in w_break]

    tol_dot_nonzero = 1e-12
    tol_output_zero = 1e-12
    tol_glue = 1e-9

    # Load u_m and coords.
    u: dict[int, list[float]] = {}
    coords: dict[int, list[float]] = {}
    for m in sorted(IN_A):
        a = load_json(IN_A[m])
        data = a.get("data") or {}
        u_key = "u_1" if m == 1 else f"u_{m}"
        coords_key = "u_1_coords_in_c1_s1" if m == 1 else f"u_{m}_coords_in_c{m}_s{m}"
        u_raw = data.get(u_key)
        coords_raw = data.get(coords_key)
        if not is_numeric_list_len(u_raw, 12):
            raise SystemExit(f"Invalid {IN_A[m].relative_to(REPO)}: data.{u_key} must be length-12 numeric list")
        if not is_numeric_list_len(coords_raw, 2):
            raise SystemExit(f"Invalid {IN_A[m].relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")
        u[m] = [float(v) for v in u_raw]
        coords[m] = [float(v) for v in coords_raw]

    # Root sign on pair1.
    dot_root = dot(w_break_by_x, u[1])
    if abs(dot_root) <= tol_dot_nonzero:
        status = "FAIL_NO_ROOT_SIGN_LIFT_FROM_W_BREAK_ON_PAIR1"
        artifact = {
            "object": "SelectorOutputSignLift_global_C_v1_rooted_transport_from_S_sel_int_w_break_strict_convention_v1",
            "stage": "F683",
            "status": status,
            "as_of": AS_OF,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "root": {"pair": "pair1", "dot_w_break_u_1": float(dot_root)},
            "no_false_pass": True,
        }
        OUT_OBJ.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJ)
        return

    s: dict[int, int] = {1: sign(dot_root)}
    u1_star = [float(s[1]) * x for x in u[1]]

    # Propagate to other charts via rooted transports.
    for m in (2, 3, 4, 5):
        O_1m = load_transport_matrix(IN_O_1M[m])
        v = matvec12(O_1m, u1_star)
        d = dot(v, u[m])
        if abs(d) <= tol_dot_nonzero:
            raise SystemExit(f"FAIL: dot(O_1{m} u1*, u_{m}) is numerically zero; cannot propagate sign")
        s[m] = sign(d)

    # Build chartwise output vectors after applying s_m to the exported coordinates.
    charts: dict[str, Any] = {}
    v_out_by_pair: dict[str, list[float]] = {}
    for m in range(1, 6):
        pair = f"pair{m}"
        y_chart = ((y_sel.get("charts") or {}).get(pair) or {})
        Y = y_chart.get("Y_sel_matrix_in_c_s_to_o")
        if not (isinstance(Y, list) and len(Y) == 2 and all(is_numeric_list_len(row, 2) for row in Y)):
            raise SystemExit(
                f"Invalid {IN_Y_SEL.relative_to(REPO)}: charts.{pair}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
            )
        Y2 = [[float(Y[0][0]), float(Y[0][1])], [float(Y[1][0]), float(Y[1][1])]]
        coords_star = [float(s[m]) * coords[m][0], float(s[m]) * coords[m][1]]
        v_out = matvec2(Y2, coords_star)
        v_out_by_pair[pair] = v_out
        charts[pair] = {
            "chart_id": pair,
            "k": m,
            "coords_in_c_s_original": coords[m],
            "sign_lift_s_m": int(s[m]),
            "coords_in_c_s_after_sign_lift": coords_star,
            "Y_sel_matrix_in_c_s_to_o": Y2,
            "v_out_in_o_plus_o_minus": v_out,
            "o_plus_component_sign": (
                None
                if abs(float(v_out[0])) <= tol_output_zero
                else ("+" if float(v_out[0]) > 0 else "-")
            ),
        }

    ref = v_out_by_pair["pair1"]
    max_abs_diff = max(
        max(abs(float(v_out_by_pair[p][0] - ref[0])), abs(float(v_out_by_pair[p][1] - ref[1])))
        for p in sorted(v_out_by_pair)
    )
    glued = bool(float(max_abs_diff) <= tol_glue)

    status = "PASS_EXPORTED" if glued else "FAIL_NOT_GLUED_WITHIN_TOLERANCE"

    artifact = {
        "object": "SelectorOutputSignLift_global_C_v1_rooted_transport_from_S_sel_int_w_break_strict_convention_v1",
        "stage": "F683",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export an explicit global sign-lift/section choice object rooted at pair1: "
            "fix sign on pair1 from the exported reflection-breaking seed weight w_break and propagate to other charts via rooted transports O_1m, "
            "to obtain chartwise output vectors that glue in the fixed output basis (o_+,o_-). "
            "This is a tracked convention/section choice; it does not claim a physical sign-sensitive orientation datum."
        ),
        "domain": {
            "configuration_space_object_ref": "fundamental_action_reconstruction/generated/c_v1_void_configuration_space_in_local_b_tilde_1_sector_v1.json",
            "configuration_space_object": "C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1",
            "charts_cover": "U_pair1 = ... = U_pair5 = C_v1 (declared global cover)",
        },
        "inputs": {
            "seed_w_break_ref": str(IN_F647.relative_to(REPO)),
            "output_channels_ref": str(IN_Y_SEL.relative_to(REPO)),
            "a_m_refs": {str(m): str(IN_A[m].relative_to(REPO)) for m in sorted(IN_A)},
            "rooted_transport_refs": {str(m): str(IN_O_1M[m].relative_to(REPO)) for m in sorted(IN_O_1M)},
            "boundary": "N512 (no operator-level groupoid promotion; axis-only edges remain projective)",
        },
        "root": {
            "pair": "pair1",
            "dot_w_break_u_1": float(dot_root),
            "root_sign_s_1": int(s[1]),
            "root_sign_source": "w_break",
        },
        "sign_lift": {
            "rule": "rooted_transport_from_w_break_on_pair1",
            "signs_by_pair": {f"pair{m}": int(s[m]) for m in sorted(s)},
            "tolerances": {"dot_nonzero_abs_tol": tol_dot_nonzero},
            "counts_as_strict_physical_orientation_datum": False,
        },
        "output_observable": {
            "basis": ["o_+", "o_-"],
            "reference_chart": "pair1",
            "output_vector_in_o_plus_o_minus": ref,
            "gluing": {"max_abs_diff_to_pair1": float(max_abs_diff), "tol": tol_glue, "glued": glued},
        },
        "charts": charts,
        "hard_limits": [
            "Convention/section choice only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Depends on rooted axis-only transport representatives on some edges; no Aut(Z_12)-invariant canonicity claim.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F683",
        "status": status,
        "object": artifact["object"],
        "root_sign_s_1": int(s[1]),
        "glued": glued,
        "max_abs_diff_to_pair1": float(max_abs_diff),
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT_OBJ.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJ)


if __name__ == "__main__":
    main()

