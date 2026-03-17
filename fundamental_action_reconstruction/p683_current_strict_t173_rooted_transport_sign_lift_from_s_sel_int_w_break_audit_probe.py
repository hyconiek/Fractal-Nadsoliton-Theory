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

OUT = GENERATED / "p683_current_strict_t173_rooted_transport_sign_lift_from_s_sel_int_w_break_audit_probe.json"
OUT_SUMMARY = GENERATED / "p683_current_strict_t173_rooted_transport_sign_lift_from_s_sel_int_w_break_audit_probe_summary.json"


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


def l2_norm(v: list[float]) -> float:
    return math.sqrt(sum(float(x) * float(x) for x in v))


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
            "stage": "P683",
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

    try:
        payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
        w_break = payload.get("w_break_by_x")
        if not is_numeric_list_len(w_break, 12):
            raise ValueError("F647.constructed_source_object.exported_payload.w_break_by_x must be length-12 numeric list")
        w_break_by_x = [float(v) for v in w_break]
    except Exception as exc:
        artifact = {
            "stage": "P683",
            "status": "INVALID_W_BREAK_INPUT_SHAPE",
            "as_of": AS_OF,
            "error": str(exc),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    tol_dot_nonzero = 1e-12
    tol_output = 1e-9
    tol_alignment = 1e-9

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

    # Root sign from w_break on pair1.
    dot_root = dot(w_break_by_x, u[1])
    dot_root_nonzero = bool(abs(dot_root) > tol_dot_nonzero)
    if not dot_root_nonzero:
        status = "FAIL_NO_ROOT_SIGN_LIFT_FROM_W_BREAK_ON_PAIR1"
        artifact = {
            "stage": "P683",
            "as_of": AS_OF,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "status": status,
            "root": {
                "pair": "pair1",
                "dot_w_break_u_1": float(dot_root),
                "dot_nonzero": dot_root_nonzero,
            },
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    s_root = sign(dot_root)
    u1_star = [float(s_root) * x for x in u[1]]
    coords_star: dict[int, list[float]] = {1: [float(s_root) * coords[1][0], float(s_root) * coords[1][1]]}

    per_pair: dict[str, Any] = {
        "pair1": {
            "pair": "pair1",
            "u_ref": str(IN_A[1].relative_to(REPO)),
            "dot_w_break_u_1": float(dot_root),
            "root_sign_s_1": int(s_root),
            "root_sign_source": "w_break",
        }
    }

    supported_pairs = ["pair1"]
    missing_pairs: list[str] = []

    # Propagate sign along O_1m edges by aligning O_1m u1* with u_m.
    for m in (2, 3, 4, 5):
        pair = f"pair{m}"
        O_1m = load_transport_matrix(IN_O_1M[m])
        v = matvec12(O_1m, u1_star)
        dot_vm_um = dot(v, u[m])
        dot_nonzero = bool(abs(dot_vm_um) > tol_dot_nonzero)
        if not dot_nonzero:
            s_m = None
            missing_pairs.append(pair)
        else:
            s_m = sign(dot_vm_um)
            supported_pairs.append(pair)
            coords_star[m] = [float(s_m) * coords[m][0], float(s_m) * coords[m][1]]

        alignment_l2 = None
        if s_m is not None:
            diff = [v[i] - float(s_m) * u[m][i] for i in range(12)]
            alignment_l2 = float(l2_norm(diff))

        per_pair[pair] = {
            "pair": pair,
            "u_ref": str(IN_A[m].relative_to(REPO)),
            "transport_ref": str(IN_O_1M[m].relative_to(REPO)),
            "dot_transport_u1_star_with_u_m": float(dot_vm_um),
            "dot_nonzero": dot_nonzero,
            "propagated_sign_s_m": s_m,
            "alignment_l2_norm_O1m_u1_star_minus_s_m_u_m": alignment_l2,
        }

    all_pairs_supported = len(missing_pairs) == 0

    # Sanity: compute output vectors after applying the rooted sign lift.
    output_signs: dict[str, str | None] = {}
    for m in sorted(IN_A):
        pair = f"pair{m}"
        y_chart = ((y_sel.get("charts") or {}).get(pair) or {})
        Y = y_chart.get("Y_sel_matrix_in_c_s_to_o")
        if not (isinstance(Y, list) and len(Y) == 2 and all(is_numeric_list_len(row, 2) for row in Y)):
            raise SystemExit(
                f"Invalid {IN_Y_SEL.relative_to(REPO)}: charts.{pair}.Y_sel_matrix_in_c_s_to_o must be a 2x2 numeric list"
            )
        Y2 = [[float(Y[0][0]), float(Y[0][1])], [float(Y[1][0]), float(Y[1][1])]]

        if m not in coords_star:
            per_pair[pair]["output_vector_o_plus_o_minus_after_rooted_sign_lift"] = None
            per_pair[pair]["output_o_plus_component_sign_after_rooted_sign_lift"] = None
            output_signs[pair] = None
            continue

        y_vec = matvec2(Y2, coords_star[m])
        sign_o_plus = None if abs(float(y_vec[0])) <= tol_output else ("+" if float(y_vec[0]) > 0 else "-")
        per_pair[pair]["output_vector_o_plus_o_minus_after_rooted_sign_lift"] = y_vec
        per_pair[pair]["output_o_plus_component_sign_after_rooted_sign_lift"] = sign_o_plus
        output_signs[pair] = sign_o_plus

    directed_output_sign_consistent = None
    if all_pairs_supported:
        directed_output_sign_consistent = len(set(output_signs[p] for p in sorted(output_signs))) == 1

    # Optional: ensure alignment residuals are tiny where defined.
    alignment_ok = True
    for m in (2, 3, 4, 5):
        pair = f"pair{m}"
        s_m = per_pair[pair]["propagated_sign_s_m"]
        if s_m is None:
            continue
        if float(per_pair[pair]["alignment_l2_norm_O1m_u1_star_minus_s_m_u_m"]) > tol_alignment:
            alignment_ok = False

    if not all_pairs_supported:
        status = "PARTIAL_ROOTED_TRANSPORT_SIGN_LIFT_SUPPORTED_ONLY_ON_SUBSET_OF_CHARTS"
    elif not bool(alignment_ok):
        status = "FAIL_TRANSPORT_ALIGNMENT_RESIDUAL_TOO_LARGE_FOR_ROOTED_SIGN_LIFT"
    elif bool(directed_output_sign_consistent):
        status = "PASS_ROOTED_TRANSPORT_SIGN_LIFT_SUPPORTED_ON_ALL_CHARTS_AND_OUTPUT_SIGN_CONSISTENT"
    else:
        status = "FAIL_OUTPUT_SIGN_MISMATCH_UNDER_ROOTED_TRANSPORT_SIGN_LIFT"

    artifact = {
        "stage": "P683",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit a rooted sign-lift: fix sign on pair1 using exported S_sel_int seed w_break, "
            "then propagate sign to pair2..pair5 by aligning exported transports O_1m u1* with exported representatives u_m, "
            "and check chartwise output-sign consistency for the exported global output channels"
        ),
        "inputs": {
            "f647_ref": str(IN_F647.relative_to(REPO)),
            "y_sel_ref": str(IN_Y_SEL.relative_to(REPO)),
            "a_m_refs": {str(m): str(IN_A[m].relative_to(REPO)) for m in sorted(IN_A)},
            "rooted_transport_refs": {str(m): str(IN_O_1M[m].relative_to(REPO)) for m in sorted(IN_O_1M)},
            "tolerances": {
                "dot_nonzero_abs_tol": tol_dot_nonzero,
                "output_component_zero_tol": tol_output,
                "alignment_l2_tol": tol_alignment,
            },
        },
        "root": {
            "pair": "pair1",
            "dot_w_break_u_1": float(dot_root),
            "root_sign_s_1": int(s_root),
            "root_sign_source": "w_break",
        },
        "per_pair": per_pair,
        "result": {
            "status": status,
            "supported_pairs": supported_pairs,
            "missing_pairs": missing_pairs,
            "all_pairs_supported": all_pairs_supported,
            "transport_alignment_ok": bool(alignment_ok),
            "directed_output_sign_consistent_if_all_pairs_supported": directed_output_sign_consistent,
            "counts_as_strict_physical_orientation_datum": False,
            "note": (
                "PASS here means only: given the exported axis-only transport representatives, one can choose a consistent vector-section sign "
                "once the root sign is fixed. Axis-only edges remain projective by construction; no physical sign/orientation claim is made."
            ),
        },
        "hard_limits": [
            "Probe-only: does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Axis-only transport edges are only projector-level by design; using them for oriented vector transport is a section/convention choice.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
            "Does not promote projector-level cocycle into an operator-level groupoid (N512 boundary).",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P683",
        "status": status,
        "root_sign_s_1": int(s_root),
        "supported_pairs": supported_pairs,
        "missing_pairs": missing_pairs,
        "all_pairs_supported": all_pairs_supported,
        "transport_alignment_ok": bool(alignment_ok),
        "directed_output_sign_consistent_if_all_pairs_supported": directed_output_sign_consistent,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

