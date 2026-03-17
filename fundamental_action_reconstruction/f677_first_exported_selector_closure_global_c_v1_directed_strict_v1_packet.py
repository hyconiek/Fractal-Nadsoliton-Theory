#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_OBJECT = GENERATED / "selector_closure_global_c_v1_directed_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f677_first_exported_selector_closure_global_c_v1_directed_strict_v1_packet_summary.json"
)

GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
GLOBAL_TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
GLOBAL_DIRECTED_STATE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"
GLOBAL_OUTPUT_OPERATOR = GENERATED / "selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dot(a: list[float], b: list[float]) -> float:
    return float(sum(float(x) * float(y) for x, y in zip(a, b, strict=True)))


def fourier_c_s_basis_vectors_on_z12(k: int) -> tuple[list[float], list[float]]:
    if k <= 0 or k >= 6:
        raise ValueError("This packet expects k in {1,2,3,4,5} for the Z_12 Fourier-degenerate pairs.")
    n = 12
    norm = math.sqrt(2.0 / n)
    c = [norm * math.cos(2.0 * math.pi * k * x / n) for x in range(n)]
    s = [norm * math.sin(2.0 * math.pi * k * x / n) for x in range(n)]
    return c, s


def is_2x2_matrix(v: Any) -> bool:
    return (
        isinstance(v, list)
        and len(v) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in v)
    )


def sign_of(x: float, tol: float) -> int:
    if abs(x) <= tol:
        return 0
    return 1 if x > 0.0 else -1


def vec_max_abs_diff(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        raise ValueError("vector length mismatch")
    return max(abs(float(a[i]) - float(b[i])) for i in range(len(a)))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    atlas = load_json(GLOBAL_ATLAS)
    transition = load_json(GLOBAL_TRANSITION)
    state = load_json(GLOBAL_DIRECTED_STATE)
    y_global = load_json(GLOBAL_OUTPUT_OPERATOR)

    u_vectors = (((state.get("outputs") or {}).get("u_vectors_directed")) or {}) if isinstance(state, dict) else {}
    charts = (y_global.get("charts") or {}) if isinstance(y_global, dict) else {}

    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    if sorted(list(charts.keys())) != sorted(expected_charts):
        raise ValueError(f"unexpected charts in global output operator object: {sorted(list(charts.keys()))}")

    tol = 1e-12

    chart_payload: dict[str, Any] = {}
    sign_lift_by_pair: dict[str, int] = {}
    v_ref: list[float] | None = None
    max_diff_to_ref = 0.0

    for k, pair in enumerate(expected_charts, start=1):
        u_key = f"u_{k}"
        u = u_vectors.get(u_key)
        if not (isinstance(u, list) and len(u) == 12 and all(isinstance(x, (int, float)) for x in u)):
            raise ValueError(f"missing_or_invalid_{u_key}_in_directed_state_object")
        u = [float(x) for x in u]

        chart = charts.get(pair) or {}
        y = chart.get("Y_sel_matrix_in_c_s_to_o")
        if not is_2x2_matrix(y):
            raise ValueError(f"missing_or_invalid_Y_sel_matrix_in_c_s_to_o_for_{pair}")
        y = [[float(y[i][j]) for j in range(2)] for i in range(2)]

        c_vec, s_vec = fourier_c_s_basis_vectors_on_z12(k)
        c_coord = dot(c_vec, u)
        s_coord = dot(s_vec, u)
        o_plus = float(y[0][0]) * c_coord + float(y[0][1]) * s_coord
        o_minus = float(y[1][0]) * c_coord + float(y[1][1]) * s_coord

        raw_vec = [o_plus, o_minus]
        sgn = sign_of(o_plus, tol=tol)
        if sgn == 0:
            raise ValueError(f"o_plus_is_zero_within_tolerance_for_{pair}_cannot_define_sign_lift")

        # Explicit per-chart sign-lift to canonicalize the directed output vector in the fixed output basis.
        s_lift = int(sgn)
        v_corr = [float(s_lift) * raw_vec[0], float(s_lift) * raw_vec[1]]

        if v_ref is None:
            v_ref = v_corr

        diff = vec_max_abs_diff(v_corr, v_ref)
        max_diff_to_ref = max(max_diff_to_ref, diff)

        sign_lift_by_pair[pair] = s_lift
        chart_payload[pair] = {
            "chart_id": pair,
            "k": k,
            "domain_basis": chart.get("domain_basis"),
            "output_basis": chart.get("output_basis"),
            "u_key": u_key,
            "coords_in_c_s": [c_coord, s_coord],
            "Y_sel_matrix_in_c_s_to_o": y,
            "v_out_raw_in_o_plus_o_minus": raw_vec,
            "raw_o_plus_sign": int(sgn),
            "sign_lift_applied": s_lift,
            "v_out_in_o_plus_o_minus": v_corr,
            "max_abs_diff_to_reference_chart_v_out": diff,
        }

    if v_ref is None:
        raise RuntimeError("internal error: missing reference output vector")

    obj: dict[str, Any] = {
        "object": "SelectorClosure_global_C_v1_directed_strict_v1",
        "stage": "F677",
        "status": "actual_exported_global_directed_selector_closure_object_on_C_v1__explicit_sign_lift__premise_based__no_false_pass",
        "as_of": "2026-03-17",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit global selector closure object on C_v1 in directed (vector-level) scope by packaging the "
            "exported global atlas/transition objects with the exported global directed selector state and the globally "
            "promoted seed-v1 selector output operator. Because the induced directed output vector differs by a chartwise sign "
            "under the fixed output bases (P674 audit), this export makes the required sign-lift/section choice explicit by "
            "exporting per-chart signs s(pair_m) ∈ {±1} and certifying that the sign-lifted output vectors glue to one global "
            "directed output vector in the fixed (o_+,o_-) basis. This export does not claim strict-core selector closure, "
            "kernel-alone/global QW-2191 discharge, or ToE closure."
        ),
        "domain": atlas.get("domain"),
        "closure_scope": {
            "level": "directed_vector_state",
            "closure_type": "global_output_vector_section_on_output_space",
            "sign_gauge": "fixed by explicit per-chart sign_lift s(pair_m) in the declared premise-based scope (no Aut(Z_12)-invariant claim)",
        },
        "inputs": {
            "global_selector_atlas": str(GLOBAL_ATLAS.relative_to(REPO)),
            "global_selector_transition": str(GLOBAL_TRANSITION.relative_to(REPO)),
            "global_directed_selector_state": str(GLOBAL_DIRECTED_STATE.relative_to(REPO)),
            "global_selector_output_operator": str(GLOBAL_OUTPUT_OPERATOR.relative_to(REPO)),
            "boundary": "N512 (no operator-level transition groupoid promotion; closure is vector-section level with explicit sign-lift)",
            "audit_note": "P674/N675 record that without an explicit sign-lift the induced directed output vector is not chart-independent under fixed output bases.",
        },
        "global_atlas_object": atlas.get("object"),
        "global_transition_object": transition.get("object"),
        "global_directed_state_object": state.get("object"),
        "global_output_operator_object": y_global.get("object"),
        "sign_lift": {
            "meaning": "Explicit per-chart sign-lift required to define a directed global closure outcome under fixed output bases.",
            "rule": "s(pair_m) := sign(o_+(pair_m)) with tolerance tol; v_out := s · v_out_raw makes o_+ >= 0 by construction.",
            "tolerance": tol,
            "signs_by_pair": sign_lift_by_pair,
            "note": "This is an explicit section choice; it is not claimed to be Aut(Z_12)-invariant (N462 boundary) and does not imply any kernel-alone QW-2191 discharge.",
        },
        "output_observable": {
            "basis": ["o_+", "o_-"],
            "output_vector_in_o_plus_o_minus": v_ref,
            "meaning": "Global chart-independent directed output vector after applying the explicit sign-lift to the chartwise induced outputs.",
        },
        "charts": chart_payload,
        "well_definedness_certificate": {
            "criterion": "max_abs_diff_across_chartwise_sign_lifted_output_vectors <= tolerance",
            "tolerance": tol,
            "max_abs_diff_to_reference_chart_v_out": max_diff_to_ref,
            "certificate_pass": bool(max_diff_to_ref <= tol),
            "note": "Vector-level equality is asserted only after applying the explicit chartwise sign-lift; this does not promote to operator-level transition groupoid identities (N512).",
        },
        "hard_limits": [
            "premise_based_directed_closure_only (depends on the exported directed state scope; no Aut(Z_12)-invariant canonicity claim)",
            "explicit_sign_lift_required (P674/N675 boundary on raw outputs under fixed output bases)",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_strict_core_selector_closure",
            "no_global_kernel_alone_QW2191_discharge",
            "no_emergent_observer_construction",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary: dict[str, Any] = {
        "stage": "F677",
        "lane": "first_exported_selector_closure_global_c_v1_directed_strict_v1_packet_only",
        "goal": "export_one_global_directed_selector_closure_object_with_explicit_sign_lift_on_C_v1",
        "status": "F677_EXECUTED_FIRST_EXPORTED_SELECTOR_CLOSURE_GLOBAL_C_V1_DIRECTED_STRICT_V1_PACKET_NO_FALSE_PASS",
        "exported_object": obj["object"],
        "exported_object_file": str(OUT_OBJECT.relative_to(REPO)),
        "chart_count": len(expected_charts),
        "certificate_pass": obj["well_definedness_certificate"]["certificate_pass"],
        "max_abs_diff_to_reference_chart_v_out": max_diff_to_ref,
        "tolerance": tol,
        "strict_core_selector_closure": False,
        "QW2191_discharge": False,
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJECT)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

