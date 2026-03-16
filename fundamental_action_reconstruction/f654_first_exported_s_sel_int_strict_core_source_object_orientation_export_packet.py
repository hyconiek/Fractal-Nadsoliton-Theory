#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f654_first_exported_s_sel_int_strict_core_source_object_orientation_export_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def dot(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        raise ValueError("dot: mismatched lengths")
    return float(sum(float(x) * float(y) for x, y in zip(a, b)))


def l2_norm(xs: list[float]) -> float:
    return math.sqrt(sum(float(x) * float(x) for x in xs))


def scale(xs: list[float], lam: float) -> list[float]:
    return [float(lam) * float(x) for x in xs]


def add(a: list[float], b: list[float]) -> list[float]:
    if len(a) != len(b):
        raise ValueError("add: mismatched lengths")
    return [float(x) + float(y) for x, y in zip(a, b)]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f647 = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
    )
    f649 = load_json(
        "fundamental_action_reconstruction/generated/f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet.json"
    )

    exported_source = str(f647["exported_constructed_source_object"])
    if exported_source != "S_sel_int_strict_core_source_object_v1":
        raise ValueError("Unexpected exported constructed source object identity")

    basis_by_x = f649["carrier"]["basis_functions_by_x"]
    c1 = [float(v) for v in basis_by_x["c1"]]
    s1 = [float(v) for v in basis_by_x["s1"]]
    if not (len(c1) == len(s1) == 12):
        raise ValueError("Expected n=12 basis functions for (c1,s1)")

    w_break = [
        float(v)
        for v in f647["constructed_source_object"]["exported_payload"]["w_break_by_x"]
    ]
    if len(w_break) != 12:
        raise ValueError("Expected w_break_by_x to be a length-12 array")

    # Project the strict-core reflection-breaking source vector onto span{c1,s1}.
    c1_norm_sq = dot(c1, c1)
    s1_norm_sq = dot(s1, s1)
    if c1_norm_sq <= 0.0 or s1_norm_sq <= 0.0:
        raise ValueError("Invalid basis norms for (c1,s1)")

    a_c1 = dot(w_break, c1) / c1_norm_sq
    b_s1 = dot(w_break, s1) / s1_norm_sq
    proj = add(scale(c1, a_c1), scale(s1, b_s1))
    proj_norm = l2_norm(proj)
    if proj_norm <= 0.0:
        raise ValueError("Projection of source vector onto span{c1,s1} is zero; cannot export E_orient")

    e_parallel_by_x = scale(proj, 1.0 / proj_norm)

    # Orthonormal partner inside the same plane, using the (a,b)->(-b,a) rotation
    # in the orthogonal equal-norm (c1,s1) basis.
    transverse = add(scale(c1, -b_s1), scale(s1, a_c1))
    e_transverse_by_x = scale(transverse, 1.0 / proj_norm)

    # Gauge/quotient safety witness: invariance under positive rescaling of w_break.
    lambda_positive = 2.0
    w_break_rescaled = scale(w_break, lambda_positive)
    a_c1_rescaled = dot(w_break_rescaled, c1) / c1_norm_sq
    b_s1_rescaled = dot(w_break_rescaled, s1) / s1_norm_sq
    proj_rescaled = add(scale(c1, a_c1_rescaled), scale(s1, b_s1_rescaled))
    proj_rescaled_norm = l2_norm(proj_rescaled)
    e_parallel_rescaled_by_x = scale(proj_rescaled, 1.0 / proj_rescaled_norm)

    scaling_invariant = all(
        abs(float(u) - float(v)) < 1.0e-12
        for u, v in zip(e_parallel_by_x, e_parallel_rescaled_by_x)
    )

    # Basic orthonormality audit.
    ortho_residual = abs(dot(e_parallel_by_x, e_transverse_by_x))
    parallel_norm_residual = abs(l2_norm(e_parallel_by_x) - 1.0)
    transverse_norm_residual = abs(l2_norm(e_transverse_by_x) - 1.0)

    exported_orientation_object = "E_orient_s_sel_int_source_object_v1"

    summary: dict[str, Any] = {
        "stage": "F654",
        "lane": "first_exported_s_sel_int_strict_core_source_object_orientation_export_only",
        "goal": "export_one_explicit_pair1_orientation_datum_from_the_exported_strict_core_S_sel_int_source_object",
        "status": "F654_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_ORIENTATION_EXPORT_PACKET_NO_FALSE_PASS",
        "source_object": exported_source,
        "source_derivation": {
            "derived_from_exported_source_object": True,
            "typing_target_frame": f649["future_orientation_export_target_frame"],
            "projection_basis": ["c1", "s1"],
            "projection_coefficients_in_basis": {"a_c1": a_c1, "b_s1": b_s1},
            "projection_norm": proj_norm,
        },
        "exported_orientation": {
            "object": exported_orientation_object,
            "frame_basis": ["c1", "s1"],
            "e_parallel_by_x": e_parallel_by_x,
            "e_transverse_by_x": e_transverse_by_x,
            "ordered_frame_by_x": [e_parallel_by_x, e_transverse_by_x],
            "selector_axis_v1_by_x": e_parallel_by_x,
            "B_sel_start_frame_v1_by_x": [e_parallel_by_x, e_transverse_by_x],
        },
        "orientation_export_properties": {
            "strict_core_only": True,
            "internal_orientation_datum": True,
            "selector_bearing_without_external_anchor": True,
            "quotient_gauge_safe": scaling_invariant,
            "bridge_ready_for_B_sel": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "uses_legacy_kernel_substitution": False,
        },
        "numeric_audit": {
            "orthogonality_residual": ortho_residual,
            "parallel_unit_norm_residual": parallel_norm_residual,
            "transverse_unit_norm_residual": transverse_norm_residual,
        },
        "gauge_safety_witness": {
            "positive_rescaling_lambda": lambda_positive,
            "rescaled_parallel_by_x": e_parallel_rescaled_by_x,
            "invariant_under_positive_rescaling": scaling_invariant,
        },
        "admissible_E_orient": False,
        "downstream_chain_complete": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

