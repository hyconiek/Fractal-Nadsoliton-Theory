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

IN_PROJECTIVE_CLOSURE_F672 = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"
IN_DIRECTED_CLOSURE_F677 = GENERATED / "selector_closure_global_c_v1_directed_strict_v1.json"
IN_DIRECTED_CLOSURE_F685 = (
    GENERATED / "selector_closure_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)
IN_DIRECTED_CLOSURE_F692 = (
    GENERATED
    / "selector_closure_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)

OUT_JSON = GENERATED / "p692_current_strict_t173_projective_vs_directed_closure_output_compatibility_audit_probe.json"
OUT_SUMMARY = (
    GENERATED
    / "p692_current_strict_t173_projective_vs_directed_closure_output_compatibility_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_2x2_numeric_matrix(obj: Any) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in obj)
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for row in obj for v in row)
    )


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


def max_abs_diff_2x2(a: list[list[float]], b: list[list[float]]) -> float:
    return max(
        abs(float(a[0][0] - b[0][0])),
        abs(float(a[0][1] - b[0][1])),
        abs(float(a[1][0] - b[1][0])),
        abs(float(a[1][1] - b[1][1])),
    )


def projector_from_vector(v: list[float]) -> list[list[float]]:
    n2 = float(v[0]) * float(v[0]) + float(v[1]) * float(v[1])
    if n2 <= 0.0 or not math.isfinite(n2):
        raise ValueError("Invalid vector norm for projector construction")
    inv = 1.0 / n2
    return [
        [float(v[0]) * float(v[0]) * inv, float(v[0]) * float(v[1]) * inv],
        [float(v[1]) * float(v[0]) * inv, float(v[1]) * float(v[1]) * inv],
    ]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_PROJECTIVE_CLOSURE_F672,
        IN_DIRECTED_CLOSURE_F677,
        IN_DIRECTED_CLOSURE_F685,
        IN_DIRECTED_CLOSURE_F692,
    ]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P692",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_JSON)
        return

    proj = load_json(IN_PROJECTIVE_CLOSURE_F672)
    b = ((proj.get("output_observable") or {}).get("output_projector_matrix_in_o_plus_o_minus"))
    if not is_2x2_numeric_matrix(b):
        artifact = {
            "stage": "P692",
            "status": "INVALID_PROJECTIVE_CLOSURE_SHAPE",
            "as_of": AS_OF,
            "error": "F672 must export output_observable.output_projector_matrix_in_o_plus_o_minus as 2x2 numeric matrix",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_JSON)
        return
    B = [[float(b[0][0]), float(b[0][1])], [float(b[1][0]), float(b[1][1])]]

    directed_inputs = {
        "F677_premise_based_directed_closure": IN_DIRECTED_CLOSURE_F677,
        "F685_w_break_rooted_directed_closure": IN_DIRECTED_CLOSURE_F685,
        "F692_sign_fixed_directed_closure": IN_DIRECTED_CLOSURE_F692,
    }

    tol_projector_match = 1e-10
    tol_invariance = 1e-10

    per_directed: dict[str, Any] = {}
    errors: list[dict[str, Any]] = []
    projectors: dict[str, list[list[float]]] = {}
    vectors: dict[str, list[float]] = {}

    for label, path in directed_inputs.items():
        try:
            obj = load_json(path)
            v = ((obj.get("output_observable") or {}).get("output_vector_in_o_plus_o_minus"))
            if not is_numeric_list_len(v, 2):
                raise ValueError("Missing output_observable.output_vector_in_o_plus_o_minus (length-2 numeric list)")
            v2 = [float(v[0]), float(v[1])]
            vectors[label] = v2

            P = projector_from_vector(v2)
            projectors[label] = P

            diff_proj = float(max_abs_diff_2x2(P, B))
            Bv = matvec2(B, v2)
            diff_invariance = float(max_abs_diff2(Bv, v2))
            norm = float(math.sqrt(v2[0] * v2[0] + v2[1] * v2[1]))

            per_directed[label] = {
                "directed_closure_ref": str(path.relative_to(REPO)),
                "directed_closure_object": str(obj.get("object") or ""),
                "output_vector_in_o_plus_o_minus": v2,
                "output_vector_norm": norm,
                "derived_output_projector_from_vector": P,
                "max_abs_diff_to_projective_output_projector": diff_proj,
                "max_abs_diff_of_Bv_minus_v": diff_invariance,
                "checks": {
                    "tol_projector_match": tol_projector_match,
                    "tol_invariance": tol_invariance,
                    "projector_matches_projective": bool(diff_proj <= tol_projector_match),
                    "vector_is_in_projective_ray": bool(diff_invariance <= tol_invariance),
                },
            }
        except Exception as exc:
            errors.append({"label": label, "path": str(path.relative_to(REPO)), "error": str(exc)})

    # Pairwise projector consistency across directed closures (sign-invariant).
    labels = sorted(projectors.keys())
    pairwise_proj_maxdiff: dict[str, float] = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a = labels[i]
            b = labels[j]
            pairwise_proj_maxdiff[f"{a}__vs__{b}"] = float(max_abs_diff_2x2(projectors[a], projectors[b]))

    all_projectors_match_projective = all(
        bool(per_directed.get(lbl, {}).get("checks", {}).get("projector_matches_projective"))
        for lbl in directed_inputs
        if lbl in per_directed
    )
    all_vectors_in_projective_ray = all(
        bool(per_directed.get(lbl, {}).get("checks", {}).get("vector_is_in_projective_ray"))
        for lbl in directed_inputs
        if lbl in per_directed
    )
    pairwise_projector_consistent = all(v <= tol_projector_match for v in pairwise_proj_maxdiff.values())

    if errors:
        status = "NOT_COMPUTABLE_ERRORS_PRESENT"
    elif all_projectors_match_projective and all_vectors_in_projective_ray and pairwise_projector_consistent:
        status = "PASS_PROJECTIVE_AND_DIRECTED_CLOSURE_OUTPUTS_DEFINE_THE_SAME_OUTPUT_RAY"
    else:
        status = "FAIL_PROJECTIVE_AND_DIRECTED_CLOSURE_OUTPUTS_NOT_COMPATIBLE_WITHIN_TOLERANCE"

    artifact = {
        "stage": "P692",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_that_exported_directed_closure_output_vectors_induce_the_same_rank_one_output_projector_as_the_exported_projective_closure_output_projector__no_false_pass",
        "inputs": {
            "projective_closure_ref": str(IN_PROJECTIVE_CLOSURE_F672.relative_to(REPO)),
            "directed_closures": {k: str(v.relative_to(REPO)) for k, v in directed_inputs.items()},
        },
        "tolerances": {
            "tol_projector_match_max_abs": tol_projector_match,
            "tol_invariance_max_abs": tol_invariance,
        },
        "projective_output_projector_in_o_plus_o_minus": B,
        "directed_outputs": per_directed,
        "pairwise_directed_projector_max_abs_diffs": pairwise_proj_maxdiff,
        "result": {
            "all_projectors_match_projective": bool(all_projectors_match_projective),
            "all_vectors_in_projective_ray": bool(all_vectors_in_projective_ray),
            "pairwise_projector_consistent_across_directed_closures": bool(pairwise_projector_consistent),
            "error_count": len(errors),
            "errors": errors,
        },
        "hard_limits": [
            "Does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
            "Projector/ray-level compatibility audit only (no operator-level groupoid promotion).",
        ],
        "no_false_pass": True,
        "status": status,
    }

    summary = {
        "stage": "P692",
        "status": status,
        "all_projectors_match_projective": bool(all_projectors_match_projective),
        "pairwise_projector_consistent": bool(pairwise_projector_consistent),
        "error_count": len(errors),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

