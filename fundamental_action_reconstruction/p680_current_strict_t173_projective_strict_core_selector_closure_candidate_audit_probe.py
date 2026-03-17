#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

N676_SUMMARY = (
    GENERATED
    / "n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json"
)
N546_SUMMARY = (
    GENERATED
    / "n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
)
F649_SUMMARY = (
    GENERATED
    / "f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet_summary.json"
)
F672_OBJECT = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"

OUT_JSON = (
    GENERATED
    / "p680_current_strict_t173_projective_strict_core_selector_closure_candidate_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p680_current_strict_t173_projective_strict_core_selector_closure_candidate_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_2x2_matrix(v: Any) -> bool:
    return (
        isinstance(v, list)
        and len(v) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in v)
        and all(isinstance(x, (int, float)) for row in v for x in row)
    )


def matmul_2x2(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]


def max_abs_2x2(a: list[list[float]]) -> float:
    return max(abs(a[0][0]), abs(a[0][1]), abs(a[1][0]), abs(a[1][1]))


def matsub_2x2(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [a[0][0] - b[0][0], a[0][1] - b[0][1]],
        [a[1][0] - b[1][0], a[1][1] - b[1][1]],
    ]


def eigvals_symmetric_2x2(m: list[list[float]]) -> tuple[float, float]:
    tr = m[0][0] + m[1][1]
    det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
    disc = tr * tr - 4.0 * det
    disc = max(disc, 0.0)
    s = math.sqrt(disc)
    return (0.5 * (tr + s), 0.5 * (tr - s))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {
        "N676_summary": N676_SUMMARY,
        "N546_summary": N546_SUMMARY,
        "F649_summary": F649_SUMMARY,
        "F672_projective_closure_object": F672_OBJECT,
    }
    missing = {k: str(v.relative_to(REPO)) for k, v in required.items() if not v.exists()}
    if missing:
        status = "P680_NOT_COMPUTABLE_MISSING_DEPENDENCY_ARTIFACTS"
        artifact = {
            "stage": "P680",
            "lane": "current_strict_t173_projective_strict_core_selector_closure_candidate_audit_probe_only",
            "status": status,
            "missing": missing,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }
        summary = {
            "stage": "P680",
            "status": status,
            "missing": missing,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    n676 = load_json(N676_SUMMARY)
    n546 = load_json(N546_SUMMARY)
    f649 = load_json(F649_SUMMARY)
    f672 = load_json(F672_OBJECT)

    checks: list[dict[str, Any]] = []
    blocking_mismatches: list[str] = []

    def check_equal(id_: str, actual: Any, expected: Any, meaning: str) -> None:
        ok = actual == expected
        checks.append(
            {"id": id_, "actual": actual, "expected": expected, "pass": ok, "meaning": meaning}
        )
        if not ok:
            blocking_mismatches.append(id_)

    def check_pred(id_: str, ok: bool, actual: Any, expected: Any, meaning: str) -> None:
        checks.append(
            {"id": id_, "actual": actual, "expected": expected, "pass": bool(ok), "meaning": meaning}
        )
        if not ok:
            blocking_mismatches.append(id_)

    check_equal(
        "s_sel_int_admissible_in_f34_sense",
        bool((n676.get("theorem_result", {}) or {}).get("admissible_S_sel_int_source_object_in_F34_sense")),
        True,
        "An admissible strict-core S_sel_int source object exists in the narrow F34 sense (N676).",
    )
    check_equal(
        "e_orient_admissible",
        bool((n546.get("theorem_result", {}) or {}).get("admissible_E_orient")),
        True,
        "An admissible strict-core orientation export exists from the S_sel_int source object (N546).",
    )

    state_support = (f649.get("state_support", {}) or {})
    dot = state_support.get("dot_w_break_with_s1")
    nonzero_s1_support = bool(state_support.get("nonzero_s1_support"))
    check_equal(
        "source_object_nonzero_s1_support",
        nonzero_s1_support,
        True,
        "Reflection-breaking witness is non-degenerate on the declared pair1 sine axis (F649).",
    )
    dot_ok = isinstance(dot, (int, float)) and abs(float(dot)) > 0.0
    check_pred(
        "dot_w_break_with_s1_nonzero",
        dot_ok,
        dot,
        "nonzero scalar",
        "The exported strict-core reflection-breaking weight has nonzero dot with s1 (prevents degenerate orientation).",
    )

    check_equal(
        "projective_closure_object_name",
        f672.get("object"),
        "SelectorClosure_global_C_v1_projective_strict_v1",
        "Global projective selector closure object has the expected typed name (F672).",
    )
    check_equal(
        "projective_closure_no_false_pass_marked",
        bool(f672.get("no_false_pass")),
        True,
        "Global projective selector closure object is explicitly marked no_false_pass.",
    )
    check_equal(
        "projective_closure_scope",
        (f672.get("closure_scope", {}) or {}).get("level"),
        "projective_ray_state",
        "Closure object is declared in projective/ray-level scope (sign gauge treated as projector-level).",
    )

    cert = (f672.get("well_definedness_certificate", {}) or {})
    check_equal(
        "closure_well_definedness_certificate_pass",
        bool(cert.get("certificate_pass")),
        True,
        "Closure object reports passing chartwise well-definedness certificate on exported data.",
    )
    tol = cert.get("tolerance")
    maxdiff = cert.get("max_abs_diff_to_reference_chart_output_projector")
    tol_ok = isinstance(tol, (int, float)) and float(tol) > 0.0
    maxdiff_ok = isinstance(maxdiff, (int, float)) and float(maxdiff) >= 0.0
    check_pred(
        "closure_certificate_tolerance_positive",
        tol_ok,
        tol,
        "> 0 scalar",
        "Certificate tolerance is a positive scalar.",
    )
    check_pred(
        "closure_certificate_maxdiff_nonnegative",
        maxdiff_ok,
        maxdiff,
        ">= 0 scalar",
        "Certificate max-diff is a nonnegative scalar.",
    )
    if tol_ok and maxdiff_ok:
        check_pred(
            "closure_certificate_maxdiff_within_tolerance",
            float(maxdiff) <= float(tol),
            {"maxdiff": float(maxdiff), "tolerance": float(tol)},
            "maxdiff <= tolerance",
            "Chartwise output projectors agree within the declared certificate tolerance.",
        )

    b = ((f672.get("output_observable", {}) or {}).get("output_projector_matrix_in_o_plus_o_minus"))
    if is_2x2_matrix(b):
        b2: list[list[float]] = [[float(b[0][0]), float(b[0][1])], [float(b[1][0]), float(b[1][1])]]
        sym_res = abs(b2[0][1] - b2[1][0])
        b_sq = matmul_2x2(b2, b2)
        idem_res = max_abs_2x2(matsub_2x2(b_sq, b2))
        tr = b2[0][0] + b2[1][1]
        det = b2[0][0] * b2[1][1] - b2[0][1] * b2[1][0]
        lam1, lam2 = eigvals_symmetric_2x2(b2)
        tol_num = 1e-10
        check_pred(
            "output_projector_symmetric",
            sym_res <= tol_num,
            sym_res,
            f"<= {tol_num}",
            "Output observable matrix is symmetric (projector sanity).",
        )
        check_pred(
            "output_projector_idempotent",
            idem_res <= tol_num,
            idem_res,
            f"<= {tol_num}",
            "Output observable matrix satisfies P^2≈P (projector sanity).",
        )
        check_pred(
            "output_projector_trace_one",
            abs(tr - 1.0) <= tol_num,
            tr,
            f"1±{tol_num}",
            "Output observable has trace≈1 (rank‑1 projector sanity).",
        )
        check_pred(
            "output_projector_det_zero",
            abs(det) <= tol_num,
            det,
            f"0±{tol_num}",
            "Output observable has det≈0 (rank‑1 projector sanity).",
        )
        check_pred(
            "output_projector_eigenvalues_near_1_and_0",
            (abs(lam1 - 1.0) <= tol_num and abs(lam2 - 0.0) <= tol_num)
            or (abs(lam2 - 1.0) <= tol_num and abs(lam1 - 0.0) <= tol_num),
            {"lambda_1": lam1, "lambda_2": lam2},
            f"{{1,0}} within {tol_num}",
            "Output observable eigenvalues are consistent with a rank‑1 projector (unique output ray in projective scope).",
        )
    else:
        check_pred(
            "output_projector_matrix_present_2x2",
            False,
            b,
            "2x2 numeric matrix",
            "Closure object exports a 2x2 output projector observable in the (o_+,o_-) basis.",
        )

    closure_candidate_supported = len(blocking_mismatches) == 0
    status = (
        "PASS_CLOSURE_CANDIDATE_SUPPORTED_PROJECTIVE_STRICT_CORE_ONLY"
        if closure_candidate_supported
        else "P680_REQUIRES_REVIEW_CLOSURE_CANDIDATE_INSUFFICIENT_OR_INCONSISTENT"
    )

    artifact = {
        "stage": "P680",
        "lane": "current_strict_t173_projective_strict_core_selector_closure_candidate_audit_probe_only",
        "goal": "audit_whether_current_exports_support_a_projective_strict_core_selector_closure_discharge_attempt_under_T173__probe_only",
        "dependencies": {
            "N676": str(N676_SUMMARY.relative_to(REPO)),
            "N546": str(N546_SUMMARY.relative_to(REPO)),
            "F649": str(F649_SUMMARY.relative_to(REPO)),
            "F672": str(F672_OBJECT.relative_to(REPO)),
        },
        "status": status,
        "closure_candidate_supported": bool(closure_candidate_supported),
        "checks": checks,
        "blocking_mismatches": blocking_mismatches,
        "hard_limits": [
            "probe_only_no_promotion_to_theorem",
            "strict_core_selector_closure_remains_false_until_theorem_exported",
            "no_directed_sign_sensitive_orientation_claim",
            "no_kernel_alone_global_QW2191_discharge",
            "no_ToE_closure",
        ],
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }
    summary = {
        "stage": "P680",
        "status": status,
        "closure_candidate_supported": bool(closure_candidate_supported),
        "blocking_mismatches": blocking_mismatches,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

