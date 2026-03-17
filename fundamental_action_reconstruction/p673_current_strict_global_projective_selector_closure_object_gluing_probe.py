#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_OBJECT = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"

OUT_JSON = GENERATED / "p673_current_strict_global_projective_selector_closure_object_gluing_probe.json"
OUT_SUMMARY = (
    GENERATED
    / "p673_current_strict_global_projective_selector_closure_object_gluing_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_2x2_matrix(v: Any) -> bool:
    return (
        isinstance(v, list)
        and len(v) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in v)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    checks: list[dict[str, Any]] = []
    blocking_mismatches: list[str] = []

    if not IN_OBJECT.exists():
        status = "P673_NOT_COMPUTABLE_MISSING_SELECTOR_CLOSURE_OBJECT_FILE"
        artifact: dict[str, Any] = {
            "stage": "P673",
            "lane": "current_strict_global_projective_selector_closure_object_gluing_probe_only",
            "status": status,
            "missing": str(IN_OBJECT.relative_to(REPO)),
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P673",
            "lane": artifact["lane"],
            "status": status,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    obj = load_json(IN_OBJECT)

    def check_equal(id_: str, actual: Any, expected: Any, meaning: str) -> None:
        ok = actual == expected
        checks.append({"id": id_, "actual": actual, "expected": expected, "pass": ok, "meaning": meaning})
        if not ok:
            blocking_mismatches.append(id_)

    check_equal(
        "object_name",
        obj.get("object"),
        "SelectorClosure_global_C_v1_projective_strict_v1",
        "Closure object has the expected typed name.",
    )
    check_equal(
        "no_false_pass_marked",
        bool(obj.get("no_false_pass")),
        True,
        "Closure object is explicitly marked no_false_pass.",
    )
    check_equal(
        "closure_scope_projective",
        (obj.get("closure_scope", {}) or {}).get("level"),
        "projective_ray_state",
        "Closure object is declared in projective/ray-level scope.",
    )

    cert = obj.get("well_definedness_certificate", {}) or {}
    tol = cert.get("tolerance")
    maxdiff = cert.get("max_abs_diff_to_reference_chart_output_projector")
    cert_pass = cert.get("certificate_pass")

    check_equal(
        "certificate_present",
        bool(cert),
        True,
        "Closure object contains an explicit well-definedness certificate payload.",
    )
    check_equal(
        "certificate_pass",
        bool(cert_pass),
        True,
        "Projector-level well-definedness certificate passes on exported data.",
    )

    # Minimal numeric sanity: certificate should provide scalar tol and maxdiff.
    tol_ok = isinstance(tol, (int, float)) and float(tol) > 0.0
    checks.append(
        {
            "id": "certificate_tolerance_scalar_positive",
            "actual": tol,
            "expected": "> 0 scalar",
            "pass": tol_ok,
            "meaning": "Certificate tolerance is a positive scalar.",
        }
    )
    if not tol_ok:
        blocking_mismatches.append("certificate_tolerance_scalar_positive")

    maxdiff_ok = isinstance(maxdiff, (int, float)) and float(maxdiff) >= 0.0
    checks.append(
        {
            "id": "certificate_maxdiff_scalar_nonnegative",
            "actual": maxdiff,
            "expected": ">= 0 scalar",
            "pass": maxdiff_ok,
            "meaning": "Certificate max-diff is a nonnegative scalar.",
        }
    )
    if not maxdiff_ok:
        blocking_mismatches.append("certificate_maxdiff_scalar_nonnegative")

    # References existence check (no implicit overlap semantics).
    refs = obj.get("inputs", {}) or {}
    ref_paths = [
        refs.get("global_selector_atlas"),
        refs.get("global_selector_transition"),
        refs.get("global_projective_selector_state"),
        refs.get("global_selector_output_operator"),
    ]
    missing_refs = [p for p in ref_paths if not isinstance(p, str) or not (REPO / p).exists()]
    checks.append(
        {
            "id": "referenced_inputs_exist",
            "actual": missing_refs,
            "expected": [],
            "pass": len(missing_refs) == 0,
            "meaning": "All referenced global inputs exist on disk (atlas/transition/state/output operator).",
        }
    )
    if missing_refs:
        blocking_mismatches.append("referenced_inputs_exist")

    # Chart payload sanity: ensure each chart has a 2x2 output projector matrix.
    charts = obj.get("charts", {}) or {}
    expected_charts = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    chart_ok = sorted(list(charts.keys())) == sorted(expected_charts)
    checks.append(
        {
            "id": "charts_cover_pair1_to_pair5",
            "actual": sorted(list(charts.keys())),
            "expected": expected_charts,
            "pass": chart_ok,
            "meaning": "Closure object carries chart payload for the declared global cover {pair1..pair5}.",
        }
    )
    if not chart_ok:
        blocking_mismatches.append("charts_cover_pair1_to_pair5")

    missing_chart_matrices = []
    for chart in expected_charts:
        entry = charts.get(chart) or {}
        b = entry.get("B_out_matrix_in_o_plus_o_minus")
        if not is_2x2_matrix(b):
            missing_chart_matrices.append(chart)

        a_ref = entry.get("A_m_ref")
        if not isinstance(a_ref, str) or not (REPO / a_ref).exists():
            missing_chart_matrices.append(f"{chart}:missing_A_m_ref")

    checks.append(
        {
            "id": "chartwise_output_projectors_present",
            "actual": missing_chart_matrices,
            "expected": [],
            "pass": len(missing_chart_matrices) == 0,
            "meaning": "Each chart has a 2x2 output projector matrix and a valid A_m reference.",
        }
    )
    if missing_chart_matrices:
        blocking_mismatches.append("chartwise_output_projectors_present")

    if blocking_mismatches:
        status = "P673_REQUIRES_REVIEW_INSUFFICIENT_OR_INVALID_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OBJECT"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OBJECT_ON_C_V1_AFTER_F672"

    artifact = {
        "stage": "P673",
        "lane": "current_strict_global_projective_selector_closure_object_gluing_probe_only",
        "status": status,
        "closure_object_file": str(IN_OBJECT.relative_to(REPO)),
        "certificate": {
            "tolerance": tol,
            "max_abs_diff_to_reference_chart_output_projector": maxdiff,
            "certificate_pass": bool(cert_pass),
        },
        "checks": checks,
        "blocking_mismatches": blocking_mismatches,
        "no_false_pass": True,
    }
    summary = {
        "stage": "P673",
        "lane": artifact["lane"],
        "status": status,
        "blocking_mismatches": blocking_mismatches,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

