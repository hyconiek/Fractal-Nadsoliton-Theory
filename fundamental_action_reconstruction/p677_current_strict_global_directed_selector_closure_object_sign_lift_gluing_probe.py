#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_OBJECT = GENERATED / "selector_closure_global_c_v1_directed_strict_v1.json"

OUT_JSON = (
    GENERATED / "p677_current_strict_global_directed_selector_closure_object_sign_lift_gluing_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p677_current_strict_global_directed_selector_closure_object_sign_lift_gluing_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_len2_vec(v: Any) -> bool:
    return isinstance(v, list) and len(v) == 2 and all(isinstance(x, (int, float)) for x in v)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    checks: list[dict[str, Any]] = []
    blocking_mismatches: list[str] = []

    if not IN_OBJECT.exists():
        status = "P677_NOT_COMPUTABLE_MISSING_DIRECTED_SELECTOR_CLOSURE_OBJECT_FILE"
        artifact: dict[str, Any] = {
            "stage": "P677",
            "lane": "current_strict_global_directed_selector_closure_object_sign_lift_gluing_probe_only",
            "status": status,
            "missing": str(IN_OBJECT.relative_to(REPO)),
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P677",
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
        "SelectorClosure_global_C_v1_directed_strict_v1",
        "Closure object has the expected typed name.",
    )
    check_equal(
        "no_false_pass_marked",
        bool(obj.get("no_false_pass")),
        True,
        "Closure object is explicitly marked no_false_pass.",
    )
    check_equal(
        "closure_scope_directed",
        (obj.get("closure_scope", {}) or {}).get("level"),
        "directed_vector_state",
        "Closure object is declared in directed/vector-level scope.",
    )

    cert = obj.get("well_definedness_certificate", {}) or {}
    tol = cert.get("tolerance")
    maxdiff = cert.get("max_abs_diff_to_reference_chart_v_out")
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
        "Directed closure well-definedness certificate passes on exported data (after explicit sign-lift).",
    )

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

    # Sign-lift payload must be explicit and cover pair1..pair5 with ±1 values.
    sign_lift = obj.get("sign_lift", {}) or {}
    signs = sign_lift.get("signs_by_pair")
    expected_pairs = [f"pair{i}" for i in range(1, 6)]
    sign_ok = (
        isinstance(signs, dict)
        and sorted(list(signs.keys())) == sorted(expected_pairs)
        and all(signs.get(p) in (-1, 1) for p in expected_pairs)
    )
    checks.append(
        {
            "id": "explicit_sign_lift_present_and_total",
            "actual": signs,
            "expected": "{pair1..pair5} -> {±1}",
            "pass": bool(sign_ok),
            "meaning": "Closure object exports an explicit per-chart sign-lift covering {pair1..pair5}.",
        }
    )
    if not sign_ok:
        blocking_mismatches.append("explicit_sign_lift_present_and_total")

    # Referenced inputs existence check.
    refs = obj.get("inputs", {}) or {}
    ref_paths = [
        refs.get("global_selector_atlas"),
        refs.get("global_selector_transition"),
        refs.get("global_directed_selector_state"),
        refs.get("global_selector_output_operator"),
    ]
    missing_refs = [p for p in ref_paths if not isinstance(p, str) or not (REPO / p).exists()]
    checks.append(
        {
            "id": "referenced_inputs_exist",
            "actual": missing_refs,
            "expected": [],
            "pass": len(missing_refs) == 0,
            "meaning": "All referenced global inputs exist on disk (atlas/transition/directed state/output operator).",
        }
    )
    if missing_refs:
        blocking_mismatches.append("referenced_inputs_exist")

    # Chart payload sanity: ensure each chart has raw and corrected output vectors.
    charts = obj.get("charts", {}) or {}
    chart_ok = sorted(list(charts.keys())) == sorted(expected_pairs)
    checks.append(
        {
            "id": "charts_cover_pair1_to_pair5",
            "actual": sorted(list(charts.keys())),
            "expected": expected_pairs,
            "pass": chart_ok,
            "meaning": "Closure object carries chart payload for the declared global cover {pair1..pair5}.",
        }
    )
    if not chart_ok:
        blocking_mismatches.append("charts_cover_pair1_to_pair5")

    missing_vectors: list[str] = []
    for pair in expected_pairs:
        entry = charts.get(pair) or {}
        if not is_len2_vec(entry.get("v_out_raw_in_o_plus_o_minus")):
            missing_vectors.append(f"{pair}:missing_raw")
        if not is_len2_vec(entry.get("v_out_in_o_plus_o_minus")):
            missing_vectors.append(f"{pair}:missing_corrected")
        if entry.get("sign_lift_applied") not in (-1, 1):
            missing_vectors.append(f"{pair}:missing_sign_lift_applied")
    checks.append(
        {
            "id": "chartwise_raw_and_corrected_vectors_present",
            "actual": missing_vectors,
            "expected": [],
            "pass": len(missing_vectors) == 0,
            "meaning": "Each chart carries raw and sign-lifted directed output vectors plus the applied sign.",
        }
    )
    if missing_vectors:
        blocking_mismatches.append("chartwise_raw_and_corrected_vectors_present")

    if blocking_mismatches:
        status = "P677_REQUIRES_REVIEW_INSUFFICIENT_OR_INVALID_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OBJECT"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OBJECT_ON_C_V1_WITH_EXPLICIT_SIGN_LIFT_AFTER_F677"

    artifact = {
        "stage": "P677",
        "lane": "current_strict_global_directed_selector_closure_object_sign_lift_gluing_probe_only",
        "status": status,
        "closure_object_file": str(IN_OBJECT.relative_to(REPO)),
        "certificate": {
            "tolerance": tol,
            "max_abs_diff_to_reference_chart_v_out": maxdiff,
            "certificate_pass": bool(cert_pass),
        },
        "checks": checks,
        "blocking_mismatches": blocking_mismatches,
        "no_false_pass": True,
    }
    summary = {
        "stage": "P677",
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

