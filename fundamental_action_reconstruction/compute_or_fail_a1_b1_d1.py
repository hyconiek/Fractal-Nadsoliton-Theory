#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

FAR = Path(__file__).resolve().parent
REPO = FAR.parent
GENERATED = FAR / "generated"
OUT_JSON = GENERATED / "compute_or_fail_a1_b1_d1.json"
OUT_SUMMARY = GENERATED / "compute_or_fail_a1_b1_d1_summary.json"

O2_PATH = GENERATED / "o2_exported_composite_a1_ext_instance.json"
O3_PATH = GENERATED / "o3_a1_ext_coefficient_evaluation_rule.json"
Q2165_PATH = REPO / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json" / "report_qw2165_l13_exhaustive_canonical_eom_gate.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_number(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        try:
            return float(text)
        except ValueError:
            return None
    return None


def find_numeric_hits(obj: Any, path: str = "$") -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            child = f"{path}.{key}"
            if key in {"a_1", "b_1", "d_1", "trace_A_1"}:
                num = parse_number(value)
                if num is not None:
                    hits.append({"path": child, "value": num})
            hits.extend(find_numeric_hits(value, child))
    elif isinstance(obj, list):
        for idx, value in enumerate(obj):
            hits.extend(find_numeric_hits(value, f"{path}[{idx}]"))
    return hits


def scan_generated_for_numeric_coeffs() -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    for path in sorted(GENERATED.glob("*.json")):
        try:
            obj = load_json(path)
        except Exception:
            continue
        for hit in find_numeric_hits(obj):
            hit["file"] = str(path.relative_to(REPO))
            hits.append(hit)
    return hits


def main() -> None:
    o2 = load_json(O2_PATH)
    o3 = load_json(O3_PATH)
    q2165 = load_json(Q2165_PATH)

    matrix = o2.get("matrix_form", [[None, None], [None, None]])
    a_raw = matrix[0][0]
    b_raw = matrix[0][1]
    d_raw = matrix[1][1]

    a_num = parse_number(a_raw)
    b_num = parse_number(b_raw)
    d_num = parse_number(d_raw)

    numeric_hits = scan_generated_for_numeric_coeffs()

    q2165_model = q2165.get("model", {})
    q2165_keys = set(q2165_model.keys())
    q2165_has_a1_export = any(key in {"A_1", "A_1_ext", "selector_block", "pair1_block"} for key in q2165_keys)
    lagrangian_density = str(q2165_model.get("lagrangian_density", ""))
    q2165_symbolic_kernel = "K_0_1" in lagrangian_density or "K_1_0" in lagrangian_density

    computable_now = all(x is not None for x in (a_num, b_num, d_num))

    if computable_now:
        trace_a1 = a_num + d_num
        delta_1 = [a_num - d_num, b_num]
        status = "COMPUTED"
        missing_inputs: list[str] = []
        computed = {
            "a_1": a_num,
            "b_1": b_num,
            "d_1": d_num,
            "trace_A_1": trace_a1,
            "Delta_1": delta_1,
        }
        reason = "persisted_A1_ext_entries_are_numeric"
    else:
        status = "NOT_COMPUTABLE_FROM_CURRENT_REPO_STATE"
        computed = {}
        missing_inputs = [
            "numeric_or_evaluable_a_1_entry_in_A_1_ext(pair1)",
            "numeric_or_evaluable_b_1_entry_in_A_1_ext(pair1)",
            "numeric_or_evaluable_d_1_entry_in_A_1_ext(pair1)",
            "actual_Route_P1_or_Route_P2_populated_witness_for_A_1_ext_entries",
            "operator_level_selector_block_export_or_equivalent_populated_extension_instance",
        ]
        reason = "persisted_A1_ext_entries_remain_symbolic_and_QW2165_does_not_export_selector_block"

    out = {
        "stage": "CF1",
        "goal": "compute_or_fail_a1_b1_d1",
        "status": status,
        "reason": reason,
        "inputs": {
            "o2_path": str(O2_PATH.relative_to(REPO)),
            "o3_path": str(O3_PATH.relative_to(REPO)),
            "q2165_path": str(Q2165_PATH.relative_to(REPO)),
        },
        "inspection": {
            "o2_matrix_form": matrix,
            "o3_current_entry_state": o3.get("current_entry_state", {}),
            "q2165_has_a1_or_selector_block_export": q2165_has_a1_export,
            "q2165_symbolic_kernel_mixing_present": q2165_symbolic_kernel,
            "numeric_coeff_hits_in_generated": numeric_hits,
        },
        "computed": computed,
        "missing_inputs": missing_inputs,
        "required_next_step": (
            "PROVIDE_POPULATED_ROUTE_P1_OR_P2_WITNESS_FOR_A1_EXT_ENTRIES"
            if status == "NOT_COMPUTABLE_FROM_CURRENT_REPO_STATE"
            else "RUN_SELECTOR_BREAKING_TEST_FROM_COMPUTED_COEFFICIENTS"
        ),
        "lane": "hypothesis_extension_only",
        "no_false_pass": True,
    }

    summary = {
        "stage": out["stage"],
        "status": out["status"],
        "reason": out["reason"],
        "computed": out["computed"],
        "required_next_step": out["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    OUT_SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
