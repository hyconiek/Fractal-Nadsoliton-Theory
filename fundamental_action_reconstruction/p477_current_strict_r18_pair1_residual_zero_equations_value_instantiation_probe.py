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

IN_R18 = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
)
IN_P459_OBJECT = GENERATED / "m_control_residual_diag_declared_value_instantiated_v1.json"

OUT_OBJECT = GENERATED / "pair1_residual_zero_equations_evaluation_under_n477_value_instance_v1.json"
OUT_JSON = GENERATED / "p477_current_strict_r18_pair1_residual_zero_equations_value_instantiation_probe.json"
OUT_SUMMARY = GENERATED / "p477_current_strict_r18_pair1_residual_zero_equations_value_instantiation_probe_summary.json"


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {"R18_reduction": IN_R18, "P459_value_instance": IN_P459_OBJECT}
    missing = [k for k, p in required.items() if not p.is_file()]
    if missing:
        payload = {
            "stage": "P477",
            "date_utc": datetime.now(timezone.utc).date().isoformat(),
            "goal": "evaluate_R18_declared_pair1_residual_zero_equations_on_the_P459_value_instantiated_residual_pullback",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {"stage": "P477", "status": payload["status"], "missing_required_inputs": missing, "no_false_pass": True},
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    r18 = _read_json(IN_R18)
    p459 = _read_json(IN_P459_OBJECT)

    reduction = r18.get("pair1_coefficient_class_reduction")
    if not isinstance(reduction, dict):
        raise ValueError("expected R18.pair1_coefficient_class_reduction dict")

    classes = reduction.get("coefficient_classes")
    if not (isinstance(classes, list) and classes):
        raise ValueError("expected R18 coefficient_classes to be a nonempty list")

    d = ((p459.get("computed") or {}).get("d_local_residual_profile_n477")) if isinstance(p459.get("computed"), dict) else None
    if not (isinstance(d, list) and len(d) == 12 and all(_is_number(x) for x in d)):
        raise ValueError("expected P459 computed d_local_residual_profile_n477 as a length-12 finite numeric list")
    d_by_slot = {f"psi{i}": float(d[i]) for i in range(12)}

    # Compute each coefficient-class Sigma sum and use it to evaluate the three independent pair1 zero equations.
    sigma_by_class: dict[str, float] = {}
    class_rows: list[dict[str, Any]] = []
    for row in classes:
        if not isinstance(row, dict):
            raise ValueError("expected each R18 coefficient class row to be a dict")
        class_symbol = row.get("class_symbol")
        carrier_slots = row.get("carrier_slots")
        signature = row.get("signature_on_pair1_entries")
        if not (isinstance(class_symbol, str) and class_symbol):
            raise ValueError("expected class_symbol string in R18 class row")
        if not (isinstance(carrier_slots, list) and all(isinstance(x, str) for x in carrier_slots)):
            raise ValueError("expected carrier_slots list[str] in R18 class row")
        if not (isinstance(signature, dict) and set(signature.keys()) == {"c1c1", "c1s1", "s1s1"}):
            raise ValueError("expected signature_on_pair1_entries dict with keys c1c1/c1s1/s1s1")
        if not all(_is_number(signature[k]) for k in ("c1c1", "c1s1", "s1s1")):
            raise ValueError("expected numeric signature_on_pair1_entries coefficients")

        sigma = 0.0
        missing_slots: list[str] = []
        for slot in carrier_slots:
            if slot not in d_by_slot:
                missing_slots.append(slot)
            else:
                sigma += float(d_by_slot[slot])
        if missing_slots:
            raise ValueError(f"unknown carrier slots in R18 class row: {missing_slots}")

        sigma_by_class[class_symbol] = float(sigma)
        class_rows.append(
            {
                "class_symbol": class_symbol,
                "carrier_slots": carrier_slots,
                "sigma_sum_value": float(sigma),
                "coeffs": {k: float(signature[k]) for k in ("c1c1", "c1s1", "s1s1")},
            }
        )

    def eval_entry(entry: str) -> float:
        acc = 0.0
        for row in classes:
            class_symbol = str(row["class_symbol"])
            coeff = float(row["signature_on_pair1_entries"][entry])
            acc += coeff * float(sigma_by_class[class_symbol])
        return float(acc)

    # Independent R18 equations.
    c1c1 = eval_entry("c1c1")
    c1s1 = eval_entry("c1s1")
    s1s1 = eval_entry("s1s1")

    # Cross-check against P459 pair1 block (direct contraction output).
    pair1_block = ((p459.get("computed") or {}).get("pair1_block")) if isinstance(p459.get("computed"), dict) else None
    if not (
        isinstance(pair1_block, list)
        and len(pair1_block) == 2
        and all(isinstance(r, list) and len(r) == 2 and all(_is_number(x) for x in r) for r in pair1_block)
    ):
        raise ValueError("expected P459 computed pair1_block to be a 2x2 finite numeric matrix")
    p459_c1c1 = float(pair1_block[0][0])
    p459_c1s1 = float(pair1_block[0][1])
    p459_s1s1 = float(pair1_block[1][1])

    # R18 coefficients are rounded (class-level), while P459 computes entries from the full R16 tensor contraction.
    # Use a looser cross-check tolerance to avoid false negatives from harmless rounding drift.
    cross_tol = 1e-9
    cross_check = {
        "tolerance": cross_tol,
        "c1c1_abs_diff": float(abs(c1c1 - p459_c1c1)),
        "c1s1_abs_diff": float(abs(c1s1 - p459_c1s1)),
        "s1s1_abs_diff": float(abs(s1s1 - p459_s1s1)),
        "pass": bool(
            abs(c1c1 - p459_c1c1) <= cross_tol
            and abs(c1s1 - p459_c1s1) <= cross_tol
            and abs(s1s1 - p459_s1s1) <= cross_tol
        ),
    }
    if not cross_check["pass"]:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P477",
                    "status": "FAIL_R18_EQUATION_EVAL_MISMATCHES_P459_PAIR1_BLOCK",
                    "cross_check": cross_check,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    zero_tol = 1e-12
    eqs = [
        {
            "entry": "c1c1",
            "value": float(c1c1),
            "abs": float(abs(c1c1)),
            "zero_tol": zero_tol,
            "is_zero": abs(c1c1) <= zero_tol,
        },
        {
            "entry": "c1s1",
            "value": float(c1s1),
            "abs": float(abs(c1s1)),
            "zero_tol": zero_tol,
            "is_zero": abs(c1s1) <= zero_tol,
        },
        {
            "entry": "s1s1",
            "value": float(s1s1),
            "abs": float(abs(s1s1)),
            "zero_tol": zero_tol,
            "is_zero": abs(s1s1) <= zero_tol,
        },
    ]

    all_zero = bool(all(e["is_zero"] for e in eqs))
    violated = [e["entry"] for e in eqs if not e["is_zero"]]

    status = "PASS_EVALUATION_ZERO_EQUATIONS_SATISFIED_UNDER_CURRENT_VALUE_INSTANCE"
    if not all_zero:
        status = "PASS_EVALUATION_ZERO_EQUATIONS_VIOLATED_UNDER_CURRENT_VALUE_INSTANCE"

    obj = {
        "object": "Pair1ResidualZeroEquationsEvaluation_under_N477_value_instance_v1",
        "status": status,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Evaluate the three independent declared pair1 residual zero equations exported by R18 on a concrete value-instantiated "
            "declared residual local-diagonal control pullback. This is an evaluation/obstruction artifact only: it does not export "
            "a strict zero witness, does not export a stationarity witness, and does not claim host matching or QW-2191 discharge."
        ),
        "inputs": {
            "r18_reduction_packet": _p(IN_R18),
            "p459_value_instance_object": _p(IN_P459_OBJECT),
        },
        "value_instance_scope": {
            "note": "P459 instantiates the residual diagonal profile using the conditional N477 Yukawa-free rewrite and a strict-derived (vpsi,g4,g6) provider; no stationarity witness is exported.",
        },
        "computed": {
            "sigma_by_class": sigma_by_class,
            "coefficient_classes": class_rows,
            "equations": eqs,
            "all_zero_equations_satisfied": bool(all_zero),
            "violated_equations": violated,
            "cross_check_against_p459_pair1_block": cross_check,
        },
        "hard_limits": [
            "Consumes P459 which uses the conditional N477 Yukawa-free diagonal residual rewrite; no stationarity witness is exported.",
            "Does not claim the declared pair1 residual equations hold in strict core or for the canonical vacuum.",
            "Does not claim any strict cancellation/host matching witness is obtained.",
            "Does not claim selector-relevant physical canonicalization within the QW-2191 O(2) family.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    payload = {
        "stage": "P477",
        "date_utc": datetime.now(timezone.utc).date().isoformat(),
        "goal": "evaluate_R18_declared_pair1_residual_zero_equations_on_the_P459_value_instantiated_residual_pullback",
        "status": status,
        "exported_object": _p(OUT_OBJECT),
        "all_zero_equations_satisfied": bool(all_zero),
        "violated_equations": violated,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P477",
        "status": status,
        "exported_object": payload["exported_object"],
        "all_zero_equations_satisfied": payload["all_zero_equations_satisfied"],
        "violated_equations": violated,
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
