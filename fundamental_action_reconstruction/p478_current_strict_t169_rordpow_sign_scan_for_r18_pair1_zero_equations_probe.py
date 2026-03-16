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

IN_R14 = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
IN_R15 = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
IN_R18 = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
)
IN_F447_PROVIDER = GENERATED / "f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json"

OUT_OBJECT = (
    GENERATED
    / "pair1_residual_zero_equations_sign_scan_under_rordpow_magnitudes_value_instance_v1.json"
)
OUT_JSON = GENERATED / "p478_current_strict_t169_rordpow_sign_scan_for_r18_pair1_zero_equations_probe.json"
OUT_SUMMARY = (
    GENERATED / "p478_current_strict_t169_rordpow_sign_scan_for_r18_pair1_zero_equations_probe_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def _extract_ktotal(r14: dict[str, Any]) -> list[list[float]]:
    rows = r14.get("host_kernel_rows")
    if not (
        isinstance(rows, list)
        and len(rows) == 12
        and all(isinstance(r, list) and len(r) == 12 and all(_is_number(x) for x in r) for r in rows)
    ):
        raise ValueError("expected R14.host_kernel_rows as a 12x12 finite numeric matrix")
    return [[float(x) for x in row] for row in rows]


def _extract_m0_sq(r15: dict[str, Any]) -> float:
    dd = r15.get("diagonal_decomposition")
    if not isinstance(dd, dict):
        raise ValueError("expected R15.diagonal_decomposition dict")
    m0_sq = dd.get("host_diagonal_floor_value")
    if not _is_number(m0_sq):
        raise ValueError("expected R15.diagonal_decomposition.host_diagonal_floor_value finite numeric")
    return float(m0_sq)


def _extract_rordpow_reference_prob_and_strict_values(f447: dict[str, Any]) -> tuple[list[float], float, float]:
    values = f447.get("values")
    if not isinstance(values, dict):
        raise ValueError("expected F447.values dict")
    rho_sq = values.get("rho_star_sq")
    lambda_psi = values.get("lambda_psi_strict")
    if not _is_number(rho_sq):
        raise ValueError("expected F447.values.rho_star_sq finite numeric")
    if not _is_number(lambda_psi):
        raise ValueError("expected F447.values.lambda_psi_strict finite numeric")

    ref = f447.get("reference")
    if not (isinstance(ref, dict) and isinstance(ref.get("r_ordpow"), str) and ref.get("r_ordpow")):
        raise ValueError("expected F447.reference.r_ordpow path string")
    ref_path = Path(str(ref["r_ordpow"]))
    if not ref_path.is_absolute():
        ref_path = REPO / ref_path
    ref_obj = _read_json(ref_path)
    q = ref_obj.get("reference_prob")
    if not (isinstance(q, list) and len(q) == 12 and all(_is_number(x) and float(x) > 0.0 for x in q)):
        raise ValueError("expected r_ordpow_z12_v1_reference_distribution.reference_prob as length-12 positive numeric list")
    qf = [float(x) for x in q]
    s = float(sum(qf))
    if not math.isfinite(s) or abs(s - 1.0) > 1e-9:
        raise ValueError("expected r_ordpow reference_prob to be normalized (sum=1)")
    return qf, float(rho_sq), float(lambda_psi)


def _extract_r18_coefficient_classes(r18: dict[str, Any]) -> list[dict[str, Any]]:
    reduction = r18.get("pair1_coefficient_class_reduction")
    if not isinstance(reduction, dict):
        raise ValueError("expected R18.pair1_coefficient_class_reduction dict")
    classes = reduction.get("coefficient_classes")
    if not (isinstance(classes, list) and len(classes) == 6):
        raise ValueError("expected R18 coefficient_classes as a length-6 list")
    out: list[dict[str, Any]] = []
    for row in classes:
        if not isinstance(row, dict):
            raise ValueError("expected each coefficient class row to be a dict")
        class_symbol = row.get("class_symbol")
        carrier_slots = row.get("carrier_slots")
        signature = row.get("signature_on_pair1_entries")
        if not (isinstance(class_symbol, str) and class_symbol):
            raise ValueError("expected class_symbol string")
        if not (isinstance(carrier_slots, list) and all(isinstance(x, str) for x in carrier_slots) and len(carrier_slots) == 2):
            raise ValueError("expected carrier_slots as length-2 list[str]")
        if not (isinstance(signature, dict) and set(signature.keys()) == {"c1c1", "c1s1", "s1s1"}):
            raise ValueError("expected signature_on_pair1_entries dict with keys c1c1/c1s1/s1s1")
        if not all(_is_number(signature[k]) for k in ("c1c1", "c1s1", "s1s1")):
            raise ValueError("expected numeric signature coefficients")
        out.append(
            {
                "class_symbol": class_symbol,
                "carrier_slots": list(carrier_slots),
                "coeffs": {k: float(signature[k]) for k in ("c1c1", "c1s1", "s1s1")},
            }
        )
    return out


def _compute_d_profile_n477(*, ktotal: list[list[float]], m0_sq: float, vpsi: list[float], g4: float) -> list[float]:
    if len(vpsi) != 12:
        raise ValueError("vpsi must be length 12")
    d: list[float] = []
    for k in range(12):
        denom = float(vpsi[k])
        if denom == 0.0:
            raise ValueError("vpsi entries must be nonzero (division premise)")
        mix = 0.0
        for j in range(12):
            if j == k:
                continue
            mix += float(ktotal[k][j]) * (float(vpsi[j]) / denom)
        d_k = -mix + 2.0 * float(g4) * (denom**2) - float(m0_sq)
        d.append(float(d_k))
    return d


def _eval_r18_pair1_entries(*, classes: list[dict[str, Any]], d: list[float]) -> dict[str, float]:
    if len(d) != 12:
        raise ValueError("d must be length 12")
    d_by_slot = {f"psi{i}": float(d[i]) for i in range(12)}
    sigma_by_class: dict[str, float] = {}
    for row in classes:
        sigma = 0.0
        for slot in row["carrier_slots"]:
            if slot not in d_by_slot:
                raise ValueError(f"unexpected slot: {slot}")
            sigma += float(d_by_slot[slot])
        sigma_by_class[str(row["class_symbol"])] = float(sigma)

    out: dict[str, float] = {}
    for entry in ("c1c1", "c1s1", "s1s1"):
        acc = 0.0
        for row in classes:
            coeff = float(row["coeffs"][entry])
            acc += coeff * float(sigma_by_class[str(row["class_symbol"])])
        out[entry] = float(acc)
    return out


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {"R14": IN_R14, "R15": IN_R15, "R18": IN_R18, "F447_provider": IN_F447_PROVIDER}
    missing = [k for k, p in required.items() if not p.is_file()]
    if missing:
        payload = {
            "stage": "P478",
            "date_utc": datetime.now(timezone.utc).date().isoformat(),
            "goal": "scan_all_sign_vectors_under_the_T169_r_ordpow_magnitude_lift_for_a_pair1_R18_declared_zero_equations_solution_under_conditional_N477",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {"stage": "P478", "status": payload["status"], "missing_required_inputs": missing, "no_false_pass": True},
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    r14 = _read_json(IN_R14)
    r15 = _read_json(IN_R15)
    r18 = _read_json(IN_R18)
    f447 = _read_json(IN_F447_PROVIDER)

    ktotal = _extract_ktotal(r14)
    m0_sq = _extract_m0_sq(r15)
    q, rho_sq, lambda_psi = _extract_rordpow_reference_prob_and_strict_values(f447)

    # T169 uniform g4 lift for the r_ordpow magnitude distribution.
    sum_q2 = float(sum(float(x) ** 2 for x in q))
    if not math.isfinite(sum_q2) or sum_q2 <= 0.0:
        raise ValueError("unexpected: sum_q2 must be finite positive")
    g4 = float(lambda_psi / sum_q2)

    # Magnitudes |vpsi_i| := sqrt(rho_*^2 * q_i).
    v_abs = [float(math.sqrt(float(rho_sq) * float(x))) for x in q]
    if any((not math.isfinite(x) or x <= 0.0) for x in v_abs):
        raise ValueError("unexpected: all magnitudes must be finite positive")

    classes = _extract_r18_coefficient_classes(r18)

    zero_tol = 1e-12
    objective = "max_abs_entry"

    scanned = 0
    any_solution = False
    solutions: list[dict[str, Any]] = []

    best: dict[str, Any] | None = None
    min_abs_by_entry = {"c1c1": float("inf"), "c1s1": float("inf"), "s1s1": float("inf")}

    # Global sign does not affect N477 ratios (s_j/s_k invariant); fix s0=+1 to halve the brute-force space.
    for mask in range(1 << 11):
        # indices 1..11 are controlled by mask bits 0..10
        s = [1]
        for i in range(1, 12):
            bit = (mask >> (i - 1)) & 1
            s.append(1 if bit == 0 else -1)

        vpsi = [float(s[i]) * float(v_abs[i]) for i in range(12)]
        d = _compute_d_profile_n477(ktotal=ktotal, m0_sq=m0_sq, vpsi=vpsi, g4=g4)
        entries = _eval_r18_pair1_entries(classes=classes, d=d)

        scanned += 1
        eqs = {
            k: {
                "value": float(entries[k]),
                "abs": float(abs(entries[k])),
                "zero_tol": float(zero_tol),
                "is_zero": bool(abs(entries[k]) <= zero_tol),
            }
            for k in ("c1c1", "c1s1", "s1s1")
        }
        for k in ("c1c1", "c1s1", "s1s1"):
            min_abs_by_entry[k] = float(min(min_abs_by_entry[k], float(eqs[k]["abs"])))

        max_abs = float(max(eqs[k]["abs"] for k in ("c1c1", "c1s1", "s1s1")))

        if best is None or max_abs < float(best["objective_value"]):
            best = {
                "objective": objective,
                "objective_value": float(max_abs),
                "sign_vector": s,
                "equations": eqs,
            }

        if all(bool(eqs[k]["is_zero"]) for k in ("c1c1", "c1s1", "s1s1")):
            any_solution = True
            if len(solutions) < 10:
                solutions.append({"sign_vector": s, "equations": eqs})

    status = "PASS_SCAN_COMPLETE_NO_SIGN_VECTOR_SATISFIES_R18_PAIR1_ZERO_EQUATIONS_UNDER_N477_ON_RORDPOW_MAGNITUDES"
    if any_solution:
        status = "PASS_SCAN_COMPLETE_AT_LEAST_ONE_SIGN_VECTOR_SATISFIES_R18_PAIR1_ZERO_EQUATIONS_UNDER_N477_ON_RORDPOW_MAGNITUDES"

    obj = {
        "object": "Pair1ResidualZeroEquationsSignScan_under_r_ordpow_magnitudes_value_instance_v1",
        "status": status,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Brute-force scan all 2^11 sign vectors (fixing s0=+1) under the fixed T169 r_ordpow magnitude lift "
            "|vpsi_i|=sqrt(rho_*^2*q_i) and the uniform g4 lift, compute the conditional N477 Yukawa-free residual profile d_k, "
            "and evaluate the R18 declared pair1 residual zero equations (c1c1,c1s1,s1s1). This is probe-only evidence: it does "
            "not export any strict stationarity witness and does not claim host matching or QW-2191 discharge."
        ),
        "inputs": {
            "r14_k_total": _p(IN_R14),
            "r15_m0_squared": _p(IN_R15),
            "r18_reduction_packet": _p(IN_R18),
            "t169_magnitudes_source_provider": _p(IN_F447_PROVIDER),
        },
        "fixed_lift_parameters": {
            "rho_star_sq": float(rho_sq),
            "lambda_psi_strict": float(lambda_psi),
            "r_ordpow_reference_prob": q,
            "g4_uniform_value": float(g4),
            "g6_assumed_zero": True,
            "global_sign_fixed": "s0=+1 (global Z2 is ratio-invariant for N477 anyway)",
        },
        "scan": {
            "scan_space_size": 2048,
            "sign_vector_masked_indices": "i=1..11",
            "zero_tolerance": float(zero_tol),
            "objective": objective,
            "any_solution_within_tolerance": bool(any_solution),
            "solutions_within_tolerance_capped": solutions,
            "best_candidate": best,
            "min_abs_by_entry_over_scan": min_abs_by_entry,
        },
        "hard_limits": [
            "Consumes the conditional N477 Yukawa-free diagonal residual rewrite; no vacuum stationarity witness is exported.",
            "Scans only sign vectors under the fixed r_ordpow magnitude lift (T169 power-law element-order constrained magnitudes) and the fixed uniform g4 lift.",
            "Does not export a strict zero witness unless a theorem-level witness is separately constructed from the scan result.",
            "Does not claim host matching, strict-core selector closure, QW-2191 discharge, or ToE closure.",
        ],
        "no_false_pass": True,
    }

    payload = {
        "stage": "P478",
        "date_utc": datetime.now(timezone.utc).date().isoformat(),
        "goal": "scan_all_sign_vectors_under_the_T169_r_ordpow_magnitude_lift_for_a_pair1_R18_declared_zero_equations_solution_under_conditional_N477",
        "status": status,
        "exported_object": _p(OUT_OBJECT),
        "scan_space_size": 2048,
        "zero_tolerance": float(zero_tol),
        "any_solution_within_tolerance": bool(any_solution),
        "best_candidate_objective": objective,
        "best_candidate_objective_value": None if best is None else float(best["objective_value"]),
        "best_candidate_equations": None if best is None else best["equations"],
        "min_abs_by_entry_over_scan": min_abs_by_entry,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P478",
        "status": status,
        "exported_object": payload["exported_object"],
        "any_solution_within_tolerance": payload["any_solution_within_tolerance"],
        "best_candidate_objective_value": payload["best_candidate_objective_value"],
        "min_abs_by_entry_over_scan": payload["min_abs_by_entry_over_scan"],
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

