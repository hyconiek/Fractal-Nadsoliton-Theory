#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

QW2122_JSON = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"
ALPHA_GEO_JSON = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_R14 = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
IN_R15 = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
IN_R18 = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
)

IN_RORDPOW = GENERATED / "r_ordpow_z12_v1_reference_distribution.json"

OUT_OBJECT = (
    GENERATED
    / "pair1_residual_zero_equations_reference_magnitude_family_scan_under_n477_v1.json"
)
OUT_JSON = GENERATED / "p479_current_strict_reference_magnitude_family_sign_scan_for_r18_pair1_zero_equations_probe.json"
OUT_SUMMARY = (
    GENERATED / "p479_current_strict_reference_magnitude_family_sign_scan_for_r18_pair1_zero_equations_probe_summary.json"
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


def _parse_alpha_geo_numeric(obj: dict[str, Any]) -> float:
    val = obj.get("value")
    if _is_number(val):
        return float(val)
    if isinstance(val, str):
        compact = val.strip().replace(" ", "")
        if compact in ("4ln(2)", "4*ln(2)", "4*log(2)"):
            return float(4.0 * math.log(2.0))
    # fallback: if computation says ln(16)
    definition = obj.get("definition")
    if isinstance(definition, dict) and "ln(16)" in str(definition.get("computation", "")):
        return float(4.0 * math.log(2.0))
    raise ValueError("could not parse alpha_geo_strict_derived_v1 numeric value")


def _extract_qw2122_strict_scalars(qw2122: dict[str, Any]) -> tuple[float, float]:
    inputs = qw2122.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("expected QW-2122 inputs dict")
    lam = inputs.get("lambda_psi_strict")
    if not _is_number(lam) or float(lam) <= 0.0:
        raise ValueError("expected QW-2122.inputs.lambda_psi_strict positive numeric")
    br = (qw2122.get("branch_results") or {}).get("broken_branch_strict") if isinstance(qw2122.get("branch_results"), dict) else None
    if not isinstance(br, dict):
        raise ValueError("expected QW-2122.branch_results.broken_branch_strict dict")
    rho_sq = br.get("rho_star_sq")
    if not _is_number(rho_sq) or float(rho_sq) <= 0.0:
        raise ValueError("expected QW-2122.branch_results.broken_branch_strict.rho_star_sq positive numeric")
    return float(rho_sq), float(lam)


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


def _extract_r18_classes(r18: dict[str, Any]) -> list[dict[str, Any]]:
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
            raise ValueError("expected signature_on_pair1_entries keys c1c1/c1s1/s1s1")
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


def _normalize_positive_weights(w: list[float]) -> list[float]:
    if not w or any((not math.isfinite(float(x)) or float(x) <= 0.0) for x in w):
        raise ValueError("weights must be finite positive")
    z = float(sum(float(x) for x in w))
    if not math.isfinite(z) or z <= 0.0:
        raise ValueError("weight sum must be finite positive")
    q = [float(x) / z for x in w]
    if abs(sum(q) - 1.0) > 1e-9:
        raise ValueError("unexpected: normalization failed (sum != 1)")
    return q


def _ord_z12_by_x() -> list[int]:
    out: list[int] = []
    for x in range(12):
        if x == 0:
            out.append(1)
        else:
            out.append(int(12 // math.gcd(x, 12)))
    return out


def _dcyc_0(x: int) -> int:
    a = int(x) % 12
    return int(min(a, 12 - a))


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


def _scan_sign_space_for_reference(
    *,
    reference_id: str,
    q: list[float],
    rho_sq: float,
    lambda_psi: float,
    ktotal: list[list[float]],
    m0_sq: float,
    classes: list[dict[str, Any]],
    zero_tol: float,
) -> dict[str, Any]:
    if not (isinstance(q, list) and len(q) == 12 and all(_is_number(x) and float(x) > 0.0 for x in q)):
        raise ValueError("q must be length-12 positive numeric list")
    qf = [float(x) for x in q]
    if abs(sum(qf) - 1.0) > 1e-9:
        raise ValueError("q must be normalized (sum=1)")

    sum_q2 = float(sum(float(x) ** 2 for x in qf))
    if not math.isfinite(sum_q2) or sum_q2 <= 0.0:
        raise ValueError("sum_q2 must be finite positive")
    g4 = float(lambda_psi / sum_q2)

    v_abs = [float(math.sqrt(float(rho_sq) * float(x))) for x in qf]

    any_solution = False
    solutions: list[dict[str, Any]] = []
    scanned = 0

    min_abs_by_entry = {"c1c1": float("inf"), "c1s1": float("inf"), "s1s1": float("inf")}
    best: dict[str, Any] | None = None

    for mask in range(1 << 11):
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
                "objective": "max_abs_entry",
                "objective_value": float(max_abs),
                "sign_vector": s,
                "equations": eqs,
            }

        if all(bool(eqs[k]["is_zero"]) for k in ("c1c1", "c1s1", "s1s1")):
            any_solution = True
            if len(solutions) < 10:
                solutions.append({"sign_vector": s, "equations": eqs})

    return {
        "reference_id": reference_id,
        "q": qf,
        "sum_q2": float(sum_q2),
        "g4_uniform_value": float(g4),
        "scan_space_size": scanned,
        "any_solution_within_tolerance": bool(any_solution),
        "solutions_within_tolerance_capped": solutions,
        "best_candidate": best,
        "min_abs_by_entry_over_scan": min_abs_by_entry,
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {
        "QW-2122": QW2122_JSON,
        "alpha_geo": ALPHA_GEO_JSON,
        "R14": IN_R14,
        "R15": IN_R15,
        "R18": IN_R18,
    }
    missing = [k for k, p in required.items() if not p.is_file()]
    if missing:
        payload = {
            "stage": "P479",
            "date_utc": datetime.now(timezone.utc).date().isoformat(),
            "goal": "scan_a_fixed_family_of_strict_reference_distributions_as_T169_magnitude_lifts_for_R18_pair1_zero_equations_under_conditional_N477",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {"stage": "P479", "status": payload["status"], "missing_required_inputs": missing, "no_false_pass": True},
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    qw2122 = _read_json(QW2122_JSON)
    alpha_geo_obj = _read_json(ALPHA_GEO_JSON)
    r14 = _read_json(IN_R14)
    r15 = _read_json(IN_R15)
    r18 = _read_json(IN_R18)

    alpha_geo = _parse_alpha_geo_numeric(alpha_geo_obj)
    rho_sq, lambda_psi = _extract_qw2122_strict_scalars(qw2122)

    ktotal = _extract_ktotal(r14)
    m0_sq = _extract_m0_sq(r15)
    classes = _extract_r18_classes(r18)

    ords = _ord_z12_by_x()
    row0 = [float(x) for x in ktotal[0]]

    def ref_uniform() -> list[float]:
        return [1.0 / 12.0 for _ in range(12)]

    def ref_ordpow() -> list[float]:
        if IN_RORDPOW.is_file():
            obj = _read_json(IN_RORDPOW)
            q = obj.get("reference_prob")
            if isinstance(q, list) and len(q) == 12 and all(_is_number(x) and float(x) > 0.0 for x in q):
                qf = [float(x) for x in q]
                if abs(sum(qf) - 1.0) <= 1e-9:
                    return qf
        w = [float(math.exp(-float(alpha_geo) * math.log(float(o)))) for o in ords]
        return _normalize_positive_weights(w)

    def ref_ord_exp() -> list[float]:
        w = [float(math.exp(-float(alpha_geo) * float(o))) for o in ords]
        return _normalize_positive_weights(w)

    def ref_dcyc_exp() -> list[float]:
        w = [float(math.exp(-float(alpha_geo) * float(_dcyc_0(i)))) for i in range(12)]
        return _normalize_positive_weights(w)

    def ref_1_plus_abs_ktotal_row0() -> list[float]:
        w = [1.0 + float(abs(x)) for x in row0]
        return _normalize_positive_weights(w)

    def ref_exp_sigma_int_alpha_geo_ktotal_row0() -> list[float]:
        sigma_int = -1.0
        w = [float(math.exp(float(sigma_int) * float(alpha_geo) * float(x))) for x in row0]
        return _normalize_positive_weights(w)

    reference_fns: list[tuple[str, str, Callable[[], list[float]]]] = [
        ("r_uniform", "Uniform reference distribution on Z_12 (translation-invariant control)", ref_uniform),
        ("r_ordpow", "r(x) ∝ ord_Z12(x)^(-alpha_geo) with alpha_geo=4 ln 2 (Aut(Z_12)-invariant; no marked direction)", ref_ordpow),
        ("r_ord_exp", "r(x) ∝ exp(-alpha_geo * ord_Z12(x)) with alpha_geo=4 ln 2 (Aut(Z_12)-invariant; no marked direction)", ref_ord_exp),
        ("r_dcyc_exp", "r(x) ∝ exp(-alpha_geo * d_cyc(0,x)) with alpha_geo=4 ln 2 (marked-site distance decay; direction-free)", ref_dcyc_exp),
        ("r_1_plus_abs_ktotal_row0", "r(x) ∝ 1 + |K_total[0,x]| (marked-site kernel row; direction-free)", ref_1_plus_abs_ktotal_row0),
        (
            "r_exp_sigma_int_alpha_geo_ktotal_row0",
            "r(x) ∝ exp(sigma_int * alpha_geo * K_total[0,x]) with sigma_int=-1, alpha_geo=4 ln 2 (marked-site kernel-shaped; direction-free)",
            ref_exp_sigma_int_alpha_geo_ktotal_row0,
        ),
    ]

    zero_tol = 1e-12
    per_reference: list[dict[str, Any]] = []
    any_solution = False

    best_overall: dict[str, Any] | None = None

    for ref_id, ref_desc, fn in reference_fns:
        q = fn()
        res = _scan_sign_space_for_reference(
            reference_id=ref_id,
            q=q,
            rho_sq=rho_sq,
            lambda_psi=lambda_psi,
            ktotal=ktotal,
            m0_sq=m0_sq,
            classes=classes,
            zero_tol=zero_tol,
        )
        res["reference_description"] = ref_desc
        per_reference.append(res)

        if bool(res.get("any_solution_within_tolerance")):
            any_solution = True

        best = res.get("best_candidate")
        if isinstance(best, dict) and _is_number(best.get("objective_value")):
            if best_overall is None or float(best["objective_value"]) < float(best_overall["objective_value"]):
                best_overall = {
                    "reference_id": ref_id,
                    "reference_description": ref_desc,
                    **best,
                }

    status = "PASS_SCAN_COMPLETE_NO_REFERENCE_HAS_SIGN_SOLUTION_FOR_R18_PAIR1_ZERO_EQUATIONS_UNDER_N477"
    if any_solution:
        status = "PASS_SCAN_COMPLETE_AT_LEAST_ONE_REFERENCE_HAS_SIGN_SOLUTION_FOR_R18_PAIR1_ZERO_EQUATIONS_UNDER_N477"

    obj = {
        "object": "Pair1ResidualZeroEquationsReferenceMagnitudeFamilyScan_under_N477_v1",
        "status": status,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Scan a fixed family of strictly-defined reference distributions to define a T169-like magnitude lift "
            "|vpsi_i|=sqrt(rho_*^2*q_i) and a uniform g4 lift g4=lambda_psi/sum q^2, then exhaustively scan all 2^11 sign vectors "
            "(fixing s0=+1) and evaluate the R18 declared pair1 residual zero equations under the conditional N477 Yukawa-free diagonal residual rewrite. "
            "This is probe-level evidence only: it does not export any vacuum stationarity witness and does not claim host matching or QW-2191 discharge."
        ),
        "inputs": {
            "qw2122_scalar_vacuum": _p(QW2122_JSON),
            "alpha_geo_strict_derived": _p(ALPHA_GEO_JSON),
            "r14_k_total": _p(IN_R14),
            "r15_m0_squared": _p(IN_R15),
            "r18_reduction_packet": _p(IN_R18),
        },
        "fixed_scalars": {"rho_star_sq": float(rho_sq), "lambda_psi_strict": float(lambda_psi), "alpha_geo": float(alpha_geo)},
        "scan": {
            "zero_tolerance": float(zero_tol),
            "sign_scan_space_size": 2048,
            "global_sign_fixed": "s0=+1 (global Z2 is ratio-invariant for N477 anyway)",
            "references_scanned": per_reference,
            "any_reference_has_solution": bool(any_solution),
            "best_overall_candidate": best_overall,
        },
        "hard_limits": [
            "Conditional on the N477 Yukawa-free diagonal residual rewrite; no vacuum stationarity witness is exported.",
            "Scans only a fixed family of reference distributions with fixed magnitudes and a fixed uniform g4 lift per reference.",
            "Does not export a strict zero witness unless a theorem-level witness is separately constructed from the scan result.",
            "Does not claim host matching, strict-core selector closure, QW-2191 discharge, or ToE closure.",
        ],
        "no_false_pass": True,
    }

    payload = {
        "stage": "P479",
        "date_utc": datetime.now(timezone.utc).date().isoformat(),
        "goal": "scan_a_fixed_family_of_strict_reference_distributions_as_T169_magnitude_lifts_for_R18_pair1_zero_equations_under_conditional_N477",
        "status": status,
        "exported_object": _p(OUT_OBJECT),
        "any_reference_has_solution": bool(any_solution),
        "best_overall_objective_value": None if best_overall is None else float(best_overall["objective_value"]),
        "best_overall_reference_id": None if best_overall is None else str(best_overall["reference_id"]),
        "no_false_pass": True,
    }

    summary = {
        "stage": "P479",
        "status": status,
        "exported_object": payload["exported_object"],
        "any_reference_has_solution": bool(any_solution),
        "best_overall_objective_value": payload["best_overall_objective_value"],
        "best_overall_reference_id": payload["best_overall_reference_id"],
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

