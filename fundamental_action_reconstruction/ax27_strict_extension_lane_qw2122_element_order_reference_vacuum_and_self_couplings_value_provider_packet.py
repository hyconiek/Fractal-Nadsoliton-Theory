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

QW2122_JSON = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"

OUT_JSON = (
    GENERATED
    / "ax27_strict_extension_lane_qw2122_element_order_reference_vacuum_and_self_couplings_value_provider_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "ax27_strict_extension_lane_qw2122_element_order_reference_vacuum_and_self_couplings_value_provider_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def zn_element_order(n: int, k: int) -> int:
    if n <= 0:
        raise ValueError("n must be positive")
    kk = k % n
    if kk == 0:
        return 1
    return n // math.gcd(kk, n)


def normalize_positive_weights(w: list[float]) -> list[float]:
    if not w:
        raise ValueError("weights must be non-empty")
    if any((not isinstance(x, (int, float))) or (not math.isfinite(float(x))) or float(x) <= 0.0 for x in w):
        raise ValueError("weights must be finite and strictly positive")
    z = float(sum(float(x) for x in w))
    if z <= 0.0:
        raise ValueError("weights must sum to positive")
    return [float(x) / z for x in w]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not QW2122_JSON.exists():
        summary = {
            "stage": "AX27",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILE",
            "missing_file": str(QW2122_JSON.relative_to(REPO)),
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    r2122 = load_json(QW2122_JSON)
    lambda_psi_strict = (r2122.get("inputs") or {}).get("lambda_psi_strict")
    rho_star_sq = ((r2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("rho_star_sq")

    missing: list[str] = []
    if not isinstance(lambda_psi_strict, (int, float)) or not math.isfinite(float(lambda_psi_strict)):
        missing.append("QW-2122.inputs.lambda_psi_strict (finite number)")
    if not isinstance(rho_star_sq, (int, float)) or not math.isfinite(float(rho_star_sq)):
        missing.append("QW-2122.branch_results.broken_branch_strict.rho_star_sq (finite number)")

    if missing:
        summary = {
            "stage": "AX27",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FIELDS",
            "missing_fields": missing,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    lambda_psi_strict_f = float(lambda_psi_strict)
    rho_star_sq_f = float(rho_star_sq)

    n = 12
    alpha_geo = 4.0 * math.log(2.0)

    # Extension-lane mapping (explicit premises):
    # - treat rho_*^2 from QW-2122 as a strict scalar norm constraint: rho_*^2 = sum_i vpsi_i^2 on the 12-slot carrier,
    # - construct a non-translation-invariant reference distribution on Z_12 from intrinsic group structure:
    #     w_i := exp(-alpha_geo * ord_Z12(i))  where ord_Z12(i) is the order of element i in Z_12,
    #     r_i := w_i / sum_j w_j,
    # - choose the per-site vacuum magnitudes by vpsi_i^2 := rho_*^2 * r_i (no zeros, no marked direction),
    # - map the scalar quartic coefficient lambda_psi to per-site self quartics by matching along the uniform ray:
    #     g4_i := 12*lambda_psi_strict,  g6_i := 0.
    orders = [zn_element_order(n, i) for i in range(n)]
    weights = [math.exp(-alpha_geo * float(o)) for o in orders]
    rprob = normalize_positive_weights(weights)

    vpsi = [math.sqrt(rho_star_sq_f * float(r)) for r in rprob]
    g4_val = 12.0 * lambda_psi_strict_f
    g4 = [float(g4_val) for _ in range(n)]
    g6 = [0.0 for _ in range(n)]

    artifact = {
        "stage": "AX27",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "STRICT_EXTENSION_ONLY_VALUE_PROVIDER_PACKET",
        "scope": "strict_extension_only",
        "goal": "provide_one_explicit_reproducible_candidate_value_provider_for_(vpsi,g4,g6)_to_enable_exploratory_P437/P434_execution_without_promoting_it_to_strict_core",
        "provenance": {
            "qw2122_source": str(QW2122_JSON.relative_to(REPO)),
            "qw2122_fields_used": {
                "lambda_psi_strict": "inputs.lambda_psi_strict",
                "rho_star_sq": "branch_results.broken_branch_strict.rho_star_sq",
            },
        },
        "explicit_mapping_premises": {
            "group_structure_reference": "w_i := exp(-alpha_geo * ord_Z12(i)) with alpha_geo=4 ln 2; r_i := w_i/sum_j w_j",
            "nonuniform_vacuum_magnitudes": "vpsi_i^2 := rho_star_sq * r_i (all positive; no marked direction)",
            "uniform_self_couplings": "g4_i := 12*lambda_psi_strict, g6_i := 0 for i=0..11",
            "warning": "This packet is NOT a strict-derived discharge of T168; it adds explicit mapping premises to run evaluation harnesses.",
        },
        "derived_reference": {
            "alpha_geo_used": alpha_geo,
            "ord_Z12": orders,
            "weights": weights,
            "reference_prob": rprob,
        },
        # Compatibility: P437 expects numeric lists at top-level keys (vpsi,g4,g6).
        "vpsi": vpsi,
        "g4": g4,
        "g6": g6,
        "values": {"vpsi": vpsi, "g4": g4, "g6": g6},
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX27",
        "status": "PASS_EXTENSION_VALUE_PROVIDER_EXPORTED",
        "export": str(OUT_JSON.relative_to(REPO)),
        "notes": "Extension-only candidate values for vpsi/g4/g6 derived from QW-2122 scalar vacuum closure + a Z_12 element-order reference distribution (non-translation-invariant; no marked direction).",
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

