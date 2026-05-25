#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2085_s1035_strict_operator_family_expansion_audit.json"
MD = GEN / "p2085_s1035_strict_operator_family_expansion_audit.md"

SCHEMA_VERSION = "p2085_s1035_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def k_legacy(d: float) -> float:
    alpha_geo = 4.0 * math.log(2.0)
    omega = math.pi / 4.0
    phi = math.pi / 6.0
    beta_tors = 0.01
    return alpha_geo * math.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: float) -> float:
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8
    return math.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def finite_diff(y: list[float], h: float) -> list[float]:
    return [(y[i + 1] - y[i - 1]) / (2.0 * h) for i in range(1, len(y) - 1)]


def finite_diff2(y: list[float], h: float) -> list[float]:
    return [(y[i + 1] - 2.0 * y[i] + y[i - 1]) / (h * h) for i in range(1, len(y) - 1)]


def rmse(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a))


def sign_changes(y: list[float]) -> int:
    return sum(1 for a, b in zip(y[:-1], y[1:]) if a != 0 and a * b < 0)


def damp_score(y: list[float]) -> float:
    total = max(1, len(y) - 1)
    return sum(1 for a, b in zip(y[:-1], y[1:]) if abs(b) <= abs(a)) / total


def g_d1(d: float) -> float:
    return d * (1.0 - d) * math.exp(-d)


def g_d2(d: float) -> float:
    return (d**2) * math.exp(-0.5 * d)


def g_damp(d: float) -> float:
    return math.exp(-d) / (1.0 + d)


def g_phase_curv(d: float) -> float:
    return math.sin(0.7 * d + 0.2) * (d**2) * math.exp(-0.8 * d)


def g_nonlinear_damp(d: float) -> float:
    return math.exp(-(d**1.2)) / (1.0 + d**1.5)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2084 = load("p2084_s1034_strict_adaptive_channel_weight_audit.json")
    ready = p2084.get("result_kind") == "PASS_ADAPTIVE_CHANNEL_WEIGHT_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    d_grid = [i / 40.0 for i in range(1, 41)]
    h = d_grid[1] - d_grid[0]
    legacy = [k_legacy(d) for d in d_grid]
    strict = [k_strict(d) for d in d_grid]

    legacy_d1 = finite_diff(legacy, h)
    legacy_d2 = finite_diff2(legacy, h)
    strict_d1 = finite_diff(strict, h)
    strict_d2 = finite_diff2(strict, h)

    strict_shape = {
        "d1_rmse": rmse(legacy_d1, strict_d1),
        "d2_rmse": rmse(legacy_d2, strict_d2),
        "zero_gap": abs(sign_changes(legacy) - sign_changes(strict)),
        "damp_gap": abs(damp_score(legacy) - damp_score(strict)),
    }

    c_grid = [-0.2, -0.1, 0.0, 0.1, 0.2]
    lam = 0.05

    operator_families = {
        "phase_curvature_plus_base": [g_d1, g_d2, g_damp, g_phase_curv],
        "nonlinear_damping_plus_base": [g_d1, g_d2, g_damp, g_nonlinear_damp],
    }

    family_results = []
    for family_name, basis in operator_families.items():
        best = None
        for c1 in c_grid:
            for c2 in c_grid:
                for c3 in c_grid:
                    for c4 in c_grid:
                        coeffs = [c1, c2, c3, c4]
                        aug = [k_strict(d) + sum(c * g(d) for c, g in zip(coeffs, basis)) for d in d_grid]
                        aug_d1 = finite_diff(aug, h)
                        aug_d2 = finite_diff2(aug, h)
                        shape = {
                            "d1_rmse": rmse(legacy_d1, aug_d1),
                            "d2_rmse": rmse(legacy_d2, aug_d2),
                            "zero_gap": abs(sign_changes(legacy) - sign_changes(aug)),
                            "damp_gap": abs(damp_score(legacy) - damp_score(aug)),
                        }
                        objective = 2.0 * shape["d1_rmse"] + 2.0 * shape["d2_rmse"] + shape["zero_gap"] + shape["damp_gap"] + lam * sum(c * c for c in coeffs)
                        improved_dims = sum(shape[k] < strict_shape[k] for k in strict_shape)
                        rec = {
                            "coefficients": {"c1": c1, "c2": c2, "c3": c3, "c4": c4},
                            "shape_error": shape,
                            "weighted_objective": objective,
                            "improved_dimensions": improved_dims,
                            "qw2191_proxy_risk": "LOW_PROXY" if max(abs(c) for c in coeffs) <= 0.2 else "ELEVATED_PROXY",
                        }
                        if best is None or rec["weighted_objective"] < best["weighted_objective"]:
                            best = rec
        assert best is not None
        best["family_name"] = family_name
        best["decision"] = "RETAIN_EXPANDED_FAMILY" if best["improved_dimensions"] >= 3 and best["qw2191_proxy_risk"] == "LOW_PROXY" else "REJECT_EXPANDED_FAMILY"
        family_results.append(best)

    robust_families = [r for r in family_results if r["decision"] == "RETAIN_EXPANDED_FAMILY"]
    selected_family = min(family_results, key=lambda r: r["weighted_objective"])

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2085",
        "stage_id": "S1035",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_OPERATOR_FAMILY_EXPANSION_AUDIT_WITH_TRACE__C3_STILL_OPEN" if ready else "OPEN_OPERATOR_FAMILY_EXPANSION_AUDIT_BLOCKED",
        "depends_on": {"p2084_present": p2084.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2084_json_sha256": file_sha256(GEN / "p2084_s1034_strict_adaptive_channel_weight_audit.json"),
        },
        "expansion_protocol": {
            "d_grid": d_grid,
            "coefficient_bound": 0.2,
            "lambda_l2": lam,
            "operator_families": list(operator_families.keys()),
            "acceptance_rule": "improved_dimensions>=3 and qw2191_proxy_risk=LOW_PROXY",
        },
        "expansion_results": {
            "strict_shape_error": strict_shape,
            "family_results": family_results,
            "robust_family_count": len(robust_families),
            "selected_best_family": selected_family,
            "functional_basis_limit_hypothesis": "OPEN" if len(robust_families) == 0 else "PARTIAL_RELIEF",
        },
        "c3_gate_update": {
            "C3_operator_family_expansion_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "family_count_expected": len(family_results) == 2,
            "all_coefficients_within_bound": all(max(abs(float(r["coefficients"][k])) for k in ("c1", "c2", "c3", "c4")) <= 0.2 for r in family_results),
            "decision_domain_ok": all(r["decision"] in {"RETAIN_EXPANDED_FAMILY", "REJECT_EXPANDED_FAMILY"} for r in family_results),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2085 S1035: operator-family expansion audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Robust family count: `{payload['expansion_results']['robust_family_count']}`",
            f"- Best family: `{selected_family['family_name']}`",
            "",
            "Expanded operator families tested to identify whether basis insufficiency blocks shape closure.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
