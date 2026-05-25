#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2083_s1033_strict_three_feature_targeted_composite_operator_audit.json"
MD = GEN / "p2083_s1033_strict_three_feature_targeted_composite_operator_audit.md"

SCHEMA_VERSION = "p2083_s1033_v1"
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
    return [(y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h) for i in range(1, len(y) - 1)]


def rmse(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a))


def sign_changes(y: list[float]) -> int:
    c = 0
    for a, b in zip(y[:-1], y[1:]):
        if a != 0 and a * b < 0:
            c += 1
    return c


def damp_score(y: list[float]) -> float:
    total = max(1, len(y) - 1)
    return sum(1 for a, b in zip(y[:-1], y[1:]) if abs(b) <= abs(a)) / total


def g_d1(d: float) -> float:
    return d * (1.0 - d) * math.exp(-d)


def g_d2(d: float) -> float:
    return (d**2) * math.exp(-0.5 * d)


def g_damp(d: float) -> float:
    return math.exp(-d) / (1.0 + d)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2082 = load("p2082_s1032_strict_two_parameter_coupled_redesign_audit.json")
    ready = p2082.get("result_kind") == "PASS_TWO_PARAMETER_COUPLED_REDESIGN_AUDIT_WITH_TRACE__C3_STILL_OPEN"

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

    c_grid = [x / 100.0 for x in range(-20, 21)]
    lam = 0.05
    # separate channel weights for operator components
    w1, w2, w3 = 1.5, 1.5, 1.0
    best = None
    for c1 in c_grid:
        for c2 in c_grid:
            for c3 in c_grid:
                aug = [k_strict(d) + w1 * c1 * g_d1(d) + w2 * c2 * g_d2(d) + w3 * c3 * g_damp(d) for d in d_grid]
                aug_d1 = finite_diff(aug, h)
                aug_d2 = finite_diff2(aug, h)
                shape = {
                    "d1_rmse": rmse(legacy_d1, aug_d1),
                    "d2_rmse": rmse(legacy_d2, aug_d2),
                    "zero_gap": abs(sign_changes(legacy) - sign_changes(aug)),
                    "damp_gap": abs(damp_score(legacy) - damp_score(aug)),
                }
                obj = 2.0 * shape["d1_rmse"] + 2.0 * shape["d2_rmse"] + shape["zero_gap"] + shape["damp_gap"] + lam * (c1**2 + c2**2 + c3**2)
                if best is None or obj < best["obj"]:
                    best = {"c1": c1, "c2": c2, "c3": c3, "shape": shape, "obj": obj}

    assert best is not None
    improved_dims = sum(best["shape"][k] < strict_shape[k] for k in strict_shape)
    decision = "RETAIN_COMPOSITE" if improved_dims >= 3 else "REJECT_COMPOSITE"

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2083",
        "stage_id": "S1033",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_THREE_FEATURE_TARGETED_COMPOSITE_OPERATOR_AUDIT_WITH_TRACE__C3_STILL_OPEN" if ready else "OPEN_THREE_FEATURE_TARGETED_COMPOSITE_OPERATOR_AUDIT_BLOCKED",
        "depends_on": {"p2082_present": p2082.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2082_json_sha256": file_sha256(GEN / "p2082_s1032_strict_two_parameter_coupled_redesign_audit.json"),
        },
        "composite_protocol": {
            "d_grid": d_grid,
            "coefficient_bound": 0.2,
            "lambda_l2": lam,
            "channel_weights": {"w_d1": w1, "w_d2": w2, "w_damp": w3},
            "retain_rule": "improved_dimensions>=3",
        },
        "composite_results": {
            "strict_shape_error": strict_shape,
            "best_c1": best["c1"],
            "best_c2": best["c2"],
            "best_c3": best["c3"],
            "best_shape_error": best["shape"],
            "weighted_objective": best["obj"],
            "improved_dimensions": improved_dims,
            "decision": decision,
        },
        "c3_gate_update": {
            "C3_three_feature_targeted_composite_operator_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "coefficients_within_bound": all(abs(float(best[k])) <= 0.2 for k in ("c1", "c2", "c3")),
            "decision_domain_ok": decision in {"RETAIN_COMPOSITE", "REJECT_COMPOSITE"},
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
            "# P2083 S1033: three-feature targeted composite operator audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- best_c1/c2/c3: `{best['c1']}/{best['c2']}/{best['c3']}`",
            f"- decision: `{decision}`",
            "",
            "Composite operator targets d1/d2/damping channels with separate weights.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
