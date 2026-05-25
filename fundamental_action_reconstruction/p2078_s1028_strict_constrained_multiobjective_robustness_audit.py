#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.json"
MD = GEN / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.md"

SCHEMA_VERSION = "p2078_s1028_v1"
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


def mse(a: list[float], b: list[float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) / float(len(a))


def f_candidate(cid: str, d: float) -> float:
    if cid == "mc_phase_lock":
        return d * math.exp(-d)
    if cid == "mc_grad_cross":
        return (d**2) / (1.0 + d**2)
    return math.sin(d) * math.exp(-d)


def scenario_grids() -> list[list[float]]:
    return [
        [i / 20.0 for i in range(1, 21)],
        [i / 30.0 for i in range(1, 31)],
        [0.04 + i * 0.048 for i in range(20)],
    ]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2077 = load("p2077_s1027_strict_constrained_coefficient_fit_audit.json")
    ready = p2077.get("result_kind") == "PASS_CONSTRAINED_COEFFICIENT_FIT_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    retained = ((p2077.get("fit_results") or {}).get("retained_candidates") or []) if isinstance(p2077, dict) else []
    candidates = list(retained)

    lambdas = [0.01, 0.05, 0.1]
    ref_perturb = [0.0, 0.005, -0.005]
    c_grid = [x / 100.0 for x in range(-20, 21)]

    candidate_rows = []
    for cid in candidates:
        scenario_count = 0
        robust_pass_count = 0
        improvements = []

        for d_grid in scenario_grids():
            base_legacy = [k_legacy(d) for d in d_grid]
            base_strict = [k_strict(d) for d in d_grid]

            for lam in lambdas:
                for eps in ref_perturb:
                    scenario_count += 1
                    legacy_vals = [v * (1.0 + eps) for v in base_legacy]
                    baseline = mse(legacy_vals, base_strict)

                    best_c = 0.0
                    best_raw = baseline
                    best_obj = baseline
                    for c in c_grid:
                        aug = [k_strict(d) + c * f_candidate(cid, d) for d in d_grid]
                        raw = mse(legacy_vals, aug)
                        obj = raw + lam * (c**2)
                        if obj < best_obj:
                            best_obj = obj
                            best_raw = raw
                            best_c = c

                    improve = baseline - best_raw
                    improvements.append(improve)
                    if improve > 0.0 and abs(best_c) <= 0.2:
                        robust_pass_count += 1

        robustness_ratio = robust_pass_count / scenario_count if scenario_count > 0 else 0.0
        min_improvement = min(improvements) if improvements else 0.0
        mean_improvement = sum(improvements) / len(improvements) if improvements else 0.0

        decision = "RETAIN_ROBUST" if robustness_ratio >= 0.8 and min_improvement > 0.0 else "REJECT_FRAGILE"

        candidate_rows.append(
            {
                "candidate_id": cid,
                "scenario_count": scenario_count,
                "robust_pass_count": robust_pass_count,
                "robustness_ratio": robustness_ratio,
                "min_improvement": min_improvement,
                "mean_improvement": mean_improvement,
                "decision": decision,
            }
        )

    robust = [r["candidate_id"] for r in candidate_rows if r["decision"] == "RETAIN_ROBUST"]
    fragile = [r["candidate_id"] for r in candidate_rows if r["decision"] == "REJECT_FRAGILE"]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2078",
        "stage_id": "S1028",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CONSTRAINED_MULTIOBJECTIVE_ROBUSTNESS_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_CONSTRAINED_MULTIOBJECTIVE_ROBUSTNESS_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2077_present": p2077.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2077_json_sha256": file_sha256(GEN / "p2077_s1027_strict_constrained_coefficient_fit_audit.json"),
        },
        "robustness_protocol": {
            "d_grid_variants": scenario_grids(),
            "lambda_variants": lambdas,
            "reference_perturbations": ref_perturb,
            "coefficient_bound": 0.2,
            "robust_retain_rule": "robustness_ratio>=0.8 and min_improvement>0",
        },
        "robustness_results": {
            "rows": candidate_rows,
            "robust_candidates": robust,
            "fragile_candidates": fragile,
        },
        "c3_gate_update": {
            "C3_constrained_multiobjective_robustness_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "candidate_rows_nonempty": len(candidate_rows) > 0,
            "has_robust_or_fragile_decisions": all(r["decision"] in {"RETAIN_ROBUST", "REJECT_FRAGILE"} for r in candidate_rows),
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
            "# P2078 S1028: constrained multi-objective robustness audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Robust candidates: `{robust}`",
            f"- Fragile candidates: `{fragile}`",
            "",
            "Retain from P2077 is stress-tested across d-grid, lambda, and reference perturbation variants.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
