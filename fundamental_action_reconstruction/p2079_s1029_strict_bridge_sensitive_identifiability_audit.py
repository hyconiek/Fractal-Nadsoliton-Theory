#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2079_s1029_strict_bridge_sensitive_identifiability_audit.json"
MD = GEN / "p2079_s1029_strict_bridge_sensitive_identifiability_audit.md"

SCHEMA_VERSION = "p2079_s1029_v1"
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


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2078 = load("p2078_s1028_strict_constrained_multiobjective_robustness_audit.json")
    ready = p2078.get("result_kind") == "PASS_CONSTRAINED_MULTIOBJECTIVE_ROBUSTNESS_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    robust_ids = ((p2078.get("robustness_results") or {}).get("robust_candidates") or []) if isinstance(p2078, dict) else []

    d_grid = [i / 20.0 for i in range(1, 21)]
    legacy_vals = [k_legacy(d) for d in d_grid]
    strict_vals = [k_strict(d) for d in d_grid]
    baseline = mse(legacy_vals, strict_vals)

    c_probe = [x / 200.0 for x in range(-40, 41)]  # [-0.2, 0.2] step 0.005
    rows = []
    for cid in robust_ids:
        err_curve = []
        for c in c_probe:
            aug = [k_strict(d) + c * f_candidate(cid, d) for d in d_grid]
            err_curve.append((c, mse(legacy_vals, aug)))

        best_c, best_err = min(err_curve, key=lambda t: t[1])

        near = [t for t in err_curve if abs(t[1] - best_err) <= 1e-4]
        c_span = max((c for c, _ in near), default=best_c) - min((c for c, _ in near), default=best_c)

        curvature_proxy = 0.0
        for i in range(1, len(err_curve) - 1):
            c0, e0 = err_curve[i - 1]
            c1, e1 = err_curve[i]
            c2, e2 = err_curve[i + 1]
            if abs(c1 - best_c) < 1e-12:
                h = c1 - c0
                if h != 0:
                    curvature_proxy = (e2 - 2 * e1 + e0) / (h**2)
                break

        bridge_signal = baseline - best_err
        identifiable = (c_span <= 0.06) and (curvature_proxy > 0.01)

        decision = "RETAIN_IDENTIFIABLE" if identifiable and bridge_signal > 0 else "REJECT_DEGENERATE"

        rows.append(
            {
                "candidate_id": cid,
                "baseline_mse": baseline,
                "best_mse": best_err,
                "best_coefficient": best_c,
                "bridge_signal_absolute": bridge_signal,
                "near_optimal_c_span": c_span,
                "curvature_proxy": curvature_proxy,
                "identifiable": identifiable,
                "decision": decision,
            }
        )

    retained = [r["candidate_id"] for r in rows if r["decision"] == "RETAIN_IDENTIFIABLE"]
    rejected = [r["candidate_id"] for r in rows if r["decision"] == "REJECT_DEGENERATE"]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2079",
        "stage_id": "S1029",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_BRIDGE_SENSITIVE_IDENTIFIABILITY_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_BRIDGE_SENSITIVE_IDENTIFIABILITY_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2078_present": p2078.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2078_json_sha256": file_sha256(GEN / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.json"),
        },
        "identifiability_protocol": {
            "coefficient_probe_grid": c_probe,
            "near_optimal_tolerance_mse": 1e-4,
            "identifiable_rule": "c_span<=0.06 and curvature_proxy>0.01 and bridge_signal>0",
        },
        "identifiability_results": {
            "rows": rows,
            "retained_identifiable_candidates": retained,
            "rejected_degenerate_candidates": rejected,
        },
        "c3_gate_update": {
            "C3_bridge_sensitive_identifiability_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "decision_domain_ok": all(r["decision"] in {"RETAIN_IDENTIFIABLE", "REJECT_DEGENERATE"} for r in rows),
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
            "# P2079 S1029: bridge-sensitive identifiability audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Retained identifiable: `{retained}`",
            f"- Rejected degenerate: `{rejected}`",
            "",
            "This stage separates bridge-signal candidates from degenerate numeric compensation.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
