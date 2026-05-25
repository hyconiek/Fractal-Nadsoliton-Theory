#!/usr/bin/env python3
"""P2048 S998: certificate falsification frontier audit.

Searches for the first epsilon_break above epsilon* where bucket ranking fails,
using adaptive scan + bisection around the break frontier.
Audit-level only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2048_s998_strict_same_scheme_certificate_falsification_frontier_audit.json"
MD = GEN / "p2048_s998_strict_same_scheme_certificate_falsification_frontier_audit.md"

SCHEMA_VERSION = "p2048_s998_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(v: Any) -> bool:
    return bool(v is True)


def linf(v: np.ndarray) -> float:
    return float(np.max(np.abs(v)))


def rank_desc(vals: list[float]) -> list[int]:
    idx = sorted(range(len(vals)), key=lambda i: (-vals[i], i))
    out = [0] * len(vals)
    for r, i in enumerate(idx):
        out[i] = r + 1
    return out


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2047 = load("p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2047 = p2047.get("gatekeeper_checks") or {}

    ready = (
        p2047.get("result_kind")
        == "PASS_LOCAL_STABILITY_CERTIFICATE_EXTRACTED_WITH_TRACE__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2047.get("epsilon_star_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    propagated = float(((p2040.get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", 0.0)))

    cert = p2047.get("certificate_extraction") or {}
    eps_star = float(cert.get("epsilon_star_max_stable", 0.0))

    scope = p2047.get("certificate_scope") or {}
    norm_condition = scope.get("norm_condition") or {}
    norm_buckets = [float(x) for x in (norm_condition.get("norm_buckets_checked") or [0.01, 0.015, 0.02, 0.025])]
    trial_seeds = [int(x) for x in (norm_condition.get("trial_seeds") or [4013, 4019, 4021])]
    vectors_per_bucket = int(norm_condition.get("vectors_per_bucket", 20))

    ranking_ref = scope.get("ranking_reference") or {}
    base_rank = ranking_ref.get("base_rank") or [4, 3, 2, 1]

    rng = np.random.default_rng(2048)
    sign_patterns = []
    for _ in range(28):
        s = np.sign(rng.normal(size=delta.shape))
        s[s == 0] = 1.0
        sign_patterns.append(s)

    def is_stable(eps: float) -> bool:
        for sgn in sign_patterns:
            M = np.eye(3, dtype=float) + delta + eps * sgn
            per_seed = []
            for sd in trial_seeds:
                rg = np.random.default_rng(sd)
                bucket_worst = []
                for n in norm_buckets:
                    worst = 0.0
                    for _ in range(vectors_per_bucket):
                        v = rg.normal(size=3)
                        nv = float(np.linalg.norm(v))
                        if nv == 0.0:
                            continue
                        a = n * (v / nv)
                        res = linf(a - M.dot(a)) + propagated
                        worst = max(worst, res)
                    bucket_worst.append(worst)
                per_seed.append(bucket_worst)
            avg = list(np.mean(np.array(per_seed), axis=0))
            if rank_desc(avg) != base_rank:
                return False
        return True

    # Adaptive scan above eps_star
    step = max(2e-7, eps_star * 0.2 + 1e-7)
    eps_lo = eps_star
    eps_hi = eps_star
    scan_trace = []
    found_break = False

    for _ in range(20):
        eps_hi += step
        st = is_stable(eps_hi)
        scan_trace.append({"epsilon": eps_hi, "stable": st})
        if not st:
            found_break = True
            break
        eps_lo = eps_hi
        step *= 1.35

    # Bisection refinement around first unstable region
    bisect_trace = []
    epsilon_break = math.inf
    epsilon_last_stable = eps_lo
    if found_break:
        lo = eps_lo
        hi = eps_hi
        for _ in range(22):
            mid = 0.5 * (lo + hi)
            st = is_stable(mid)
            bisect_trace.append({"epsilon": mid, "stable": st})
            if st:
                lo = mid
            else:
                hi = mid
        epsilon_last_stable = lo
        epsilon_break = hi

    frontier_gap = (epsilon_break - epsilon_last_stable) if math.isfinite(epsilon_break) else math.inf

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2048",
        "stage_id": "S998",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CERTIFICATE_FALSIFICATION_FRONTIER_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready and math.isfinite(epsilon_last_stable)
            else "OPEN_CERTIFICATE_FALSIFICATION_FRONTIER_AUDIT_BLOCKED"
        ),
        "route": "strict_only_certificate_frontier_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2047_present": p2047.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2047_json_sha256": file_sha256(GEN / "p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json"),
        },
        "frontier_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "norm_buckets": norm_buckets,
            "trial_seeds": trial_seeds,
            "vectors_per_bucket": vectors_per_bucket,
            "sign_patterns_checked": len(sign_patterns),
        },
        "certificate_frontier": {
            "epsilon_star_input": eps_star,
            "epsilon_last_stable": epsilon_last_stable,
            "epsilon_break_first_detected": epsilon_break,
            "frontier_gap": frontier_gap,
            "scan_trace": scan_trace,
            "bisect_trace": bisect_trace,
            "break_found": found_break,
            "statement": (
                "Frontier audit locates the first detected instability epsilon_break above epsilon* "
                "in the present local adversarial model."
            ),
        },
        "c3_gate_update": {
            "C3_certificate_falsification_frontier_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "scan_trace_nonempty": len(scan_trace) > 0,
            "epsilon_last_stable_finite": math.isfinite(epsilon_last_stable),
            "break_or_censored": (found_break and math.isfinite(epsilon_break)) or (not found_break),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2048 S998: certificate falsification frontier audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported adaptive scan+bisection frontier search for first epsilon_break",
        "above epsilon* where ranking stability fails (or censoring if not found).",
        "",
        "## Gate update",
        "",
        "- `C3`: frontier audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
