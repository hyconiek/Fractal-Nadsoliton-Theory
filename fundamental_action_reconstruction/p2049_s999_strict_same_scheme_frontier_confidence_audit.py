#!/usr/bin/env python3
"""P2049 S999: frontier confidence audit.

Repeats P2048-style frontier search across multiple trial-seed/sign-pattern
ensembles and exports confidence bounds for epsilon_break (or right-censored
lower bounds when break is not detected).

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
OUT = GEN / "p2049_s999_strict_same_scheme_frontier_confidence_audit.json"
MD = GEN / "p2049_s999_strict_same_scheme_frontier_confidence_audit.md"

SCHEMA_VERSION = "p2049_s999_v1"
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
    p2048 = load("p2048_s998_strict_same_scheme_certificate_falsification_frontier_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2048 = p2048.get("gatekeeper_checks") or {}

    ready = (
        p2048.get("result_kind")
        == "PASS_CERTIFICATE_FALSIFICATION_FRONTIER_AUDIT_WITH_TRACE__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2048.get("epsilon_last_stable_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    propagated = float(((p2040.get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", 0.0)))

    frontier = p2048.get("certificate_frontier") or {}
    eps_star = float(frontier.get("epsilon_star_input", 0.0))

    fscope = p2048.get("frontier_scope") or {}
    norm_buckets = [float(x) for x in (fscope.get("norm_buckets") or [0.01, 0.015, 0.02, 0.025])]
    vectors_per_bucket = int(fscope.get("vectors_per_bucket", 20))

    cert_scope = (load("p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json").get("certificate_scope") or {})
    ranking_ref = cert_scope.get("ranking_reference") or {}
    base_rank = ranking_ref.get("base_rank") or [4, 3, 2, 1]

    # Ensembles for confidence audit
    ensemble_count = 12
    ensemble_rng = np.random.default_rng(2049)

    def is_stable(eps: float, trial_seeds: list[int], sign_patterns: list[np.ndarray]) -> bool:
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

    ensemble_rows = []
    epsilon_break_samples = []
    censored_lower_bounds = []

    for eidx in range(ensemble_count):
        trial_seeds = [int(x) for x in ensemble_rng.integers(5000, 9000, size=3)]
        sign_patterns = []
        for _ in range(18):
            s = np.sign(ensemble_rng.normal(size=delta.shape))
            s[s == 0] = 1.0
            sign_patterns.append(s)

        step = max(2e-7, eps_star * 0.2 + 1e-7)
        eps_lo = eps_star
        eps_hi = eps_star
        found_break = False

        for _ in range(18):
            eps_hi += step
            st = is_stable(eps_hi, trial_seeds, sign_patterns)
            if not st:
                found_break = True
                break
            eps_lo = eps_hi
            step *= 1.35

        eps_last_stable = eps_lo
        eps_break = math.inf

        if found_break:
            lo = eps_lo
            hi = eps_hi
            for _ in range(20):
                mid = 0.5 * (lo + hi)
                st = is_stable(mid, trial_seeds, sign_patterns)
                if st:
                    lo = mid
                else:
                    hi = mid
            eps_last_stable = lo
            eps_break = hi
            epsilon_break_samples.append(eps_break)
        else:
            censored_lower_bounds.append(eps_last_stable)

        ensemble_rows.append(
            {
                "ensemble_id": eidx,
                "trial_seeds": trial_seeds,
                "sign_pattern_count": len(sign_patterns),
                "break_found": found_break,
                "epsilon_last_stable": eps_last_stable,
                "epsilon_break": eps_break,
            }
        )

    if epsilon_break_samples:
        arr = np.array(epsilon_break_samples, dtype=float)
        ci_low = float(np.quantile(arr, 0.025))
        ci_high = float(np.quantile(arr, 0.975))
        center = float(np.median(arr))
    else:
        ci_low = math.inf
        ci_high = math.inf
        center = math.inf

    right_censored_lower = float(np.max(censored_lower_bounds)) if censored_lower_bounds else -math.inf

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2049",
        "stage_id": "S999",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_FRONTIER_CONFIDENCE_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready and (math.isfinite(center) or bool(censored_lower_bounds))
            else "OPEN_FRONTIER_CONFIDENCE_AUDIT_BLOCKED"
        ),
        "route": "strict_only_frontier_confidence_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2048_present": p2048.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2048_json_sha256": file_sha256(GEN / "p2048_s998_strict_same_scheme_certificate_falsification_frontier_audit.json"),
        },
        "audit_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "ensemble_count": ensemble_count,
            "norm_buckets": norm_buckets,
            "vectors_per_bucket": vectors_per_bucket,
            "epsilon_star_input": eps_star,
        },
        "ensemble_results": ensemble_rows,
        "frontier_confidence": {
            "break_samples_count": len(epsilon_break_samples),
            "right_censored_count": len(censored_lower_bounds),
            "epsilon_break_median": center,
            "epsilon_break_ci_level": 0.95,
            "epsilon_break_ci_low": ci_low,
            "epsilon_break_ci_high": ci_high,
            "right_censored_lower_bound": right_censored_lower,
            "interpretation": (
                "If break_samples_count>0, CI summarizes detected epsilon_break variation. "
                "If right_censored_count>0, right_censored_lower_bound records unresolved lower frontier."
            ),
        },
        "c3_gate_update": {
            "C3_frontier_confidence_audit": "COMPUTED",
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
            "ensemble_nonempty": len(ensemble_rows) > 0,
            "break_or_censored_present": (len(epsilon_break_samples) > 0) or (len(censored_lower_bounds) > 0),
            "confidence_payload_consistent": (len(epsilon_break_samples) > 0 and math.isfinite(center)) or (len(censored_lower_bounds) > 0),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2049 S999: frontier confidence audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Repeated P2048 frontier runs across multiple ensembles and exported CI",
        "for epsilon_break (or right-censored lower bound where applicable).",
        "",
        "## Gate update",
        "",
        "- `C3`: frontier confidence audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
