#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2126 = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
OUT = GEN / "p2127_s1077_strict_cmp2_bootstrap_backend_evidence_stresstest.json"
MD = GEN / "p2127_s1077_strict_cmp2_bootstrap_backend_evidence_stresstest.md"

SCHEMA_VERSION = "p2127_s1077_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2127
N_BOOT = 400


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def rank_of_max(weights: dict[str, float]) -> str:
    return max(weights.items(), key=lambda kv: kv[1])[0]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2126 = load(IN_2126)

    pre_ready = p2126.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_WITH_TRACE"
    rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    n = len(rows)

    rng = np.random.default_rng(RNG_SEED)

    if n > 0:
        coverage = np.array([1.0 if r.get("posterior_predictive_covered", False) else 0.0 for r in rows], dtype=float)
        top_model = [rank_of_max((r.get("posterior_weights_backend_evidence", {}) or {})) for r in rows]

        boot_cov = []
        boot_top_counts = {"M1_nn": [], "M2_monotone": [], "M3_nn_tiebreak": [], "M4_monotone_penalized": []}

        for _ in range(N_BOOT):
            idx = rng.integers(0, n, size=n)
            sample_cov = coverage[idx]
            boot_cov.append(float(np.mean(sample_cov)))

            sample_top = [top_model[i] for i in idx]
            for k in boot_top_counts:
                boot_top_counts[k].append(sample_top.count(k) / n)

        cov_ci = [float(np.quantile(boot_cov, 0.025)), float(np.quantile(boot_cov, 0.975))]

        rank_stability = {
            k: {
                "mean_frequency": float(np.mean(v)),
                "ci95": [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))],
            }
            for k, v in boot_top_counts.items()
        }
    else:
        cov_ci = [None, None]
        rank_stability = {}

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2127",
        "stage_id": "S1077",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BOOTSTRAP_BACKEND_EVIDENCE_STRESSTEST_WITH_TRACE"
            if pre_ready and n > 0
            else "OPEN_STRICT_CMP2_BOOTSTRAP_BACKEND_EVIDENCE_STRESSTEST_BLOCKED"
        ),
        "depends_on": {
            "p2126_present": p2126.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "bootstrap_backend_evidence_stresstest": {
            "n_rows": n,
            "n_bootstrap": N_BOOT,
            "rng_seed": RNG_SEED,
            "coverage_rate_point": float((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("coverage_rate", 0.0)) if n > 0 else None,
            "coverage_rate_ci95": cov_ci,
            "posterior_top_model_rank_stability": rank_stability,
            "scope_limit": "operational bootstrap stress-test only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_BOOTSTRAP_STRESSTEST_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2128_candidate",
            "goal": "run block-bootstrap by s-clusters and compare CI inflation vs iid bootstrap",
        },
        "c3_gate_update": {
            "C3_cmp2_backend_evidence_weighted_posterior_calibration": "COMPUTED",
            "C3_cmp2_bootstrap_stresstest": "COMPUTED" if n > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "bootstrap_executed": n > 0,
            "coverage_ci_exported": n > 0,
            "rank_stability_exported": n > 0,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2127 S1077: strict CMP2 bootstrap backend-evidence stress-test",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Bootstrap samples: `{N_BOOT}`",
            f"- Coverage CI95: `{cov_ci}`",
            "",
            "This stage stress-tests backend-evidence weighting with bootstrap resampling and exports uncertainty intervals for coverage and posterior rank stability.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
