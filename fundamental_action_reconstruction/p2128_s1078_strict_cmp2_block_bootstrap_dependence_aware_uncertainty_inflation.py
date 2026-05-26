#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2126 = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
IN_2127 = GEN / "p2127_s1077_strict_cmp2_bootstrap_backend_evidence_stresstest.json"
OUT = GEN / "p2128_s1078_strict_cmp2_block_bootstrap_dependence_aware_uncertainty_inflation.json"
MD = GEN / "p2128_s1078_strict_cmp2_block_bootstrap_dependence_aware_uncertainty_inflation.md"

SCHEMA_VERSION = "p2128_s1078_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2128
N_BOOT = 400


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def rank_of_max(weights: dict[str, float]) -> str:
    return max(weights.items(), key=lambda kv: kv[1])[0]


def ci_width(ci: list[float | None]) -> float | None:
    if ci[0] is None or ci[1] is None:
        return None
    return float(ci[1] - ci[0])


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2126 = load(IN_2126)
    p2127 = load(IN_2127)

    rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    pre_ready = p2126.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_WITH_TRACE"
    iid_ready = p2127.get("result_kind") == "PASS_STRICT_CMP2_BOOTSTRAP_BACKEND_EVIDENCE_STRESSTEST_WITH_TRACE"

    rng = np.random.default_rng(RNG_SEED)
    n = len(rows)
    coverage_ci_block = [None, None]
    rank_stability = {}

    if pre_ready and n > 0:
        by_cluster: dict[float, list[dict[str, Any]]] = {}
        for row in rows:
            by_cluster.setdefault(float(row.get("backend_s", 0.0)), []).append(row)
        cluster_keys = sorted(by_cluster.keys())
        n_clusters = len(cluster_keys)

        boot_cov = []
        boot_top_counts = {"M1_nn": [], "M2_monotone": [], "M3_nn_tiebreak": [], "M4_monotone_penalized": []}

        for _ in range(N_BOOT):
            sampled_clusters = [cluster_keys[i] for i in rng.integers(0, n_clusters, size=n_clusters)]
            sample_rows = [r for c in sampled_clusters for r in by_cluster[c]]
            m = len(sample_rows)
            coverage = [1.0 if r.get("posterior_predictive_covered", False) else 0.0 for r in sample_rows]
            boot_cov.append(float(np.mean(coverage)) if m else 0.0)

            top = [rank_of_max((r.get("posterior_weights_backend_evidence", {}) or {})) for r in sample_rows]
            for k in boot_top_counts:
                boot_top_counts[k].append(top.count(k) / m if m else 0.0)

        coverage_ci_block = [float(np.quantile(boot_cov, 0.025)), float(np.quantile(boot_cov, 0.975))]
        rank_stability = {
            k: {
                "mean_frequency": float(np.mean(v)),
                "ci95": [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))],
            }
            for k, v in boot_top_counts.items()
        }
    else:
        n_clusters = 0

    iid_obj = p2127.get("bootstrap_backend_evidence_stresstest", {}) or {}
    iid_cov_ci = iid_obj.get("coverage_rate_ci95", [None, None])
    iid_rank = iid_obj.get("posterior_top_model_rank_stability", {}) or {}

    coverage_inflation = None
    w_iid = ci_width(iid_cov_ci)
    w_blk = ci_width(coverage_ci_block)
    if w_iid is not None and w_blk is not None and w_iid > 0:
        coverage_inflation = float(w_blk / w_iid)

    rank_ci_width_inflation = {}
    for model, blk in rank_stability.items():
        iid_ci = (iid_rank.get(model, {}) or {}).get("ci95", [None, None])
        blk_ci = blk.get("ci95", [None, None])
        w_i = ci_width(iid_ci)
        w_b = ci_width(blk_ci)
        rank_ci_width_inflation[model] = float(w_b / w_i) if (w_i is not None and w_i > 0 and w_b is not None) else None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2128",
        "stage_id": "S1078",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_BLOCK_BOOTSTRAP_DEPENDENCE_AWARE_UNCERTAINTY_INFLATION_WITH_TRACE" if pre_ready and iid_ready and n > 0 else "OPEN_STRICT_CMP2_BLOCK_BOOTSTRAP_DEPENDENCE_AWARE_UNCERTAINTY_INFLATION_BLOCKED",
        "depends_on": {
            "p2126_present": p2126.get("_missing") is None,
            "p2127_present": p2127.get("_missing") is None,
            "preconditions_ready": pre_ready and iid_ready,
        },
        "dependence_aware_uncertainty_inflation_object": {
            "n_rows": n,
            "n_clusters": n_clusters,
            "cluster_key": "backend_s",
            "n_bootstrap": N_BOOT,
            "rng_seed": RNG_SEED,
            "iid_coverage_rate_ci95": iid_cov_ci,
            "block_coverage_rate_ci95": coverage_ci_block,
            "coverage_ci_width_inflation_block_over_iid": coverage_inflation,
            "iid_posterior_top_model_rank_stability": iid_rank,
            "block_posterior_top_model_rank_stability": rank_stability,
            "rank_ci_width_inflation_block_over_iid": rank_ci_width_inflation,
            "scope_limit": "dependence-aware uncertainty quantification only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2129_candidate",
            "goal": "run sensitivity sweep over block definitions and report robustness of uncertainty inflation factors",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready and iid_ready,
            "block_bootstrap_executed": pre_ready and iid_ready and n > 0,
            "inflation_object_exported": pre_ready and iid_ready and n > 0,
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
            "# P2128 S1078: strict CMP2 block-bootstrap dependence-aware uncertainty inflation",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Bootstrap samples: `{N_BOOT}`",
            f"- Block coverage CI95: `{coverage_ci_block}`",
            f"- Coverage CI inflation (block/iid): `{coverage_inflation}`",
            "",
            "This stage adds dependence-aware block bootstrap by backend_s clusters and compares uncertainty widths against iid bootstrap.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
