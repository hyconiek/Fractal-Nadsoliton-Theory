#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2126 = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
IN_2133 = GEN / "p2133_s1083_strict_cmp2_real_extension_merge_contract.json"
EXT = GEN / "cmp2_backend_rows_extension_v1.json"
OUT = GEN / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.json"
MD = GEN / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.md"

SCHEMA_VERSION = "p2134_s1084_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2134
N_BOOT = 400
MODELS = ["M1_nn", "M2_monotone", "M3_nn_tiebreak", "M4_monotone_penalized"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def ci_width(ci: list[float | None]) -> float | None:
    if ci[0] is None or ci[1] is None:
        return None
    return float(ci[1] - ci[0])


def top_model(row: dict[str, Any]) -> str:
    w = row.get("posterior_weights_backend_evidence", {}) or {}
    return max(w.items(), key=lambda kv: kv[1])[0]


def iid_bootstrap(rows: list[dict[str, Any]], rng: np.random.Generator) -> dict[str, Any]:
    n = len(rows)
    cov = np.array([1.0 if r.get("posterior_predictive_covered", False) else 0.0 for r in rows], dtype=float)
    tops = [top_model(r) for r in rows]
    boot_cov = []
    boot_top = {m: [] for m in MODELS}
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        boot_cov.append(float(np.mean(cov[idx])))
        s = [tops[i] for i in idx]
        for m in MODELS:
            boot_top[m].append(s.count(m) / n)
    return {
        "coverage_rate_ci95": [float(np.quantile(boot_cov, 0.025)), float(np.quantile(boot_cov, 0.975))],
        "posterior_top_model_rank_stability": {
            m: {"ci95": [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))]}
            for m, v in boot_top.items()
        },
    }


def block_bootstrap(rows: list[dict[str, Any]], rng: np.random.Generator) -> dict[str, Any]:
    by_cluster: dict[float, list[dict[str, Any]]] = {}
    for r in rows:
        by_cluster.setdefault(float(r.get("backend_s", 0.0)), []).append(r)
    keys = sorted(by_cluster.keys())

    boot_cov = []
    boot_top = {m: [] for m in MODELS}
    for _ in range(N_BOOT):
        sampled = [keys[i] for i in rng.integers(0, len(keys), size=len(keys))]
        sample = [rr for k in sampled for rr in by_cluster[k]]
        m = len(sample)
        cov = [1.0 if r.get("posterior_predictive_covered", False) else 0.0 for r in sample]
        boot_cov.append(float(np.mean(cov)) if m else 0.0)
        tops = [top_model(r) for r in sample]
        for model in MODELS:
            boot_top[model].append(tops.count(model) / m if m else 0.0)

    return {
        "coverage_rate_ci95": [float(np.quantile(boot_cov, 0.025)), float(np.quantile(boot_cov, 0.975))],
        "posterior_top_model_rank_stability": {
            m: {"ci95": [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))]}
            for m, v in boot_top.items()
        },
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2126 = load(IN_2126)
    p2133 = load(IN_2133)
    ext = load(EXT) if EXT.exists() else {"rows": []}

    base_rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    new_rows = ext.get("rows") or []
    merge_ready = bool(((p2133.get("real_extension_merge_contract", {}) or {}).get("merge_ready", False)))

    required = {"cmp_bin_index", "backend_s", "posterior_weights_backend_evidence", "posterior_predictive_covered"}
    valid_new_rows = [r for r in new_rows if isinstance(r, dict) and required.issubset(r.keys())]

    rerun = {}
    if merge_ready and base_rows and len(valid_new_rows) == len(new_rows) and new_rows:
        merged = base_rows + valid_new_rows
        rng = np.random.default_rng(RNG_SEED)

        iid_base = iid_bootstrap(base_rows, rng)
        blk_base = block_bootstrap(base_rows, rng)
        iid_merged = iid_bootstrap(merged, rng)
        blk_merged = block_bootstrap(merged, rng)

        infl_base = None
        bwb = ci_width(blk_base["coverage_rate_ci95"])
        iwb = ci_width(iid_base["coverage_rate_ci95"])
        if bwb is not None and iwb is not None and iwb > 0:
            infl_base = float(bwb / iwb)

        infl_merged = None
        bwm = ci_width(blk_merged["coverage_rate_ci95"])
        iwm = ci_width(iid_merged["coverage_rate_ci95"])
        if bwm is not None and iwm is not None and iwm > 0:
            infl_merged = float(bwm / iwm)

        rerun = {
            "n_rows_base": len(base_rows),
            "n_rows_new": len(valid_new_rows),
            "n_rows_merged": len(merged),
            "base_iid_coverage_rate_ci95": iid_base["coverage_rate_ci95"],
            "base_block_coverage_rate_ci95": blk_base["coverage_rate_ci95"],
            "base_coverage_ci_width_inflation_block_over_iid": infl_base,
            "merged_iid_coverage_rate_ci95": iid_merged["coverage_rate_ci95"],
            "merged_block_coverage_rate_ci95": blk_merged["coverage_rate_ci95"],
            "merged_coverage_ci_width_inflation_block_over_iid": infl_merged,
            "inflation_delta_merged_minus_base": (float(infl_merged - infl_base) if infl_base is not None and infl_merged is not None else None),
        }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2134",
        "stage_id": "S1084",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_NONSYNTHETIC_RERUN_COMPARISON_AUDIT_WITH_TRACE" if rerun else "OPEN_STRICT_CMP2_NONSYNTHETIC_RERUN_COMPARISON_AUDIT_BLOCKED",
        "depends_on": {
            "p2133_merge_ready": merge_ready,
            "extension_present": EXT.exists(),
        },
        "nonsynthetic_rerun_comparison_audit": {
            "n_bootstrap": N_BOOT,
            "rng_seed": RNG_SEED,
            "required_row_keys": sorted(required),
            "rerun_result": rerun,
            "scope_limit": "non-synthetic rerun comparison only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2135_candidate",
            "goal": "if P2134 passes with finite inflation metrics, extend comparison to full P2129 block-definition variants on merged real dataset",
        },
        "gatekeeper_checks": {
            "merge_ready_required": merge_ready,
            "nonsynthetic_rerun_executed": bool(rerun),
            "finite_inflation_metrics_available": bool(rerun and rerun.get("merged_coverage_ci_width_inflation_block_over_iid") is not None),
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
            "# P2134 S1084: strict CMP2 non-synthetic rerun comparison audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Merge ready: `{merge_ready}`",
            "",
            "This stage reruns IID/block bootstrap comparison on merged real (non-synthetic) extension when available.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
