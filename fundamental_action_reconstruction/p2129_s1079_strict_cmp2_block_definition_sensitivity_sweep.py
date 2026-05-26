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
OUT = GEN / "p2129_s1079_strict_cmp2_block_definition_sensitivity_sweep.json"
MD = GEN / "p2129_s1079_strict_cmp2_block_definition_sensitivity_sweep.md"

SCHEMA_VERSION = "p2129_s1079_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2129
N_BOOT = 400
MODELS = ["M1_nn", "M2_monotone", "M3_nn_tiebreak", "M4_monotone_penalized"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def rank_of_max(weights: dict[str, float]) -> str:
    return max(weights.items(), key=lambda kv: kv[1])[0]


def ci_width(ci: list[float | None]) -> float | None:
    lo, hi = ci
    if lo is None or hi is None:
        return None
    return float(hi - lo)


def quantized_block_ids(rows: list[dict[str, Any]], n_quantiles: int) -> list[int]:
    s = np.array([float(r.get("backend_s", 0.0)) for r in rows], dtype=float)
    if len(s) == 0:
        return []
    edges = np.quantile(s, np.linspace(0.0, 1.0, n_quantiles + 1))
    out = []
    for v in s:
        idx = int(np.searchsorted(edges, v, side="right") - 1)
        idx = max(0, min(idx, n_quantiles - 1))
        out.append(idx)
    return out


def adjacency_block_ids(rows: list[dict[str, Any]], block_size: int) -> list[int]:
    return [i // block_size for i in range(len(rows))]


def block_bootstrap(rows: list[dict[str, Any]], block_ids: list[int], rng: np.random.Generator) -> dict[str, Any]:
    by_block: dict[int, list[dict[str, Any]]] = {}
    for row, bid in zip(rows, block_ids):
        by_block.setdefault(int(bid), []).append(row)
    blocks = sorted(by_block.keys())
    n_blocks = len(blocks)

    boot_cov = []
    boot_top = {m: [] for m in MODELS}
    for _ in range(N_BOOT):
        sampled = [blocks[i] for i in rng.integers(0, n_blocks, size=n_blocks)]
        sample_rows = [r for b in sampled for r in by_block[b]]
        m = len(sample_rows)
        cov = [1.0 if r.get("posterior_predictive_covered", False) else 0.0 for r in sample_rows]
        boot_cov.append(float(np.mean(cov)) if m else 0.0)
        tops = [rank_of_max((r.get("posterior_weights_backend_evidence", {}) or {})) for r in sample_rows]
        for model in MODELS:
            boot_top[model].append(tops.count(model) / m if m else 0.0)

    rank_stability = {
        model: {
            "mean_frequency": float(np.mean(vals)),
            "ci95": [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))],
        }
        for model, vals in boot_top.items()
    }
    return {
        "n_blocks": n_blocks,
        "coverage_rate_ci95": [float(np.quantile(boot_cov, 0.025)), float(np.quantile(boot_cov, 0.975))],
        "posterior_top_model_rank_stability": rank_stability,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2126 = load(IN_2126)
    p2127 = load(IN_2127)
    rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])

    pre_ready = p2126.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_WITH_TRACE"
    iid_ready = p2127.get("result_kind") == "PASS_STRICT_CMP2_BOOTSTRAP_BACKEND_EVIDENCE_STRESSTEST_WITH_TRACE"

    iid_obj = p2127.get("bootstrap_backend_evidence_stresstest", {}) or {}
    iid_cov_ci = iid_obj.get("coverage_rate_ci95", [None, None])
    iid_rank = iid_obj.get("posterior_top_model_rank_stability", {}) or {}

    variants = {}
    if pre_ready and iid_ready and rows:
        rng = np.random.default_rng(RNG_SEED)
        specs = {
            "backend_s_unique": [int(float(r.get("backend_s", 0.0)) * 1_000_000) for r in rows],
            "backend_s_quantile_q2": quantized_block_ids(rows, 2),
            "backend_s_quantile_q4": quantized_block_ids(rows, 4),
            "bin_adjacency_k2": adjacency_block_ids(rows, 2),
        }

        for name, block_ids in specs.items():
            res = block_bootstrap(rows, block_ids, rng)
            cov_w = ci_width(res["coverage_rate_ci95"])
            iid_w = ci_width(iid_cov_ci)
            cov_infl = float(cov_w / iid_w) if (cov_w is not None and iid_w is not None and iid_w > 0) else None

            rank_infl = {}
            for model in MODELS:
                rw = ci_width((res["posterior_top_model_rank_stability"].get(model, {}) or {}).get("ci95", [None, None]))
                iw = ci_width((iid_rank.get(model, {}) or {}).get("ci95", [None, None]))
                rank_infl[model] = float(rw / iw) if (rw is not None and iw is not None and iw > 0) else None

            variants[name] = {
                "block_definition": name,
                "n_blocks": res["n_blocks"],
                "coverage_rate_ci95": res["coverage_rate_ci95"],
                "coverage_ci_width_inflation_block_over_iid": cov_infl,
                "posterior_top_model_rank_stability": res["posterior_top_model_rank_stability"],
                "rank_ci_width_inflation_block_over_iid": rank_infl,
            }

    finite_cov_infl = [v["coverage_ci_width_inflation_block_over_iid"] for v in variants.values() if v.get("coverage_ci_width_inflation_block_over_iid") is not None]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2129",
        "stage_id": "S1079",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_BLOCK_DEFINITION_SENSITIVITY_SWEEP_WITH_TRACE" if (pre_ready and iid_ready and len(rows) > 0 and variants) else "OPEN_STRICT_CMP2_BLOCK_DEFINITION_SENSITIVITY_SWEEP_BLOCKED",
        "depends_on": {
            "p2126_present": p2126.get("_missing") is None,
            "p2127_present": p2127.get("_missing") is None,
            "preconditions_ready": pre_ready and iid_ready,
        },
        "block_definition_sensitivity_sweep": {
            "n_rows": len(rows),
            "n_bootstrap": N_BOOT,
            "rng_seed": RNG_SEED,
            "iid_reference": {
                "coverage_rate_ci95": iid_cov_ci,
                "posterior_top_model_rank_stability": iid_rank,
            },
            "variants": variants,
            "coverage_inflation_stability_summary": {
                "n_finite": len(finite_cov_infl),
                "min": float(min(finite_cov_infl)) if finite_cov_infl else None,
                "max": float(max(finite_cov_infl)) if finite_cov_infl else None,
                "spread": float(max(finite_cov_infl) - min(finite_cov_infl)) if len(finite_cov_infl) >= 2 else None,
            },
            "scope_limit": "sensitivity sweep on block definitions only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2130_candidate",
            "goal": "increase sample support (more bins/backends) before interpreting CI inflation as stable structural evidence",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready and iid_ready,
            "sweep_executed": bool(variants),
            "inflation_stability_compared": bool(variants),
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
            "# P2129 S1079: strict CMP2 block-definition sensitivity sweep",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Sweep variants: `{len(variants)}`",
            f"- Coverage inflation spread: `{payload['block_definition_sensitivity_sweep']['coverage_inflation_stability_summary']['spread']}`",
            "",
            "This stage compares uncertainty inflation stability across multiple dependence-aware block definitions.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
