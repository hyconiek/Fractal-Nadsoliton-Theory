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
IN_2134 = GEN / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.json"
EXT = GEN / "cmp2_backend_rows_extension_v1.json"
OUT = GEN / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json"
MD = GEN / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.md"

SCHEMA_VERSION = "p2135_s1085_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2135
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


def block_ids_quantile(rows: list[dict[str, Any]], q: int) -> list[int]:
    s = np.array([float(r.get("backend_s", 0.0)) for r in rows], dtype=float)
    edges = np.quantile(s, np.linspace(0.0, 1.0, q + 1))
    out = []
    for v in s:
        idx = int(np.searchsorted(edges, v, side="right") - 1)
        out.append(max(0, min(idx, q - 1)))
    return out


def block_ids_adj(rows: list[dict[str, Any]], k: int) -> list[int]:
    return [i // k for i in range(len(rows))]


def block_ids_unique(rows: list[dict[str, Any]]) -> list[int]:
    return [int(float(r.get("backend_s", 0.0)) * 1_000_000) for r in rows]


def block_bootstrap(rows: list[dict[str, Any]], block_ids: list[int], rng: np.random.Generator) -> dict[str, Any]:
    by = {}
    for r, bid in zip(rows, block_ids):
        by.setdefault(int(bid), []).append(r)
    keys = sorted(by.keys())
    boot_cov = []
    boot_top = {m: [] for m in MODELS}
    for _ in range(N_BOOT):
        sampled = [keys[i] for i in rng.integers(0, len(keys), size=len(keys))]
        sample = [rr for k in sampled for rr in by[k]]
        n = len(sample)
        cov = [1.0 if r.get("posterior_predictive_covered", False) else 0.0 for r in sample]
        boot_cov.append(float(np.mean(cov)) if n else 0.0)
        tops = [top_model(r) for r in sample]
        for m in MODELS:
            boot_top[m].append(tops.count(m) / n if n else 0.0)
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
    p2134 = load(IN_2134)
    ext = load(EXT) if EXT.exists() else {"rows": []}

    base_rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    new_rows = ext.get("rows") or []
    merge_ready = bool(((p2133.get("real_extension_merge_contract", {}) or {}).get("merge_ready", False)))
    p2134_pass = p2134.get("result_kind") == "PASS_STRICT_CMP2_NONSYNTHETIC_RERUN_COMPARISON_AUDIT_WITH_TRACE"

    required = {"cmp_bin_index", "backend_s", "posterior_weights_backend_evidence", "posterior_predictive_covered"}
    valid_new_rows = [r for r in new_rows if isinstance(r, dict) and required.issubset(r.keys())]

    variants = {}
    spread = None
    if merge_ready and p2134_pass and base_rows and new_rows and len(valid_new_rows) == len(new_rows):
        merged = base_rows + valid_new_rows
        rng = np.random.default_rng(RNG_SEED)
        iid = iid_bootstrap(merged, rng)
        iw = ci_width(iid["coverage_rate_ci95"])

        specs = {
            "backend_s_unique": block_ids_unique(merged),
            "backend_s_quantile_q2": block_ids_quantile(merged, 2),
            "backend_s_quantile_q4": block_ids_quantile(merged, 4),
            "bin_adjacency_k2": block_ids_adj(merged, 2),
        }

        infls = []
        for name, bids in specs.items():
            blk = block_bootstrap(merged, bids, rng)
            bw = ci_width(blk["coverage_rate_ci95"])
            infl = float(bw / iw) if (bw is not None and iw is not None and iw > 0) else None
            if infl is not None:
                infls.append(infl)
            variants[name] = {
                "coverage_rate_ci95": blk["coverage_rate_ci95"],
                "coverage_ci_width_inflation_block_over_iid": infl,
                "posterior_top_model_rank_stability": blk["posterior_top_model_rank_stability"],
            }

        if len(infls) >= 2:
            spread = float(max(infls) - min(infls))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2135",
        "stage_id": "S1085",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_MERGED_REAL_BLOCK_VARIANT_STABILITY_AUDIT_WITH_TRACE" if variants else "OPEN_STRICT_CMP2_MERGED_REAL_BLOCK_VARIANT_STABILITY_AUDIT_BLOCKED",
        "depends_on": {
            "p2133_merge_ready": merge_ready,
            "p2134_pass": p2134_pass,
            "extension_present": EXT.exists(),
        },
        "merged_real_block_variant_stability_audit": {
            "n_bootstrap": N_BOOT,
            "rng_seed": RNG_SEED,
            "variants": variants,
            "inflation_spread_across_variants": spread,
            "scope_limit": "merged-real block-variant stability audit only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2136_candidate",
            "goal": "if P2135 passes with finite spread, export final comparison memo and prioritize real backend growth for tighter CI stability",
        },
        "gatekeeper_checks": {
            "preconditions_ready": merge_ready and p2134_pass,
            "variant_audit_executed": bool(variants),
            "finite_spread_available": spread is not None,
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
            "# P2135 S1085: strict CMP2 merged-real block-variant stability audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Variants computed: `{len(variants)}`",
            "",
            "This stage extends P2134 to full P2129 block variants on merged real dataset when prerequisites are satisfied.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
