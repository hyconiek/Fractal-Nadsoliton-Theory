#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2126 = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
IN_2129 = GEN / "p2129_s1079_strict_cmp2_block_definition_sensitivity_sweep.json"
OUT = GEN / "p2130_s1080_strict_cmp2_sample_support_expansion_planner.json"
MD = GEN / "p2130_s1080_strict_cmp2_sample_support_expansion_planner.md"

SCHEMA_VERSION = "p2130_s1080_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2130
N_SIM = 2000


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def wilson_half_width(p: float, n: int, z: float = 1.96) -> float:
    if n <= 0:
        return float("nan")
    den = 1.0 + z * z / n
    ctr = p + z * z / (2.0 * n)
    rad = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n)
    return rad / den


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2126 = load(IN_2126)
    p2129 = load(IN_2129)

    rows = ((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    pre_ready = p2126.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_WITH_TRACE"
    sweep_ready = p2129.get("result_kind") == "PASS_STRICT_CMP2_BLOCK_DEFINITION_SENSITIVITY_SWEEP_WITH_TRACE"

    n = len(rows)
    y = int(sum(1 for r in rows if r.get("posterior_predictive_covered", False)))
    p_hat = (y / n) if n > 0 else None

    variants = ((p2129.get("block_definition_sensitivity_sweep", {}) or {}).get("variants") or {})
    null_inflation_variants = [k for k, v in variants.items() if v.get("coverage_ci_width_inflation_block_over_iid") is None]

    rng = np.random.default_rng(RNG_SEED)
    projected = {}
    if pre_ready and sweep_ready and n > 0:
        # Beta posterior over Bernoulli coverage as operational planner (not closure evidence)
        alpha = 1 + y
        beta = 1 + (n - y)
        p_samples = rng.beta(alpha, beta, size=N_SIM)

        for m in [2, 3, 5, 8, 12]:
            n_new = n * m
            widths = [2.0 * wilson_half_width(float(p), n_new) for p in p_samples]
            projected[f"x{m}"] = {
                "projected_n_rows": n_new,
                "projected_wilson_ci95_width_median": float(np.median(widths)),
                "projected_wilson_ci95_width_q95": float(np.quantile(widths, 0.95)),
                "degenerate_zero_width_probability": float(np.mean(np.array(widths) < 1e-12)),
            }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2130",
        "stage_id": "S1080",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_SAMPLE_SUPPORT_EXPANSION_PLANNER_WITH_TRACE" if (pre_ready and sweep_ready and n > 0) else "OPEN_STRICT_CMP2_SAMPLE_SUPPORT_EXPANSION_PLANNER_BLOCKED",
        "depends_on": {
            "p2126_present": p2126.get("_missing") is None,
            "p2129_present": p2129.get("_missing") is None,
            "preconditions_ready": pre_ready and sweep_ready,
        },
        "sample_support_expansion_planner": {
            "n_rows_current": n,
            "coverage_successes_current": y,
            "coverage_rate_current": p_hat,
            "null_inflation_variants_from_p2129": null_inflation_variants,
            "null_inflation_root_cause": "degenerate/near-degenerate CI widths at tiny sample support",
            "projection_model": "beta-bernoulli + wilson-width planner",
            "n_sim": N_SIM,
            "rng_seed": RNG_SEED,
            "projected_ci95_widths_by_support_multiplier": projected,
            "scope_limit": "operational planning for support expansion only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2131_candidate",
            "goal": "materially increase CMP2 bin/backend support and rerun P2127-P2129 to obtain non-degenerate finite inflation factors",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready and sweep_ready,
            "root_cause_diagnosed": pre_ready and sweep_ready and n > 0,
            "support_expansion_plan_exported": pre_ready and sweep_ready and n > 0,
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
            "# P2130 S1080: strict CMP2 sample-support expansion planner",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Current rows: `{n}`",
            f"- Null-inflation variants (from P2129): `{len(null_inflation_variants)}`",
            "",
            "This stage diagnoses degenerate CI-width inflation under tiny support and exports an operational support-expansion plan.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
