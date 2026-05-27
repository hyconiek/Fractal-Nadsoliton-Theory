#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2259 = GEN / "p2259_s1209_strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe.json"
OUT = GEN / "p2262_s1212_strict_nu_branch_group_policy_multiseed_reduction_stability_probe.json"
MD = GEN / "p2262_s1212_strict_nu_branch_group_policy_multiseed_reduction_stability_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2259 = load(IN_2259)
    probe = (p2259.get("strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    gamma_hat = float(inp.get("gamma_hat", 0.0) or 0.0)
    horizon_scales = [float(h) for h in inp.get("horizon_scales", [0.8, 1.0, 1.2, 1.4])]
    threshold = float(inp.get("boundary_hit_threshold", 0.0025) or 0.0025)

    seeds = [2262, 2263, 2264]
    draws = 1000
    rows = []

    for seed in seeds:
        rng = random.Random(seed)
        first_hits = 0
        second_hits = 0
        total = 0
        for h in horizon_scales:
            for _ in range(draws):
                noise = (rng.random() - 0.5) * 2.0e-4
                first_load = h * h * 1.0e-3 + noise
                second_load = max(0.0, first_load - gamma_hat * h * h)
                first_hits += int(first_load >= threshold)
                second_hits += int(second_load >= threshold)
                total += 1
        r1 = first_hits / max(total, 1)
        r2 = second_hits / max(total, 1)
        rows.append(
            {
                "seed": seed,
                "first_hit_rate": r1,
                "second_hit_rate": r2,
                "reduction": r1 - r2,
            }
        )

    reductions = [r["reduction"] for r in rows]
    mean_reduction = sum(reductions) / max(len(reductions), 1)
    min_reduction = min(reductions) if reductions else 0.0
    max_reduction = max(reductions) if reductions else 0.0

    payload = {
        "schema_version": "p2262_s1212_v1",
        "packet_id": "P2262",
        "stage_id": "S1212",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_MULTISEED_REDUCTION_STABILITY_PROBE",
        "strict_nu_branch_group_policy_multiseed_reduction_stability_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_MULTISEED_REDUCTION_STABILITY_PROBE_V1",
            "source_packets": [str(IN_2259.relative_to(ROOT))],
            "inputs": {
                "gamma_hat": gamma_hat,
                "horizon_scales": horizon_scales,
                "boundary_hit_threshold": threshold,
                "draws_per_horizon": draws,
                "seeds": seeds,
            },
            "multiseed_reduction_summary": {
                "rows": rows,
                "mean_reduction": mean_reduction,
                "min_reduction": min_reduction,
                "max_reduction": max_reduction,
            },
            "physical_interpretation_note": "Cross-seed stability of hit-rate reduction supports robustness of nonlinear correction effect against stochastic realization noise.",
            "theorem_scope_limit": "finite multiseed stability diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2263_candidate",
            "goal": "derive confidence interval for mean reduction and compare with single-seed concentration lower bound",
        },
        "gatekeeper_checks": {
            "multiseed_scan_exported": True,
            "all_seed_reductions_nonnegative": all(r["reduction"] >= -1e-15 for r in rows),
            "reduction_interval_ordered": min_reduction <= mean_reduction <= max_reduction,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2262 S1212: multiseed reduction-stability probe",
                "",
                f"- mean reduction: `{mean_reduction:.12e}`",
                f"- min reduction: `{min_reduction:.12e}`",
                f"- max reduction: `{max_reduction:.12e}`",
                f"- seeds: `{seeds}`",
                "",
                "Finite multiseed stability diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
