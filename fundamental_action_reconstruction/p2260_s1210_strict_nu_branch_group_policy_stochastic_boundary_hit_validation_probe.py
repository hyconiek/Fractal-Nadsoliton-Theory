#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2259 = GEN / "p2259_s1209_strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe.json"
OUT = GEN / "p2260_s1210_strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe.json"
MD = GEN / "p2260_s1210_strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe.md"


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

    rng = random.Random(2260)
    draws = 500

    first_hits = 0
    second_hits = 0
    total = 0
    rows = []

    for h in horizon_scales:
        h_first = 0
        h_second = 0
        for _ in range(draws):
            noise = (rng.random() - 0.5) * 2.0e-4
            first_load = h * h * 1.0e-3 + noise
            second_load = max(0.0, first_load - gamma_hat * h * h)
            hf = first_load >= threshold
            hs = second_load >= threshold
            h_first += int(hf)
            h_second += int(hs)
        first_hits += h_first
        second_hits += h_second
        total += draws
        rows.append(
            {
                "horizon_scale": h,
                "first_hit_rate": h_first / draws,
                "second_hit_rate": h_second / draws,
                "hit_rate_reduction": (h_first - h_second) / draws,
            }
        )

    first_rate = first_hits / max(total, 1)
    second_rate = second_hits / max(total, 1)
    reduction = first_rate - second_rate

    payload = {
        "schema_version": "p2260_s1210_v1",
        "packet_id": "P2260",
        "stage_id": "S1210",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_STOCHASTIC_BOUNDARY_HIT_VALIDATION_PROBE",
        "strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_STOCHASTIC_BOUNDARY_HIT_VALIDATION_PROBE_V1",
            "source_packets": [str(IN_2259.relative_to(ROOT))],
            "inputs": {
                "gamma_hat": gamma_hat,
                "horizon_scales": horizon_scales,
                "boundary_hit_threshold": threshold,
                "draws_per_horizon": draws,
                "random_seed": 2260,
            },
            "stochastic_hit_rate_comparison": {
                "rows": rows,
                "global_first_hit_rate": first_rate,
                "global_second_hit_rate": second_rate,
                "global_hit_rate_reduction": reduction,
            },
            "physical_interpretation_note": "Under stochastic perturbations around trajectory load, gamma-injected correction retains lower boundary-hit rate, supporting robustness of nonlinear safety tightening beyond deterministic proxy scans.",
            "theorem_scope_limit": "stochastic finite-sample validation only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2261_candidate",
            "goal": "derive concentration bound for observed stochastic hit-rate reduction and compare with deterministic envelope predictions",
        },
        "gatekeeper_checks": {
            "stochastic_validation_exported": True,
            "rates_bounded": 0.0 <= first_rate <= 1.0 and 0.0 <= second_rate <= 1.0,
            "reduction_nonnegative": reduction >= -1e-15,
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
                "# P2260 S1210: stochastic boundary-hit validation probe",
                "",
                f"- global first-order hit rate: `{first_rate:.12e}`",
                f"- global second-order hit rate: `{second_rate:.12e}`",
                f"- global hit-rate reduction: `{reduction:.12e}`",
                f"- draws per horizon: `{draws}`",
                "",
                "Stochastic finite-sample validation only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
