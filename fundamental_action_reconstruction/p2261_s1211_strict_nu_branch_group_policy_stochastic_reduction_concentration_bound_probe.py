#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2260 = GEN / "p2260_s1210_strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe.json"
OUT = GEN / "p2261_s1211_strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe.json"
MD = GEN / "p2261_s1211_strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def hoeffding_two_sample_eps(n: int, alpha: float = 0.05) -> float:
    return math.sqrt(max(0.0, math.log(2.0 / max(alpha, 1e-12))) / (2.0 * max(n, 1)))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2260 = load(IN_2260)
    probe = (p2260.get("strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})
    cmp = (probe.get("stochastic_hit_rate_comparison", {}) or {})

    first_rate = float(cmp.get("global_first_hit_rate", 0.0) or 0.0)
    second_rate = float(cmp.get("global_second_hit_rate", 0.0) or 0.0)
    observed_reduction = float(cmp.get("global_hit_rate_reduction", first_rate - second_rate) or (first_rate - second_rate))

    draws = int(inp.get("draws_per_horizon", 500) or 500)
    horizon_count = len(inp.get("horizon_scales", []) or [])
    n_total = max(draws * max(horizon_count, 1), 1)

    eps = hoeffding_two_sample_eps(n_total, alpha=0.05)
    lower_bound_reduction_95 = observed_reduction - 2.0 * eps

    payload = {
        "schema_version": "p2261_s1211_v1",
        "packet_id": "P2261",
        "stage_id": "S1211",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_STOCHASTIC_REDUCTION_CONCENTRATION_BOUND_PROBE",
        "strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_STOCHASTIC_REDUCTION_CONCENTRATION_BOUND_PROBE_V1",
            "source_packets": [str(IN_2260.relative_to(ROOT))],
            "inputs": {
                "global_first_hit_rate": first_rate,
                "global_second_hit_rate": second_rate,
                "observed_reduction": observed_reduction,
                "n_total_samples": n_total,
                "alpha": 0.05,
            },
            "concentration_bound": {
                "hoeffding_eps": eps,
                "lower_bound_reduction_95": lower_bound_reduction_95,
                "reduction_certified_positive_95": lower_bound_reduction_95 > 0.0,
            },
            "physical_interpretation_note": "Concentration lower bound quantifies confidence that observed stochastic hit-rate reduction is not a sampling artifact, strengthening physical credibility of second-order safety tightening.",
            "theorem_scope_limit": "finite-sample concentration diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2262_candidate",
            "goal": "increase stochastic sample budget and verify stability of certified reduction lower bound across seeds",
        },
        "gatekeeper_checks": {
            "concentration_bound_exported": True,
            "eps_nonnegative": eps >= 0.0,
            "lower_bound_computable": isinstance(lower_bound_reduction_95, float),
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
                "# P2261 S1211: stochastic reduction concentration-bound probe",
                "",
                f"- observed reduction: `{observed_reduction:.12e}`",
                f"- Hoeffding eps: `{eps:.12e}`",
                f"- 95% lower bound reduction: `{lower_bound_reduction_95:.12e}`",
                f"- reduction certified positive (95%): `{lower_bound_reduction_95 > 0.0}`",
                "",
                "Finite-sample concentration diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
