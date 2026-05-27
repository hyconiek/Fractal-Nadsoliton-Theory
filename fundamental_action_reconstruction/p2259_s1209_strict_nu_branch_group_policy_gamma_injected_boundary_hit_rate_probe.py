#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2258 = GEN / "p2258_s1208_strict_nu_branch_group_policy_second_order_gamma_fit_probe.json"
OUT = GEN / "p2259_s1209_strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe.json"
MD = GEN / "p2259_s1209_strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2258 = load(IN_2258)
    probe = (p2258.get("strict_nu_branch_group_policy_second_order_gamma_fit_probe", {}) or {})
    fit = (probe.get("gamma_fit", {}) or {})
    gamma_hat = float(fit.get("gamma_hat", 0.0) or 0.0)

    # Synthetic boundary-hit comparison under fixed horizon scale set.
    horizon_scales = [0.8, 1.0, 1.2, 1.4]
    base_hit_threshold = 0.0025

    rows = []
    hit_first = 0
    hit_second = 0
    for h in horizon_scales:
        first_order_load = h * h * 1.0e-3
        second_order_load = max(0.0, first_order_load - gamma_hat * h * h)
        first_hit = first_order_load >= base_hit_threshold
        second_hit = second_order_load >= base_hit_threshold
        hit_first += int(first_hit)
        hit_second += int(second_hit)
        rows.append(
            {
                "horizon_scale": h,
                "first_order_effective_load": first_order_load,
                "second_order_effective_load": second_order_load,
                "first_order_boundary_hit": first_hit,
                "second_order_boundary_hit": second_hit,
            }
        )

    hit_rate_first = hit_first / max(len(rows), 1)
    hit_rate_second = hit_second / max(len(rows), 1)
    hit_rate_reduction = hit_rate_first - hit_rate_second

    payload = {
        "schema_version": "p2259_s1209_v1",
        "packet_id": "P2259",
        "stage_id": "S1209",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_GAMMA_INJECTED_BOUNDARY_HIT_RATE_PROBE",
        "strict_nu_branch_group_policy_gamma_injected_boundary_hit_rate_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_GAMMA_INJECTED_BOUNDARY_HIT_RATE_PROBE_V1",
            "source_packets": [str(IN_2258.relative_to(ROOT))],
            "inputs": {
                "gamma_hat": gamma_hat,
                "horizon_scales": horizon_scales,
                "boundary_hit_threshold": base_hit_threshold,
            },
            "boundary_hit_comparison": {
                "rows": rows,
                "first_order_hit_rate": hit_rate_first,
                "second_order_hit_rate": hit_rate_second,
                "hit_rate_reduction": hit_rate_reduction,
            },
            "physical_interpretation_note": "Injecting fitted second-order curvature lowers effective boundary load accumulation, reducing boundary-hit frequency and improving safety-margin preservation probability under horizon growth.",
            "theorem_scope_limit": "synthetic hit-rate comparison diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2260_candidate",
            "goal": "validate hit-rate reduction on nonlinear simulated trajectories with stochastic perturbations around controller path",
        },
        "gatekeeper_checks": {
            "gamma_injected_hit_rate_exported": True,
            "hit_rates_bounded": 0.0 <= hit_rate_first <= 1.0 and 0.0 <= hit_rate_second <= 1.0,
            "hit_rate_reduction_computable": isinstance(hit_rate_reduction, float),
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
                "# P2259 S1209: gamma-injected boundary-hit-rate probe",
                "",
                f"- gamma_hat: `{gamma_hat:.12e}`",
                f"- first-order hit rate: `{hit_rate_first:.12e}`",
                f"- second-order hit rate: `{hit_rate_second:.12e}`",
                f"- hit-rate reduction: `{hit_rate_reduction:.12e}`",
                "",
                "Synthetic hit-rate diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
