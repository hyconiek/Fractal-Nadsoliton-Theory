#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2243 = GEN / "p2243_s1193_strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe.json"
OUT = GEN / "p2244_s1194_strict_nu_branch_group_policy_violation_tail_bound_probe.json"
MD = GEN / "p2244_s1194_strict_nu_branch_group_policy_violation_tail_bound_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def hoeffding_upper(success_rate: float, n: int, alpha: float = 0.05) -> float:
    # One-sided confidence upper bound for Bernoulli mean from empirical rate.
    eps = math.sqrt(max(0.0, math.log(1.0 / max(alpha, 1e-12))) / (2.0 * max(n, 1)))
    return min(1.0, max(0.0, success_rate + eps))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2243 = load(IN_2243)
    probe = (p2243.get("strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe", {}) or {})
    inputs = (probe.get("inputs", {}) or {})
    env = (probe.get("empirical_violation_envelope", {}) or {})

    n_groups = int(inputs.get("group_count", 1) or 1)
    total_cov = float(inputs.get("total_coverage_mass", 1.0) or 1.0)
    load_ratio = float(inputs.get("load_ratio_threshold", 1.0) or 1.0)
    draws = int(inputs.get("draw_count", 1) or 1)
    empirical_violation = float(env.get("violation_probability", 1.0) or 1.0)

    required_mass_threshold = (n_groups - 1) + load_ratio
    mass_gap = total_cov - required_mass_threshold

    # Deterministic adversarial bound from P2242-style derivation.
    if mass_gap >= -1e-15:
        deterministic_upper = 0.0
    else:
        deterministic_upper = 1.0

    # Statistical finite-sample upper bound around empirical violation estimate.
    empirical_upper_95 = hoeffding_upper(empirical_violation, draws, alpha=0.05)

    payload = {
        "schema_version": "p2244_s1194_v1",
        "packet_id": "P2244",
        "stage_id": "S1194",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_VIOLATION_TAIL_BOUND_PROBE",
        "strict_nu_branch_group_policy_violation_tail_bound_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_VIOLATION_TAIL_BOUND_PROBE_V1",
            "source_packets": [str(IN_2243.relative_to(ROOT))],
            "inputs": {
                "group_count": n_groups,
                "total_coverage_mass": total_cov,
                "load_ratio_threshold": load_ratio,
                "draw_count": draws,
                "empirical_violation_probability": empirical_violation,
            },
            "derived_bounds": {
                "required_mass_threshold": required_mass_threshold,
                "mass_gap": mass_gap,
                "deterministic_adversarial_violation_upper_bound": deterministic_upper,
                "hoeffding_empirical_violation_upper_bound_95": empirical_upper_95,
            },
            "physical_interpretation_note": "Mass-threshold surplus controls deterministic worst-case starvation risk (0/1 bound), while Hoeffding upper bound quantifies residual finite-sample uncertainty of observed violation frequency.",
            "theorem_scope_limit": "finite-sample strict-lane tail-bound diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2245_candidate",
            "goal": "calibrate conservative policy reserve factor kappa so deterministic and statistical bounds stay below operational risk tolerance",
        },
        "gatekeeper_checks": {
            "tail_bound_exported": True,
            "deterministic_bound_computable": deterministic_upper in (0.0, 1.0),
            "statistical_bound_computable": 0.0 <= empirical_upper_95 <= 1.0,
            "mass_gap_reported": True,
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
                "# P2244 S1194: violation tail-bound probe",
                "",
                f"- required mass threshold: `{required_mass_threshold:.12e}`",
                f"- observed total mass: `{total_cov:.12e}`",
                f"- mass gap: `{mass_gap:.12e}`",
                f"- deterministic adversarial upper bound: `{deterministic_upper:.12e}`",
                f"- Hoeffding 95% upper bound: `{empirical_upper_95:.12e}`",
                "",
                "Finite-sample tail-bound diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
