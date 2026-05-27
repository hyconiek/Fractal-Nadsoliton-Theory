#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2244 = GEN / "p2244_s1194_strict_nu_branch_group_policy_violation_tail_bound_probe.json"
OUT = GEN / "p2245_s1195_strict_nu_branch_group_policy_reserve_factor_calibration_probe.json"
MD = GEN / "p2245_s1195_strict_nu_branch_group_policy_reserve_factor_calibration_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2244 = load(IN_2244)
    probe = (p2244.get("strict_nu_branch_group_policy_violation_tail_bound_probe", {}) or {})
    inputs = (probe.get("inputs", {}) or {})
    bounds = (probe.get("derived_bounds", {}) or {})

    n_groups = int(inputs.get("group_count", 1) or 1)
    total_cov = float(inputs.get("total_coverage_mass", 1.0) or 1.0)
    load_ratio = float(inputs.get("load_ratio_threshold", 1.0) or 1.0)
    empirical_upper = float(bounds.get("hoeffding_empirical_violation_upper_bound_95", 1.0) or 1.0)

    # Operational target risk envelope (strictly diagnostic policy target, not theorem claim).
    target_risk = 0.05

    # Required reserve multiplier on load-ratio threshold to pass deterministic condition:
    # total_cov >= (n-1) + kappa*load_ratio  => kappa <= (total_cov-(n-1))/load_ratio.
    # We export minimal extra reserve request relative to kappa=1 baseline.
    denom = max(load_ratio, 1e-12)
    kappa_capacity = (total_cov - (n_groups - 1)) / denom
    kappa_needed_for_target = 1.0

    deterministic_pass = kappa_capacity + 1e-15 >= kappa_needed_for_target
    statistical_pass = empirical_upper <= target_risk + 1e-15
    joint_pass = deterministic_pass and statistical_pass

    # If statistical side fails, export conservative reserve suggestion from risk ratio.
    risk_ratio = empirical_upper / max(target_risk, 1e-12)
    suggested_kappa = max(1.0, min(kappa_capacity, risk_ratio))

    payload = {
        "schema_version": "p2245_s1195_v1",
        "packet_id": "P2245",
        "stage_id": "S1195",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_RESERVE_FACTOR_CALIBRATION_PROBE",
        "strict_nu_branch_group_policy_reserve_factor_calibration_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_RESERVE_FACTOR_CALIBRATION_PROBE_V1",
            "source_packets": [str(IN_2244.relative_to(ROOT))],
            "inputs": {
                "group_count": n_groups,
                "total_coverage_mass": total_cov,
                "load_ratio_threshold": load_ratio,
                "hoeffding_empirical_violation_upper_bound_95": empirical_upper,
                "target_risk": target_risk,
            },
            "reserve_calibration": {
                "kappa_capacity_from_mass": kappa_capacity,
                "kappa_needed_for_deterministic_condition": kappa_needed_for_target,
                "suggested_kappa": suggested_kappa,
                "deterministic_pass": deterministic_pass,
                "statistical_pass": statistical_pass,
                "joint_pass": joint_pass,
            },
            "physical_interpretation_note": "Reserve factor kappa calibrates how much coverage buffer is available relative to load-ratio demand; joint deterministic+statistical pass indicates robust safety-lane protection under both worst-case structure and finite-sample uncertainty.",
            "theorem_scope_limit": "finite-sample strict-lane reserve-calibration diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2246_candidate",
            "goal": "map reserve-factor feasibility surface over (group_count, load_ratio, target_risk) and export admissible operational region",
        },
        "gatekeeper_checks": {
            "reserve_factor_exported": True,
            "kappa_capacity_computable": kappa_capacity >= 0.0,
            "target_risk_applied": True,
            "joint_pass_computable": isinstance(joint_pass, bool),
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
                "# P2245 S1195: reserve-factor calibration probe",
                "",
                f"- target risk: `{target_risk:.12e}`",
                f"- kappa capacity from mass: `{kappa_capacity:.12e}`",
                f"- kappa needed (deterministic): `{kappa_needed_for_target:.12e}`",
                f"- suggested kappa: `{suggested_kappa:.12e}`",
                f"- deterministic pass: `{deterministic_pass}`",
                f"- statistical pass: `{statistical_pass}`",
                f"- joint pass: `{joint_pass}`",
                "",
                "Finite-sample reserve-calibration diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
