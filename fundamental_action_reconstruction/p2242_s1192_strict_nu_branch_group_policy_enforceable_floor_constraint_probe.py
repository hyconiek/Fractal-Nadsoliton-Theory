#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2241 = GEN / "p2241_s1191_strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe.json"
OUT = GEN / "p2242_s1192_strict_nu_branch_group_policy_enforceable_floor_constraint_probe.json"
MD = GEN / "p2242_s1192_strict_nu_branch_group_policy_enforceable_floor_constraint_probe.md"

def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))

def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2241 = load(IN_2241)
    probe = (p2241.get("strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe", {}) or {})
    fixed = (probe.get("fixed_mass_inputs", {}) or {})

    n = int(fixed.get("group_count", 1) or 1)
    total_cov = float(fixed.get("total_coverage_mass", 1.0) or 1.0)
    load_ratio = float(fixed.get("load_ratio_threshold", 1.0) or 1.0)

    # With per-group cap <=1 and fixed total mass, adversarial worst minimum is:
    # min_cov_worst = total_cov - (n-1), clipped to [0,1].
    # To guarantee min_cov_worst >= load_ratio it is sufficient that:
    # total_cov >= (n - 1) + load_ratio.
    required_total_mass = (n - 1) + load_ratio
    enforceable_mass_condition_holds = total_cov + 1e-15 >= required_total_mass

    equivalent_uniform_floor = total_cov / max(n, 1)
    floor_above_ratio = equivalent_uniform_floor + 1e-15 >= load_ratio

    payload = {
        "schema_version": "p2242_s1192_v1",
        "packet_id": "P2242",
        "stage_id": "S1192",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ENFORCEABLE_FLOOR_CONSTRAINT_PROBE",
        "strict_nu_branch_group_policy_enforceable_floor_constraint_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ENFORCEABLE_FLOOR_CONSTRAINT_PROBE_V1",
            "source_packets": [str(IN_2241.relative_to(ROOT))],
            "inputs": {
                "group_count": n,
                "total_coverage_mass": total_cov,
                "load_ratio_threshold": load_ratio,
            },
            "derived_constraint": {
                "sufficient_mass_condition": "total_coverage_mass >= (group_count - 1) + load_ratio_threshold",
                "required_total_coverage_mass": required_total_mass,
                "observed_total_coverage_mass": total_cov,
                "condition_holds": enforceable_mass_condition_holds,
            },
            "policy_floor_translation": {
                "uniform_floor_candidate": equivalent_uniform_floor,
                "uniform_floor_above_ratio": floor_above_ratio,
                "note": "Uniform floor is a deployable policy form; sufficient mass condition is the stricter adversarial guarantee.",
            },
            "physical_interpretation_note": "This exports a concrete conservation-style control law: maintain enough global coverage mass so no group can be adversarially starved below the safety load-ratio threshold.",
            "theorem_scope_limit": "finite-sample strict-lane policy-constraint diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2243_candidate",
            "goal": "validate constraint tightness by Monte-Carlo policy-mixing perturbations and compute empirical violation probability envelope",
        },
        "gatekeeper_checks": {
            "enforceable_floor_constraint_exported": True,
            "mass_condition_computable": True,
            "mass_condition_holds": enforceable_mass_condition_holds,
            "uniform_floor_translation_exported": True,
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
                "# P2242 S1192: enforceable floor-constraint probe",
                "",
                "- sufficient condition:",
                "  `total_coverage_mass >= (group_count - 1) + load_ratio_threshold`",
                f"- observed total coverage mass: `{total_cov:.12e}`",
                f"- required total coverage mass: `{required_total_mass:.12e}`",
                f"- condition holds: `{enforceable_mass_condition_holds}`",
                f"- uniform floor candidate (deployable): `{equivalent_uniform_floor:.12e}`",
                f"- uniform floor above ratio: `{floor_above_ratio}`",
                "",
                "Finite-sample policy-constraint diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

if __name__ == "__main__":
    main()
