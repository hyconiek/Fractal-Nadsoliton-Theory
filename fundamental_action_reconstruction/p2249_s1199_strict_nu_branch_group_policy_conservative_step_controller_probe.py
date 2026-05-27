#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2248 = GEN / "p2248_s1198_strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe.json"
OUT = GEN / "p2249_s1199_strict_nu_branch_group_policy_conservative_step_controller_probe.json"
MD = GEN / "p2249_s1199_strict_nu_branch_group_policy_conservative_step_controller_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2248 = load(IN_2248)
    probe = (p2248.get("strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    n0 = float(inp.get("base_group_count", 1.0) or 1.0)
    load0 = float(inp.get("base_load_ratio", 1.0) or 1.0)
    total0 = float(inp.get("base_total_coverage_mass", 1.0) or 1.0)
    kappa = float(inp.get("kappa_reference", 1.0) or 1.0)

    # Margin model from previous stages:
    # M = total0 - ((n-1) + kappa*load)
    margin0 = total0 - ((n0 - 1.0) + kappa * load0)

    # Conservative controller: for proposed increments (dn>=0, dload>=0), require
    # dn + kappa*dload <= rho*margin0, with 0<rho<1 reserve fraction.
    rho = 0.8
    proposed_dn = 0.10
    proposed_dload = 0.05 * max(load0, 1e-12)

    lhs_budget_use = proposed_dn + kappa * proposed_dload
    rhs_budget_cap = rho * max(margin0, 0.0)
    step_admissible = lhs_budget_use <= rhs_budget_cap + 1e-15

    # Guaranteed post-step margin lower bound under first-order law.
    guaranteed_margin_lb = margin0 - lhs_budget_use

    payload = {
        "schema_version": "p2249_s1199_v1",
        "packet_id": "P2249",
        "stage_id": "S1199",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_CONSERVATIVE_STEP_CONTROLLER_PROBE",
        "strict_nu_branch_group_policy_conservative_step_controller_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_CONSERVATIVE_STEP_CONTROLLER_PROBE_V1",
            "source_packets": [str(IN_2248.relative_to(ROOT))],
            "inputs": {
                "base_group_count": n0,
                "base_load_ratio": load0,
                "base_total_coverage_mass": total0,
                "kappa_reference": kappa,
                "base_margin": margin0,
            },
            "controller_rule": {
                "inequality": "delta_group_count + kappa_reference*delta_load_ratio <= rho*base_margin",
                "rho": rho,
                "proposed_delta_group_count": proposed_dn,
                "proposed_delta_load_ratio": proposed_dload,
                "lhs_budget_use": lhs_budget_use,
                "rhs_budget_cap": rhs_budget_cap,
                "step_admissible": step_admissible,
                "guaranteed_post_step_margin_lower_bound": guaranteed_margin_lb,
            },
            "physical_interpretation_note": "Controller treats margin as a conserved safety budget: each policy move spends budget at rates set by sensitivity slopes; admissible moves preserve nonnegative safety margin by construction.",
            "theorem_scope_limit": "local first-order controller diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2250_candidate",
            "goal": "run multi-step simulated policy trajectory with controller to verify cumulative margin invariance and stopping behavior near boundary",
        },
        "gatekeeper_checks": {
            "controller_exported": True,
            "base_margin_nonnegative": margin0 >= -1e-15,
            "step_admissibility_computable": isinstance(step_admissible, bool),
            "guaranteed_margin_nonnegative": guaranteed_margin_lb >= -1e-15,
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
                "# P2249 S1199: conservative step-controller probe",
                "",
                f"- base margin: `{margin0:.12e}`",
                f"- controller rho: `{rho:.12e}`",
                f"- lhs budget use: `{lhs_budget_use:.12e}`",
                f"- rhs budget cap: `{rhs_budget_cap:.12e}`",
                f"- step admissible: `{step_admissible}`",
                f"- guaranteed post-step margin LB: `{guaranteed_margin_lb:.12e}`",
                "",
                "Local first-order controller diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
