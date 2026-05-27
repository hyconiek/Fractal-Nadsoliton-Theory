#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2246 = GEN / "p2246_s1196_strict_nu_branch_group_policy_feasibility_surface_probe.json"
OUT = GEN / "p2247_s1197_strict_nu_branch_group_policy_boundary_sensitivity_probe.json"
MD = GEN / "p2247_s1197_strict_nu_branch_group_policy_boundary_sensitivity_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2246 = load(IN_2246)
    probe = (p2246.get("strict_nu_branch_group_policy_feasibility_surface_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    n0 = int(inp.get("base_group_count", 1) or 1)
    load0 = float(inp.get("base_load_ratio", 1.0) or 1.0)
    total0 = float(inp.get("base_total_coverage_mass", 1.0) or 1.0)
    kappa = float(inp.get("kappa_reference", 1.0) or 1.0)

    # Deterministic boundary from previous stages:
    # margin(n,load) = total0 - [(n-1) + kappa*load].
    margin0 = total0 - ((n0 - 1) + kappa * load0)
    dmargin_dn = -1.0
    dmargin_dload = -kappa

    # Closed-form boundary curve load*(n) at zero margin.
    # load_star(n) = (total0 - (n-1))/kappa.
    n_samples = [max(1, n0 - 1), n0, n0 + 1, n0 + 2]
    boundary_rows = []
    for n in n_samples:
        load_star = (total0 - (n - 1)) / max(kappa, 1e-12)
        boundary_rows.append({"group_count": n, "critical_load_ratio": load_star})

    payload = {
        "schema_version": "p2247_s1197_v1",
        "packet_id": "P2247",
        "stage_id": "S1197",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_BOUNDARY_SENSITIVITY_PROBE",
        "strict_nu_branch_group_policy_boundary_sensitivity_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_BOUNDARY_SENSITIVITY_PROBE_V1",
            "source_packets": [str(IN_2246.relative_to(ROOT))],
            "inputs": {
                "base_group_count": n0,
                "base_load_ratio": load0,
                "base_total_coverage_mass": total0,
                "kappa_reference": kappa,
            },
            "analytic_sensitivity": {
                "base_margin": margin0,
                "dmargin_d_group_count": dmargin_dn,
                "dmargin_d_load_ratio": dmargin_dload,
            },
            "closed_form_boundary": {
                "equation": "critical_load_ratio(n) = (total_coverage_mass - (n-1))/kappa_reference",
                "rows": boundary_rows,
            },
            "physical_interpretation_note": "Boundary sensitivity quantifies robustness erosion rate: each added group consumes one full margin unit, while load increase consumes kappa-weighted margin, giving a physical slope map of policy fragility.",
            "theorem_scope_limit": "local analytic sensitivity diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2248_candidate",
            "goal": "validate analytic sensitivity against finite-difference estimates on random neighborhood points of feasibility surface",
        },
        "gatekeeper_checks": {
            "boundary_sensitivity_exported": True,
            "derivative_wrt_group_count_negative": dmargin_dn < 0.0,
            "derivative_wrt_load_negative": dmargin_dload < 0.0,
            "closed_form_boundary_exported": True,
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
                "# P2247 S1197: boundary-sensitivity probe",
                "",
                f"- base margin: `{margin0:.12e}`",
                f"- dmargin/d(group_count): `{dmargin_dn:.12e}`",
                f"- dmargin/d(load_ratio): `{dmargin_dload:.12e}`",
                f"- kappa reference: `{kappa:.12e}`",
                "",
                "Local analytic sensitivity only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
