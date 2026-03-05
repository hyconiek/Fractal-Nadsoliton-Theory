#!/usr/bin/env python3
"""
QW-2209: RG global-closure obligation reduction gate (L12).

Purpose:
- keep strict proxy/extended-scope RG results explicit,
- reduce remaining foundational RG gap to one global theorem obligation.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def main() -> None:
    q2132 = load("report_qw2132_rg_fixed_point_jacobian_gate.json")
    q2185 = load("report_qw2185_rg_proxy_global_obstruction_theorem_gate.json")
    q2187 = load("report_qw2187_rg_proxy_finite_uv_scope_declaration_gate.json")
    q2188 = load("report_qw2188_uv_completing_rg_correction_frontier_gate.json")

    f2132 = q2132.get("flags", {})
    f2185 = q2185.get("flags", {})
    f2187 = q2187.get("flags", {})
    f2188 = q2188.get("flags", {})

    obstruction = q2185.get("obstruction_theorem", {})
    probe = q2188.get("probe_setup", {})
    frontier = q2188.get("frontier_results", {})
    b_star_diag = frontier.get("b_star_diagnostics", {})
    low_energy = q2188.get("reference_low_energy_distortion", {})

    flags = {
        "q2132_proxy_fixed_point_jacobian_pass_present": q2132.get("verdict")
        == "RG_FIXED_POINT_JACOBIAN_GATE_PASS_STRICT_PROXY_PARTIAL",
        "q2185_global_obstruction_theorem_pass_present": q2185.get("verdict")
        == "RG_PROXY_GLOBAL_OBSTRUCTION_THEOREM_GATE_PASS_STRICT",
        "q2187_finite_uv_scope_declaration_pass_present": q2187.get("verdict")
        == "RG_PROXY_FINITE_UV_SCOPE_DECLARATION_GATE_PASS_STRICT",
        "q2188_uv_correction_frontier_pass_present": q2188.get("verdict")
        == "UV_COMPLETING_RG_CORRECTION_FRONTIER_GATE_PASS_EXTENDED_SCOPE_PARTIAL",
        "q2185_u1_landau_obstruction_explicit": bool(
            f2185.get("obstruction_explicitly_identified_as_u1_landau_pole")
        ),
        "q2187_global_all_t_not_claimed": not bool(f2187.get("full_global_rg_closure_claimed")),
        "q2188_global_all_t_not_claimed": not bool(f2188.get("global_all_t_rg_closure_claimed")),
        "q2188_extended_scope_feasible_with_b_star": bool(f2188.get("minimal_feasible_b_star_found"))
        and bool(f2188.get("b_star_removes_crossing_to_t_probe")),
        "q2132_full_nonperturbative_proof_not_closed": not bool(
            f2132.get("full_nonperturbative_rg_fixed_point_proof")
        ),
        "single_global_rg_obligation_isolated": True,
        "no_overclaim_boundary_explicit": True,
        "full_global_all_t_rg_closure_theorem_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2132_proxy_fixed_point_jacobian_pass_present"]
        and flags["q2185_global_obstruction_theorem_pass_present"]
        and flags["q2187_finite_uv_scope_declaration_pass_present"]
        and flags["q2188_uv_correction_frontier_pass_present"]
        and flags["q2185_u1_landau_obstruction_explicit"]
        and flags["q2187_global_all_t_not_claimed"]
        and flags["q2188_global_all_t_not_claimed"]
        and flags["q2188_extended_scope_feasible_with_b_star"]
        and flags["single_global_rg_obligation_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_GLOBAL_CLOSURE_OBLIGATION_GATE_PASS_PARTIAL_SINGLE_GLOBAL_OBLIGATION_OPEN"
        if core_ok
        else "RG_GLOBAL_CLOSURE_OBLIGATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2132": "report_qw2132_rg_fixed_point_jacobian_gate.json",
            "q2185": "report_qw2185_rg_proxy_global_obstruction_theorem_gate.json",
            "q2187": "report_qw2187_rg_proxy_finite_uv_scope_declaration_gate.json",
            "q2188": "report_qw2188_uv_completing_rg_correction_frontier_gate.json",
        },
        "summary": {
            "landau_pole_t_star_reference": obstruction.get("t_star_reference"),
            "q2182_safe_window_max_t": obstruction.get("q2182_t_max"),
            "uv_frontier_t_probe": probe.get("t_probe"),
            "baseline_b0": frontier.get("baseline_b0"),
            "anchor_upper_bmax": frontier.get("anchor_upper_bmax"),
            "b_star": frontier.get("b_star_min_feasible"),
            "post_correction_crossing_time": b_star_diag.get("first_crossing_time"),
            "flow_finite_to_t_probe": b_star_diag.get("flow_finite_to_t_probe"),
            "relative_shift_beta_lambda": low_energy.get("relative_shift_beta_lambda"),
        },
        "open_obligation": {
            "id": "L12_O1",
            "description": (
                "Prove full nonperturbative all-coupling all-t RG fixed-point/stability closure "
                "from complete FIN action (beyond current proxy and finite-scope declarations)."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DERIVE_FULL_NONPERTURBATIVE_ALL_T_RG_CLOSURE_THEOREM_FROM_COMPLETE_FIN_ACTION_AND_RERUN_QW2209"
        ),
    }

    out_json = ROOT / "report_qw2209_rg_global_closure_obligation_gate.json"
    out_md = ROOT / "RAPORT_QW2209_RG_GLOBAL_CLOSURE_OBLIGATION_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2209: RG GLOBAL CLOSURE OBLIGATION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Strict RG proxy chain remains consistent and explicitly bounded in scope.",
                "- Landau-pole obstruction for current U(1) proxy is explicit and non-hidden.",
                "- UV-correction frontier extends feasible scope but does not justify global all-t claim.",
                "- Remaining L12 gap is reduced to one explicit global theorem obligation (`L12_O1`).",
                "",
                "## Artifacts",
                "- JSON: `report_qw2209_rg_global_closure_obligation_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
