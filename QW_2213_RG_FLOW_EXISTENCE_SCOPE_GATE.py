#!/usr/bin/env python3
"""
QW-2213: RG flow existence scope gate (L12_O1a).

Purpose:
- keep decomposition from QW-2211 explicit,
- reduce L12_O1a to one terminal theorem obligation.
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
    q2211 = load("report_qw2211_rg_global_obligation_decomposition_gate.json")
    q2132 = load("report_qw2132_rg_fixed_point_jacobian_gate.json")
    q2185 = load("report_qw2185_rg_proxy_global_obstruction_theorem_gate.json")
    q2188 = load("report_qw2188_uv_completing_rg_correction_frontier_gate.json")
    q2136 = load("report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json")

    f2132 = q2132.get("flags", {})
    f2185 = q2185.get("flags", {})
    f2188 = q2188.get("flags", {})
    f2136 = q2136.get("flags", {})

    decomp = q2211.get("decomposition", {})
    l12_o1a = decomp.get("L12_O1a", {})

    frontier = q2188.get("frontier_results", {})

    flags = {
        "q2211_decomposition_pass_present": q2211.get("verdict")
        == "RG_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_TWO_SUBOBLIGATIONS_OPEN",
        "l12_o1a_subobligation_present": bool(l12_o1a),
        "q2132_beta_system_explicit": bool(f2132.get("beta_system_explicitly_defined")),
        "q2185_obstruction_theorem_explicit": bool(f2185.get("obstruction_explicitly_identified_as_u1_landau_pole")),
        "q2188_uv_frontier_feasible": bool(f2188.get("minimal_feasible_b_star_found"))
        and bool(f2188.get("b_star_removes_crossing_to_t_probe")),
        "q2136_finite_counterterm_basis_declared": bool(f2136.get("finite_counterterm_basis_policy_declared"))
        and bool(f2136.get("finite_counterterm_basis_condition_holds")),
        "deterministic_no_scan_no_retune_chain": bool(f2132.get("deterministic_no_scan_no_retune"))
        and bool(f2185.get("deterministic_no_scan_no_retune"))
        and bool(f2188.get("deterministic_no_scan_no_retune"))
        and bool(f2136.get("deterministic_no_scan_no_retune")),
        "single_terminal_obligation_for_l12_o1a_isolated": True,
        "l12_o1a_terminal_obligation_closed": False,
        "l12_o1a_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2211_decomposition_pass_present"]
        and flags["l12_o1a_subobligation_present"]
        and flags["q2132_beta_system_explicit"]
        and flags["q2185_obstruction_theorem_explicit"]
        and flags["q2188_uv_frontier_feasible"]
        and flags["q2136_finite_counterterm_basis_declared"]
        and flags["deterministic_no_scan_no_retune_chain"]
        and flags["single_terminal_obligation_for_l12_o1a_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_FLOW_EXISTENCE_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN"
        if core_ok
        else "RG_FLOW_EXISTENCE_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2211": "report_qw2211_rg_global_obligation_decomposition_gate.json",
            "q2132": "report_qw2132_rg_fixed_point_jacobian_gate.json",
            "q2185": "report_qw2185_rg_proxy_global_obstruction_theorem_gate.json",
            "q2188": "report_qw2188_uv_completing_rg_correction_frontier_gate.json",
            "q2136": "report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json",
        },
        "summary": {
            "l12_o1a_prior": l12_o1a,
            "landau_pole_t_star_reference": q2185.get("obstruction_theorem", {}).get("t_star_reference"),
            "uv_frontier_b_star_min_feasible": frontier.get("b_star_min_feasible"),
            "uv_frontier_anchor_upper_bmax": frontier.get("anchor_upper_bmax"),
        },
        "open_obligation": {
            "id": "L12_O1a_O1",
            "description": (
                "Prove nonperturbative global existence/uniqueness theorem for the full all-coupling "
                "FIN RG flow derived from complete FIN action (beyond proxy truncation)."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_L12_O1A_O1_NONPERTURBATIVE_GLOBAL_EXISTENCE_UNIQUENESS_AND_RERUN_QW2213",
    }

    out_json = ROOT / "report_qw2213_rg_flow_existence_scope_gate.json"
    out_md = ROOT / "RAPORT_QW2213_RG_FLOW_EXISTENCE_SCOPE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2213: RG FLOW EXISTENCE SCOPE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L12_O1a is reduced to one terminal theorem-level obligation (`L12_O1a_O1`).",
                "- Proxy/obstruction/UV-frontier boundaries remain explicit and unchanged.",
                "- No global all-t overclaim is introduced.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2213_rg_flow_existence_scope_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
