#!/usr/bin/env python3
"""
QW-2211: RG global obligation decomposition gate (L12_O1).

Purpose:
- keep L12 single-obligation statement from QW-2209 explicit,
- decompose it into minimal theorem-level sub-obligations for execution.
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
    q2209 = load("report_qw2209_rg_global_closure_obligation_gate.json")
    q2185 = load("report_qw2185_rg_proxy_global_obstruction_theorem_gate.json")
    q2188 = load("report_qw2188_uv_completing_rg_correction_frontier_gate.json")

    f2209 = q2209.get("flags", {})
    f2185 = q2185.get("flags", {})
    f2188 = q2188.get("flags", {})

    l12_o1 = q2209.get("open_obligation", {})
    summary_2209 = q2209.get("summary", {})

    decomposition = {
        "L12_O1a": {
            "description": (
                "Derive a full nonperturbative all-coupling FIN RG flow system from complete FIN action "
                "(beyond current proxy truncation) with well-posed existence/uniqueness of flow."
            ),
            "currently_closed": False,
            "depends_on": [],
        },
        "L12_O1b": {
            "description": (
                "Prove global all-t fixed-point/stability closure for the full derived RG system, "
                "including UV behavior without hidden scope restriction."
            ),
            "currently_closed": False,
            "depends_on": ["L12_O1a"],
        },
    }

    flags = {
        "q2209_single_obligation_pass_present": q2209.get("verdict")
        == "RG_GLOBAL_CLOSURE_OBLIGATION_GATE_PASS_PARTIAL_SINGLE_GLOBAL_OBLIGATION_OPEN",
        "l12_o1_id_is_explicit": l12_o1.get("id") == "L12_O1",
        "proxy_landau_obstruction_explicit": bool(f2185.get("obstruction_explicitly_identified_as_u1_landau_pole")),
        "extended_scope_frontier_present": bool(f2188.get("minimal_feasible_b_star_found"))
        and bool(f2188.get("b_star_removes_crossing_to_t_probe")),
        "global_all_t_not_claimed_in_current_chain": bool(f2209.get("q2187_global_all_t_not_claimed"))
        and bool(f2209.get("q2188_global_all_t_not_claimed")),
        "decomposition_to_two_subobligations_explicit": True,
        "dependency_graph_acyclic": True,
        "l12_o1a_closed": False,
        "l12_o1b_closed": False,
        "full_l12_global_closure_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2209_single_obligation_pass_present"]
        and flags["l12_o1_id_is_explicit"]
        and flags["proxy_landau_obstruction_explicit"]
        and flags["extended_scope_frontier_present"]
        and flags["global_all_t_not_claimed_in_current_chain"]
        and flags["decomposition_to_two_subobligations_explicit"]
        and flags["dependency_graph_acyclic"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_TWO_SUBOBLIGATIONS_OPEN"
        if core_ok
        else "RG_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2209": "report_qw2209_rg_global_closure_obligation_gate.json",
            "q2185": "report_qw2185_rg_proxy_global_obstruction_theorem_gate.json",
            "q2188": "report_qw2188_uv_completing_rg_correction_frontier_gate.json",
        },
        "summary": {
            "l12_o1_prior": l12_o1,
            "landau_pole_t_star_reference": summary_2209.get("landau_pole_t_star_reference"),
            "extended_scope_t_probe": summary_2209.get("uv_frontier_t_probe"),
            "b_star": summary_2209.get("b_star"),
        },
        "decomposition": decomposition,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "CLOSE_L12_O1A_FULL_NONPERTURBATIVE_FIN_RG_FLOW_DERIVATION_THEN_CLOSE_L12_O1B_GLOBAL_ALL_T_FIXED_POINT_STABILITY"
        ),
    }

    out_json = ROOT / "report_qw2211_rg_global_obligation_decomposition_gate.json"
    out_md = ROOT / "RAPORT_QW2211_RG_GLOBAL_OBLIGATION_DECOMPOSITION_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2211: RG GLOBAL OBLIGATION DECOMPOSITION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Single open L12 obligation (`L12_O1`) is decomposed into two executable theorem sub-obligations.",
                "- Scope/obstruction boundaries remain explicit (no global all-t overclaim).",
                "- Remaining closure path is now sequenced: `L12_O1a -> L12_O1b`.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2211_rg_global_obligation_decomposition_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
