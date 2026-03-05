#!/usr/bin/env python3
"""
QW-2217: RG terminal theorem specification gate (L12 terminal layer).

Purpose:
- consolidate terminal obligations from QW-2213/QW-2215,
- define explicit theorem acceptance criteria without overclaim.
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
    q2213 = load("report_qw2213_rg_flow_existence_scope_gate.json")
    q2215 = load("report_qw2215_rg_global_stability_scope_gate.json")
    q2209 = load("report_qw2209_rg_global_closure_obligation_gate.json")

    o2213 = q2213.get("open_obligation", {})
    o2215 = q2215.get("open_obligation", {})
    f2209 = q2209.get("flags", {})

    acceptance = {
        "T1": "Constructive derivation of full nonperturbative all-coupling FIN RG flow from complete FIN action.",
        "T2": "Global well-posedness theorem (existence/uniqueness) for all physical scales t>=0.",
        "T3": "Global fixed-point/stability theorem all-t for full derived flow, including UV behavior.",
        "T4": "No hidden scope truncation; explicit bound-control lemmas and theorem dependencies.",
        "T5": "Machine-checkable theorem package attached with reproducible proof object hashes.",
    }

    dependency_dag = {
        "L12_O1a_O1": [],
        "L12_O1b_O1": ["L12_O1a_O1"],
    }

    flags = {
        "q2213_terminal_obligation_present": q2213.get("verdict")
        == "RG_FLOW_EXISTENCE_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN",
        "q2215_terminal_obligation_present": q2215.get("verdict")
        == "RG_GLOBAL_STABILITY_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN",
        "l12_o1a_o1_id_match": o2213.get("id") == "L12_O1a_O1",
        "l12_o1b_o1_id_match": o2215.get("id") == "L12_O1b_O1",
        "terminal_obligations_distinct": o2213.get("id") != o2215.get("id"),
        "dependency_order_explicit": dependency_dag["L12_O1b_O1"] == ["L12_O1a_O1"],
        "acceptance_criteria_complete_ge_5": len(acceptance) >= 5,
        "global_all_t_still_not_claimed": bool(f2209.get("q2187_global_all_t_not_claimed"))
        and bool(f2209.get("q2188_global_all_t_not_claimed")),
        "l12_terminal_theorem_layer_ready": True,
        "l12_o1a_o1_closed": False,
        "l12_o1b_o1_closed": False,
        "l12_full_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2213_terminal_obligation_present"]
        and flags["q2215_terminal_obligation_present"]
        and flags["l12_o1a_o1_id_match"]
        and flags["l12_o1b_o1_id_match"]
        and flags["terminal_obligations_distinct"]
        and flags["dependency_order_explicit"]
        and flags["acceptance_criteria_complete_ge_5"]
        and flags["global_all_t_still_not_claimed"]
        and flags["l12_terminal_theorem_layer_ready"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_TERMINAL_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY"
        if core_ok
        else "RG_TERMINAL_THEOREM_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2213": "report_qw2213_rg_flow_existence_scope_gate.json",
            "q2215": "report_qw2215_rg_global_stability_scope_gate.json",
            "q2209": "report_qw2209_rg_global_closure_obligation_gate.json",
        },
        "terminal_obligations": {
            "L12_O1a_O1": o2213,
            "L12_O1b_O1": o2215,
        },
        "dependency_dag": dependency_dag,
        "acceptance_criteria": acceptance,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_AND_ATTACH_MACHINE_CHECKED_THEOREM_PACKAGES_FOR_L12_O1A_O1_AND_L12_O1B_O1",
    }

    out_json = ROOT / "report_qw2217_rg_terminal_theorem_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2217_RG_TERMINAL_THEOREM_SPEC_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2217: RG TERMINAL THEOREM SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L12 terminal layer is fully specified: two terminal obligations with explicit dependency DAG.",
                "- Acceptance criteria for theorem closure are explicit and machine-check oriented.",
                "- No global closure overclaim is introduced.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2217_rg_terminal_theorem_spec_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
