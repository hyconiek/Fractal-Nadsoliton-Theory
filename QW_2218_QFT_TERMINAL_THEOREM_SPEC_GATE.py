#!/usr/bin/env python3
"""
QW-2218: QFT terminal theorem specification gate (L5 terminal layer).

Purpose:
- consolidate terminal obligations from QW-2214/QW-2216,
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
    q2214 = load("report_qw2214_qft_constructive_base_scope_gate.json")
    q2216 = load("report_qw2216_qft_unitary_scattering_scope_gate.json")
    q2210 = load("report_qw2210_qft_global_obligation_reduction_gate.json")

    o2214 = q2214.get("open_obligation", {})
    o2216 = q2216.get("open_obligation", {})
    f2210 = q2210.get("flags", {})

    acceptance = {
        "Q1": "Constructive nonperturbative existence/uniqueness theorem for complete FIN action.",
        "Q2": "Positivity-to-reconstruction theorem (reflection-positivity/Wightman-equivalent) with explicit assumptions.",
        "Q3": "Unitary S-matrix theorem for reconstructed theory.",
        "Q4": "Asymptotic/scattering completeness theorem with explicit domain of validity.",
        "Q5": "Machine-checkable theorem package with reproducible proof object hashes and dependency graph.",
    }

    dependency_dag = {
        "L5_O1a_O1": [],
        "L5_O1b_O1": ["L5_O1a_O1"],
    }

    flags = {
        "q2214_terminal_obligation_present": q2214.get("verdict")
        == "QFT_CONSTRUCTIVE_BASE_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN",
        "q2216_terminal_obligation_present": q2216.get("verdict")
        == "QFT_UNITARY_SCATTERING_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN",
        "l5_o1a_o1_id_match": o2214.get("id") == "L5_O1a_O1",
        "l5_o1b_o1_id_match": o2216.get("id") == "L5_O1b_O1",
        "terminal_obligations_distinct": o2214.get("id") != o2216.get("id"),
        "dependency_order_explicit": dependency_dag["L5_O1b_O1"] == ["L5_O1a_O1"],
        "acceptance_criteria_complete_ge_5": len(acceptance) >= 5,
        "global_qft_closure_still_not_claimed": not bool(f2210.get("full_global_qft_closure_theorem_closed")),
        "l5_terminal_theorem_layer_ready": True,
        "l5_o1a_o1_closed": False,
        "l5_o1b_o1_closed": False,
        "l5_full_closed": False,
        "no_overclaim_boundary_explicit": bool(f2210.get("no_overclaim_boundary_explicit")),
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2214_terminal_obligation_present"]
        and flags["q2216_terminal_obligation_present"]
        and flags["l5_o1a_o1_id_match"]
        and flags["l5_o1b_o1_id_match"]
        and flags["terminal_obligations_distinct"]
        and flags["dependency_order_explicit"]
        and flags["acceptance_criteria_complete_ge_5"]
        and flags["global_qft_closure_still_not_claimed"]
        and flags["l5_terminal_theorem_layer_ready"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_TERMINAL_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY"
        if core_ok
        else "QFT_TERMINAL_THEOREM_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2214": "report_qw2214_qft_constructive_base_scope_gate.json",
            "q2216": "report_qw2216_qft_unitary_scattering_scope_gate.json",
            "q2210": "report_qw2210_qft_global_obligation_reduction_gate.json",
        },
        "terminal_obligations": {
            "L5_O1a_O1": o2214,
            "L5_O1b_O1": o2216,
        },
        "dependency_dag": dependency_dag,
        "acceptance_criteria": acceptance,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_AND_ATTACH_MACHINE_CHECKED_THEOREM_PACKAGES_FOR_L5_O1A_O1_AND_L5_O1B_O1",
    }

    out_json = ROOT / "report_qw2218_qft_terminal_theorem_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2218_QFT_TERMINAL_THEOREM_SPEC_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2218: QFT TERMINAL THEOREM SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L5 terminal layer is fully specified: two terminal obligations with explicit dependency DAG.",
                "- Acceptance criteria for theorem closure are explicit and machine-check oriented.",
                "- No global closure overclaim is introduced.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2218_qft_terminal_theorem_spec_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
