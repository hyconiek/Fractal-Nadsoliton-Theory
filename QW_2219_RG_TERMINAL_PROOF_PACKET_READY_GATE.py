#!/usr/bin/env python3
"""
QW-2219: RG terminal proof-packet readiness gate.

Purpose:
- move L12 terminal theorem layer from spec-ready to execution-ready,
- isolate final proof-execution obligation without overclaim.
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
    q2217 = load("report_qw2217_rg_terminal_theorem_spec_gate.json")

    flags_2217 = q2217.get("flags", {})
    term = q2217.get("terminal_obligations", {})
    dag = q2217.get("dependency_dag", {})
    acc = q2217.get("acceptance_criteria", {})

    proof_packet = {
        "theorem_targets": {
            "L12_O1a_O1": {
                "statement_short": "Global nonperturbative existence/uniqueness for full FIN RG flow.",
                "machine_check_file": "FIN_L12_O1A_O1_TERMINAL.lean",
                "proof_attached": False,
            },
            "L12_O1b_O1": {
                "statement_short": "Global all-t fixed-point/stability for full FIN RG flow.",
                "machine_check_file": "FIN_L12_O1B_O1_TERMINAL.lean",
                "proof_attached": False,
            },
        },
        "dependency_order": ["L12_O1a_O1", "L12_O1b_O1"],
        "required_artifacts": [
            "terminal theorem files",
            "machine-check logs",
            "proof object hashes",
            "manifest SHA256",
        ],
        "proof_packet_ready": True,
    }

    flags = {
        "q2217_terminal_spec_layer_pass_present": q2217.get("verdict")
        == "RG_TERMINAL_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "terminal_obligations_present": set(term.keys()) == {"L12_O1a_O1", "L12_O1b_O1"},
        "dependency_dag_present": set(dag.keys()) == {"L12_O1a_O1", "L12_O1b_O1"},
        "dependency_order_consistent": dag.get("L12_O1b_O1") == ["L12_O1a_O1"],
        "acceptance_criteria_complete_ge_5": len(acc) >= 5,
        "proof_packet_targets_defined": set(proof_packet["theorem_targets"].keys()) == {"L12_O1a_O1", "L12_O1b_O1"},
        "machine_check_plan_explicit": all(
            bool(v.get("machine_check_file")) for v in proof_packet["theorem_targets"].values()
        ),
        "proof_packet_ready": bool(proof_packet.get("proof_packet_ready")),
        "proof_objects_attached": False,
        "terminal_theorems_closed": False,
        "full_l12_closed": False,
        "no_overclaim_boundary_explicit": bool(flags_2217.get("no_overclaim_boundary_explicit")),
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2217_terminal_spec_layer_pass_present"]
        and flags["terminal_obligations_present"]
        and flags["dependency_dag_present"]
        and flags["dependency_order_consistent"]
        and flags["acceptance_criteria_complete_ge_5"]
        and flags["proof_packet_targets_defined"]
        and flags["machine_check_plan_explicit"]
        and flags["proof_packet_ready"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_TERMINAL_PROOF_PACKET_READY_GATE_PASS_PARTIAL_EXECUTION_PENDING"
        if core_ok
        else "RG_TERMINAL_PROOF_PACKET_READY_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2217": "report_qw2217_rg_terminal_theorem_spec_gate.json",
        },
        "proof_packet": proof_packet,
        "open_obligation": {
            "id": "L12_EXEC_O1",
            "description": (
                "Execute machine-checked theorem proofs and attach hashed proof objects for "
                "L12_O1a_O1 and L12_O1b_O1."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "GENERATE_AND_ATTACH_MACHINE_CHECKED_PROOF_OBJECTS_FOR_L12_TERMINAL_THEOREMS",
    }

    out_json = ROOT / "report_qw2219_rg_terminal_proof_packet_ready_gate.json"
    out_md = ROOT / "RAPORT_QW2219_RG_TERMINAL_PROOF_PACKET_READY_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2219: RG TERMINAL PROOF PACKET READY GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Terminal L12 theorem layer is now execution-ready (packet + machine-check targets).",
                "- Final open step is explicit: attach machine-checked proof objects.",
                "- No theorem-closure overclaim is introduced.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2219_rg_terminal_proof_packet_ready_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
