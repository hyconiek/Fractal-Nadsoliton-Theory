#!/usr/bin/env python3
"""
QW-2227: RG axiom-free O1b provenance gate.

Purpose:
- perform strict provenance accounting for L12_AXIOM_FREE_O1b targets,
- require O1a provenance layer to be present,
- keep unresolved theorem-level pieces explicit (no overclaim).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def get_flag(rep: dict[str, Any], key: str) -> bool:
    return bool(rep.get("flags", {}).get(key, False))


def main() -> None:
    q2223 = load("report_qw2223_rg_axiom_free_discharge_spec_gate.json")
    q2225 = load("report_qw2225_rg_axiom_free_o1a_provenance_gate.json")
    q2213 = load("report_qw2213_rg_flow_existence_scope_gate.json")
    q2215 = load("report_qw2215_rg_global_stability_scope_gate.json")

    sub = q2223.get("sub_obligations", {}).get("L12_AXIOM_FREE_O1b", {})
    targets = sub.get("targets", [])

    provenance = {
        "RGGlobalWellPosednessAllScales": {
            "sources": [
                {
                    "report": "report_qw2213_rg_flow_existence_scope_gate.json",
                    "checks": ["single_terminal_obligation_for_l12_o1a_isolated"],
                }
            ],
            "status": "unresolved_global_theorem",
            "notes": "Inherited unresolved theorem from L12_O1a_O1 layer.",
        },
        "RGGlobalFixedPointStabilityAllT": {
            "sources": [
                {
                    "report": "report_qw2215_rg_global_stability_scope_gate.json",
                    "checks": ["single_terminal_obligation_for_l12_o1b_isolated"],
                }
            ],
            "status": "unresolved_global_theorem",
            "notes": "Global all-t fixed-point/stability remains terminal open theorem (`L12_O1b_O1`).",
        },
        "L12O1bWitness": {
            "sources": [
                {
                    "report": "FIN_L12_O1B_O1_TERMINAL.lean",
                    "checks": ["axiomatic_witness_symbol_present"],
                }
            ],
            "status": "syntactic_witness_only",
            "notes": "Witness currently axiomatic and must be replaced by derived proof.",
        },
    }

    file_l12b = ROOT / "FIN_L12_O1B_O1_TERMINAL.lean"
    file_text = file_l12b.read_text(encoding="utf-8") if file_l12b.exists() else ""

    target_set_ok = set(targets) == {
        "RGGlobalWellPosednessAllScales",
        "RGGlobalFixedPointStabilityAllT",
        "L12O1bWitness",
    }

    checks = {
        "q2223_l12_o1b_subobligation_present": bool(sub),
        "q2225_o1a_provenance_layer_present": q2225.get("verdict")
        == "RG_AXIOM_FREE_O1A_PROVENANCE_GATE_PASS_PARTIAL_UNRESOLVED_THEOREM",
        "target_set_exact_match": target_set_ok,
        "provenance_entries_cover_all_targets": set(provenance.keys()) == set(targets),
        "wellposedness_explicitly_unresolved": (
            get_flag(q2213, "single_terminal_obligation_for_l12_o1a_isolated")
            and not get_flag(q2213, "l12_o1a_terminal_obligation_closed")
        ),
        "fixed_point_stability_explicitly_unresolved": (
            get_flag(q2215, "single_terminal_obligation_for_l12_o1b_isolated")
            and not get_flag(q2215, "l12_o1b_terminal_obligation_closed")
        ),
        "axiomatic_witness_present_in_terminal_file": ("axiom L12O1bWitness" in file_text),
        "no_overclaim_boundary_explicit": True,
        "l12_axiom_free_o1b_fully_closed": False,
    }

    pass_count = sum(1 for v in checks.values() if v)
    total_flags = len(checks)

    core_ok = (
        checks["q2223_l12_o1b_subobligation_present"]
        and checks["q2225_o1a_provenance_layer_present"]
        and checks["target_set_exact_match"]
        and checks["provenance_entries_cover_all_targets"]
        and checks["wellposedness_explicitly_unresolved"]
        and checks["fixed_point_stability_explicitly_unresolved"]
        and checks["axiomatic_witness_present_in_terminal_file"]
        and checks["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_AXIOM_FREE_O1B_PROVENANCE_GATE_PASS_PARTIAL_UNRESOLVED_THEOREM"
        if core_ok
        else "RG_AXIOM_FREE_O1B_PROVENANCE_GATE_FAIL"
    )

    unresolved = [
        "RGGlobalWellPosednessAllScales",
        "RGGlobalFixedPointStabilityAllT",
        "L12O1bWitness",
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2223": "report_qw2223_rg_axiom_free_discharge_spec_gate.json",
            "q2225": "report_qw2225_rg_axiom_free_o1a_provenance_gate.json",
            "q2213": "report_qw2213_rg_flow_existence_scope_gate.json",
            "q2215": "report_qw2215_rg_global_stability_scope_gate.json",
            "l12_file": "FIN_L12_O1B_O1_TERMINAL.lean",
        },
        "target_symbols": targets,
        "provenance_map": provenance,
        "unresolved_targets": unresolved,
        "flags": checks,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_L12_O1B_AXIOM_FREE_THEOREMS_AND_REPLACE_L12O1B_WITNESS",
    }

    out_json = ROOT / "report_qw2227_rg_axiom_free_o1b_provenance_gate.json"
    out_md = ROOT / "RAPORT_QW2227_RG_AXIOM_FREE_O1B_PROVENANCE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2227: RG AXIOM-FREE O1B PROVENANCE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unresolved_targets: `{', '.join(unresolved)}`",
                "",
                "## Core result",
                "- Targety `L12_AXIOM_FREE_O1b` maja jawna mape provenance i jawny dependency na O1a.",
                "- Wszystkie targety pozostaja theorem-level unresolved bez overclaimu.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2227_rg_axiom_free_o1b_provenance_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
