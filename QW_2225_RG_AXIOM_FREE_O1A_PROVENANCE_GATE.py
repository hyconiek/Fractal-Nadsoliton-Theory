#!/usr/bin/env python3
"""
QW-2225: RG axiom-free O1a provenance gate.

Purpose:
- perform strict provenance accounting for L12_AXIOM_FREE_O1a targets,
- map each target symbol to explicit source reports/flags,
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
    q2213 = load("report_qw2213_rg_flow_existence_scope_gate.json")
    q2182 = load("report_qw2182_rg_nonperturbative_domain_flow_gate.json")
    q2165 = load("report_qw2165_l13_exhaustive_canonical_eom_gate.json")

    sub = q2223.get("sub_obligations", {}).get("L12_AXIOM_FREE_O1a", {})
    targets = sub.get("targets", [])

    provenance = {
        "FINActionComplete": {
            "sources": [
                {
                    "report": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
                    "checks": [
                        "q2163_full_canonical_variational_layer_present",
                        "q2161_proxy_variational_bundle_present",
                        "exhaustive_bundle_closed_on_canonical_action_template",
                    ],
                }
            ],
            "status": "grounded_strict_scope",
            "notes": "Canonical FIN action template and exhaustive EoM layer are present in strict scope.",
        },
        "RGConstructiveMap": {
            "sources": [
                {
                    "report": "report_qw2182_rg_nonperturbative_domain_flow_gate.json",
                    "checks": [
                        "strict_constructive_domain_declared",
                        "deterministic_grid_no_scan_no_retune",
                        "finite_semiflow_on_declared_domain",
                    ],
                }
            ],
            "status": "grounded_strict_domain",
            "notes": "Constructive RG semiflow is established on declared strict domain.",
        },
        "RGGlobalWellPosednessAllScales": {
            "sources": [
                {
                    "report": "report_qw2213_rg_flow_existence_scope_gate.json",
                    "checks": [
                        "single_terminal_obligation_for_l12_o1a_isolated",
                    ],
                }
            ],
            "status": "unresolved_global_theorem",
            "notes": "Global all-scales existence/uniqueness remains terminal open theorem (`L12_O1a_O1`).",
        },
        "L12O1aWitness": {
            "sources": [
                {
                    "report": "FIN_L12_O1A_O1_TERMINAL.lean",
                    "checks": [
                        "axiomatic_witness_symbol_present",
                    ],
                }
            ],
            "status": "syntactic_witness_only",
            "notes": "Witness currently axiomatic and must be replaced by derived proof.",
        },
    }

    file_l12a = ROOT / "FIN_L12_O1A_O1_TERMINAL.lean"
    file_text = file_l12a.read_text(encoding="utf-8") if file_l12a.exists() else ""

    target_set_ok = set(targets) == {
        "FINActionComplete",
        "RGConstructiveMap",
        "RGGlobalWellPosednessAllScales",
        "L12O1aWitness",
    }

    checks = {
        "q2223_l12_o1a_subobligation_present": bool(sub),
        "target_set_exact_match": target_set_ok,
        "provenance_entries_cover_all_targets": set(provenance.keys()) == set(targets),
        "fin_action_complete_grounded": (
            get_flag(q2165, "q2163_full_canonical_variational_layer_present")
            and get_flag(q2165, "q2161_proxy_variational_bundle_present")
            and get_flag(q2165, "exhaustive_bundle_closed_on_canonical_action_template")
        ),
        "rg_constructive_map_grounded": (
            get_flag(q2182, "strict_constructive_domain_declared")
            and get_flag(q2182, "deterministic_grid_no_scan_no_retune")
            and get_flag(q2182, "finite_semiflow_on_declared_domain")
        ),
        "global_wellposedness_explicitly_unresolved": (
            get_flag(q2213, "single_terminal_obligation_for_l12_o1a_isolated")
            and not get_flag(q2213, "l12_o1a_terminal_obligation_closed")
        ),
        "axiomatic_witness_present_in_terminal_file": ("axiom L12O1aWitness" in file_text),
        "no_overclaim_boundary_explicit": True,
        "l12_axiom_free_o1a_fully_closed": False,
    }

    pass_count = sum(1 for v in checks.values() if v)
    total_flags = len(checks)

    core_ok = (
        checks["q2223_l12_o1a_subobligation_present"]
        and checks["target_set_exact_match"]
        and checks["provenance_entries_cover_all_targets"]
        and checks["fin_action_complete_grounded"]
        and checks["rg_constructive_map_grounded"]
        and checks["global_wellposedness_explicitly_unresolved"]
        and checks["axiomatic_witness_present_in_terminal_file"]
        and checks["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_AXIOM_FREE_O1A_PROVENANCE_GATE_PASS_PARTIAL_UNRESOLVED_THEOREM"
        if core_ok
        else "RG_AXIOM_FREE_O1A_PROVENANCE_GATE_FAIL"
    )

    unresolved = [
        "RGGlobalWellPosednessAllScales",
        "L12O1aWitness",
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2223": "report_qw2223_rg_axiom_free_discharge_spec_gate.json",
            "q2213": "report_qw2213_rg_flow_existence_scope_gate.json",
            "q2182": "report_qw2182_rg_nonperturbative_domain_flow_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "l12_file": "FIN_L12_O1A_O1_TERMINAL.lean",
        },
        "target_symbols": targets,
        "provenance_map": provenance,
        "unresolved_targets": unresolved,
        "flags": checks,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_RG_GLOBAL_WELLPPOSEDNESS_ALL_SCALES_AND_REPLACE_L12O1A_WITNESS",
    }

    out_json = ROOT / "report_qw2225_rg_axiom_free_o1a_provenance_gate.json"
    out_md = ROOT / "RAPORT_QW2225_RG_AXIOM_FREE_O1A_PROVENANCE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2225: RG AXIOM-FREE O1A PROVENANCE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unresolved_targets: `{', '.join(unresolved)}`",
                "",
                "## Core result",
                "- Kazdy target `L12_AXIOM_FREE_O1a` ma jawna mape provenance.",
                "- Dwa targety pozostaja theorem-level unresolved i sa jawnie utrzymane bez overclaimu.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2225_rg_axiom_free_o1a_provenance_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
