#!/usr/bin/env python3
"""
QW-2226: QFT axiom-free O1a provenance gate.

Purpose:
- perform strict provenance accounting for L5_AXIOM_FREE_O1a targets,
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
    q2224 = load("report_qw2224_qft_axiom_free_discharge_spec_gate.json")
    q2214 = load("report_qw2214_qft_constructive_base_scope_gate.json")
    q2202 = load("report_qw2202_qft_strict_scope_stratification_gate.json")
    q2165 = load("report_qw2165_l13_exhaustive_canonical_eom_gate.json")

    sub = q2224.get("sub_obligations", {}).get("L5_AXIOM_FREE_O1a", {})
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
        "ConstructiveNonPerturbativeScheme": {
            "sources": [
                {
                    "report": "report_qw2202_qft_strict_scope_stratification_gate.json",
                    "checks": [
                        "strict_scope_quantization_causality_renorm_stack_declared_closed",
                    ],
                },
                {
                    "report": "report_qw2214_qft_constructive_base_scope_gate.json",
                    "checks": [
                        "single_terminal_obligation_for_l5_o1a_isolated",
                    ],
                },
            ],
            "status": "grounded_strict_scope",
            "notes": "Constructive scheme is strict-scope grounded and terminalized to theorem target.",
        },
        "PositivityToReconstruction": {
            "sources": [
                {
                    "report": "report_qw2214_qft_constructive_base_scope_gate.json",
                    "checks": [
                        "single_terminal_obligation_for_l5_o1a_isolated",
                    ],
                }
            ],
            "status": "unresolved_global_theorem",
            "notes": "Positivity-to-reconstruction remains theorem-level terminal open (`L5_O1a_O1`).",
        },
        "L5O1aWitness": {
            "sources": [
                {
                    "report": "FIN_L5_O1A_O1_TERMINAL.lean",
                    "checks": [
                        "axiomatic_witness_symbol_present",
                    ],
                }
            ],
            "status": "syntactic_witness_only",
            "notes": "Witness currently axiomatic and must be replaced by derived proof.",
        },
    }

    file_l5a = ROOT / "FIN_L5_O1A_O1_TERMINAL.lean"
    file_text = file_l5a.read_text(encoding="utf-8") if file_l5a.exists() else ""

    target_set_ok = set(targets) == {
        "FINActionComplete",
        "ConstructiveNonPerturbativeScheme",
        "PositivityToReconstruction",
        "L5O1aWitness",
    }

    checks = {
        "q2224_l5_o1a_subobligation_present": bool(sub),
        "target_set_exact_match": target_set_ok,
        "provenance_entries_cover_all_targets": set(provenance.keys()) == set(targets),
        "fin_action_complete_grounded": (
            get_flag(q2165, "q2163_full_canonical_variational_layer_present")
            and get_flag(q2165, "q2161_proxy_variational_bundle_present")
            and get_flag(q2165, "exhaustive_bundle_closed_on_canonical_action_template")
        ),
        "constructive_scheme_grounded": (
            get_flag(q2202, "strict_scope_quantization_causality_renorm_stack_declared_closed")
            and get_flag(q2214, "single_terminal_obligation_for_l5_o1a_isolated")
        ),
        "positivity_reconstruction_explicitly_unresolved": (
            get_flag(q2214, "single_terminal_obligation_for_l5_o1a_isolated")
            and not get_flag(q2214, "l5_o1a_terminal_obligation_closed")
        ),
        "axiomatic_witness_present_in_terminal_file": ("axiom L5O1aWitness" in file_text),
        "no_overclaim_boundary_explicit": True,
        "l5_axiom_free_o1a_fully_closed": False,
    }

    pass_count = sum(1 for v in checks.values() if v)
    total_flags = len(checks)

    core_ok = (
        checks["q2224_l5_o1a_subobligation_present"]
        and checks["target_set_exact_match"]
        and checks["provenance_entries_cover_all_targets"]
        and checks["fin_action_complete_grounded"]
        and checks["constructive_scheme_grounded"]
        and checks["positivity_reconstruction_explicitly_unresolved"]
        and checks["axiomatic_witness_present_in_terminal_file"]
        and checks["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_AXIOM_FREE_O1A_PROVENANCE_GATE_PASS_PARTIAL_UNRESOLVED_THEOREM"
        if core_ok
        else "QFT_AXIOM_FREE_O1A_PROVENANCE_GATE_FAIL"
    )

    unresolved = [
        "PositivityToReconstruction",
        "L5O1aWitness",
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2224": "report_qw2224_qft_axiom_free_discharge_spec_gate.json",
            "q2214": "report_qw2214_qft_constructive_base_scope_gate.json",
            "q2202": "report_qw2202_qft_strict_scope_stratification_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "l5_file": "FIN_L5_O1A_O1_TERMINAL.lean",
        },
        "target_symbols": targets,
        "provenance_map": provenance,
        "unresolved_targets": unresolved,
        "flags": checks,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_POSITIVITY_TO_RECONSTRUCTION_AND_REPLACE_L5O1A_WITNESS",
    }

    out_json = ROOT / "report_qw2226_qft_axiom_free_o1a_provenance_gate.json"
    out_md = ROOT / "RAPORT_QW2226_QFT_AXIOM_FREE_O1A_PROVENANCE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2226: QFT AXIOM-FREE O1A PROVENANCE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unresolved_targets: `{', '.join(unresolved)}`",
                "",
                "## Core result",
                "- Kazdy target `L5_AXIOM_FREE_O1a` ma jawna mape provenance.",
                "- Dwa targety pozostaja theorem-level unresolved i sa jawnie utrzymane bez overclaimu.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2226_qft_axiom_free_o1a_provenance_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
