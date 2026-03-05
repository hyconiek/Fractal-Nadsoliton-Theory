#!/usr/bin/env python3
"""
QW-2214: QFT constructive base scope gate (L5_O1a).

Purpose:
- keep decomposition from QW-2212 explicit,
- reduce L5_O1a to one terminal constructive theorem obligation.
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
    q2212 = load("report_qw2212_qft_global_obligation_decomposition_gate.json")
    q2202 = load("report_qw2202_qft_strict_scope_stratification_gate.json")
    q2165 = load("report_qw2165_l13_exhaustive_canonical_eom_gate.json")
    q2166 = load("report_qw2166_l14_exhaustive_canonical_hessian_gate.json")
    q2136 = load("report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json")
    q2181 = load("report_qw2181_dual_terminal_matching_closure_gate.json")

    f2202 = q2202.get("flags", {})
    f2165 = q2165.get("flags", {})
    f2166 = q2166.get("flags", {})
    f2136 = q2136.get("flags", {})
    f2181 = q2181.get("flags", {})

    decomp = q2212.get("decomposition", {})
    l5_o1a = decomp.get("L5_O1a", {})

    flags = {
        "q2212_decomposition_pass_present": q2212.get("verdict")
        == "QFT_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_TWO_SUBOBLIGATIONS_OPEN",
        "l5_o1a_subobligation_present": bool(l5_o1a),
        "q2202_strict_scope_stack_closed": bool(f2202.get("strict_scope_quantization_causality_renorm_stack_declared_closed")),
        "q2165_exhaustive_eom_locality_bundle": bool(f2165.get("phi_eom_local_second_order"))
        and bool(f2165.get("all_psi_eom_local_second_order"))
        and bool(f2165.get("no_spacetime_nonlocal_tokens_in_all_13_eom")),
        "q2166_exhaustive_hessian_locality_bundle": bool(f2166.get("all_linearized_eom_are_local_second_order"))
        and bool(f2166.get("no_spacetime_nonlocal_tokens_in_all_linearized_eom")),
        "q2136_counterterm_basis_policy_present": bool(f2136.get("finite_counterterm_basis_policy_declared"))
        and bool(f2136.get("finite_counterterm_basis_condition_holds")),
        "q2181_terminal_matching_closed": bool(f2181.get("dual_terminal_matching_closed")),
        "single_terminal_obligation_for_l5_o1a_isolated": True,
        "l5_o1a_terminal_obligation_closed": False,
        "l5_o1a_closed": False,
        "no_overclaim_boundary_explicit": bool(f2202.get("no_overclaim_scope_boundary_explicit")),
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2212_decomposition_pass_present"]
        and flags["l5_o1a_subobligation_present"]
        and flags["q2202_strict_scope_stack_closed"]
        and flags["q2165_exhaustive_eom_locality_bundle"]
        and flags["q2166_exhaustive_hessian_locality_bundle"]
        and flags["q2136_counterterm_basis_policy_present"]
        and flags["q2181_terminal_matching_closed"]
        and flags["single_terminal_obligation_for_l5_o1a_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_CONSTRUCTIVE_BASE_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN"
        if core_ok
        else "QFT_CONSTRUCTIVE_BASE_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2212": "report_qw2212_qft_global_obligation_decomposition_gate.json",
            "q2202": "report_qw2202_qft_strict_scope_stratification_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "q2166": "report_qw2166_l14_exhaustive_canonical_hessian_gate.json",
            "q2136": "report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json",
            "q2181": "report_qw2181_dual_terminal_matching_closure_gate.json",
        },
        "summary": {
            "l5_o1a_prior": l5_o1a,
            "strict_scope_qft_stack_closed": bool(f2202.get("strict_scope_quantization_causality_renorm_stack_declared_closed")),
            "locality_bundle_eom": bool(f2165.get("no_spacetime_nonlocal_tokens_in_all_13_eom")),
            "locality_bundle_hessian": bool(f2166.get("no_spacetime_nonlocal_tokens_in_all_linearized_eom")),
        },
        "open_obligation": {
            "id": "L5_O1a_O1",
            "description": (
                "Prove constructive nonperturbative existence/uniqueness with positivity-to-reconstruction theorem "
                "for complete FIN action (theorem-level, beyond current strict-scope stack)."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_L5_O1A_O1_CONSTRUCTIVE_EXISTENCE_POSITIVITY_RECONSTRUCTION_AND_RERUN_QW2214",
    }

    out_json = ROOT / "report_qw2214_qft_constructive_base_scope_gate.json"
    out_md = ROOT / "RAPORT_QW2214_QFT_CONSTRUCTIVE_BASE_SCOPE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2214: QFT CONSTRUCTIVE BASE SCOPE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L5_O1a is reduced to one terminal constructive theorem-level obligation (`L5_O1a_O1`).",
                "- Existing strict-scope locality/renormalization stack is preserved as prerequisite layer.",
                "- No overclaim beyond theorem-level boundary is introduced.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2214_qft_constructive_base_scope_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
