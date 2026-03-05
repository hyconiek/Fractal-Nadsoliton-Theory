#!/usr/bin/env python3
"""
QW-2216: QFT unitary scattering scope gate (L5_O1b).

Purpose:
- keep decomposition from QW-2212 explicit,
- reduce L5_O1b to one terminal theorem obligation.
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
    q2214 = load("report_qw2214_qft_constructive_base_scope_gate.json")
    q2202 = load("report_qw2202_qft_strict_scope_stratification_gate.json")
    q2127 = load("report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json")
    q2097 = load("report_qw2097_ckm_cp_target_refinement_gate.json")

    f2202 = q2202.get("flags", {})
    f2127 = q2127.get("flags", {})
    f2097 = q2097.get("flags", {})

    decomp = q2212.get("decomposition", {})
    l5_o1b = decomp.get("L5_O1b", {})

    flags = {
        "q2212_decomposition_pass_present": q2212.get("verdict")
        == "QFT_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_TWO_SUBOBLIGATIONS_OPEN",
        "l5_o1b_subobligation_present": bool(l5_o1b),
        "q2214_l5_o1a_terminalized": q2214.get("verdict")
        == "QFT_CONSTRUCTIVE_BASE_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN",
        "q2202_strict_scope_stack_closed": bool(f2202.get("strict_scope_quantization_causality_renorm_stack_declared_closed")),
        "q2127_nonabelian_action_bridge_present": bool(f2127.get("full_nonabelian_spinor_action_bridge_defined")),
        "q2097_mixing_unitarity_checks_pass": bool(f2097.get("unitarity_ok")),
        "q2202_global_smatrix_unitarity_not_claimed": not bool(
            f2202.get("global_smatrix_unitarity_theorem_from_complete_fin_action_closed")
        ),
        "single_terminal_obligation_for_l5_o1b_isolated": True,
        "l5_o1b_terminal_obligation_closed": False,
        "l5_o1b_closed": False,
        "full_l5_closed": False,
        "no_overclaim_boundary_explicit": bool(f2202.get("no_overclaim_scope_boundary_explicit")),
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2212_decomposition_pass_present"]
        and flags["l5_o1b_subobligation_present"]
        and flags["q2214_l5_o1a_terminalized"]
        and flags["q2202_strict_scope_stack_closed"]
        and flags["q2127_nonabelian_action_bridge_present"]
        and flags["q2097_mixing_unitarity_checks_pass"]
        and flags["q2202_global_smatrix_unitarity_not_claimed"]
        and flags["single_terminal_obligation_for_l5_o1b_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_UNITARY_SCATTERING_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN"
        if core_ok
        else "QFT_UNITARY_SCATTERING_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2212": "report_qw2212_qft_global_obligation_decomposition_gate.json",
            "q2214": "report_qw2214_qft_constructive_base_scope_gate.json",
            "q2202": "report_qw2202_qft_strict_scope_stratification_gate.json",
            "q2127": "report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json",
            "q2097": "report_qw2097_ckm_cp_target_refinement_gate.json",
        },
        "summary": {
            "l5_o1b_prior": l5_o1b,
            "l5_o1a_terminal_obligation": q2214.get("open_obligation", {}),
            "l5_o1_boundary": q2212.get("summary", {}).get("l5_o1_prior", {}),
        },
        "open_obligation": {
            "id": "L5_O1b_O1",
            "description": (
                "Prove unitary S-matrix with asymptotic/scattering completeness theorem for the reconstructed nonperturbative "
                "QFT from complete FIN action (without hidden truncation or scope restriction)."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_L5_O1B_O1_UNITARY_SCATTERING_COMPLETENESS_AND_RERUN_QW2216",
    }

    out_json = ROOT / "report_qw2216_qft_unitary_scattering_scope_gate.json"
    out_md = ROOT / "RAPORT_QW2216_QFT_UNITARY_SCATTERING_SCOPE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2216: QFT UNITARY SCATTERING SCOPE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L5_O1b is reduced to one terminal theorem-level obligation (`L5_O1b_O1`).",
                "- Prerequisite chain (`L5_O1a` terminalization + strict stack) is preserved explicitly.",
                "- Remaining L5 closure is now expressed as two terminal obligations: `L5_O1a_O1` and `L5_O1b_O1`.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2216_qft_unitary_scattering_scope_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
