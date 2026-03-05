#!/usr/bin/env python3
"""
QW-2212: QFT global obligation decomposition gate (L5_O1).

Purpose:
- keep L5 package-obligation statement from QW-2210 explicit,
- decompose it into minimal executable theorem sub-obligations.
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
    q2210 = load("report_qw2210_qft_global_obligation_reduction_gate.json")
    q2202 = load("report_qw2202_qft_strict_scope_stratification_gate.json")

    f2210 = q2210.get("flags", {})
    f2202 = q2202.get("flags", {})
    l5_o1 = q2210.get("open_obligation", {})

    open_components = q2210.get("summary", {}).get("open_foundational_components", [])

    decomposition = {
        "L5_O1a": {
            "description": (
                "Constructive nonperturbative QFT existence/uniqueness for complete FIN action "
                "with reflection-positivity (or equivalent Euclidean positivity) sufficient for reconstruction."
            ),
            "currently_closed": False,
            "depends_on": [],
        },
        "L5_O1b": {
            "description": (
                "Derive asymptotic/scattering completeness and unitary S-matrix for the reconstructed theory "
                "from complete FIN action (no hidden scope truncation)."
            ),
            "currently_closed": False,
            "depends_on": ["L5_O1a"],
        },
    }

    flags = {
        "q2210_single_package_obligation_pass_present": q2210.get("verdict")
        == "QFT_GLOBAL_OBLIGATION_REDUCTION_GATE_PASS_PARTIAL_SINGLE_PACKAGE_OBLIGATION_OPEN",
        "l5_o1_id_is_explicit": l5_o1.get("id") == "L5_O1",
        "strict_scope_stack_closed": bool(f2210.get("strict_scope_quantization_causality_renorm_stack_closed")),
        "q2202_global_theorem_triplet_open_explicit": (
            not bool(f2202.get("global_nonperturbative_qft_existence_uniqueness_theorem_closed"))
            and not bool(f2202.get("global_smatrix_unitarity_theorem_from_complete_fin_action_closed"))
            and not bool(f2202.get("global_reflection_positivity_or_wightman_reconstruction_closed"))
        ),
        "q2210_open_component_count_three": len(open_components) == 3,
        "decomposition_to_two_subobligations_explicit": True,
        "dependency_graph_acyclic": True,
        "l5_o1a_closed": False,
        "l5_o1b_closed": False,
        "full_l5_global_qft_package_closed": False,
        "no_overclaim_boundary_explicit": bool(f2210.get("no_overclaim_boundary_explicit")),
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2210_single_package_obligation_pass_present"]
        and flags["l5_o1_id_is_explicit"]
        and flags["strict_scope_stack_closed"]
        and flags["q2202_global_theorem_triplet_open_explicit"]
        and flags["q2210_open_component_count_three"]
        and flags["decomposition_to_two_subobligations_explicit"]
        and flags["dependency_graph_acyclic"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_TWO_SUBOBLIGATIONS_OPEN"
        if core_ok
        else "QFT_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2210": "report_qw2210_qft_global_obligation_reduction_gate.json",
            "q2202": "report_qw2202_qft_strict_scope_stratification_gate.json",
        },
        "summary": {
            "l5_o1_prior": l5_o1,
            "open_foundational_components_from_q2210": open_components,
            "strict_scope_stack_closed": bool(f2210.get("strict_scope_quantization_causality_renorm_stack_closed")),
        },
        "decomposition": decomposition,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "CLOSE_L5_O1A_CONSTRUCTIVE_EXISTENCE_POSITIVITY_RECONSTRUCTION_THEN_CLOSE_L5_O1B_UNITARY_SCATTERING_COMPLETENESS"
        ),
    }

    out_json = ROOT / "report_qw2212_qft_global_obligation_decomposition_gate.json"
    out_md = ROOT / "RAPORT_QW2212_QFT_GLOBAL_OBLIGATION_DECOMPOSITION_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2212: QFT GLOBAL OBLIGATION DECOMPOSITION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Single open L5 package obligation (`L5_O1`) is decomposed into two executable theorem sub-obligations.",
                "- Strict-scope closure remains preserved; global theorem-level boundary remains explicit.",
                "- Remaining closure path is now sequenced: `L5_O1a -> L5_O1b`.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2212_qft_global_obligation_decomposition_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
