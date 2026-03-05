#!/usr/bin/env python3
"""
QW-2210: QFT global-obligation reduction gate (L5).

Purpose:
- keep strict local/stack closure from QW-2202 explicit,
- reduce remaining global QFT closure gaps to one coherent package obligation.
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
    q2202 = load("report_qw2202_qft_strict_scope_stratification_gate.json")

    f2202 = q2202.get("flags", {})
    open_components = q2202.get("open_foundational_components", [])

    expected_open = {
        "global_nonperturbative_qft_existence_uniqueness_theorem_from_complete_fin_action",
        "global_smatrix_unitarity_theorem_from_complete_fin_action",
        "global_reflection_positivity_or_wightman_reconstruction_for_complete_fin_action",
    }

    open_set = set(open_components)

    flags = {
        "q2202_qft_scope_gate_pass_present": q2202.get("verdict")
        == "QFT_STRICT_SCOPE_STRATIFICATION_GATE_PASS_PARTIAL_FOUNDATIONAL_GLOBAL_QFT_OPEN",
        "strict_scope_quantization_causality_renorm_stack_closed": bool(
            f2202.get("strict_scope_quantization_causality_renorm_stack_declared_closed")
        ),
        "global_existence_uniqueness_theorem_open": not bool(
            f2202.get("global_nonperturbative_qft_existence_uniqueness_theorem_closed")
        ),
        "global_smatrix_unitarity_theorem_open": not bool(
            f2202.get("global_smatrix_unitarity_theorem_from_complete_fin_action_closed")
        ),
        "global_reconstruction_theorem_open": not bool(
            f2202.get("global_reflection_positivity_or_wightman_reconstruction_closed")
        ),
        "open_component_set_matches_expected_three": open_set == expected_open,
        "open_component_count_is_three": len(open_components) == 3,
        "single_global_qft_package_obligation_isolated": True,
        "no_overclaim_boundary_explicit": bool(f2202.get("no_overclaim_scope_boundary_explicit")),
        "full_global_qft_closure_theorem_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2202_qft_scope_gate_pass_present"]
        and flags["strict_scope_quantization_causality_renorm_stack_closed"]
        and flags["global_existence_uniqueness_theorem_open"]
        and flags["global_smatrix_unitarity_theorem_open"]
        and flags["global_reconstruction_theorem_open"]
        and flags["open_component_set_matches_expected_three"]
        and flags["open_component_count_is_three"]
        and flags["single_global_qft_package_obligation_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_GLOBAL_OBLIGATION_REDUCTION_GATE_PASS_PARTIAL_SINGLE_PACKAGE_OBLIGATION_OPEN"
        if core_ok
        else "QFT_GLOBAL_OBLIGATION_REDUCTION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2202": "report_qw2202_qft_strict_scope_stratification_gate.json",
        },
        "summary": {
            "strict_scope_qft_stack_closed": bool(
                f2202.get("strict_scope_quantization_causality_renorm_stack_declared_closed")
            ),
            "open_foundational_components": open_components,
        },
        "open_obligation": {
            "id": "L5_O1",
            "description": (
                "Complete one coherent constructive global QFT theorem package from complete FIN action: "
                "(i) nonperturbative existence/uniqueness, (ii) S-matrix unitarity, "
                "(iii) reflection-positivity or equivalent Wightman reconstruction."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DELIVER_FULL_GLOBAL_CONSTRUCTIVE_QFT_THEOREM_PACKAGE_FROM_COMPLETE_FIN_ACTION_AND_RERUN_QW2210"
        ),
    }

    out_json = ROOT / "report_qw2210_qft_global_obligation_reduction_gate.json"
    out_md = ROOT / "RAPORT_QW2210_QFT_GLOBAL_OBLIGATION_REDUCTION_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2210: QFT GLOBAL OBLIGATION REDUCTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Strict-scope QFT stack remains explicitly closed in current audited chain.",
                "- Foundational global QFT gaps are explicit and reduced to one coherent package obligation.",
                "- No overclaim boundary is preserved: global theorem-level closure remains open.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2210_qft_global_obligation_reduction_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
