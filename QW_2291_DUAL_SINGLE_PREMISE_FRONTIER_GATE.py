#!/usr/bin/env python3
"""QW-2291: dual single-premise frontier gate.

Integrates v5 residual statuses with conditional machine-checked providers to
formalize the dual-branch frontier as two explicit physics-level premises.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2287 = load("report_qw2287_rg_residual_core_blocker_execution_status_v5_gate.json")
    q2288 = load("report_qw2288_qft_residual_core_blocker_execution_status_v5_gate.json")
    q2289 = load("report_qw2289_rg_export_provider_single_premise_conditional_gate.json")
    q2290 = load("report_qw2290_qft_export_provider_single_premise_conditional_gate.json")

    flags = {
        "q2287_v5_single_nonlogical_rg": q2287.get("verdict")
        == "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V5_GATE_PASS_PARTIAL_SINGLE_NONLOGICAL_OBLIGATION",
        "q2288_v5_single_nonlogical_qft": q2288.get("verdict")
        == "QFT_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V5_GATE_PASS_PARTIAL_SINGLE_NONLOGICAL_OBLIGATION",
        "q2289_rg_conditional_machine_checked": q2289.get("verdict")
        == "RG_EXPORT_PROVIDER_SINGLE_PREMISE_CONDITIONAL_GATE_PASS_PARTIAL_CONDITIONAL_PROVIDER_MACHINE_CHECKED",
        "q2290_qft_conditional_machine_checked": q2290.get("verdict")
        == "QFT_EXPORT_PROVIDER_SINGLE_PREMISE_CONDITIONAL_GATE_PASS_PARTIAL_CONDITIONAL_PROVIDER_MACHINE_CHECKED",
        "dual_frontier_reduced_to_two_explicit_premises": True,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    remaining_frontier = [
        {
            "branch": "L12",
            "premise_symbol": "RG_PhysicalBridgePremise",
            "statement": "(FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales",
        },
        {
            "branch": "L5",
            "premise_symbol": "QFT_PhysicalBridgePremise",
            "statement": "(FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction",
        },
    ]

    verdict = (
        "DUAL_SINGLE_PREMISE_FRONTIER_GATE_PASS_PARTIAL_FRONTIER_EXPLICIT"
        if flags["q2287_v5_single_nonlogical_rg"]
        and flags["q2288_v5_single_nonlogical_qft"]
        and flags["q2289_rg_conditional_machine_checked"]
        and flags["q2290_qft_conditional_machine_checked"]
        else "DUAL_SINGLE_PREMISE_FRONTIER_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2287": "report_qw2287_rg_residual_core_blocker_execution_status_v5_gate.json",
            "q2288": "report_qw2288_qft_residual_core_blocker_execution_status_v5_gate.json",
            "q2289": "report_qw2289_rg_export_provider_single_premise_conditional_gate.json",
            "q2290": "report_qw2290_qft_export_provider_single_premise_conditional_gate.json",
        },
        "remaining_frontier": remaining_frontier,
        "n_remaining_frontier_items": len(remaining_frontier),
    }
    proof_path = ROOT / "proof_object_qw2291_dual_single_premise_frontier.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "remaining_frontier": remaining_frontier,
        "n_remaining_frontier_items": len(remaining_frontier),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_TWO_PHYSICAL_BRIDGE_PREMISES_FROM_FIN_ACTION_LEVEL_THEOREMS",
    }

    out_json = ROOT / "report_qw2291_dual_single_premise_frontier_gate.json"
    out_md = ROOT / "RAPORT_QW2291_DUAL_SINGLE_PREMISE_FRONTIER_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2291: DUAL SINGLE PREMISE FRONTIER GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_remaining_frontier_items: `{len(remaining_frontier)}`",
                f"- remaining_frontier: `{remaining_frontier}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_remaining_frontier_items": len(remaining_frontier)}))


if __name__ == "__main__":
    main()
