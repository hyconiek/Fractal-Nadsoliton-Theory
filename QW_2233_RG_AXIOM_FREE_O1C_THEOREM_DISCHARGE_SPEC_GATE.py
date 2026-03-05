#!/usr/bin/env python3
"""
QW-2233: RG axiom-free O1c theorem-discharge spec gate.

Purpose:
- turn the remaining C1 theorem-discharge into explicit proof obligations,
- export a strict dependency DAG and acceptance criteria,
- keep closure boundary explicit (no overclaim).
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
    q2231 = load("report_qw2231_rg_axiom_free_o1c_execution_step_gate.json")
    q2229 = load("report_qw2229_rg_axiom_free_o1c_attachment_spec_gate.json")

    cand_a = ROOT / "FIN_L12_O1A_O1_O1C_STEP.lean"
    cand_b = ROOT / "FIN_L12_O1B_O1_O1C_STEP.lean"

    text_a = cand_a.read_text(encoding="utf-8") if cand_a.exists() else ""
    text_b = cand_b.read_text(encoding="utf-8") if cand_b.exists() else ""

    derived_or_pending_present = ("DerivedOrPending" in text_a) and ("DerivedOrPending" in text_b)

    obligations = {
        "scope": "L12_AXIOM_FREE_O1C",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "obligations": [
            {
                "id": "RG_C1_1",
                "statement": "(FINActionComplete AND RGConstructiveMap) -> RGGlobalWellPosednessAllScales",
                "source_symbol": "RGGlobalWellPosednessAllScales_DerivedOrPending",
                "depends_on": [],
            },
            {
                "id": "RG_C1_2",
                "statement": "RGGlobalWellPosednessAllScales -> RGGlobalFixedPointStabilityAllT",
                "source_symbol": "RGGlobalFixedPointStabilityAllT_DerivedOrPending",
                "depends_on": ["RG_C1_1"],
            },
            {
                "id": "RG_C1_3",
                "statement": "compose RG_C1_1 + RG_C1_2 into final O1c theorem attachment",
                "source_symbol": "L12_O1A_O1_O1C_STEP + L12_O1B_O1_O1C_STEP",
                "depends_on": ["RG_C1_1", "RG_C1_2"],
            },
        ],
        "acceptance_criteria": {
            "D1": "all obligations have machine-checkable theorem statements",
            "D2": "dependency DAG is acyclic",
            "D3": "all DerivedOrPending symbols have mapped replacement obligations",
            "D4": "final proof object hash is attached after discharge",
        },
    }

    obligations_path = ROOT / "spec_qw2233_rg_o1c_theorem_discharge_obligations.json"
    obligations_path.write_text(json.dumps(obligations, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_object = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": {
            "q2231": "report_qw2231_rg_axiom_free_o1c_execution_step_gate.json",
            "q2229": "report_qw2229_rg_axiom_free_o1c_attachment_spec_gate.json",
        },
        "candidate_files": {
            cand_a.name: sha256_file(cand_a) if cand_a.exists() else None,
            cand_b.name: sha256_file(cand_b) if cand_b.exists() else None,
        },
        "obligations_file": obligations_path.name,
        "obligations_file_sha256": sha256_file(obligations_path),
        "scope_boundary": {
            "theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2233_rg_o1c_theorem_discharge_spec.json"
    proof_path.write_text(json.dumps(proof_object, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2231_execution_step_pass_present": q2231.get("verdict")
        == "RG_AXIOM_FREE_O1C_EXECUTION_STEP_GATE_PASS_PARTIAL_WITNESS_REMOVED_THEOREM_PENDING",
        "q2229_attachment_spec_pass_present": q2229.get("verdict")
        == "RG_AXIOM_FREE_O1C_ATTACHMENT_SPEC_GATE_PASS_PARTIAL_DISCHARGE_PENDING",
        "candidate_files_present": cand_a.exists() and cand_b.exists(),
        "derived_or_pending_markers_present": bool(derived_or_pending_present),
        "obligations_file_written": obligations_path.exists(),
        "obligation_count_ge_3": len(obligations["obligations"]) >= 3,
        "dependency_dag_acyclic_declared": True,
        "acceptance_criteria_defined": True,
        "proof_object_generated": proof_path.exists(),
        "proof_object_hash_recorded": True,
        "c1_theorem_discharge_completed": False,
        "o1c_fully_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2231_execution_step_pass_present"]
        and flags["q2229_attachment_spec_pass_present"]
        and flags["candidate_files_present"]
        and flags["derived_or_pending_markers_present"]
        and flags["obligations_file_written"]
        and flags["obligation_count_ge_3"]
        and flags["dependency_dag_acyclic_declared"]
        and flags["acceptance_criteria_defined"]
        and flags["proof_object_generated"]
        and flags["proof_object_hash_recorded"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_AXIOM_FREE_O1C_THEOREM_DISCHARGE_SPEC_GATE_PASS_PARTIAL_PROOFS_PENDING"
        if core_ok
        else "RG_AXIOM_FREE_O1C_THEOREM_DISCHARGE_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2231": "report_qw2231_rg_axiom_free_o1c_execution_step_gate.json",
            "q2229": "report_qw2229_rg_axiom_free_o1c_attachment_spec_gate.json",
        },
        "obligations_file": obligations_path.name,
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DISCHARGE_RG_C1_1_RG_C1_2_RG_C1_3_AND_REPLACE_DERIVED_OR_PENDING_SYMBOLS",
    }

    out_json = ROOT / "report_qw2233_rg_axiom_free_o1c_theorem_discharge_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2233_RG_AXIOM_FREE_O1C_THEOREM_DISCHARGE_SPEC_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2233: RG AXIOM-FREE O1C THEOREM-DISCHARGE SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Sformalizowano C1 theorem-discharge do jawnych obligacji `RG_C1_1..RG_C1_3`.",
                "- Granica closure pozostaje jawna: proof obligations sa zdefiniowane, ale nieudowodnione.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2233_rg_axiom_free_o1c_theorem_discharge_spec_gate.json`",
                "- Obligations: `spec_qw2233_rg_o1c_theorem_discharge_obligations.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
