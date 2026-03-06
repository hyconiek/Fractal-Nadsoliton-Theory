#!/usr/bin/env python3
"""QW-2285: RG export-provider logical nonderivability gate.

Formalizes that the propositional skeleton
  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales
is not a tautology, so it cannot be derived from empty logical context alone.
"""

from __future__ import annotations

import hashlib
import itertools
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


def formula(a: bool, b: bool, c: bool) -> bool:
    return (not (a and b)) or c


def main() -> None:
    q2283 = load("report_qw2283_rg_residual_core_blocker_execution_status_v4_gate.json")

    table = []
    for a, b, c in itertools.product([False, True], repeat=3):
        value = formula(a, b, c)
        table.append(
            {
                "FINActionComplete": a,
                "RGConstructiveMap": b,
                "RGGlobalWellPosednessAllScales": c,
                "formula_value": value,
            }
        )

    countermodels = [row for row in table if not row["formula_value"]]
    tautology = len(countermodels) == 0
    principal_countermodel = countermodels[0] if countermodels else None

    flags = {
        "q2283_v4_present": q2283.get("verdict")
        == "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V4_GATE_PASS_PARTIAL_SINGLE_SYMBOL_MINIMAL_OBSTRUCTION",
        "truth_table_complete": len(table) == 8,
        "countermodel_exists": len(countermodels) > 0,
        "tautology_rejected": not tautology,
        "principal_countermodel_is_a1_b1_c0": principal_countermodel
        == {
            "FINActionComplete": True,
            "RGConstructiveMap": True,
            "RGGlobalWellPosednessAllScales": False,
            "formula_value": False,
        },
        "logical_empty_context_derivation_impossible": not tautology,
        "residual_obligation_remains_nonlogical": True,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_EXPORT_PROVIDER_LOGICAL_NONDERIVABILITY_GATE_PASS_OBSTRUCTION_FORMALLY_PROVED"
        if flags["q2283_v4_present"] and flags["tautology_rejected"]
        else "RG_EXPORT_PROVIDER_LOGICAL_NONDERIVABILITY_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2283_rg_residual_core_blocker_execution_status_v4_gate.json",
        "formula": "(FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales",
        "truth_table": table,
        "countermodels": countermodels,
        "principal_countermodel": principal_countermodel,
    }
    proof_path = ROOT / "proof_object_qw2285_rg_export_provider_logical_nonderivability.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": proof_obj["source"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "formula": proof_obj["formula"],
        "n_countermodels": len(countermodels),
        "principal_countermodel": principal_countermodel,
        "tautology": tautology,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "ADD_NONLOGICAL_PHYSICAL_DERIVATION_PREMISE_FOR_RG_EXPORT_PROVIDER",
    }

    out_json = ROOT / "report_qw2285_rg_export_provider_logical_nonderivability_gate.json"
    out_md = ROOT / "RAPORT_QW2285_RG_EXPORT_PROVIDER_LOGICAL_NONDERIVABILITY_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2285: RG EXPORT PROVIDER LOGICAL NONDERIVABILITY GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_countermodels: `{len(countermodels)}`",
                f"- principal_countermodel: `{principal_countermodel}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_countermodels": len(countermodels), "tautology": tautology}))


if __name__ == "__main__":
    main()
