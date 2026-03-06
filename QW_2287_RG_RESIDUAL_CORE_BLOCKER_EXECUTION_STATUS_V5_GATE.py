#!/usr/bin/env python3
"""QW-2287: RG residual core-blocker execution status v5 gate.

Adds logical-nonderivability evidence on top of v4 to classify remaining
obligation as a single nonlogical physics-level derivation requirement.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
EXPECTED_EXPORT = "RG_CanonicalAction_to_WellPosedness_EXPORT"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2283 = load("report_qw2283_rg_residual_core_blocker_execution_status_v4_gate.json")
    q2285 = load("report_qw2285_rg_export_provider_logical_nonderivability_gate.json")

    n_total = int(q2283.get("n_obligations_total", 0))
    n_satisfied_v4 = int(q2283.get("n_obligations_satisfied_strict_v4", 0))
    isolated = q2283.get("isolated_unknown_identifiers", [])
    nonderivable = q2285.get("tautology") is False

    n_satisfied_v5 = 1 if n_satisfied_v4 == 1 else 0

    flags = {
        "q2283_v4_present": q2283.get("verdict")
        == "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V4_GATE_PASS_PARTIAL_SINGLE_SYMBOL_MINIMAL_OBSTRUCTION",
        "q2285_logical_nonderivability_present": q2285.get("verdict")
        == "RG_EXPORT_PROVIDER_LOGICAL_NONDERIVABILITY_GATE_PASS_OBSTRUCTION_FORMALLY_PROVED",
        "single_symbol_isolated": isolated == [EXPECTED_EXPORT],
        "logical_empty_context_nonderivable": nonderivable,
        "remaining_obligation_classified_nonlogical_single_symbol": isolated == [EXPECTED_EXPORT] and nonderivable,
        "all_obligations_satisfied_strict_v5": n_satisfied_v5 == n_total and n_total > 0,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V5_GATE_PASS_PARTIAL_SINGLE_NONLOGICAL_OBLIGATION"
        if flags["q2283_v4_present"] and flags["q2285_logical_nonderivability_present"]
        else "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V5_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2283": "report_qw2283_rg_residual_core_blocker_execution_status_v4_gate.json",
            "q2285": "report_qw2285_rg_export_provider_logical_nonderivability_gate.json",
        },
        "n_obligations_total": n_total,
        "n_obligations_satisfied_strict_v4": n_satisfied_v4,
        "n_obligations_satisfied_strict_v5": n_satisfied_v5,
        "isolated_unknown_identifiers": isolated,
        "tautology": q2285.get("tautology"),
    }
    proof_path = ROOT / "proof_object_qw2287_rg_residual_core_blocker_execution_status_v5.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_obligations_total": n_total,
        "n_obligations_satisfied_strict_v5": n_satisfied_v5,
        "isolated_unknown_identifiers": isolated,
        "tautology": q2285.get("tautology"),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "CONSTRUCT_NONLOGICAL_PHYSICS_DERIVATION_FOR_RG_EXPORT_PROVIDER",
    }

    out_json = ROOT / "report_qw2287_rg_residual_core_blocker_execution_status_v5_gate.json"
    out_md = ROOT / "RAPORT_QW2287_RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V5_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2287: RG RESIDUAL CORE BLOCKER EXECUTION STATUS V5 GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- obligations satisfied v5: `{n_satisfied_v5}/{n_total}`",
                f"- isolated_unknown_identifiers: `{isolated}`",
                f"- tautology: `{q2285.get('tautology')}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations_satisfied_strict_v5": n_satisfied_v5, "n_obligations_total": n_total}))


if __name__ == "__main__":
    main()
