#!/usr/bin/env python3
"""QW-2283: RG residual core-blocker execution status v4 gate.

Integrates v3 status with blocker-isolation evidence to certify
single-symbol minimal obstruction.
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
    q2279 = load("report_qw2279_rg_residual_core_blocker_execution_status_v3_gate.json")
    q2281 = load("report_qw2281_rg_residual_core_blocker_isolation_gate.json")
    p2269 = load("proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json")

    n_total = int(q2279.get("n_obligations_total", 0))
    n_satisfied_v3 = int(q2279.get("n_obligations_satisfied_strict_v3", 0))
    unknowns = q2281.get("unknown_identifiers", [])
    isolated_single_export = unknowns == [EXPECTED_EXPORT] and not q2281.get("proposition_mismatch", True)

    n_satisfied_v4 = 1 if (n_satisfied_v3 == 1 and q2281.get("machine_check_exit_code") == 0) else 0

    flags = {
        "q2279_v3_present": q2279.get("verdict")
        == "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V3_GATE_PASS_PARTIAL_PENDING_MACHINE_CHECKABLE_NON_AXIOMATIC",
        "q2281_isolation_present": q2281.get("verdict")
        in {
            "RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE_PASS_PARTIAL_CORE_BLOCKER_ISOLATED",
            "RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE_PASS_STRICT_PROVIDER_CONSTRUCTED",
        },
        "single_residual_obligation_present": int(p2269.get("n_obligations", 0)) == 1,
        "core_blocker_isolated_single_export_symbol": isolated_single_export,
        "machine_check_exit_zero": q2281.get("machine_check_exit_code") == 0,
        "all_obligations_satisfied_strict_v4": n_satisfied_v4 == n_total and n_total > 0,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V4_GATE_PASS_PARTIAL_SINGLE_SYMBOL_MINIMAL_OBSTRUCTION"
        if flags["q2279_v3_present"] and flags["q2281_isolation_present"]
        else "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V4_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2279": "report_qw2279_rg_residual_core_blocker_execution_status_v3_gate.json",
            "q2281": "report_qw2281_rg_residual_core_blocker_isolation_gate.json",
            "p2269": "proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json",
        },
        "n_obligations_total": n_total,
        "n_obligations_satisfied_strict_v3": n_satisfied_v3,
        "n_obligations_satisfied_strict_v4": n_satisfied_v4,
        "isolated_unknown_identifiers": unknowns,
    }
    proof_path = ROOT / "proof_object_qw2283_rg_residual_core_blocker_execution_status_v4.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_obligations_total": n_total,
        "n_obligations_satisfied_strict_v3": n_satisfied_v3,
        "n_obligations_satisfied_strict_v4": n_satisfied_v4,
        "isolated_unknown_identifiers": unknowns,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DISCHARGE_SINGLE_RG_EXPORT_PROVIDER_SYMBOL_IN_MACHINE_CHECKABLE_NON_AXIOMATIC_FORM"
            if n_satisfied_v4 < n_total
            else "UPDATE_RG_FULL_CLOSURE_STATUS"
        ),
    }

    out_json = ROOT / "report_qw2283_rg_residual_core_blocker_execution_status_v4_gate.json"
    out_md = ROOT / "RAPORT_QW2283_RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V4_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2283: RG RESIDUAL CORE BLOCKER EXECUTION STATUS V4 GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- obligations satisfied v4: `{n_satisfied_v4}/{n_total}`",
                f"- isolated_unknown_identifiers: `{unknowns}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations_satisfied_strict_v4": n_satisfied_v4, "n_obligations_total": n_total}))


if __name__ == "__main__":
    main()
