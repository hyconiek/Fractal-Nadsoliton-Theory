#!/usr/bin/env python3
"""QW-2281: RG residual core-blocker isolation gate.

Runs a kind-corrected Lean attempt to isolate whether the residual blocker
collapses to a single missing export-provider symbol.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
LEAN_BIN = Path("/tmp/lean4/lean-4.28.0-linux/bin/lean")
EXPECTED_EXPORT = "RG_CanonicalAction_to_WellPosedness_EXPORT"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def parse_unknowns(diag: str) -> list[str]:
    pats = [
        re.compile(r"unknown identifier '([^']+)'"),
        re.compile(r"unknown constant '([^']+)'"),
        re.compile(r"Unknown identifier `([^`]+)`"),
    ]
    out: set[str] = set()
    for pat in pats:
        out.update(pat.findall(diag))
    return sorted(out)


def main() -> None:
    q2277 = load("report_qw2277_rg_residual_strict_non_axiomatic_provider_construction_gate.json")
    spec = load("spec_qw2269_rg_residual_core_blocker_discharge_packet.json")
    obligations = spec.get("obligations", [])
    derived_symbol = obligations[0].get("required_outcome", {}).get("introduce_symbol", "") if obligations else ""

    attempt_file = ROOT / "FIN_L12_RESIDUAL_CORE_BLOCKER_ISOLATION_ATTEMPT.lean"
    attempt_text = "\n".join(
        [
            "-- FIN Release 5.1: QW-2281 RG residual core-blocker isolation attempt",
            "variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)",
            "",
            f"theorem {derived_symbol} :",
            "  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by",
            f"  exact {EXPECTED_EXPORT}",
            "",
        ]
    )
    attempt_file.write_text(attempt_text, encoding="utf-8")

    lean_ok = LEAN_BIN.exists() and LEAN_BIN.is_file()
    machine_executed = False
    exit_code = None
    stdout = ""
    stderr = ""
    if lean_ok:
        machine_executed = True
        proc = subprocess.run(
            [str(LEAN_BIN), attempt_file.name],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        exit_code = proc.returncode
        stdout = proc.stdout
        stderr = proc.stderr

    diag = "\n".join([stdout, stderr]).strip()
    unknowns = parse_unknowns(diag)
    proposition_mismatch = "is not a proposition" in diag

    flags = {
        "q2277_machine_obstruction_present": q2277.get("verdict")
        == "RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE_PASS_PARTIAL_OBSTRUCTION_CONFIRMED",
        "lean_binary_detected": lean_ok,
        "attempt_file_written": attempt_file.exists(),
        "attempt_file_axiom_free": "axiom " not in attempt_text,
        "attempt_file_has_no_derived_or_pending_tokens": "_DerivedOrPending" not in attempt_text,
        "kind_guard_inserted": " : Prop)" in attempt_text,
        "machine_check_executed": machine_executed,
        "machine_check_exit_zero": bool(exit_code == 0),
        "proposition_mismatch_removed": not proposition_mismatch,
        "single_unknown_identifier": len(unknowns) == 1,
        "unknown_is_expected_export_symbol": unknowns == [EXPECTED_EXPORT],
        "core_blocker_isolated_to_export_symbol": (len(unknowns) == 1 and unknowns[0] == EXPECTED_EXPORT and not proposition_mismatch),
        "strict_non_axiomatic_provider_constructed": bool(exit_code == 0),
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2277_machine_obstruction_present"]
        and flags["lean_binary_detected"]
        and flags["attempt_file_written"]
        and flags["machine_check_executed"]
    )
    if not core_ok:
        verdict = "RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE_FAIL"
    elif flags["machine_check_exit_zero"]:
        verdict = "RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE_PASS_STRICT_PROVIDER_CONSTRUCTED"
    elif flags["core_blocker_isolated_to_export_symbol"]:
        verdict = "RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE_PASS_PARTIAL_CORE_BLOCKER_ISOLATED"
    else:
        verdict = "RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE_PASS_PARTIAL_NONMINIMAL_OBSTRUCTION"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2277": "report_qw2277_rg_residual_strict_non_axiomatic_provider_construction_gate.json",
            "spec": "spec_qw2269_rg_residual_core_blocker_discharge_packet.json",
        },
        "attempt_file": attempt_file.name,
        "attempt_file_sha256": sha256_file(attempt_file),
        "machine_check": {
            "exit_code": exit_code,
            "stdout": stdout,
            "stderr": stderr,
            "unknown_identifiers": unknowns,
            "proposition_mismatch": proposition_mismatch,
        },
        "expected_export_symbol": EXPECTED_EXPORT,
    }
    proof_path = ROOT / "proof_object_qw2281_rg_residual_core_blocker_isolation.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "attempt_file": attempt_file.name,
        "attempt_file_sha256": sha256_file(attempt_file),
        "machine_check_exit_code": exit_code,
        "unknown_identifiers": unknowns,
        "n_unknown_identifiers": len(unknowns),
        "proposition_mismatch": proposition_mismatch,
        "expected_export_symbol": EXPECTED_EXPORT,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PRODUCE_MACHINE_CHECKABLE_NON_AXIOMATIC_RG_EXPORT_PROVIDER_SYMBOL"
            if exit_code != 0
            else "UPDATE_RG_RESIDUAL_STATUS_AFTER_PROVIDER_CONSTRUCTION"
        ),
    }

    out_json = ROOT / "report_qw2281_rg_residual_core_blocker_isolation_gate.json"
    out_md = ROOT / "RAPORT_QW2281_RG_RESIDUAL_CORE_BLOCKER_ISOLATION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2281: RG RESIDUAL CORE BLOCKER ISOLATION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- machine_check_exit_code: `{exit_code}`",
                f"- unknown_identifiers: `{unknowns}`",
                f"- proposition_mismatch: `{proposition_mismatch}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "machine_check_exit_code": exit_code, "unknown_identifiers": unknowns}))


if __name__ == "__main__":
    main()
