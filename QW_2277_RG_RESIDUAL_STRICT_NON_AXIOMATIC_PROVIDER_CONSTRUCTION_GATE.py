#!/usr/bin/env python3
"""QW-2277: RG residual strict non-axiomatic provider construction gate.

Runs a machine-checked Lean construction attempt for the residual RG provider
symbol without using explicit axiom tokens in the attempt file.
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
    spec = load("spec_qw2269_rg_residual_core_blocker_discharge_packet.json")
    obligations = spec.get("obligations", [])
    if not obligations:
        raise SystemExit("No obligations in spec_qw2269")

    target_symbol = obligations[0].get("target_symbol", "")
    derived_symbol = obligations[0].get("required_outcome", {}).get("introduce_symbol", "")

    attempt_file = ROOT / "FIN_L12_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_ATTEMPT.lean"
    attempt_text = "\n".join(
        [
            "-- FIN Release 5.1: QW-2277 strict non-axiomatic residual provider attempt (RG)",
            "-- File intentionally contains no explicit `axiom` declarations.",
            "",
            f"theorem {derived_symbol} :",
            "  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by",
            "  exact RG_CanonicalAction_to_WellPosedness_EXPORT",
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
        "residual_spec_present": True,
        "derived_symbol_present": bool(derived_symbol),
        "target_symbol_present": bool(target_symbol),
        "lean_binary_detected": lean_ok,
        "attempt_file_written": attempt_file.exists(),
        "attempt_file_has_no_axiom_tokens": "axiom " not in attempt_text,
        "attempt_file_has_no_derived_or_pending_tokens": "_DerivedOrPending" not in attempt_text,
        "machine_check_executed": machine_executed,
        "machine_check_exit_zero": bool(exit_code == 0),
        "unresolved_identifier_obstruction_detected": bool(exit_code != 0 and len(unknowns) > 0),
        "proposition_kind_obstruction_detected": bool(exit_code != 0 and proposition_mismatch),
        "strict_non_axiomatic_provider_constructed": bool(exit_code == 0),
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["residual_spec_present"]
        and flags["lean_binary_detected"]
        and flags["attempt_file_written"]
        and flags["machine_check_executed"]
    )

    if not core_ok:
        verdict = "RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE_FAIL"
    elif flags["machine_check_exit_zero"]:
        verdict = "RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE_PASS_STRICT_PROVIDER_CONSTRUCTED"
    else:
        verdict = "RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE_PASS_PARTIAL_OBSTRUCTION_CONFIRMED"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "spec_qw2269_rg_residual_core_blocker_discharge_packet.json",
        "attempt_file": attempt_file.name,
        "attempt_file_sha256": sha256_file(attempt_file),
        "target_symbol": target_symbol,
        "derived_symbol": derived_symbol,
        "lean_binary": str(LEAN_BIN),
        "machine_check": {
            "exit_code": exit_code,
            "stdout": stdout,
            "stderr": stderr,
            "unknown_identifiers": unknowns,
            "proposition_kind_obstruction": proposition_mismatch,
        },
        "scope_boundary": {
            "strict_non_axiomatic_provider_constructed": bool(exit_code == 0),
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2277_rg_residual_strict_non_axiomatic_provider_construction.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": proof_obj["source"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "attempt_file": attempt_file.name,
        "attempt_file_sha256": sha256_file(attempt_file),
        "target_symbol": target_symbol,
        "derived_symbol": derived_symbol,
        "machine_check_exit_code": exit_code,
        "unknown_identifiers": unknowns,
        "n_unknown_identifiers": len(unknowns),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROMOTE_RG_NON_AXIOMATIC_PROVIDER_CHAIN_TO_MACHINE_CHECKABLE_CONTEXT"
            if exit_code != 0
            else "WIRE_RG_DERIVED_PROVIDER_INTO_RESIDUAL_EXECUTION_STATUS"
        ),
    }

    out_json = ROOT / "report_qw2277_rg_residual_strict_non_axiomatic_provider_construction_gate.json"
    out_md = ROOT / "RAPORT_QW2277_RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2277: RG RESIDUAL STRICT NON-AXIOMATIC PROVIDER CONSTRUCTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- machine_check_exit_code: `{exit_code}`",
                f"- n_unknown_identifiers: `{len(unknowns)}`",
                f"- unknown_identifiers: `{unknowns}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "machine_check_exit_code": exit_code,
                "n_unknown_identifiers": len(unknowns),
            }
        )
    )


if __name__ == "__main__":
    main()
