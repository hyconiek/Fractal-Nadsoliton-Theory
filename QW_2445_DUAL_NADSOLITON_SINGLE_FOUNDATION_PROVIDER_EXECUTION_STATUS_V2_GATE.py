#!/usr/bin/env python3
"""QW-2445: dual Nadsoliton single-foundation provider execution-status V2 gate.

Runs provider attempts using runtime discovered in QW-2444.
Keeps strict boundaries: no fake full-closure claim.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

L12_ATTEMPT = ROOT / "FIN_L12_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_ATTEMPT_V2.lean"
L5_ATTEMPT = ROOT / "FIN_L5_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_ATTEMPT_V2.lean"

L12_EXPECTED = "RG_CanonicalAction_to_WellPosedness_EXPORT"
L5_EXPECTED = "QFT_CanonicalAction_to_Positivity_EXPORT"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def parse_unknowns(text: str) -> list[str]:
    pats = [
        re.compile(r"unknown identifier '([^']+)'"),
        re.compile(r"unknown constant '([^']+)'"),
        re.compile(r"Unknown identifier `([^`]+)`"),
    ]
    out: set[str] = set()
    for pat in pats:
        out.update(pat.findall(text))
    return sorted(out)


def run_lean(lean_bin: str, path: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["ELAN_HOME"] = str(ROOT / ".elan")
    env["HOME"] = str(ROOT / ".home_lean")
    env["PATH"] = f"{ROOT / '.elan/bin'}:{env.get('PATH', '')}"
    proc = subprocess.run([lean_bin, str(path.name)], cwd=ROOT, capture_output=True, text=True, check=False, env=env)
    return {
        "file": path.name,
        "lean_bin": lean_bin,
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "unknown_identifiers": parse_unknowns(proc.stdout + "\n" + proc.stderr),
    }


def main() -> None:
    q2441 = load("report_qw2441_dual_nadsoliton_single_foundation_export_provider_packet_gate.json")
    q2444 = load("report_qw2444_lean_runtime_discovery_gate.json")

    lean_bin = q2444.get("selected_runtime")

    l12_code = "\n".join(
        [
            "-- QW-2445 strict attempt V2: single-foundation -> RG canonical export",
            "set_option autoImplicit false",
            "variable (NadsolitonSingleFoundation RG_WellPosedness_Target : Prop)",
            "",
            "theorem RG_NadsolitonSingleFoundationToWellPosedness_EXPORT_V2 :",
            "  NadsolitonSingleFoundation ->",
            "  RG_WellPosedness_Target := by",
            "  intro hFoundation",
            "  exact RG_CanonicalAction_to_WellPosedness_EXPORT",
            "",
        ]
    )
    L12_ATTEMPT.write_text(l12_code, encoding="utf-8")

    l5_code = "\n".join(
        [
            "-- QW-2445 strict attempt V2: single-foundation -> QFT canonical export",
            "set_option autoImplicit false",
            "variable (NadsolitonSingleFoundation QFT_Positivity_Target : Prop)",
            "",
            "theorem QFT_NadsolitonSingleFoundationToPositivity_EXPORT_V2 :",
            "  NadsolitonSingleFoundation ->",
            "  QFT_Positivity_Target := by",
            "  intro hFoundation",
            "  exact QFT_CanonicalAction_to_Positivity_EXPORT",
            "",
        ]
    )
    L5_ATTEMPT.write_text(l5_code, encoding="utf-8")

    l12_run: dict[str, Any] | None = None
    l5_run: dict[str, Any] | None = None

    if lean_bin:
        l12_run = run_lean(str(lean_bin), L12_ATTEMPT)
        l5_run = run_lean(str(lean_bin), L5_ATTEMPT)

    flags = {
        "q2441_packet_ready": q2441.get("verdict")
        == "DUAL_NADSOLITON_SINGLE_FOUNDATION_EXPORT_PROVIDER_PACKET_GATE_PASS_PACKET_READY",
        "q2444_runtime_available": q2444.get("verdict") == "LEAN_RUNTIME_DISCOVERY_GATE_PASS_RUNTIME_AVAILABLE",
        "execution_attempt_performed": l12_run is not None and l5_run is not None,
        "l12_expected_unknown_symbol_isolated": bool(l12_run and L12_EXPECTED in set(l12_run.get("unknown_identifiers", []))),
        "l5_expected_unknown_symbol_isolated": bool(l5_run and L5_EXPECTED in set(l5_run.get("unknown_identifiers", []))),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if not flags["q2444_runtime_available"]:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"
        required_next_step = "ATTACH_OR_INSTALL_LEAN_RUNTIME_AND_RERUN_QW2445"
    elif flags["execution_attempt_performed"] and flags["l12_expected_unknown_symbol_isolated"] and flags["l5_expected_unknown_symbol_isolated"]:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_CANONICAL_EXPORT_SYMBOLS"
        required_next_step = "EXTRACT_DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT"
    else:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE_FAIL"
        required_next_step = "REPAIR_QW2445_EXECUTION_PIPELINE"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2441": "report_qw2441_dual_nadsoliton_single_foundation_export_provider_packet_gate.json",
            "q2444": "report_qw2444_lean_runtime_discovery_gate.json",
            "l12_attempt": L12_ATTEMPT.name,
            "l5_attempt": L5_ATTEMPT.name,
        },
        "runs": {
            "l12": l12_run,
            "l5": l5_run,
        },
        "scope_boundary": {
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "unknown_identifiers": {
            "l12": (l12_run or {}).get("unknown_identifiers", []),
            "l5": (l5_run or {}).get("unknown_identifiers", []),
        },
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json"
    out_md = ROOT / "RAPORT_QW2445_DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2445: DUAL NADSOLITON SINGLE FOUNDATION PROVIDER EXECUTION STATUS V2 GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unknown_l12: `{out['unknown_identifiers']['l12']}`",
                f"- unknown_l5: `{out['unknown_identifiers']['l5']}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "unknown_identifiers": out["unknown_identifiers"]}))


if __name__ == "__main__":
    main()
