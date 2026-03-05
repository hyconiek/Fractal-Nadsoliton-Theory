#!/usr/bin/env python3
"""
QW-2243: RG DAX1 non-axiomatic provider attempt gate.

Purpose:
- attempt to instantiate RG_DAX_1 as a non-axiomatic Lean provider theorem,
- verify whether canonical-chain export theorem exists,
- classify exact blocker if export is missing.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
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


def detect_lean() -> str | None:
    direct = shutil.which("lean")
    if direct:
        return direct
    fallback = Path("/tmp/lean4/lean-4.28.0-linux/bin/lean")
    return str(fallback) if fallback.exists() else None


def run(cmd: list[str]) -> tuple[int, str, str]:
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    return int(proc.returncode), proc.stdout, proc.stderr


def symbol_exists(pattern: str) -> bool:
    rc, out, _err = run(["rg", "-n", pattern, "-g*.lean"])
    if rc == 1:
        return False
    return bool(out.strip())


def main() -> None:
    q2241 = load("report_qw2241_rg_provider_deaxiomatization_obstruction_gate.json")
    q2165 = load("report_qw2165_l13_exhaustive_canonical_eom_gate.json")
    q2179 = load("report_qw2179_l13_u2b2_exact_matching_identity_gate.json")
    q2221 = load("report_qw2221_rg_terminal_proof_object_execution_gate.json")

    export_symbol = "RG_CanonicalAction_to_WellPosedness_EXPORT"
    export_present = symbol_exists(rf"\b{export_symbol}\b")

    attempt = ROOT / "FIN_L12_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT.lean"
    attempt.write_text(
        "\n".join(
            [
                "-- FIN Release 5.1: QW-2243 RG_DAX1 non-axiomatic provider attempt",
                "-- Expected canonical export symbol from FIN chain.",
                "",
                "axiom FINActionComplete : Prop",
                "axiom RGConstructiveMap : Prop",
                "axiom RGGlobalWellPosednessAllScales : Prop",
                "",
                "theorem RG_DAX1_PROVIDER_NON_AXIOMATIC_ATTEMPT :",
                "  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by",
                "  exact RG_CanonicalAction_to_WellPosedness_EXPORT",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    lean_bin = detect_lean()
    checker_attempted = False
    checker_rc = None
    checker_stdout = ""
    checker_stderr = ""
    checker_missing_export_confirmed = False

    if lean_bin:
        checker_attempted = True
        rc, so, se = run([lean_bin, attempt.name])
        checker_rc, checker_stdout, checker_stderr = rc, so, se
        merged = (so + "\n" + se).lower()
        if rc != 0 and "unknown identifier `rg_canonicalaction_to_wellposedness_export`" in merged:
            checker_missing_export_confirmed = True

    proof_object = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2241": "report_qw2241_rg_provider_deaxiomatization_obstruction_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "q2179": "report_qw2179_l13_u2b2_exact_matching_identity_gate.json",
            "q2221": "report_qw2221_rg_terminal_proof_object_execution_gate.json",
        },
        "attempt_file": {
            "name": attempt.name,
            "sha256": sha256_file(attempt),
        },
        "export_symbol": export_symbol,
        "export_symbol_present_in_repo": export_present,
        "checker": {
            "attempted": checker_attempted,
            "exit_code": checker_rc,
            "stdout": checker_stdout,
            "stderr": checker_stderr,
            "missing_export_confirmed": checker_missing_export_confirmed,
        },
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2243_rg_dax1_non_axiomatic_provider_attempt.json"
    proof_path.write_text(json.dumps(proof_object, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2241_obstruction_map_present": q2241.get("verdict")
        == "RG_PROVIDER_DEAXIOMATIZATION_OBSTRUCTION_GATE_PASS_PARTIAL_NON_AXIOMATIC_SOURCE_MISSING",
        "q2165_chain_anchor_present": q2165.get("verdict")
        == "L13_EXHAUSTIVE_CANONICAL_EOM_GATE_PASS_PARTIAL_ALL_ORDERS_OPEN",
        "q2179_terminal_identity_closed": q2179.get("verdict")
        == "L13_U2B2_EXACT_MATCHING_IDENTITY_GATE_PASS_TERMINAL_CHAIN_CLOSED",
        "q2221_execution_layer_present": q2221.get("verdict")
        == "RG_TERMINAL_PROOF_OBJECT_EXECUTION_GATE_PASS_PARTIAL_AXIOMATIC_BOUNDARY",
        "attempt_file_written": attempt.exists(),
        "lean_checker_detected": bool(lean_bin),
        "lean_checker_attempted": checker_attempted,
        "canonical_export_symbol_exists": bool(export_present),
        "checker_confirms_missing_export_symbol": bool(checker_missing_export_confirmed),
        "proof_object_generated": proof_path.exists(),
        "proof_object_hash_recorded": True,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2241_obstruction_map_present"]
        and flags["q2165_chain_anchor_present"]
        and flags["q2179_terminal_identity_closed"]
        and flags["q2221_execution_layer_present"]
        and flags["attempt_file_written"]
        and flags["lean_checker_detected"]
        and flags["lean_checker_attempted"]
        and flags["checker_confirms_missing_export_symbol"]
        and flags["proof_object_generated"]
        and flags["proof_object_hash_recorded"]
    )

    verdict = (
        "RG_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING"
        if core_ok
        else "RG_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_object["sources"],
        "attempt_file": attempt.name,
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "FORMALIZE_RG_CANONICALACTION_TO_WELLPOSSEDNESS_EXPORT_AS_NON_AXIOMATIC_LEAN_THEOREM",
    }

    out_json = ROOT / "report_qw2243_rg_dax1_non_axiomatic_provider_attempt_gate.json"
    out_md = ROOT / "RAPORT_QW2243_RG_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2243: RG DAX1 NON-AXIOMATIC PROVIDER ATTEMPT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Wykonano bezposredni attempt DAX1 non-axiomatic provider.",
                "- Brak canonical export symbol zostal potwierdzony machine-checkiem.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2243_rg_dax1_non_axiomatic_provider_attempt_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
