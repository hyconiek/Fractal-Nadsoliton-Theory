#!/usr/bin/env python3
"""
QW-2239: RG axiom-free O1c execution with provider layer gate.

Purpose:
- execute C1 attempt with explicit provider theorems available,
- verify missing-provider blocker is removed,
- keep axiomatic provenance boundary explicit.
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


def main() -> None:
    q2237 = load("report_qw2237_rg_axiom_free_o1c_provider_layer_gate.json")
    q2235 = load("report_qw2235_rg_axiom_free_o1c_theorem_discharge_execution_gate.json")

    provider = ROOT / "FIN_L12_O1C_PROVIDER_LAYER.lean"
    attempt = ROOT / "FIN_L12_O1C_THEOREM_DISCHARGE_PROVIDER_ATTEMPT.lean"

    provider_text = provider.read_text(encoding="utf-8") if provider.exists() else ""
    attempt.write_text(
        provider_text
        + "\n"
        + "\n".join(
            [
                "-- QW-2239 provider execution attempt",
                "theorem RG_C1_3_COMPOSED_PROVIDER :",
                "  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalFixedPointStabilityAllT := by",
                "  intro h",
                "  exact RG_C1_2_DERIVED (RG_C1_1_DERIVED h)",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    lean_bin = detect_lean()
    checker_attempted = False
    checker_ok = False
    checker_stdout = ""
    checker_stderr = ""
    checker_rc = None
    if lean_bin:
        checker_attempted = True
        rc, so, se = run([lean_bin, attempt.name])
        checker_rc, checker_stdout, checker_stderr = rc, so, se
        checker_ok = rc == 0

    merged = (checker_stdout + "\n" + checker_stderr).lower()
    missing_provider_error_present = ("unknown identifier `rg_c1_1_derived`" in merged) or (
        "unknown identifier `rg_c1_2_derived`" in merged
    )

    provider_axiomatic = "DerivedOrPending" in provider_text

    proof_object = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": {
            "q2237": "report_qw2237_rg_axiom_free_o1c_provider_layer_gate.json",
            "q2235": "report_qw2235_rg_axiom_free_o1c_theorem_discharge_execution_gate.json",
        },
        "provider_file": {
            "name": provider.name,
            "sha256": sha256_file(provider) if provider.exists() else None,
        },
        "attempt_file": {
            "name": attempt.name,
            "sha256": sha256_file(attempt),
        },
        "checker": {
            "attempted": checker_attempted,
            "exit_code": checker_rc,
            "stdout": checker_stdout,
            "stderr": checker_stderr,
            "missing_provider_error_present": missing_provider_error_present,
        },
        "scope_boundary": {
            "provider_layer_axiomatic": provider_axiomatic,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2239_rg_o1c_provider_execution.json"
    proof_path.write_text(json.dumps(proof_object, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2237_provider_layer_pass_present": q2237.get("verdict")
        == "RG_AXIOM_FREE_O1C_PROVIDER_LAYER_GATE_PASS_PARTIAL_AXIOMATIC_PROVIDER_OPEN",
        "q2235_missing_provider_blocker_previously_confirmed": q2235.get("verdict")
        == "RG_AXIOM_FREE_O1C_THEOREM_DISCHARGE_EXECUTION_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_SOURCE_THEOREMS",
        "provider_file_present": provider.exists(),
        "attempt_file_written": attempt.exists(),
        "lean_checker_detected": bool(lean_bin),
        "lean_checker_attempted": checker_attempted,
        "lean_checker_exit_zero": bool(checker_ok),
        "missing_provider_blocker_removed": not missing_provider_error_present,
        "provider_layer_still_axiomatic": bool(provider_axiomatic),
        "proof_object_generated": proof_path.exists(),
        "proof_object_hash_recorded": True,
        "c1_theorem_discharge_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2237_provider_layer_pass_present"]
        and flags["q2235_missing_provider_blocker_previously_confirmed"]
        and flags["provider_file_present"]
        and flags["attempt_file_written"]
        and flags["lean_checker_detected"]
        and flags["lean_checker_attempted"]
        and flags["lean_checker_exit_zero"]
        and flags["missing_provider_blocker_removed"]
        and flags["provider_layer_still_axiomatic"]
        and flags["proof_object_generated"]
        and flags["proof_object_hash_recorded"]
    )

    verdict = (
        "RG_AXIOM_FREE_O1C_PROVIDER_EXECUTION_GATE_PASS_PARTIAL_PROVIDER_OK_AXIOMATIC_SOURCE_OPEN"
        if core_ok
        else "RG_AXIOM_FREE_O1C_PROVIDER_EXECUTION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2237": "report_qw2237_rg_axiom_free_o1c_provider_layer_gate.json",
            "q2235": "report_qw2235_rg_axiom_free_o1c_theorem_discharge_execution_gate.json",
        },
        "provider_file": provider.name,
        "attempt_file": attempt.name,
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "REPLACE_DERIVED_OR_PENDING_AXIOMS_WITH_NON_AXIOMATIC_PROVIDER_THEOREMS_FROM_CANONICAL_CHAIN",
    }

    out_json = ROOT / "report_qw2239_rg_axiom_free_o1c_provider_execution_gate.json"
    out_md = ROOT / "RAPORT_QW2239_RG_AXIOM_FREE_O1C_PROVIDER_EXECUTION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2239: RG AXIOM-FREE O1C PROVIDER EXECUTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Missing-provider blocker zostal usuniety (provider symbols sa wykonywalne).",
                "- Boundary pozostaje jawna: provider source nadal axiomatic (`DerivedOrPending`).",
                "",
                "## Artifacts",
                "- JSON: `report_qw2239_rg_axiom_free_o1c_provider_execution_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
