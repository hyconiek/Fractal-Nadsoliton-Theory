#!/usr/bin/env python3
"""QW-2332: dual kernel-identity-regularity execution status gate.

Runs machine-check attempts for both identity-regularity symbols and isolates
next blocker-cut at kernel-identity-conservation layer.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
LEAN_FALLBACK = Path("/tmp/lean4/lean-4.28.0-linux/bin/lean")

L12_ATTEMPT_FILE = ROOT / "FIN_L12_KERNEL_IDENTITY_REGULARITY_DISCHARGE_ATTEMPT.lean"
L5_ATTEMPT_FILE = ROOT / "FIN_L5_KERNEL_IDENTITY_REGULARITY_DISCHARGE_ATTEMPT.lean"

L12_TARGET = "RG_KernelIdentityRegularityToWellPosedness_Theorem"
L5_TARGET = "QFT_KernelIdentityRegularityToPositivity_Theorem"
L12_BLOCKER = "RG_KernelIdentityConservationToWellPosedness_Theorem"
L5_BLOCKER = "QFT_KernelIdentityConservationToPositivity_Theorem"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def theorem_defined(text: str, theorem_name: str) -> bool:
    return bool(re.search(rf"^theorem\s+{re.escape(theorem_name)}\s*:", text, flags=re.M))


def detect_lean() -> Path | None:
    direct = shutil.which("lean")
    if direct:
        return Path(direct)
    if LEAN_FALLBACK.exists() and LEAN_FALLBACK.is_file():
        return LEAN_FALLBACK
    return None


def run_lean(lean_bin: Path, target: Path) -> tuple[int, str, str]:
    proc = subprocess.run(
        [str(lean_bin), target.name],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout, proc.stderr


def write_attempt_files() -> None:
    L12_ATTEMPT_FILE.write_text(
        "\n".join(
            [
                "-- FIN Release 5.1: QW-2332 L12 kernel-identity-regularity discharge attempt",
                "-- Scope: derive identity-regularity theorem from identity-conservation theorem.",
                "",
                "variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)",
                "",
                "theorem RG_KernelIdentityRegularityToWellPosedness_Theorem :",
                "  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by",
                "  exact RG_KernelIdentityConservationToWellPosedness_Theorem",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    L5_ATTEMPT_FILE.write_text(
        "\n".join(
            [
                "-- FIN Release 5.1: QW-2332 L5 kernel-identity-regularity discharge attempt",
                "-- Scope: derive identity-regularity theorem from identity-conservation theorem.",
                "",
                "variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)",
                "",
                "theorem QFT_KernelIdentityRegularityToPositivity_Theorem :",
                "  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by",
                "  exact QFT_KernelIdentityConservationToPositivity_Theorem",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    q2331 = load("report_qw2331_dual_kernel_identity_regularity_discharge_packet_gate.json")
    write_attempt_files()

    l12_text = L12_ATTEMPT_FILE.read_text(encoding="utf-8")
    l5_text = L5_ATTEMPT_FILE.read_text(encoding="utf-8")

    lean_bin = detect_lean()
    l12_rc = None
    l12_stdout = ""
    l12_stderr = ""
    l5_rc = None
    l5_stdout = ""
    l5_stderr = ""

    if lean_bin:
        l12_rc, l12_stdout, l12_stderr = run_lean(lean_bin, L12_ATTEMPT_FILE)
        l5_rc, l5_stdout, l5_stderr = run_lean(lean_bin, L5_ATTEMPT_FILE)

    l12_merged = (l12_stdout + "\n" + l12_stderr).lower()
    l5_merged = (l5_stdout + "\n" + l5_stderr).lower()

    l12_blocker_confirmed = l12_rc != 0 and f"unknown identifier `{L12_BLOCKER.lower()}`" in l12_merged
    l5_blocker_confirmed = l5_rc != 0 and f"unknown identifier `{L5_BLOCKER.lower()}`" in l5_merged

    flags = {
        "q2331_packet_ready": q2331.get("verdict")
        == "DUAL_KERNEL_IDENTITY_REGULARITY_DISCHARGE_PACKET_GATE_PASS_PACKET_READY",
        "l12_attempt_file_written": L12_ATTEMPT_FILE.exists(),
        "l5_attempt_file_written": L5_ATTEMPT_FILE.exists(),
        "l12_target_theorem_defined": theorem_defined(l12_text, L12_TARGET),
        "l5_target_theorem_defined": theorem_defined(l5_text, L5_TARGET),
        "l12_contains_no_axiom_tokens": "axiom " not in l12_text,
        "l5_contains_no_axiom_tokens": "axiom " not in l5_text,
        "l12_contains_no_derived_or_pending_tokens": "_DerivedOrPending" not in l12_text,
        "l5_contains_no_derived_or_pending_tokens": "_DerivedOrPending" not in l5_text,
        "lean_binary_detected": bool(lean_bin),
        "machine_check_executed_both": l12_rc is not None and l5_rc is not None,
        "l12_blocker_confirmed_missing_identity_conservation_symbol": l12_blocker_confirmed,
        "l5_blocker_confirmed_missing_identity_conservation_symbol": l5_blocker_confirmed,
        "traceability_to_identity_conservation_symbols_explicit": L12_BLOCKER in l12_text and L5_BLOCKER in l5_text,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2331_packet_ready"]
        and flags["l12_attempt_file_written"]
        and flags["l5_attempt_file_written"]
        and flags["l12_target_theorem_defined"]
        and flags["l5_target_theorem_defined"]
        and flags["l12_contains_no_axiom_tokens"]
        and flags["l5_contains_no_axiom_tokens"]
        and flags["l12_contains_no_derived_or_pending_tokens"]
        and flags["l5_contains_no_derived_or_pending_tokens"]
        and flags["lean_binary_detected"]
        and flags["machine_check_executed_both"]
        and flags["l12_blocker_confirmed_missing_identity_conservation_symbol"]
        and flags["l5_blocker_confirmed_missing_identity_conservation_symbol"]
        and flags["traceability_to_identity_conservation_symbols_explicit"]
    )

    verdict = (
        "DUAL_KERNEL_IDENTITY_REGULARITY_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_IDENTITY_CONSERVATION_THEOREMS"
        if core_ok
        else "DUAL_KERNEL_IDENTITY_REGULARITY_EXECUTION_STATUS_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2331_dual_kernel_identity_regularity_discharge_packet_gate.json",
        "attempt_files": {
            "l12": {
                "name": L12_ATTEMPT_FILE.name,
                "sha256": sha256_file(L12_ATTEMPT_FILE),
                "target_theorem": L12_TARGET,
                "missing_identity_conservation_symbol": L12_BLOCKER,
            },
            "l5": {
                "name": L5_ATTEMPT_FILE.name,
                "sha256": sha256_file(L5_ATTEMPT_FILE),
                "target_theorem": L5_TARGET,
                "missing_identity_conservation_symbol": L5_BLOCKER,
            },
        },
        "machine_check": {
            "lean_binary": str(lean_bin) if lean_bin else None,
            "l12": {"exit_code": l12_rc, "stdout": l12_stdout, "stderr": l12_stderr},
            "l5": {"exit_code": l5_rc, "stdout": l5_stdout, "stderr": l5_stderr},
        },
        "blocker_cut": [
            {"branch": "L12", "symbol": L12_BLOCKER},
            {"branch": "L5", "symbol": L5_BLOCKER},
        ],
        "scope_boundary": {
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2332_dual_kernel_identity_regularity_execution_status.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2331_dual_kernel_identity_regularity_discharge_packet_gate.json",
        "attempt_files": {"l12": L12_ATTEMPT_FILE.name, "l5": L5_ATTEMPT_FILE.name},
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "machine_check_exit_codes": {"l12": l12_rc, "l5": l5_rc},
        "blocker_cut": proof_obj["blocker_cut"],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "EXTRACT_MINIMAL_KERNEL_IDENTITY_CONSERVATION_BLOCKER_CUT",
    }

    out_json = ROOT / "report_qw2332_dual_kernel_identity_regularity_execution_status_gate.json"
    out_md = ROOT / "RAPORT_QW2332_DUAL_KERNEL_IDENTITY_REGULARITY_EXECUTION_STATUS_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2332: DUAL KERNEL IDENTITY REGULARITY EXECUTION STATUS GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- machine_check_exit_codes: `{out['machine_check_exit_codes']}`",
                f"- blocker_cut: `{out['blocker_cut']}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "machine_check_exit_codes": out["machine_check_exit_codes"],
                "blocker_cut_size": len(out["blocker_cut"]),
            }
        )
    )


if __name__ == "__main__":
    main()
