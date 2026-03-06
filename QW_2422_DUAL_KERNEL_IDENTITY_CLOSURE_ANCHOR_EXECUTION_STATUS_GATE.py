#!/usr/bin/env python3
"""QW-2422: dual kernel-identity-closure anchor execution status gate.

Runs machine-check for kernel-identity-closure noncyclic anchor candidates and
isolates kernel-identity-locality blockers.
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

L12_FILE = ROOT / "FIN_L12_KERNEL_IDENTITY_CLOSURE_NONCYCLIC_ANCHOR_ATTEMPT.lean"
L5_FILE = ROOT / "FIN_L5_KERNEL_IDENTITY_CLOSURE_NONCYCLIC_ANCHOR_ATTEMPT.lean"

L12_BLOCKER = "RG_KernelIdentityLocalityToWellPosedness_Theorem"
L5_BLOCKER = "QFT_KernelIdentityLocalityToPositivity_Theorem"

UNKNOWN_RE = re.compile(r"unknown identifier `([^`]+)`", flags=re.I)


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def detect_lean() -> Path | None:
    direct = shutil.which("lean")
    if direct:
        return Path(direct)
    if LEAN_FALLBACK.exists() and LEAN_FALLBACK.is_file():
        return LEAN_FALLBACK
    return None


def run_lean(lean_bin: Path, target: Path) -> tuple[int, str, str]:
    proc = subprocess.run([str(lean_bin), target.name], cwd=ROOT, check=False, capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def hard_hygiene(path: Path) -> bool:
    txt = path.read_text(encoding="utf-8")
    return "axiom " not in txt and "_DerivedOrPending" not in txt


def extract_unknown(stdout: str, stderr: str) -> str | None:
    m = UNKNOWN_RE.search(f"{stdout}\n{stderr}")
    return m.group(1) if m else None


def main() -> None:
    q2420 = load("report_qw2420_dual_kernel_identity_closure_chain_reuse_admission_gate.json")
    q2421 = load("report_qw2421_dual_kernel_identity_closure_noncyclic_anchor_obligation_packet_gate.json")

    lean_bin = detect_lean()

    l12_rc = l5_rc = None
    l12_out = l12_err = l5_out = l5_err = ""

    if lean_bin and L12_FILE.exists() and L5_FILE.exists():
        l12_rc, l12_out, l12_err = run_lean(lean_bin, L12_FILE)
        l5_rc, l5_out, l5_err = run_lean(lean_bin, L5_FILE)

    l12_unknown = extract_unknown(l12_out, l12_err)
    l5_unknown = extract_unknown(l5_out, l5_err)

    flags = {
        "q2420_admission_allowed": bool(q2420.get("flags", {}).get("admission_allowed", False)),
        "q2421_packet_ready": q2421.get("verdict")
        == "DUAL_KERNEL_IDENTITY_CLOSURE_NONCYCLIC_ANCHOR_OBLIGATION_PACKET_GATE_PASS_PACKET_READY",
        "l12_candidate_exists": L12_FILE.exists(),
        "l5_candidate_exists": L5_FILE.exists(),
        "l12_candidate_hard_hygiene": hard_hygiene(L12_FILE) if L12_FILE.exists() else False,
        "l5_candidate_hard_hygiene": hard_hygiene(L5_FILE) if L5_FILE.exists() else False,
        "lean_binary_detected": bool(lean_bin),
        "machine_check_executed_both": l12_rc is not None and l5_rc is not None,
        "l12_blocker_confirmed_missing_identity_locality_symbol": l12_rc != 0 and l12_unknown == L12_BLOCKER,
        "l5_blocker_confirmed_missing_identity_locality_symbol": l5_rc != 0 and l5_unknown == L5_BLOCKER,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if (
        flags["q2420_admission_allowed"]
        and flags["q2421_packet_ready"]
        and flags["l12_candidate_exists"]
        and flags["l5_candidate_exists"]
        and flags["l12_candidate_hard_hygiene"]
        and flags["l5_candidate_hard_hygiene"]
        and flags["lean_binary_detected"]
        and flags["machine_check_executed_both"]
        and flags["l12_blocker_confirmed_missing_identity_locality_symbol"]
        and flags["l5_blocker_confirmed_missing_identity_locality_symbol"]
    ):
        verdict = "DUAL_KERNEL_IDENTITY_CLOSURE_ANCHOR_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_IDENTITY_LOCALITY_THEOREMS"
        required_next_step = "CHECK_IDENTITY_LOCALITY_FRONTIER_ALIGNMENT_WITH_QW2320"
        blocker_cut = [
            {"branch": "L12", "symbol": L12_BLOCKER},
            {"branch": "L5", "symbol": L5_BLOCKER},
        ]
    elif flags["q2420_admission_allowed"] and flags["machine_check_executed_both"] and l12_rc == 0 and l5_rc == 0:
        verdict = "DUAL_KERNEL_IDENTITY_CLOSURE_ANCHOR_EXECUTION_STATUS_GATE_PASS_PARTIAL_MACHINE_CHECKED"
        required_next_step = "RUN_DUAL_KERNEL_IDENTITY_CLOSURE_ANCHOR_INTEGRATION_GATE"
        blocker_cut = []
    else:
        verdict = "DUAL_KERNEL_IDENTITY_CLOSURE_ANCHOR_EXECUTION_STATUS_GATE_FAIL"
        required_next_step = "REPAIR_DUAL_KERNEL_IDENTITY_CLOSURE_ANCHOR_EXECUTION_PIPELINE"
        blocker_cut = []

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2420": "report_qw2420_dual_kernel_identity_closure_chain_reuse_admission_gate.json",
            "q2421": "report_qw2421_dual_kernel_identity_closure_noncyclic_anchor_obligation_packet_gate.json",
        },
        "attempt_files": {
            "l12": {"name": L12_FILE.name if L12_FILE.exists() else None, "sha256": sha256_file(L12_FILE) if L12_FILE.exists() else None},
            "l5": {"name": L5_FILE.name if L5_FILE.exists() else None, "sha256": sha256_file(L5_FILE) if L5_FILE.exists() else None},
        },
        "machine_check": {
            "lean_binary": str(lean_bin) if lean_bin else None,
            "l12": {"exit_code": l12_rc, "stdout": l12_out, "stderr": l12_err, "unknown_identifier": l12_unknown},
            "l5": {"exit_code": l5_rc, "stdout": l5_out, "stderr": l5_err, "unknown_identifier": l5_unknown},
        },
        "blocker_cut": blocker_cut,
        "scope_boundary": {
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2422_dual_kernel_identity_closure_anchor_execution_status.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "attempt_files": {"l12": L12_FILE.name, "l5": L5_FILE.name},
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "machine_check_exit_codes": {"l12": l12_rc, "l5": l5_rc},
        "blocker_cut": blocker_cut,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2422_dual_kernel_identity_closure_anchor_execution_status_gate.json"
    out_md = ROOT / "RAPORT_QW2422_DUAL_KERNEL_IDENTITY_CLOSURE_ANCHOR_EXECUTION_STATUS_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2422: DUAL KERNEL IDENTITY CLOSURE ANCHOR EXECUTION STATUS GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- machine_check_exit_codes: `{out['machine_check_exit_codes']}`",
                f"- blocker_cut: `{blocker_cut}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "exit": out["machine_check_exit_codes"]}))


if __name__ == "__main__":
    main()
