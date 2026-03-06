#!/usr/bin/env python3
"""QW-2387: dual kernel-identity anchor execution status gate.

Runs machine-check for admitted noncircular anchor candidates and classifies
result without overclaiming closure.
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

L12_FILE = ROOT / "FIN_L12_KERNEL_IDENTITY_LOCALITY_ANCHOR_ATTEMPT.lean"
L5_FILE = ROOT / "FIN_L5_KERNEL_IDENTITY_LOCALITY_ANCHOR_ATTEMPT.lean"

L12_EXPECTED_PROVIDER = "RG_ActionLevel_PhysicalBridge_Derivation"
L5_EXPECTED_PROVIDER = "QFT_ActionLevel_PhysicalBridge_Derivation"

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
    proc = subprocess.run(
        [str(lean_bin), target.name],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout, proc.stderr


def extract_unknown(stderr: str, stdout: str) -> str | None:
    merged = f"{stdout}\n{stderr}"
    m = UNKNOWN_RE.search(merged)
    if not m:
        return None
    return m.group(1)


def hard_hygiene(path: Path) -> bool:
    txt = path.read_text(encoding="utf-8")
    return "axiom " not in txt and "_DerivedOrPending" not in txt


def main() -> None:
    q2386 = load("report_qw2386_dual_kernel_identity_anchor_evidence_admission_gate.json")

    lean_bin = detect_lean()

    l12_rc = l5_rc = None
    l12_out = l12_err = l5_out = l5_err = ""

    if lean_bin and L12_FILE.exists() and L5_FILE.exists():
        l12_rc, l12_out, l12_err = run_lean(lean_bin, L12_FILE)
        l5_rc, l5_out, l5_err = run_lean(lean_bin, L5_FILE)

    l12_unknown = extract_unknown(l12_err, l12_out)
    l5_unknown = extract_unknown(l5_err, l5_out)

    l12_expected_missing = l12_unknown == L12_EXPECTED_PROVIDER
    l5_expected_missing = l5_unknown == L5_EXPECTED_PROVIDER

    flags = {
        "q2386_admission_allowed": bool(q2386.get("flags", {}).get("admission_allowed", False)),
        "l12_candidate_exists": L12_FILE.exists(),
        "l5_candidate_exists": L5_FILE.exists(),
        "l12_candidate_hard_hygiene": hard_hygiene(L12_FILE) if L12_FILE.exists() else False,
        "l5_candidate_hard_hygiene": hard_hygiene(L5_FILE) if L5_FILE.exists() else False,
        "lean_binary_detected": bool(lean_bin),
        "machine_check_executed_both": l12_rc is not None and l5_rc is not None,
        "l12_expected_provider_missing_blocker_confirmed": l12_rc != 0 and l12_expected_missing,
        "l5_expected_provider_missing_blocker_confirmed": l5_rc != 0 and l5_expected_missing,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    both_success = l12_rc == 0 and l5_rc == 0
    both_expected_blocked = (
        flags["l12_expected_provider_missing_blocker_confirmed"]
        and flags["l5_expected_provider_missing_blocker_confirmed"]
    )

    if (
        flags["q2386_admission_allowed"]
        and flags["l12_candidate_exists"]
        and flags["l5_candidate_exists"]
        and flags["l12_candidate_hard_hygiene"]
        and flags["l5_candidate_hard_hygiene"]
        and flags["lean_binary_detected"]
        and flags["machine_check_executed_both"]
        and both_expected_blocked
    ):
        verdict = "DUAL_KERNEL_IDENTITY_ANCHOR_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_ACTION_LEVEL_ANCHOR_PROVIDERS"
        required_next_step = "EXTRACT_DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_BLOCKER_CUT"
        blocker_cut = [
            {"branch": "L12", "symbol": L12_EXPECTED_PROVIDER},
            {"branch": "L5", "symbol": L5_EXPECTED_PROVIDER},
        ]
    elif (
        flags["q2386_admission_allowed"]
        and flags["machine_check_executed_both"]
        and both_success
    ):
        verdict = "DUAL_KERNEL_IDENTITY_ANCHOR_EXECUTION_STATUS_GATE_PASS_PARTIAL_ANCHOR_MACHINE_CHECKED"
        required_next_step = "RUN_DUAL_KERNEL_IDENTITY_ANCHOR_INTEGRATION_GATE"
        blocker_cut = []
    else:
        verdict = "DUAL_KERNEL_IDENTITY_ANCHOR_EXECUTION_STATUS_GATE_FAIL"
        required_next_step = "REPAIR_DUAL_KERNEL_IDENTITY_ANCHOR_EXECUTION_PIPELINE"
        blocker_cut = []

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2386_dual_kernel_identity_anchor_evidence_admission_gate.json",
        "attempt_files": {
            "l12": {"name": L12_FILE.name, "sha256": sha256_file(L12_FILE) if L12_FILE.exists() else None},
            "l5": {"name": L5_FILE.name, "sha256": sha256_file(L5_FILE) if L5_FILE.exists() else None},
        },
        "machine_check": {
            "lean_binary": str(lean_bin) if lean_bin else None,
            "l12": {
                "exit_code": l12_rc,
                "stdout": l12_out,
                "stderr": l12_err,
                "unknown_identifier": l12_unknown,
                "expected_provider": L12_EXPECTED_PROVIDER,
            },
            "l5": {
                "exit_code": l5_rc,
                "stdout": l5_out,
                "stderr": l5_err,
                "unknown_identifier": l5_unknown,
                "expected_provider": L5_EXPECTED_PROVIDER,
            },
        },
        "blocker_cut": blocker_cut,
        "scope_boundary": {
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2387_dual_kernel_identity_anchor_execution_status.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2386_dual_kernel_identity_anchor_evidence_admission_gate.json",
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

    out_json = ROOT / "report_qw2387_dual_kernel_identity_anchor_execution_status_gate.json"
    out_md = ROOT / "RAPORT_QW2387_DUAL_KERNEL_IDENTITY_ANCHOR_EXECUTION_STATUS_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2387: DUAL KERNEL IDENTITY ANCHOR EXECUTION STATUS GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- exit_codes: `{out['machine_check_exit_codes']}`",
                f"- blocker_cut: `{blocker_cut}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "exit": out["machine_check_exit_codes"]}))


if __name__ == "__main__":
    main()
