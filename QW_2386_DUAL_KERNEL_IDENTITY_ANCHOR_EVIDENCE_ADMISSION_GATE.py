#!/usr/bin/env python3
"""QW-2386: dual anchor-evidence admission gate.

Admits execution only if noncircular anchor-candidate proof files exist for both
branches and satisfy hard hygiene prechecks.
"""

from __future__ import annotations

import hashlib
import json
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


def file_hard_hygiene_ok(path: Path) -> bool:
    txt = path.read_text(encoding="utf-8")
    return "axiom " not in txt and "_DerivedOrPending" not in txt


def choose_candidate(pattern: str) -> Path | None:
    files = sorted(ROOT.glob(pattern))
    return files[0] if files else None


def main() -> None:
    q2385 = load("report_qw2385_dual_kernel_identity_anchor_obligation_packet_gate.json")

    l12_file = choose_candidate("FIN_L12_*ANCHOR*_ATTEMPT.lean")
    l5_file = choose_candidate("FIN_L5_*ANCHOR*_ATTEMPT.lean")

    l12_exists = l12_file is not None
    l5_exists = l5_file is not None

    l12_hygiene = file_hard_hygiene_ok(l12_file) if l12_file else False
    l5_hygiene = file_hard_hygiene_ok(l5_file) if l5_file else False

    both_present = l12_exists and l5_exists
    both_hygiene = l12_hygiene and l5_hygiene
    admission_allowed = both_present and both_hygiene

    flags = {
        "q2385_packet_ready": q2385.get("verdict")
        == "DUAL_KERNEL_IDENTITY_ANCHOR_OBLIGATION_PACKET_GATE_PASS_PACKET_READY",
        "l12_anchor_candidate_file_present": l12_exists,
        "l5_anchor_candidate_file_present": l5_exists,
        "l12_anchor_candidate_hard_hygiene": l12_hygiene,
        "l5_anchor_candidate_hard_hygiene": l5_hygiene,
        "admission_allowed": admission_allowed,
        "admission_denied": not admission_allowed,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if flags["q2385_packet_ready"] and not admission_allowed:
        verdict = "DUAL_KERNEL_IDENTITY_ANCHOR_EVIDENCE_ADMISSION_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_ANCHOR_CANDIDATES"
        required_next_step = "AUTHOR_NONCIRCULAR_ANCHOR_PROOF_CANDIDATE_FILES"
    elif flags["q2385_packet_ready"] and admission_allowed:
        verdict = "DUAL_KERNEL_IDENTITY_ANCHOR_EVIDENCE_ADMISSION_GATE_PASS_ADMITTED"
        required_next_step = "RUN_DUAL_KERNEL_IDENTITY_ANCHOR_EXECUTION_STATUS_GATE"
    else:
        verdict = "DUAL_KERNEL_IDENTITY_ANCHOR_EVIDENCE_ADMISSION_GATE_FAIL"
        required_next_step = "REPAIR_ANCHOR_OBLIGATION_PACKET"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2385": "report_qw2385_dual_kernel_identity_anchor_obligation_packet_gate.json",
        },
        "candidate_files": {
            "l12": str(l12_file.name) if l12_file else None,
            "l5": str(l5_file.name) if l5_file else None,
        },
        "candidate_hard_hygiene": {
            "l12": l12_hygiene,
            "l5": l5_hygiene,
        },
        "scope_boundary": {
            "admission_allowed": admission_allowed,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2386_dual_kernel_identity_anchor_evidence_admission.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "candidate_files": proof_obj["candidate_files"],
        "candidate_hard_hygiene": proof_obj["candidate_hard_hygiene"],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2386_dual_kernel_identity_anchor_evidence_admission_gate.json"
    out_md = ROOT / "RAPORT_QW2386_DUAL_KERNEL_IDENTITY_ANCHOR_EVIDENCE_ADMISSION_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2386: DUAL KERNEL IDENTITY ANCHOR EVIDENCE ADMISSION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- l12_anchor_candidate: `{proof_obj['candidate_files']['l12']}`",
                f"- l5_anchor_candidate: `{proof_obj['candidate_files']['l5']}`",
                f"- admission_allowed: `{admission_allowed}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "admission_allowed": admission_allowed}))


if __name__ == "__main__":
    main()
