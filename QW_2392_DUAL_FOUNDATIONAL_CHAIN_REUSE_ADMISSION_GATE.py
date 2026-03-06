#!/usr/bin/env python3
"""QW-2392: dual foundational-chain reuse admission gate.

Decides whether reusing historical foundational chain is admissible without
new noncyclic evidence.
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


def choose_candidate(pattern: str) -> Path | None:
    files = sorted(ROOT.glob(pattern))
    return files[0] if files else None


def hygiene(path: Path) -> bool:
    txt = path.read_text(encoding="utf-8")
    return "axiom " not in txt and "_DerivedOrPending" not in txt


def main() -> None:
    q2391 = load("report_qw2391_dual_action_level_anchor_frontier_alignment_gate.json")
    q2381 = load("report_qw2381_dual_kernel_cycle_recurrence_gate.json")

    l12_cand = choose_candidate("FIN_L12_FOUNDATIONAL_*ANCHOR*_ATTEMPT.lean")
    l5_cand = choose_candidate("FIN_L5_FOUNDATIONAL_*ANCHOR*_ATTEMPT.lean")

    l12_present = l12_cand is not None
    l5_present = l5_cand is not None
    l12_ok = hygiene(l12_cand) if l12_cand else False
    l5_ok = hygiene(l5_cand) if l5_cand else False

    aligned = bool(q2391.get("flags", {}).get("foundational_frontier_aligned", False))
    prior_cycle = q2381.get("verdict") == "DUAL_KERNEL_CYCLE_RECURRENCE_GATE_PASS_BLOCKER_LOOP_CONFIRMED"

    can_reuse_without_new_evidence = False
    admission_allowed = (l12_present and l5_present and l12_ok and l5_ok) and not can_reuse_without_new_evidence

    flags = {
        "q2391_frontier_alignment_confirmed": aligned,
        "q2381_cycle_recurrence_confirmed": prior_cycle,
        "l12_foundational_anchor_candidate_present": l12_present,
        "l5_foundational_anchor_candidate_present": l5_present,
        "l12_foundational_anchor_candidate_hard_hygiene": l12_ok,
        "l5_foundational_anchor_candidate_hard_hygiene": l5_ok,
        "reuse_without_new_evidence_forbidden": aligned and prior_cycle,
        "admission_allowed": admission_allowed,
        "admission_denied": not admission_allowed,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if flags["reuse_without_new_evidence_forbidden"] and not admission_allowed:
        verdict = "DUAL_FOUNDATIONAL_CHAIN_REUSE_ADMISSION_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_FOUNDATIONAL_ANCHOR_CANDIDATES"
        required_next_step = "BUILD_DUAL_FOUNDATIONAL_NONCYCLIC_ANCHOR_OBLIGATION_PACKET"
    elif admission_allowed:
        verdict = "DUAL_FOUNDATIONAL_CHAIN_REUSE_ADMISSION_GATE_PASS_ADMITTED"
        required_next_step = "RUN_DUAL_FOUNDATIONAL_ANCHOR_EXECUTION_STATUS_GATE"
    else:
        verdict = "DUAL_FOUNDATIONAL_CHAIN_REUSE_ADMISSION_GATE_FAIL"
        required_next_step = "REPAIR_FOUNDATIONAL_REUSE_ADMISSION_PIPELINE"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2391": "report_qw2391_dual_action_level_anchor_frontier_alignment_gate.json",
            "q2381": "report_qw2381_dual_kernel_cycle_recurrence_gate.json",
        },
        "candidate_files": {
            "l12": l12_cand.name if l12_cand else None,
            "l5": l5_cand.name if l5_cand else None,
        },
        "scope_boundary": {
            "admission_allowed": admission_allowed,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2392_dual_foundational_chain_reuse_admission.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "candidate_files": proof_obj["candidate_files"],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2392_dual_foundational_chain_reuse_admission_gate.json"
    out_md = ROOT / "RAPORT_QW2392_DUAL_FOUNDATIONAL_CHAIN_REUSE_ADMISSION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2392: DUAL FOUNDATIONAL CHAIN REUSE ADMISSION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- admission_allowed: `{admission_allowed}`",
                f"- candidate_files: `{proof_obj['candidate_files']}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "admission_allowed": admission_allowed}))


if __name__ == "__main__":
    main()
