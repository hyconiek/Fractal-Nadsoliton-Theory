#!/usr/bin/env python3
"""QW-2383: dual noncyclic step admission gate.

Evaluates whether the immediate next-step proposal is admissible under
QW-2382 anti-recurrence constraints.
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


def normalize_cut(rows: list[dict[str, Any]]) -> list[tuple[str, str]]:
    return sorted((str(r.get("branch", "")), str(r.get("symbol", ""))) for r in rows)


def main() -> None:
    q2382 = load("report_qw2382_dual_noncyclic_strategy_packet_gate.json")
    q2320 = load("report_qw2320_dual_kernel_identity_closure_execution_status_gate.json")
    q2380 = load("report_qw2380_dual_kernel_identity_closure_execution_status_gate.json")

    candidate_step = str(q2380.get("required_next_step", ""))
    historical_step = str(q2320.get("required_next_step", ""))

    candidate_cut = normalize_cut(q2380.get("blocker_cut", []))
    historical_cut = normalize_cut(q2320.get("blocker_cut", []))

    repeat_step = candidate_step == historical_step
    repeat_cut = candidate_cut == historical_cut

    violations = {
        "NC1_REPEAT_STEP_FORBIDDEN": repeat_step,
        "NC2_BLOCKER_NOVELTY_REQUIRED": repeat_cut,
        "NC3_DUAL_BRANCH_PROGRESS_GUARD": repeat_cut and repeat_step,
    }

    flags = {
        "q2382_strategy_packet_ready": q2382.get("verdict") == "DUAL_NONCYCLIC_STRATEGY_PACKET_GATE_PASS_PACKET_READY",
        "candidate_step_detected": bool(candidate_step),
        "historical_step_detected": bool(historical_step),
        "repeat_step_detected": repeat_step,
        "repeat_blocker_cut_detected": repeat_cut,
        "hard_violation_present": any(violations.values()),
        "admission_denied": True,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_NONCYCLIC_STEP_ADMISSION_GATE_PASS_REPEAT_STEP_REJECTED"
        if flags["q2382_strategy_packet_ready"]
        and flags["candidate_step_detected"]
        and flags["historical_step_detected"]
        and flags["hard_violation_present"]
        and flags["admission_denied"]
        else "DUAL_NONCYCLIC_STEP_ADMISSION_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2382": "report_qw2382_dual_noncyclic_strategy_packet_gate.json",
            "q2320": "report_qw2320_dual_kernel_identity_closure_execution_status_gate.json",
            "q2380": "report_qw2380_dual_kernel_identity_closure_execution_status_gate.json",
        },
        "candidate_step": candidate_step,
        "historical_step": historical_step,
        "candidate_blocker_cut": q2380.get("blocker_cut", []),
        "historical_blocker_cut": q2320.get("blocker_cut", []),
        "violations": violations,
        "scope_boundary": {
            "admission_denied": True,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2383_dual_noncyclic_step_admission.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "candidate_step": candidate_step,
        "historical_step": historical_step,
        "violations": violations,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DESIGN_NONREPEATING_EXECUTION_PACKET",
    }

    out_json = ROOT / "report_qw2383_dual_noncyclic_step_admission_gate.json"
    out_md = ROOT / "RAPORT_QW2383_DUAL_NONCYCLIC_STEP_ADMISSION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2383: DUAL NONCYCLIC STEP ADMISSION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- candidate_step: `{candidate_step}`",
                f"- historical_step: `{historical_step}`",
                f"- violations: `{violations}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "hard_violation_present": flags["hard_violation_present"]}))


if __name__ == "__main__":
    main()
