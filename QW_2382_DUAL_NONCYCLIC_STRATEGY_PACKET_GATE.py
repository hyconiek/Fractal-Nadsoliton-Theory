#!/usr/bin/env python3
"""QW-2382: dual noncyclic strategy packet gate.

Builds strict anti-recurrence strategy packet after confirmed blocker-cut loop.
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
    q2381 = load("report_qw2381_dual_kernel_cycle_recurrence_gate.json")
    q2320 = load("report_qw2320_dual_kernel_identity_closure_execution_status_gate.json")
    q2380 = load("report_qw2380_dual_kernel_identity_closure_execution_status_gate.json")

    base_cut = q2320.get("blocker_cut", [])
    curr_cut = q2380.get("blocker_cut", [])

    constraints = [
        {
            "id": "NC1_REPEAT_STEP_FORBIDDEN",
            "rule": "required_next_step must differ from last historical step for identical blocker-cut",
            "severity": "hard",
        },
        {
            "id": "NC2_BLOCKER_NOVELTY_REQUIRED",
            "rule": "new execution attempt must include at least one blocker symbol outside baseline recurrence set",
            "severity": "hard",
        },
        {
            "id": "NC3_DUAL_BRANCH_PROGRESS_GUARD",
            "rule": "both L12 and L5 branches must show non-repeating transition evidence",
            "severity": "hard",
        },
        {
            "id": "NC4_NO_FALSE_CLOSURE_CLAIM",
            "rule": "all_strict_obligations_fully_closed remains false until theorem-level discharge exists",
            "severity": "hard",
        },
    ]

    flags = {
        "q2381_recurrence_confirmed": q2381.get("verdict")
        == "DUAL_KERNEL_CYCLE_RECURRENCE_GATE_PASS_BLOCKER_LOOP_CONFIRMED",
        "baseline_and_current_cuts_identical": normalize_cut(base_cut) == normalize_cut(curr_cut),
        "historical_and_current_required_next_step_identical": q2320.get("required_next_step")
        == q2380.get("required_next_step"),
        "noncyclic_constraints_count_ge_four": len(constraints) >= 4,
        "hard_constraints_only": all(c.get("severity") == "hard" for c in constraints),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_NONCYCLIC_STRATEGY_PACKET_GATE_PASS_PACKET_READY"
        if flags["q2381_recurrence_confirmed"]
        and flags["baseline_and_current_cuts_identical"]
        and flags["historical_and_current_required_next_step_identical"]
        and flags["noncyclic_constraints_count_ge_four"]
        and flags["hard_constraints_only"]
        else "DUAL_NONCYCLIC_STRATEGY_PACKET_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2381": "report_qw2381_dual_kernel_cycle_recurrence_gate.json",
            "q2320": "report_qw2320_dual_kernel_identity_closure_execution_status_gate.json",
            "q2380": "report_qw2380_dual_kernel_identity_closure_execution_status_gate.json",
        },
        "constraints": constraints,
        "baseline_blocker_cut": base_cut,
        "current_blocker_cut": curr_cut,
        "historical_required_next_step": q2320.get("required_next_step"),
        "current_required_next_step": q2380.get("required_next_step"),
        "scope_boundary": {
            "strategy_packet_ready": True,
            "execution_admitted": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    spec_path = ROOT / "spec_qw2382_dual_noncyclic_strategy_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": spec["sources"] | {"spec": spec_path.name},
        "constraints": constraints,
        "baseline_cut_normalized": normalize_cut(base_cut),
        "current_cut_normalized": normalize_cut(curr_cut),
    }
    proof_path = ROOT / "proof_object_qw2382_dual_noncyclic_strategy_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "constraints": constraints,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "RUN_NONCYCLIC_STEP_ADMISSION_GATE",
    }

    out_json = ROOT / "report_qw2382_dual_noncyclic_strategy_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2382_DUAL_NONCYCLIC_STRATEGY_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2382: DUAL NONCYCLIC STRATEGY PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- constraints_count: `{len(constraints)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "constraints_count": len(constraints)}))


if __name__ == "__main__":
    main()
