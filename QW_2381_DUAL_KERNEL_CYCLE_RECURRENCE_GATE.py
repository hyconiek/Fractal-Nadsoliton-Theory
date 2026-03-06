#!/usr/bin/env python3
"""QW-2381: dual-kernel cycle recurrence gate.

Confirms whether the current dual blocker-cut reproduces an earlier blocker-cut,
which indicates cycle recurrence rather than net theorem-level progress.
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


def normalize_blocker_cut(rows: list[dict[str, Any]]) -> list[tuple[str, str]]:
    pairs = []
    for row in rows:
        pairs.append((str(row.get("branch", "")), str(row.get("symbol", ""))))
    return sorted(pairs)


def main() -> None:
    baseline = load("report_qw2320_dual_kernel_identity_closure_execution_status_gate.json")
    current = load("report_qw2380_dual_kernel_identity_closure_execution_status_gate.json")

    base_cut = baseline.get("blocker_cut", [])
    curr_cut = current.get("blocker_cut", [])

    base_norm = normalize_blocker_cut(base_cut)
    curr_norm = normalize_blocker_cut(curr_cut)

    flags = {
        "baseline_and_current_have_two_blockers": len(base_cut) == 2 and len(curr_cut) == 2,
        "baseline_has_l12_l5_branches": sorted(r.get("branch") for r in base_cut) == ["L12", "L5"],
        "current_has_l12_l5_branches": sorted(r.get("branch") for r in curr_cut) == ["L12", "L5"],
        "blocker_cuts_identical_as_sets": base_norm == curr_norm,
        "current_verdict_is_partial_blocked": current.get("verdict")
        == "DUAL_KERNEL_IDENTITY_CLOSURE_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_IDENTITY_LOCALITY_THEOREMS",
        "execution_completed_still_false": current.get("flags", {}).get("execution_completed") is False,
        "all_strict_obligations_fully_closed_still_false": current.get("flags", {}).get("all_strict_obligations_fully_closed") is False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_KERNEL_CYCLE_RECURRENCE_GATE_PASS_BLOCKER_LOOP_CONFIRMED"
        if all(flags.values())
        else "DUAL_KERNEL_CYCLE_RECURRENCE_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "baseline": "report_qw2320_dual_kernel_identity_closure_execution_status_gate.json",
        "current": "report_qw2380_dual_kernel_identity_closure_execution_status_gate.json",
        "baseline_blocker_cut": base_cut,
        "current_blocker_cut": curr_cut,
        "scope_boundary": {
            "cycle_recurrence_assessed": True,
            "theorem_level_progress_assessed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2381_dual_kernel_cycle_recurrence.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "baseline": "report_qw2320_dual_kernel_identity_closure_execution_status_gate.json",
            "current": "report_qw2380_dual_kernel_identity_closure_execution_status_gate.json",
            "spec": spec_path.name,
        },
        "baseline_blocker_cut_normalized": base_norm,
        "current_blocker_cut_normalized": curr_norm,
    }
    proof_path = ROOT / "proof_object_qw2381_dual_kernel_cycle_recurrence.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "baseline_blocker_cut": base_cut,
        "current_blocker_cut": curr_cut,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DESIGN_NONCYCLIC_STRATEGY_PACKET",
    }

    out_json = ROOT / "report_qw2381_dual_kernel_cycle_recurrence_gate.json"
    out_md = ROOT / "RAPORT_QW2381_DUAL_KERNEL_CYCLE_RECURRENCE_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2381: DUAL KERNEL CYCLE RECURRENCE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- baseline_blocker_cut: `{base_cut}`",
                f"- current_blocker_cut: `{curr_cut}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": f"{pass_count}/{total_flags}"}))


if __name__ == "__main__":
    main()
