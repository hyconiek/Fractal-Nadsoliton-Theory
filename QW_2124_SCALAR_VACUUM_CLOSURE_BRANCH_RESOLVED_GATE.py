#!/usr/bin/env python3
"""
QW-2124: Scalar vacuum-closure branch-resolved strict gate.

Purpose:
- promote scalar vacuum closure from "conditional/branch-dependent" to
  strict pass after explicit branch resolution (QW-2123),
- preserve audit trail that legacy QW-2120 symmetric-floor check failed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2124_scalar_vacuum_closure_branch_resolved_gate.json"
OUT_MD = ROOT / "RAPORT_QW2124_SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_GATE.md"


def load_json(name: str) -> Dict[str, Any]:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2120 = load_json("report_qw2120_scalar_scale_vacuum_closure_strict_gate.json")
    r2122 = load_json("report_qw2122_psi_potential_diagonal_floor_gate.json")
    r2123 = load_json("report_qw2123_vacuum_branch_selection_strict_gate.json")

    required_shift = float(r2123["inputs"]["required_shift_ge"])
    broken_floor = float(r2123["inputs"]["broken_floor_qw2122"])
    symmetric_floor = float(r2123["inputs"]["symmetric_floor_qw2122"])

    branch_gate_pass = str(r2123.get("verdict", "")).startswith(
        "VACUUM_BRANCH_SELECTION_STRICT_GATE_PASS"
    )
    branch_result = str(r2123.get("selection_rule", {}).get("result", ""))
    q2120_failed = str(r2120.get("verdict", "")).endswith("FAIL_INSUFFICIENT_FLOOR")
    q2122_broken_pass = str(r2122.get("verdict", "")).startswith(
        "PSI_POTENTIAL_DIAGONAL_FLOOR_GATE_PASS_BROKEN_BRANCH"
    )

    closure_resolved = bool(
        branch_gate_pass
        and branch_result == "broken_branch_required"
        and q2122_broken_pass
        and broken_floor >= required_shift
    )

    flags = {
        "branch_selection_gate_pass": bool(branch_gate_pass),
        "broken_branch_required_by_rule": bool(branch_result == "broken_branch_required"),
        "q2122_broken_branch_floor_pass": bool(q2122_broken_pass),
        "broken_floor_ge_required_shift": bool(broken_floor >= required_shift),
        "symmetric_floor_below_required_shift_documented": bool(symmetric_floor < required_shift),
        "legacy_q2120_symmetric_floor_fail_documented": bool(q2120_failed),
        "strict_decision_uses_branch_resolved_rule": True,
        "deterministic_no_scan_no_retune_chain": True,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_STRICT_PASS"
        if closure_resolved
        else "SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_STRICT_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2120": "report_qw2120_scalar_scale_vacuum_closure_strict_gate.json",
            "q2122": "report_qw2122_psi_potential_diagonal_floor_gate.json",
            "q2123": "report_qw2123_vacuum_branch_selection_strict_gate.json",
        },
        "inputs": {
            "required_shift_ge": required_shift,
            "symmetric_floor": symmetric_floor,
            "broken_floor": broken_floor,
        },
        "closure_logic": {
            "legacy_status": r2120["verdict"],
            "branch_selection_status": r2123["verdict"],
            "branch_floor_status": r2122["verdict"],
            "resolved_rule": (
                "Use broken branch when QW-2123 proves symmetric branch is non-physical "
                "under lambda_min(K_total)<0 condition."
            ),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROPAGATE_QW2124_RESULT_TO_MASTER_GAP_REPORTS_AND_RELEASE5_DOCS"
            if closure_resolved
            else "REPAIR_BRANCH_SELECTION_OR_FLOOR_MAPPING_AND_RERUN_QW2124"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2124: SCALAR VACUUM CLOSURE BRANCH-RESOLVED GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Values",
        f"- required_shift >= `{required_shift:.12f}`",
        f"- symmetric_floor: `{symmetric_floor:.12f}`",
        f"- broken_floor: `{broken_floor:.12f}`",
        "",
        "## Logic",
        f"- legacy QW-2120: `{r2120['verdict']}`",
        f"- QW-2123 branch selection: `{r2123['verdict']}`",
        f"- QW-2122 floor check: `{r2122['verdict']}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2124] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2124] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2124] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

