#!/usr/bin/env python3
"""
QW-2123: Vacuum-branch selection strict gate.

Purpose:
- resolve QW-2122 branch ambiguity with an explicit strict rule:
  if K_total has a negative mode (QW-2118), symmetric branch is not
  the physical strict closure branch and broken branch must be used.
- keep verdict tied only to strict-chain artifacts (no exploratory channel).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2123_vacuum_branch_selection_strict_gate.json"
OUT_MD = ROOT / "RAPORT_QW2123_VACUUM_BRANCH_SELECTION_STRICT_GATE.md"


def load_json(name: str) -> Dict[str, Any]:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def _bool(x: Any) -> bool:
    return bool(x)


def main() -> None:
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")
    r2120 = load_json("report_qw2120_scalar_scale_vacuum_closure_strict_gate.json")
    r2122 = load_json("report_qw2122_psi_potential_diagonal_floor_gate.json")

    eig = r2118["eigen_spectrum"]
    lambda_min = float(eig["lambda_min"])
    required_shift = float(r2118["vacuum_shift_condition"]["required_uniform_mass_shift_ge"])

    strict_floor_2120 = float(r2120["strict_floor_components"]["strict_floor_used_for_verdict"])
    symmetric_floor = float(r2122["branch_results"]["symmetric_branch_reference"]["diag_floor"])
    broken_floor = float(r2122["branch_results"]["broken_branch_strict"]["diag_floor"])
    rho_star_sq = float(r2122["branch_results"]["broken_branch_strict"]["rho_star_sq"])

    branch_rule_satisfied = bool(
        lambda_min < 0.0
        and symmetric_floor < required_shift
        and broken_floor >= required_shift
        and rho_star_sq > 0.0
    )

    flags = {
        "q2118_conditional_vacuum_shift_loaded": _bool(
            str(r2118.get("verdict", "")).endswith("CONDITIONAL_VACUUM_SHIFT")
        ),
        "k_total_has_negative_mode": _bool(lambda_min < 0.0),
        "required_shift_positive": _bool(required_shift > 0.0),
        "symmetric_floor_matches_qw2120_floor": _bool(abs(symmetric_floor - strict_floor_2120) <= 1e-12),
        "symmetric_floor_below_required_shift": _bool(symmetric_floor < required_shift),
        "broken_floor_ge_required_shift": _bool(broken_floor >= required_shift),
        "broken_branch_vacuum_defined": _bool(rho_star_sq > 0.0),
        "q2122_deterministic_no_scan_no_retune": _bool(
            r2122["flags"].get("deterministic_no_scan_no_retune", False)
        ),
        "branch_selection_rule_satisfied": _bool(branch_rule_satisfied),
        "selection_independent_of_exploratory_channel": True,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "VACUUM_BRANCH_SELECTION_STRICT_GATE_PASS_BROKEN_BRANCH_REQUIRED"
        if branch_rule_satisfied
        else "VACUUM_BRANCH_SELECTION_STRICT_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2118": "report_qw2118_ktotal_spectral_tripartition_gate.json",
            "q2120": "report_qw2120_scalar_scale_vacuum_closure_strict_gate.json",
            "q2122": "report_qw2122_psi_potential_diagonal_floor_gate.json",
        },
        "inputs": {
            "lambda_min_k_total": lambda_min,
            "required_shift_ge": required_shift,
            "strict_floor_qw2120": strict_floor_2120,
            "symmetric_floor_qw2122": symmetric_floor,
            "broken_floor_qw2122": broken_floor,
            "broken_branch_rho_star_sq": rho_star_sq,
        },
        "selection_rule": {
            "definition": (
                "If lambda_min(K_total)<0 and symmetric_floor<required_shift<=broken_floor, "
                "select broken branch as strict physical closure branch."
            ),
            "result": "broken_branch_required" if branch_rule_satisfied else "not_resolved",
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "RUN_QW2124_SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_GATE"
            if branch_rule_satisfied
            else "REPAIR_BRANCH_SELECTION_RULE_INPUTS_AND_RERUN_QW2123"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2123: VACUUM BRANCH SELECTION STRICT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Inputs",
        f"- lambda_min(K_total): `{lambda_min:.12f}`",
        f"- required_shift >= `{required_shift:.12f}`",
        f"- symmetric_floor: `{symmetric_floor:.12f}`",
        f"- broken_floor: `{broken_floor:.12f}`",
        "",
        "## Selection rule",
        f"- {out['selection_rule']['definition']}",
        f"- result: `{out['selection_rule']['result']}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2123] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2123] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2123] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

