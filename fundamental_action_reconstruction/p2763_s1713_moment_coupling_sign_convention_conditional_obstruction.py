#!/usr/bin/env python3
"""P2763/S1713: moment-to-coupling sign-convention conditional obstruction.

P2762 left a narrow admissible follow-up: test the sign-convention atom, but
only conditionally on a future canonical reference-cell/action-density theorem.
This audit asks whether current artifacts already select the signs by which
P1562's signed moments become positive Lagrangian coupling candidates.  They do
not: M2 and M3 are negative while kappa_gr_eff and epsilon_mix_eff are positive,
so a sign-rectification convention is being used, but no exported theorem
chooses that convention over the equally formal opposite-sign branches.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
P1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
P1866 = GEN / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json"
P2761 = GEN / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.json"
P2762 = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
OUT = GEN / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.json"
MD = GEN / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2762_conditional_next": r"P2762|sign-convention theorem|canonical reference-cell/action-density|P2697-P2762",
    "p2761_sign_gap": r"P2761|sign convention theorem|unit/reference|variational-normalization",
    "p1562_signed_moments": r"P1562|M2|M3|kappa_gr_eff|epsilon_mix_eff|derived_lagrangian_coefficients",
    "later_eom_block": r"P1563|P1866|nonproxy EOM|OPEN_OBSTRUCTION_WITH_TRACE|missing_theorems",
    "closure_boundary": r"L_total|ToE closure|role transfer|bridge closure|selector closure",
}

NEGATIVE_EXPORT_FLAGS = [
    "sign_convention_theorem_exported",
    "unique_coupling_sign_branch_selected",
    "reference_cell_action_density_closure_imported",
    "physical_coupling_provenance_theorem_exported",
    "nonproxy_variational_insertion_theorem_exported",
    "role_bearing_ltotal_promoted",
    "toe_closure_exported",
]

SIGN_BRANCH_TARGETS = ["kappa_gr_eff", "epsilon_mix_eff"]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {"row_count": len(rows), "rows": rows, "hit_counts": {r["lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def sign_branch_scan(p1562: dict[str, Any]) -> dict[str, Any]:
    moments = {k: float(v) for k, v in p1562.get("moments", {}).items()}
    coeffs = {k: float(v) for k, v in p1562.get("derived_lagrangian_coefficients", {}).items()}
    branch_rows = []
    for s_kappa in (-1, 1):
        for s_epsilon in (-1, 1):
            branch_rows.append({
                "branch": {"kappa_gr_eff_sign": s_kappa, "epsilon_mix_eff_sign": s_epsilon},
                "lambda_sm_eff": coeffs["lambda_sm_eff"],
                "kappa_gr_eff": s_kappa * abs(coeffs["kappa_gr_eff"]),
                "epsilon_mix_eff": s_epsilon * abs(coeffs["epsilon_mix_eff"]),
                "preserves_numeric_magnitudes": True,
            })
    rows = [
        {
            "object": "lambda_sm_eff",
            "source_moment_sign": 1 if moments["M1"] > 0 else -1,
            "candidate_coupling_sign": 1 if coeffs["lambda_sm_eff"] > 0 else -1,
            "sign_rectification_needed": False,
            "unique_sign_theorem_exported": False,
            "blocker": "The dimensionless row is sign-aligned, but P2762 still blocks physical normalization/provenance.",
        },
        {
            "object": "kappa_gr_eff",
            "source_moment_sign": 1 if moments["M2"] > 0 else -1,
            "candidate_coupling_sign": 1 if coeffs["kappa_gr_eff"] > 0 else -1,
            "sign_rectification_needed": moments["M2"] * coeffs["kappa_gr_eff"] < 0,
            "unique_sign_theorem_exported": False,
            "blocker": "P1562 maps a negative moment-side datum to a positive coupling candidate, but current artifacts do not export the sign convention theorem choosing this branch.",
        },
        {
            "object": "epsilon_mix_eff",
            "source_moment_sign": 1 if moments["M3"] > 0 else -1,
            "candidate_coupling_sign": 1 if coeffs["epsilon_mix_eff"] > 0 else -1,
            "sign_rectification_needed": moments["M3"] * coeffs["epsilon_mix_eff"] < 0,
            "unique_sign_theorem_exported": False,
            "blocker": "P1562 maps a negative moment-side datum to a positive coupling candidate, but current artifacts do not export the sign convention theorem choosing this branch.",
        },
    ]
    return {
        "moments": moments,
        "derived_lagrangian_coefficients": coeffs,
        "rows": rows,
        "row_count": len(rows),
        "rectified_row_count": sum(1 for row in rows if row["sign_rectification_needed"]),
        "branch_family": branch_rows,
        "branch_family_count": len(branch_rows),
        "all_branches_preserve_magnitudes": all(row["preserves_numeric_magnitudes"] for row in branch_rows),
        "unique_branch_selected_by_current_artifacts": False,
    }


def conditional_acceptance_matrix(scan: dict[str, Any], sign_matrix: dict[str, Any], p2762: dict[str, Any], p1563: dict[str, Any], p1866: dict[str, Any]) -> dict[str, Any]:
    p2762_blocks_reference = p2762.get("acceptance_matrix", {}).get("accepted_as_reference_cell_action_density_theorem") is False
    later_eom_block = p1563.get("toe_closed") is False and p1866.get("status") == "OPEN_OBSTRUCTION_WITH_TRACE"
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "signed_moments_and_coefficients_present": bool(sign_matrix["moments"]) and bool(sign_matrix["derived_lagrangian_coefficients"]),
        "negative_moment_to_positive_coupling_rectification_detected": sign_matrix["rectified_row_count"] == 2,
        "finite_opposite_sign_branch_family_present": sign_matrix["branch_family_count"] == 4 and sign_matrix["all_branches_preserve_magnitudes"],
        "p2762_reference_cell_action_density_still_open": p2762_blocks_reference,
        "later_nonproxy_eom_closure_still_open": later_eom_block,
        "unique_sign_branch_selected_by_current_artifacts": False,
        "sign_convention_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_sign_convention_theorem": False,
        "accepted_as_conditional_theorem": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "A finite four-branch sign family preserves the P1562 coupling magnitudes, while current artifacts do not export a rule selecting the positive kappa/epsilon branch. Because P2762 leaves reference-cell/action-density normalization open, even a future sign rule would have to be explicitly conditional rather than physical-coupling closure.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    s = payload["sign_convention_matrix"]
    lines = [
        "# P2763/S1713 moment-to-coupling sign-convention conditional obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Sign rows",
    ]
    for row in s["rows"]:
        lines.append(f"- {row['object']}: source_sign={row['source_moment_sign']}; coupling_sign={row['candidate_coupling_sign']}; rectification={row['sign_rectification_needed']}; theorem_exported={row['unique_sign_theorem_exported']}")
    lines.extend([
        "",
        "## Finite branch family",
        f"- branch_family_count={s['branch_family_count']}",
        f"- all_branches_preserve_magnitudes={s['all_branches_preserve_magnitudes']}",
        f"- unique_branch_selected_by_current_artifacts={s['unique_branch_selected_by_current_artifacts']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p1563 = read_json(P1563)
    p1866 = read_json(P1866)
    p2761 = read_json(P2761)
    p2762 = read_json(P2762)
    scan = evidence_scan()
    sign_matrix = sign_branch_scan(p1562)
    acceptance = conditional_acceptance_matrix(scan, sign_matrix, p2762, p1563, p1866)
    payload = {
        "status": "P2763_MOMENT_COUPLING_SIGN_CONVENTION_CONDITIONAL_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P1563": sha(P1563), "P1866": sha(P1866), "P2761": sha(P2761), "P2762": sha(P2762)},
        "input_statuses": {"P1562": p1562.get("status"), "P1563": p1563.get("status"), "P1866": p1866.get("status"), "P2761": p2761.get("status"), "P2762": p2762.get("status")},
        "audited_atom": "moment-to-coupling sign-convention theorem, explicitly conditional on future reference-cell/action-density closure",
        "content_evidence_scan": scan,
        "sign_convention_matrix": sign_matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the positive P1562 kappa/epsilon signs to physical coupling signs.  The next honest proof-grade move should target the remaining provenance atom that is not merely sign choice: field/curvature normalization compatibility with the still-open canonical reference-cell/action-density theorem.  If that object cannot be supplied, preserve the P2697-P2763 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2763/S1713 moment-coupling sign-convention conditional obstruction", "## P2763/S1713 moment-coupling sign-convention conditional obstruction\n\n`P2763/S1713` audits the sign-convention atom left after P2762.  The P1562 moment side has `M2<0` and `M3<0`, while the derived candidates `kappa_gr_eff` and `epsilon_mix_eff` are positive; a finite four-branch sign family preserves the same magnitudes but flips the two dimensionful coupling signs.  Current artifacts do not export a theorem selecting the positive branch, and P2762 still leaves canonical reference-cell/action-density normalization open.  No sign-convention theorem, physical-coupling provenance theorem, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2763/S1713 sign-convention Ltotal guard", "## P2763/S1713 sign-convention Ltotal guard\n\n`P2763/S1713` adds no variational source term.  It records that positive `kappa_gr_eff` and `epsilon_mix_eff` are sign-branch candidates, not physical Lagrangian sign theorems, because no branch-selecting sign convention is exported and P2762 leaves reference-cell/action-density normalization open.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current moment-coupling sign-convention conditional obstruction guardrail (P2763/S1713, 2026-06-15)", "## Current moment-coupling sign-convention conditional obstruction guardrail (P2763/S1713, 2026-06-15)\n\n- P2763 audits the P2762-recommended sign-convention atom only as a conditional provenance test.\n- P1562 has negative moment-side rows `M2` and `M3` but positive `kappa_gr_eff` and `epsilon_mix_eff`; the finite four-branch sign family preserves the same magnitudes, and current artifacts do not export a theorem selecting the positive branch.\n- Do not promote the signs or coefficients to physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  The next admissible move is field/curvature normalization compatibility with the still-open canonical reference-cell/action-density theorem, or preservation of the P2697-P2763 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
