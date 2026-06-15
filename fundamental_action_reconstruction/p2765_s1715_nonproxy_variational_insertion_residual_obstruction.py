#!/usr/bin/env python3
"""P2765/S1715: nonproxy variational-insertion residual obstruction.

P2764 left one conditional provenance atom: a nonproxy EOM residual insertion
test for P1562-derived coupling candidates, explicitly conditional on the still
missing reference-cell/action-density, sign, and field/curvature normalization
theorems.  This audit checks whether current artifacts already supply the
nonproxy residual support needed for that insertion.  They do not: P1563 exports
formal/skeleton EOM strings, while P1866 keeps only proxy 1D EOM exports and
lists missing 4D covariant Einstein/Euler-Lagrange witness tables.
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
P2762 = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
P2763 = GEN / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.json"
P2764 = GEN / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.json"
OUT = GEN / "p2765_s1715_nonproxy_variational_insertion_residual_obstruction.json"
MD = GEN / "p2765_s1715_nonproxy_variational_insertion_residual_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2764_next_atom": r"P2764|variational-insertion|nonproxy EOM residual|P2697-P2764",
    "p1563_formal_eom": r"P1563|euler_lagrange_equations|lagrangian_density|E_full_variational_proof_log_machine_checkable",
    "p1866_proxy_block": r"P1866|OPEN_OBSTRUCTION_WITH_TRACE|proxy_1d|4D covariant|Einstein witness",
    "provenance_prerequisites": r"reference-cell/action-density|sign-convention|field/curvature normalization|physical-coupling provenance",
    "closure_boundary": r"role-bearing L_total|ToE closure|role transfer|bridge closure|selector closure",
}

NEGATIVE_EXPORT_FLAGS = [
    "nonproxy_variational_insertion_theorem_exported",
    "four_d_covariant_eom_residual_table_exported",
    "metric_einstein_residual_closed",
    "mixed_sector_residual_closed",
    "reference_cell_action_density_closure_imported",
    "sign_convention_theorem_exported",
    "field_curvature_normalization_theorem_exported",
    "physical_coupling_provenance_theorem_exported",
    "role_bearing_ltotal_promoted",
    "toe_closure_exported",
]

INSERTION_TARGETS = {
    "lambda_sm_eff": {
        "candidate_term": "lambda_sm_eff * phi^4 / scalar potential sector",
        "required_nonproxy_residuals": ["scalar_4d_euler_lagrange_residual"],
        "accepted_proxy_keys": ["eom_phi_proxy_1d"],
    },
    "kappa_gr_eff": {
        "candidate_term": "kappa_gr_eff * R / Einstein-Hilbert metric sector",
        "required_nonproxy_residuals": ["metric_4d_einstein_residual"],
        "accepted_proxy_keys": [],
    },
    "epsilon_mix_eff": {
        "candidate_term": "epsilon_mix_eff * phi^2*R or psi*R mixed sector",
        "required_nonproxy_residuals": ["mixed_scalar_metric_4d_residual", "metric_backreaction_residual"],
        "accepted_proxy_keys": [],
    },
}


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


def residual_support_matrix(p1563: dict[str, Any], p1866: dict[str, Any], p2762: dict[str, Any], p2763: dict[str, Any], p2764: dict[str, Any]) -> dict[str, Any]:
    p1563_eom_keys = sorted(p1563.get("euler_lagrange_equations", {}).keys())
    p1866_eom_keys = sorted(p1866.get("eom_export", {}).keys())
    blockers = p1866.get("strict_core_closure_blockers", [])
    rows = []
    for coupling, target in INSERTION_TARGETS.items():
        proxy_hits = [key for key in target["accepted_proxy_keys"] if key in p1866_eom_keys]
        nonproxy_exported = False
        rows.append({
            "coupling": coupling,
            "candidate_term": target["candidate_term"],
            "required_nonproxy_residuals": target["required_nonproxy_residuals"],
            "p1563_formal_eom_keys": p1563_eom_keys,
            "p1866_proxy_eom_hits": proxy_hits,
            "nonproxy_residual_table_exported": nonproxy_exported,
            "accepted_variational_insertion": False,
            "blocker": "Only formal/proxy EOM support is present; the required 4D nonproxy residual row is not exported under the still-open provenance prerequisites.",
        })
    missing_nonproxy_rows = sorted({req for target in INSERTION_TARGETS.values() for req in target["required_nonproxy_residuals"]})
    return {
        "p1563_formal_eom_keys": p1563_eom_keys,
        "p1866_eom_export_keys": p1866_eom_keys,
        "p1866_strict_core_closure_blockers": blockers,
        "rows": rows,
        "row_count": len(rows),
        "accepted_variational_insertion_count": sum(1 for row in rows if row["accepted_variational_insertion"]),
        "missing_nonproxy_residual_rows": missing_nonproxy_rows,
        "missing_nonproxy_residual_row_count": len(missing_nonproxy_rows),
        "has_proxy_eom_exports": bool(p1866_eom_keys),
        "has_4d_covariant_blocker": any("4D covariant" in str(blocker) for blocker in blockers),
        "p2762_reference_cell_action_density_still_open": p2762.get("acceptance_matrix", {}).get("accepted_as_reference_cell_action_density_theorem") is False,
        "p2763_sign_convention_still_open": p2763.get("acceptance_matrix", {}).get("accepted_as_sign_convention_theorem") is False,
        "p2764_field_curvature_normalization_still_open": p2764.get("acceptance_matrix", {}).get("accepted_as_field_curvature_normalization_theorem") is False,
    }


def acceptance_matrix(scan: dict[str, Any], matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "formal_eom_exports_present": bool(matrix["p1563_formal_eom_keys"]),
        "proxy_eom_exports_present": matrix["has_proxy_eom_exports"],
        "p1866_has_4d_covariant_blocker": matrix["has_4d_covariant_blocker"],
        "no_variational_insertion_accepted_without_nonproxy_rows": matrix["accepted_variational_insertion_count"] == 0,
        "reference_cell_action_density_still_open": matrix["p2762_reference_cell_action_density_still_open"],
        "sign_convention_still_open": matrix["p2763_sign_convention_still_open"],
        "field_curvature_normalization_still_open": matrix["p2764_field_curvature_normalization_still_open"],
        "four_d_covariant_eom_residual_table_exported": False,
        "nonproxy_variational_insertion_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_nonproxy_variational_insertion_theorem": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "P1563/P1866 provide formal/proxy EOM support, but the required 4D nonproxy scalar, metric, and mixed residual rows are absent; P2762-P2764 provenance prerequisites also remain open.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["nonproxy_variational_insertion_matrix"]
    lines = [
        "# P2765/S1715 nonproxy variational-insertion residual obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Residual support",
        f"- p1563_formal_eom_keys={m['p1563_formal_eom_keys']}",
        f"- p1866_eom_export_keys={m['p1866_eom_export_keys']}",
        f"- missing_nonproxy_residual_row_count={m['missing_nonproxy_residual_row_count']}",
        "",
        "## Rows",
    ]
    for row in m["rows"]:
        lines.append(f"- {row['coupling']}: accepted={row['accepted_variational_insertion']}; required={row['required_nonproxy_residuals']}; proxy_hits={row['p1866_proxy_eom_hits']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p1563 = read_json(P1563)
    p1866 = read_json(P1866)
    p2762 = read_json(P2762)
    p2763 = read_json(P2763)
    p2764 = read_json(P2764)
    scan = evidence_scan()
    matrix = residual_support_matrix(p1563, p1866, p2762, p2763, p2764)
    acceptance = acceptance_matrix(scan, matrix)
    payload = {
        "status": "P2765_NONPROXY_VARIATIONAL_INSERTION_RESIDUAL_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P1563": sha(P1563), "P1866": sha(P1866), "P2762": sha(P2762), "P2763": sha(P2763), "P2764": sha(P2764)},
        "input_statuses": {"P1562": p1562.get("status"), "P1563": p1563.get("status"), "P1866": p1866.get("status"), "P2762": p2762.get("status"), "P2763": p2763.get("status"), "P2764": p2764.get("status")},
        "audited_atom": "nonproxy EOM residual variational insertion, conditional on future provenance theorems",
        "content_evidence_scan": scan,
        "nonproxy_variational_insertion_matrix": matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote P1563/P1866 formal or proxy EOM exports to a nonproxy variational-insertion theorem.  The honest next move is a post-P2761-to-P2765 provenance-state reconciliation: either supply a genuinely new theorem fixing at least one of reference-cell/action-density, sign, field/curvature normalization, or 4D nonproxy residual rows, or preserve the P2697-P2765 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2765/S1715 nonproxy variational-insertion residual obstruction", "## P2765/S1715 nonproxy variational-insertion residual obstruction\n\n`P2765/S1715` audits the variational-insertion atom left after P2764.  P1563 exports formal/skeleton EOM strings and P1866 exports only proxy 1D EOM keys (`eom_phi_proxy_1d`, `eom_A_proxy_1d`) while explicitly listing missing 4D covariant Euler-Lagrange/Einstein witness tables.  The required nonproxy scalar, metric, and mixed residual rows for `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff` are absent, and P2762-P2764 provenance prerequisites remain open.  No nonproxy variational-insertion theorem, physical-coupling provenance theorem, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2765/S1715 nonproxy variational-insertion Ltotal guard", "## P2765/S1715 nonproxy variational-insertion Ltotal guard\n\n`P2765/S1715` adds no variational source term.  It records that formal/proxy EOM exports cannot become a nonproxy variational-insertion theorem while 4D covariant scalar, metric, and mixed residual rows are missing and P2762-P2764 provenance prerequisites remain open.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current nonproxy variational-insertion residual obstruction guardrail (P2765/S1715, 2026-06-15)", "## Current nonproxy variational-insertion residual obstruction guardrail (P2765/S1715, 2026-06-15)\n\n- P2765 audits the nonproxy variational-insertion residual atom left after P2764.\n- P1563 formal/skeleton EOM strings and P1866 proxy 1D EOM exports do not supply the required 4D nonproxy scalar, metric, and mixed residual rows; P2762 reference-cell/action-density, P2763 sign, and P2764 field/curvature normalization prerequisites remain open.\n- Do not promote this to physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  The next admissible move is a post-P2761-to-P2765 provenance-state reconciliation or a genuinely new theorem fixing one named provenance atom; otherwise preserve the P2697-P2765 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
