#!/usr/bin/env python3
"""P2761/S1711: kernel-moment -> physical-coupling provenance obstruction.

P2760 selected G5 as the most honest next target.  This audit tests whether the
P1562 strict-kernel moments already suffice to promote lambda_sm_eff,
kappa_gr_eff, and epsilon_mix_eff into physical Lagrangian couplings.  The
answer is no: the numerical moment map is real, but units/reference scales,
sign conventions, and variational insertion normalizations are not exported.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
P1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
P1866 = GEN / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json"
P2760 = GEN / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.json"
OUT = GEN / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.json"
MD = GEN / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2760_g5_obligation": r"G5|kernel-moments-to-physical-couplings|physical-couplings provenance|units/sign/variational",
    "p1562_moment_map": r"P1562|lambda_sm_eff|kappa_gr_eff|epsilon_mix_eff|derived_lagrangian_coefficients",
    "p1563_missing_exports": r"P1563|missing_exports|missing_theorems|qw2191_closed|toe_closed",
    "p1866_open_obstruction": r"P1866|OPEN_OBSTRUCTION_WITH_TRACE|reverse_requirement|strict_core_closure_blockers",
    "ltotal_boundary": r"L_total|role-bearing L_total|ToE closure|selector closure|role transfer",
}

NEGATIVE_EXPORT_FLAGS = [
    "physical_coupling_provenance_theorem_exported",
    "unit_reference_source_exported",
    "sign_convention_theorem_exported",
    "variational_normalization_exported",
    "p1562_stale_flags_revalidated",
    "role_bearing_ltotal_promoted",
    "toe_closure_exported",
]

COUPLING_TARGETS = {
    "lambda_sm_eff": {
        "expected_mass_dimension_4d": 0,
        "variational_role": "quartic/scalar-sector dimensionless effective coupling",
    },
    "kappa_gr_eff": {
        "expected_mass_dimension_4d": 2,
        "variational_role": "Einstein-Hilbert R coefficient / inverse gravitational coupling scale",
    },
    "epsilon_mix_eff": {
        "expected_mass_dimension_4d": 1,
        "variational_role": "schematic psi*R or mixed-sector coefficient requiring field and curvature normalization",
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


def infer_ratio_candidates(moments: dict[str, float], coeffs: dict[str, float]) -> dict[str, Any]:
    m0 = moments["M0"]
    ratios = {
        "M1_over_M0": moments["M1"] / m0,
        "M2_over_M0": moments["M2"] / m0,
        "M3_over_M0": moments["M3"] / m0,
        "abs_M2_over_M0": abs(moments["M2"] / m0),
        "abs_M3_over_M0": abs(moments["M3"] / m0),
        "sqrt_abs_M2_over_M0": math.sqrt(abs(moments["M2"] / m0)),
        "sqrt_abs_M3_over_M0": math.sqrt(abs(moments["M3"] / m0)),
    }
    matches: dict[str, list[dict[str, Any]]] = {name: [] for name in coeffs}
    for cname, cval in coeffs.items():
        for rname, rval in ratios.items():
            if abs(cval - rval) < 1e-10:
                matches[cname].append({"ratio": rname, "value": rval, "transform": "identity"})
            if abs(cval + rval) < 1e-10:
                matches[cname].append({"ratio": rname, "value": rval, "transform": "minus"})
    return {"ratios": ratios, "exact_simple_ratio_matches": matches}


def dimension_obligation_matrix(p1562: dict[str, Any], p1563: dict[str, Any], p1866: dict[str, Any]) -> dict[str, Any]:
    moments = {k: float(v) for k, v in p1562.get("moments", {}).items()}
    coeffs = {k: float(v) for k, v in p1562.get("derived_lagrangian_coefficients", {}).items()}
    ratio_audit = infer_ratio_candidates(moments, coeffs)
    rows = []
    missing_global = [
        "canonical physical length/reference cell mapping for moment powers",
        "action-density normalization and hbar/c unit convention",
        "field normalization Z_phi/Z_A/Z_psi and curvature normalization",
        "sign convention theorem for negative moments and positive couplings",
        "variational insertion theorem tying coefficients to nonproxy EOM residuals",
    ]
    for name, meta in COUPLING_TARGETS.items():
        rows.append(
            {
                "coupling": name,
                "p1562_numeric_value": coeffs.get(name),
                "expected_mass_dimension_4d": meta["expected_mass_dimension_4d"],
                "variational_role": meta["variational_role"],
                "simple_ratio_matches": ratio_audit["exact_simple_ratio_matches"].get(name, []),
                "unit_source_exported": False,
                "sign_convention_exported": False,
                "variational_normalization_exported": False,
                "accepted_as_physical_coupling": False,
                "blocker": "P1562 supplies a number, but current artifacts do not export the reference unit, sign convention, and variational normalization needed to interpret it as this physical coupling.",
            }
        )
    stale_flag_detected = bool(p1562.get("qw2191_closed")) or bool(p1562.get("toe_closed"))
    later_blocks = p1563.get("qw2191_closed") is False and p1563.get("toe_closed") is False and p1866.get("status") == "OPEN_OBSTRUCTION_WITH_TRACE"
    return {
        "moments": moments,
        "derived_lagrangian_coefficients": coeffs,
        "ratio_audit": ratio_audit,
        "rows": rows,
        "row_count": len(rows),
        "accepted_physical_coupling_count": sum(1 for row in rows if row["accepted_as_physical_coupling"]),
        "missing_global_provenance_atoms": missing_global,
        "missing_global_provenance_atom_count": len(missing_global),
        "p1562_stale_closure_flag_detected": stale_flag_detected,
        "later_artifacts_block_closure": later_blocks,
        "stale_flags_quarantined": stale_flag_detected and later_blocks,
    }


def acceptance_matrix(scan: dict[str, Any], matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "g5_content_evidence_present": scan["all_patterns_have_hits"],
        "p1562_moment_numbers_present": bool(matrix["moments"]) and bool(matrix["derived_lagrangian_coefficients"]),
        "all_three_couplings_audited": matrix["row_count"] == 3,
        "no_coupling_accepted_without_provenance": matrix["accepted_physical_coupling_count"] == 0,
        "stale_p1562_flags_quarantined": matrix["stale_flags_quarantined"],
        "unit_reference_source_exported": False,
        "sign_convention_theorem_exported": False,
        "variational_normalization_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_g5_physical_coupling_provenance_theorem": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The moment-to-coupling numerics are present, but every target coupling lacks exported unit/reference, sign, and variational-normalization provenance; stale P1562 closure flags remain quarantined by later artifacts.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["dimension_obligation_matrix"]
    lines = [
        "# P2761/S1711 kernel-moment to physical-coupling provenance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Moment/coupling audit",
        f"- moments={m['moments']}",
        f"- derived_lagrangian_coefficients={m['derived_lagrangian_coefficients']}",
        f"- row_count={m['row_count']}",
        f"- accepted_physical_coupling_count={m['accepted_physical_coupling_count']}",
        f"- missing_global_provenance_atom_count={m['missing_global_provenance_atom_count']}",
        f"- stale_flags_quarantined={m['stale_flags_quarantined']}",
        "",
        "## Coupling rows",
    ]
    for row in m["rows"]:
        lines.append(f"- {row['coupling']}: accepted={row['accepted_as_physical_coupling']}; expected_mass_dimension_4d={row['expected_mass_dimension_4d']}; blocker={row['blocker']}")
    lines.extend(["", "## Missing global provenance atoms"])
    for atom in m["missing_global_provenance_atoms"]:
        lines.append(f"- {atom}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p1563 = read_json(P1563)
    p1866 = read_json(P1866)
    p2760 = read_json(P2760)
    scan = evidence_scan()
    matrix = dimension_obligation_matrix(p1562, p1563, p1866)
    acceptance = acceptance_matrix(scan, matrix)
    payload = {
        "status": "P2761_KERNEL_MOMENT_COUPLING_PROVENANCE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P1563": sha(P1563), "P1866": sha(P1866), "P2760": sha(P2760)},
        "input_statuses": {"P1562": p1562.get("status"), "P1563": p1563.get("status"), "P1866": p1866.get("status"), "P2760": p2760.get("status")},
        "audited_gap": "P2760 G5 kernel-moments-to-physical-couplings provenance",
        "content_evidence_scan": scan,
        "dimension_obligation_matrix": matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote P1562 moment-derived numbers into physical couplings yet.  The next proof-grade move should target exactly one missing global provenance atom: construct a canonical physical length/reference-cell and action-density normalization theorem for the strict moment map; only after that rerun the sign and variational-normalization rows.  Without that theorem, preserve the P2697-P2761 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2761/S1711 kernel-moment physical-coupling provenance obstruction", "## P2761/S1711 kernel-moment physical-coupling provenance obstruction\n\n`P2761/S1711` attacks the P2760 `G5` target directly.  It audits P1562 moments and derived coefficients `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff` against the minimum provenance needed for physical Lagrangian couplings: canonical reference units, action-density normalization, field/curvature normalization, sign conventions, and variational insertion into nonproxy EOM residuals.  The numeric moment map is present, but all three coupling rows fail provenance acceptance; stale P1562 `qw2191_closed/toe_closed` flags remain quarantined by P1563/P1866/current guardrails.  No physical-coupling provenance theorem, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2761/S1711 kernel-moment coupling Ltotal guard", "## P2761/S1711 kernel-moment coupling Ltotal guard\n\n`P2761/S1711` adds no variational source term.  It shows that P1562 moment-derived coefficients remain numeric candidates rather than accepted physical Lagrangian couplings because unit/reference, sign, and variational-normalization provenance is missing.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current kernel-moment physical-coupling provenance obstruction guardrail (P2761/S1711, 2026-06-15)", "## Current kernel-moment physical-coupling provenance obstruction guardrail (P2761/S1711, 2026-06-15)\n\n- P2761 attacks P2760 `G5` and audits P1562 moment-derived coefficients `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff` as physical-coupling candidates.\n- The numeric moment map is real, but no canonical reference unit/action-density normalization, sign convention theorem, field/curvature normalization, or nonproxy variational insertion theorem is exported; stale P1562 closure flags remain quarantined.\n- Do not promote these coefficients to role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  The next admissible move is one missing provenance atom, preferably a canonical physical length/reference-cell and action-density normalization theorem for the strict moment map.\n")
    return payload


if __name__ == "__main__":
    main()
