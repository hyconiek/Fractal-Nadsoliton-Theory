#!/usr/bin/env python3
"""P2764/S1714: field/curvature normalization compatibility obstruction.

P2763 left a narrow provenance atom: field/curvature normalization compatibility
with the still-open canonical reference-cell/action-density theorem.  This audit
checks whether current artifacts already fix the scalar, curvature, and mixed
field normalizations needed to interpret P1562 coefficients as physical
Lagrangian couplings.  They do not: independent positive field/curvature
renormalizations preserve the formal coefficient magnitudes but change the
physical insertion coefficients unless a common normalization theorem is added.
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
OUT = GEN / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.json"
MD = GEN / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2763_next_atom": r"P2763|field/curvature normalization compatibility|P2697-P2763|sign-convention",
    "p2762_reference_open": r"P2762|reference-cell/action-density|canonical reference cell|scale-orbit",
    "p1562_candidate_terms": r"P1562|lambda_sm_eff|kappa_gr_eff|epsilon_mix_eff|L_total",
    "field_curvature_language": r"field/curvature normalization|Z_phi|Z_A|Z_psi|curvature normalization|nonproxy EOM",
    "closure_boundary": r"role-bearing L_total|ToE closure|role transfer|bridge closure|selector closure",
}

NEGATIVE_EXPORT_FLAGS = [
    "field_curvature_normalization_theorem_exported",
    "common_normalization_compatibility_theorem_exported",
    "reference_cell_action_density_closure_imported",
    "sign_convention_theorem_exported",
    "physical_coupling_provenance_theorem_exported",
    "nonproxy_variational_insertion_theorem_exported",
    "role_bearing_ltotal_promoted",
    "toe_closure_exported",
]

NORMALIZATION_SECTORS = {
    "lambda_sm_eff": {
        "formal_term": "lambda_sm_eff * phi^4",
        "normalization_factors": {"Z_phi": -4},
        "blocker": "A scalar field normalization theorem is needed before a quartic coefficient becomes physical.",
    },
    "kappa_gr_eff": {
        "formal_term": "kappa_gr_eff * R",
        "normalization_factors": {"Z_R": -1},
        "blocker": "A curvature/metric normalization theorem is needed before an Einstein-Hilbert coefficient becomes physical.",
    },
    "epsilon_mix_eff": {
        "formal_term": "epsilon_mix_eff * phi^2 * R",
        "normalization_factors": {"Z_phi": -2, "Z_R": -1},
        "blocker": "A common scalar-curvature normalization compatibility theorem is needed before the mixed term becomes physical.",
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


def apply_factor(value: float, factors: dict[str, int], z: dict[str, float]) -> float:
    out = value
    for name, power in factors.items():
        out *= z[name] ** power
    return out


def finite_normalization_orbit(coeffs: dict[str, float]) -> dict[str, Any]:
    samples = [
        {"Z_phi": 0.5, "Z_R": 0.5},
        {"Z_phi": 0.5, "Z_R": 2.0},
        {"Z_phi": 1.0, "Z_R": 1.0},
        {"Z_phi": 2.0, "Z_R": 0.5},
        {"Z_phi": 2.0, "Z_R": 2.0},
    ]
    rows = []
    for z in samples:
        normalized = {}
        for name, sector in NORMALIZATION_SECTORS.items():
            normalized[name] = apply_factor(coeffs[name], sector["normalization_factors"], z)
        rows.append({"normalization": z, "normalized_coefficients": normalized})
    distinct = len({tuple(round(v, 12) for v in row["normalized_coefficients"].values()) for row in rows})
    per_coeff_varies = {
        name: len({round(row["normalized_coefficients"][name], 12) for row in rows}) > 1
        for name in NORMALIZATION_SECTORS
    }
    return {
        "sampled_positive_normalizations": samples,
        "rows": rows,
        "distinct_normalized_coefficient_triples": distinct,
        "per_coefficient_varies_under_normalization": per_coeff_varies,
        "all_coefficients_normalization_sensitive": all(per_coeff_varies.values()),
    }


def compatibility_matrix(p1562: dict[str, Any], p2762: dict[str, Any], p2763: dict[str, Any], p1563: dict[str, Any], p1866: dict[str, Any]) -> dict[str, Any]:
    coeffs = {k: float(v) for k, v in p1562.get("derived_lagrangian_coefficients", {}).items()}
    orbit = finite_normalization_orbit(coeffs)
    rows = []
    for name, sector in NORMALIZATION_SECTORS.items():
        rows.append({
            "coupling": name,
            "formal_term": sector["formal_term"],
            "normalization_factors": sector["normalization_factors"],
            "normalization_sensitive": orbit["per_coefficient_varies_under_normalization"][name],
            "field_or_curvature_normalization_exported": False,
            "compatible_with_reference_cell_action_density_theorem": False,
            "accepted_as_physical_normalization": False,
            "blocker": sector["blocker"],
        })
    return {
        "derived_lagrangian_coefficients": coeffs,
        "finite_normalization_orbit": orbit,
        "rows": rows,
        "row_count": len(rows),
        "accepted_normalization_row_count": sum(1 for row in rows if row["accepted_as_physical_normalization"]),
        "p2762_reference_cell_action_density_still_open": p2762.get("acceptance_matrix", {}).get("accepted_as_reference_cell_action_density_theorem") is False,
        "p2763_sign_convention_still_open": p2763.get("acceptance_matrix", {}).get("accepted_as_sign_convention_theorem") is False,
        "later_nonproxy_eom_closure_still_open": p1563.get("toe_closed") is False and p1866.get("status") == "OPEN_OBSTRUCTION_WITH_TRACE",
    }


def acceptance_matrix(scan: dict[str, Any], matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "p1562_candidate_coefficients_present": bool(matrix["derived_lagrangian_coefficients"]),
        "finite_normalization_orbit_exhibits_nonuniqueness": matrix["finite_normalization_orbit"]["distinct_normalized_coefficient_triples"] > 1,
        "all_rows_normalization_sensitive": matrix["finite_normalization_orbit"]["all_coefficients_normalization_sensitive"],
        "p2762_reference_cell_action_density_still_open": matrix["p2762_reference_cell_action_density_still_open"],
        "p2763_sign_convention_still_open": matrix["p2763_sign_convention_still_open"],
        "later_nonproxy_eom_closure_still_open": matrix["later_nonproxy_eom_closure_still_open"],
        "field_curvature_normalization_theorem_exported": False,
        "common_normalization_compatibility_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_field_curvature_normalization_theorem": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The finite normalization orbit shows that positive scalar/curvature normalizations change all three candidate coefficients. Current artifacts do not export a common field-curvature normalization theorem, and P2762/P2763 keep reference-cell and sign provenance open.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["field_curvature_normalization_matrix"]
    lines = [
        "# P2764/S1714 field-curvature normalization compatibility obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite normalization orbit",
        f"- distinct_normalized_coefficient_triples={m['finite_normalization_orbit']['distinct_normalized_coefficient_triples']}",
        f"- all_coefficients_normalization_sensitive={m['finite_normalization_orbit']['all_coefficients_normalization_sensitive']}",
        "",
        "## Rows",
    ]
    for row in m["rows"]:
        lines.append(f"- {row['coupling']}: term={row['formal_term']}; sensitive={row['normalization_sensitive']}; accepted={row['accepted_as_physical_normalization']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p1563 = read_json(P1563)
    p1866 = read_json(P1866)
    p2762 = read_json(P2762)
    p2763 = read_json(P2763)
    scan = evidence_scan()
    matrix = compatibility_matrix(p1562, p2762, p2763, p1563, p1866)
    acceptance = acceptance_matrix(scan, matrix)
    payload = {
        "status": "P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P1563": sha(P1563), "P1866": sha(P1866), "P2762": sha(P2762), "P2763": sha(P2763)},
        "input_statuses": {"P1562": p1562.get("status"), "P1563": p1563.get("status"), "P1866": p1866.get("status"), "P2762": p2762.get("status"), "P2763": p2763.get("status")},
        "audited_atom": "field/curvature normalization compatibility with still-open canonical reference-cell/action-density theorem",
        "content_evidence_scan": scan,
        "field_curvature_normalization_matrix": matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote field or curvature rescaling conventions to physical coupling provenance.  The next honest proof-grade move is the remaining variational-insertion atom: a nonproxy EOM residual insertion test explicitly conditional on future reference-cell/action-density, sign, and field/curvature normalization theorems.  If that cannot be supplied, preserve the P2697-P2764 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2764/S1714 field-curvature normalization compatibility obstruction", "## P2764/S1714 field-curvature normalization compatibility obstruction\n\n`P2764/S1714` audits the field/curvature normalization atom left after P2763.  A finite positive-normalization orbit over scalar and curvature factors changes all three candidate coefficient rows (`lambda_sm_eff`, `kappa_gr_eff`, `epsilon_mix_eff`), showing that current artifacts do not export a common normalization compatibility theorem for the formal Lagrangian terms.  P2762 reference-cell/action-density closure and P2763 sign-convention closure remain open.  No field/curvature normalization theorem, physical-coupling provenance theorem, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2764/S1714 field-curvature normalization Ltotal guard", "## P2764/S1714 field-curvature normalization Ltotal guard\n\n`P2764/S1714` adds no variational source term.  It records that scalar and curvature rescalings change the formal P1562 coefficient insertions unless a common field/curvature normalization theorem is exported together with the still-missing reference-cell/action-density and sign-convention provenance.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current field-curvature normalization compatibility obstruction guardrail (P2764/S1714, 2026-06-15)", "## Current field-curvature normalization compatibility obstruction guardrail (P2764/S1714, 2026-06-15)\n\n- P2764 audits field/curvature normalization compatibility with the still-open canonical reference-cell/action-density theorem.\n- A finite positive-normalization orbit shows that scalar and curvature rescalings change all three candidate coefficient rows (`lambda_sm_eff`, `kappa_gr_eff`, `epsilon_mix_eff`); no common normalization compatibility theorem is exported, and P2762/P2763 remain open.\n- Do not promote these conventions to physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  The next admissible move is a nonproxy variational-insertion residual test explicitly conditional on future provenance theorems, or preservation of the P2697-P2764 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
