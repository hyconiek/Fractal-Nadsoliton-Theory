#!/usr/bin/env python3
"""P2762/S1712: reference-cell/action-density normalization obstruction.

P2761 left one preferred atom: a canonical physical length/reference-cell plus
action-density normalization theorem for the strict moment map.  This audit
checks whether current artifacts already determine such a theorem.  They do
not: the moment numerics are present, but the physical length scale, action
cell volume, action-density unit, and field/curvature normalizations remain a
continuous gauge family rather than a unique exported normalization.
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
P2689 = GEN / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.json"
P2760 = GEN / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.json"
P2761 = GEN / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.json"
OUT = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
MD = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2761_preferred_atom": r"P2761|canonical physical length/reference-cell|action-density normalization theorem|strict moment map",
    "p2689_uv_unit_block": r"P2689|entropy reference-cell|bit-to-length|canonical length/UV unit|scale-orbit quotient",
    "p1562_moment_coefficients": r"P1562|M0|M1|M2|M3|lambda_sm_eff|kappa_gr_eff|epsilon_mix_eff",
    "normalization_boundary": r"action-density normalization|field/curvature normalization|variational insertion|nonproxy EOM",
    "closure_boundary": r"L_total|ToE closure|role transfer|bridge closure|selector closure",
}

NEGATIVE_EXPORT_FLAGS = [
    "canonical_reference_cell_theorem_exported",
    "action_density_normalization_theorem_exported",
    "unique_physical_length_scale_exported",
    "moment_map_physical_units_fixed",
    "physical_coupling_provenance_theorem_exported",
    "role_bearing_ltotal_promoted",
    "toe_closure_exported",
]

COUPLING_ROWS = {
    "lambda_sm_eff": {"expected_mass_dimension_4d": 0, "minimal_unknowns": ["ell", "A_action", "Z_phi"]},
    "kappa_gr_eff": {"expected_mass_dimension_4d": 2, "minimal_unknowns": ["ell", "A_action", "Z_g"]},
    "epsilon_mix_eff": {"expected_mass_dimension_4d": 1, "minimal_unknowns": ["ell", "A_action", "Z_psi", "Z_R"]},
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


def finite_scale_orbit_witness(coeffs: dict[str, float]) -> dict[str, Any]:
    """Show numerically that a reference length changes dimensionalized rows.

    This is not a physics theorem; it is a finite witness that the exported
    numbers alone do not choose one physical unit.  We use the expected 4d mass
    dimensions from P2761: dimensionless, mass^2, mass^1.  With length ell,
    mass^d corresponds to ell^-d, so every positive ell gives a distinct but
    internally consistent dimensionalization unless a canonical ell theorem is
    supplied.
    """
    scales = [0.25, 0.5, 1.0, 2.0, 4.0]
    rows = []
    for ell in scales:
        dimmed = {}
        for name, meta in COUPLING_ROWS.items():
            d = meta["expected_mass_dimension_4d"]
            dimmed[name] = coeffs[name] * (ell ** (-d))
        rows.append({"ell": ell, "dimensionalized_coefficients": dimmed})
    kappa_values = [r["dimensionalized_coefficients"]["kappa_gr_eff"] for r in rows]
    epsilon_values = [r["dimensionalized_coefficients"]["epsilon_mix_eff"] for r in rows]
    return {
        "tested_positive_reference_lengths": scales,
        "rows": rows,
        "lambda_invariant_under_ell": all(abs(r["dimensionalized_coefficients"]["lambda_sm_eff"] - rows[0]["dimensionalized_coefficients"]["lambda_sm_eff"]) < 1e-12 for r in rows),
        "kappa_changes_with_ell": max(kappa_values) - min(kappa_values) > 1e-12,
        "epsilon_changes_with_ell": max(epsilon_values) - min(epsilon_values) > 1e-12,
        "distinct_dimensionalizations": len({tuple(round(v, 12) for v in r["dimensionalized_coefficients"].values()) for r in rows}),
    }


def normalization_matrix(p1562: dict[str, Any], p2689: dict[str, Any], p2761: dict[str, Any]) -> dict[str, Any]:
    coeffs = {k: float(v) for k, v in p1562.get("derived_lagrangian_coefficients", {}).items()}
    moments = {k: float(v) for k, v in p1562.get("moments", {}).items()}
    witness = finite_scale_orbit_witness(coeffs)
    rows = []
    for name, meta in COUPLING_ROWS.items():
        rows.append({
            "coupling": name,
            "numeric_candidate": coeffs.get(name),
            "expected_mass_dimension_4d": meta["expected_mass_dimension_4d"],
            "minimal_normalization_unknowns": meta["minimal_unknowns"],
            "reference_cell_fixed": False,
            "action_density_unit_fixed": False,
            "field_or_curvature_norm_fixed": False,
            "accepted_physical_normalization": False,
            "blocker": "A positive reference length/action-density normalization family remains; current artifacts do not export the canonical cell and normalization theorem needed to choose one member.",
        })
    p2689_decision = p2689.get("decision", {})
    p2761_acceptance = p2761.get("acceptance_matrix", {}).get("facts", {})
    missing_atoms = [
        "canonical selector-free physical length/reference cell",
        "action-density unit convention for converting dimensionless sums into integral density",
        "cell volume / measure normalization for strict kernel moments M0..M3",
        "field and curvature normalization compatible with the same cell",
        "variational insertion check after fixing the above normalizations",
    ]
    return {
        "moments": moments,
        "derived_lagrangian_coefficients": coeffs,
        "finite_scale_orbit_witness": witness,
        "rows": rows,
        "row_count": len(rows),
        "accepted_normalized_coupling_count": sum(1 for row in rows if row["accepted_physical_normalization"]),
        "missing_normalization_atoms": missing_atoms,
        "missing_normalization_atom_count": len(missing_atoms),
        "p2689_already_blocks_unconditional_uv_unit": p2689.get("status", "").endswith("NO_FALSE_PASS") and not p2689_decision.get("uv_unit_exported", False),
        "p2761_already_blocks_unit_reference_source": p2761_acceptance.get("unit_reference_source_exported") is False,
        "continuous_scale_gauge_unfixed": witness["distinct_dimensionalizations"] > 1,
    }


def acceptance_matrix(scan: dict[str, Any], matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "p1562_coefficients_present": bool(matrix["derived_lagrangian_coefficients"]),
        "finite_scale_orbit_exhibits_nonuniqueness": matrix["continuous_scale_gauge_unfixed"],
        "p2689_blocks_unconditional_uv_unit": matrix["p2689_already_blocks_unconditional_uv_unit"],
        "p2761_blocks_unit_reference_source": matrix["p2761_already_blocks_unit_reference_source"],
        "canonical_reference_cell_theorem_exported": False,
        "action_density_normalization_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_reference_cell_action_density_theorem": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The audit finds a real scale-orbit/nonuniqueness witness: changing the positive reference length changes dimensionful kappa/epsilon while leaving the numeric moment map intact. P2689 and P2761 already block an unconditional unit/reference source, and no new canonical reference-cell/action-density theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["normalization_obligation_matrix"]
    lines = [
        "# P2762/S1712 reference-cell and action-density normalization obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite scale-orbit witness",
        f"- tested_positive_reference_lengths={m['finite_scale_orbit_witness']['tested_positive_reference_lengths']}",
        f"- distinct_dimensionalizations={m['finite_scale_orbit_witness']['distinct_dimensionalizations']}",
        f"- lambda_invariant_under_ell={m['finite_scale_orbit_witness']['lambda_invariant_under_ell']}",
        f"- kappa_changes_with_ell={m['finite_scale_orbit_witness']['kappa_changes_with_ell']}",
        f"- epsilon_changes_with_ell={m['finite_scale_orbit_witness']['epsilon_changes_with_ell']}",
        "",
        "## Normalization rows",
    ]
    for row in m["rows"]:
        lines.append(f"- {row['coupling']}: accepted={row['accepted_physical_normalization']}; expected_mass_dimension_4d={row['expected_mass_dimension_4d']}; unknowns={row['minimal_normalization_unknowns']}")
    lines.extend(["", "## Missing normalization atoms"])
    for atom in m["missing_normalization_atoms"]:
        lines.append(f"- {atom}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p2689 = read_json(P2689)
    p2760 = read_json(P2760)
    p2761 = read_json(P2761)
    scan = evidence_scan()
    matrix = normalization_matrix(p1562, p2689, p2761)
    acceptance = acceptance_matrix(scan, matrix)
    payload = {
        "status": "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P1562": sha(P1562), "P2689": sha(P2689), "P2760": sha(P2760), "P2761": sha(P2761)},
        "input_statuses": {"P1562": p1562.get("status"), "P2689": p2689.get("status"), "P2760": p2760.get("status"), "P2761": p2761.get("status")},
        "audited_atom": "canonical physical length/reference-cell and action-density normalization theorem for the strict moment map",
        "content_evidence_scan": scan,
        "normalization_obligation_matrix": matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the strict moment-map coefficients to physical Lagrangian couplings through an arbitrary unit choice.  The next honest proof-grade move should target exactly one remaining provenance atom after P2762: a sign-convention theorem for the moment-to-coupling map, but only if it is explicitly conditional on a future canonical reference-cell/action-density theorem; otherwise preserve the P2697-P2762 no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2762/S1712 reference-cell action-density normalization obstruction", "## P2762/S1712 reference-cell action-density normalization obstruction\n\n`P2762/S1712` attacks the P2761 preferred provenance atom: canonical physical length/reference-cell plus action-density normalization for the strict moment map.  It gives a finite scale-orbit witness over positive reference lengths: the dimensionless `lambda_sm_eff` row is unchanged, while dimensionful `kappa_gr_eff` and `epsilon_mix_eff` vary with the chosen length scale.  Together with the P2689 UV-unit block and P2761 unit-reference block, this shows that current artifacts do not export a unique reference cell, action-density unit, or physical moment normalization.  No physical-coupling provenance theorem, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2762/S1712 reference-cell action-density Ltotal guard", "## P2762/S1712 reference-cell action-density Ltotal guard\n\n`P2762/S1712` adds no variational source term.  It records that the strict moment coefficients remain gauge-normalization candidates until a canonical reference cell and action-density theorem fixes the positive length/cell scale and compatible field/curvature normalizations.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current reference-cell action-density normalization obstruction guardrail (P2762/S1712, 2026-06-15)", "## Current reference-cell action-density normalization obstruction guardrail (P2762/S1712, 2026-06-15)\n\n- P2762 attacks the P2761 preferred atom: canonical physical length/reference-cell and action-density normalization for the strict moment map.\n- A finite positive-length scale-orbit witness shows that the numeric moment coefficients remain nonunique as physical dimensionful couplings: `lambda_sm_eff` is scale-invariant, but `kappa_gr_eff` and `epsilon_mix_eff` change with the reference length unless a canonical cell/action-density theorem is exported.\n- Do not promote arbitrary unit choices to physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  The next admissible move is one remaining provenance atom, preferably a sign-convention theorem explicitly conditional on future canonical reference-cell/action-density closure, or preservation of the P2697-P2762 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
