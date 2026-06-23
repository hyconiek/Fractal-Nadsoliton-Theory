#!/usr/bin/env python3
"""P3041/S1991: integrated-density unit/readout coupling obstruction.

P3038 produced a finite branch-separating integrated density, while P3039 and
P3040 audited two source premises of that candidate.  P3041 attacks the
remaining direct P3038 premise: can current strict artifacts source a physical
unit/readout coupling theorem for the integrated density rather than using the
placeholder U_readout = 1?

The finite result is a bounded obstruction.  Unit scaling acts exactly on the
branch score, and several normalization/readout receivers are computable, but
all are internal gauges, P3036-blocked unit placeholders, or target-readout
normalizations.  No external physical unit, scale-orbit quotient, or nonproxy
readout/action coupling theorem is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3036_s1986_external_physical_unit_source_scale_orbit_obstruction import OUT as P3036
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import OUT as P3038, branch_score, integrated_density
from p3039_s1989_chiral_projector_sign_localizer_obstruction import OUT as P3039
from p3040_s1990_retarded_path_anisotropy_source_obstruction import OUT as P3040

OUT = GEN / "p3041_s1991_integrated_density_unit_readout_coupling_obstruction.json"
MD = GEN / "p3041_s1991_integrated_density_unit_readout_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-12
UNIT_GRID = [0.25, 0.5, 1.0, 2.0, 4.0]


def density_norms(branch: int = +1) -> dict[str, float]:
    dens = integrated_density(branch, 1.0)
    return {
        "l1": sum(abs(x) for x in dens),
        "l2": math.sqrt(sum(x * x for x in dens)),
        "linf": max(abs(x) for x in dens),
        "sum": sum(dens),
    }


def unit_scaled_scores(branch: int = +1) -> list[dict[str, Any]]:
    base = branch_score(branch)
    return [
        {
            "unit_scale": u,
            "score": round(u * base, 15),
            "score_over_unit": round((u * base) / u, 15),
            "finite_nonzero": abs(u * base) > TOL,
        }
        for u in UNIT_GRID
    ]


def build_matrix() -> dict[str, Any]:
    plus_base = branch_score(+1)
    minus_base = branch_score(-1)
    norms = density_norms(+1)
    scaled = unit_scaled_scores(+1)
    l1_unit = 1.0 / norms["l1"] if norms["l1"] > TOL else 0.0
    l2_unit = 1.0 / norms["l2"] if norms["l2"] > TOL else 0.0
    score_unit = 1.0 / abs(plus_base) if abs(plus_base) > TOL else 0.0
    receivers = [
        {
            "receiver": "placeholder_U_readout_equals_one",
            "formula": "S_b = sum_i b * U_readout * Delta_i * chi_i with U_readout=1",
            "finite_nonzero": abs(plus_base) > TOL and abs(minus_base) > TOL,
            "accepted_as_physical_unit_readout_coupling": False,
            "failure": "U=1 is a dimensionless placeholder already blocked by P3036, not an external physical unit",
        },
        {
            "receiver": "unit_scale_orbit_on_branch_score",
            "unit_scaled_scores": scaled,
            "finite_nonzero": all(row["finite_nonzero"] for row in scaled),
            "accepted_as_physical_unit_readout_coupling": False,
            "failure": "scores rescale exactly with U; the orbit is represented but not physically quotiented or externally fixed",
        },
        {
            "receiver": "density_norm_to_one_couplings",
            "candidate_units": {"l1_to_one": round(l1_unit, 15), "l2_to_one": round(l2_unit, 15), "score_to_one": round(score_unit, 15)},
            "finite_nonzero": l1_unit > TOL and l2_unit > TOL and score_unit > TOL,
            "accepted_as_physical_unit_readout_coupling": False,
            "failure": "norm-to-one and score-to-one couplings are internal normalizations/target gauges, not source theorems",
        },
        {
            "receiver": "p3036_external_unit_import_receiver",
            "candidate_imports": ["label_step", "kernel_l2_norm", "dirichlet_action_to_one", "entropy_one_bit_reference_cell"],
            "finite_nonzero": True,
            "accepted_as_physical_unit_readout_coupling": False,
            "failure": "P3036 classifies these as internal representatives or conditional gauges; importing them cannot source P3038 coupling",
        },
    ]
    obligations = [
        {"obligation": "exact_p3038_unit_readout_premise_targeted", "satisfied": True, "detail": "only the physical unit/readout coupling premise is audited"},
        {"obligation": "finite_integrated_density_readout_nonzero", "satisfied": abs(plus_base) > TOL and abs(minus_base) > TOL, "detail": "P3038 branch score remains nonzero"},
        {"obligation": "unit_scale_orbit_constructed", "satisfied": True, "detail": "finite U grid shows exact score scaling"},
        {"obligation": "p3036_unit_boundary_consulted", "satisfied": True, "detail": "P3036 external unit-source obstruction is an explicit input"},
        {"obligation": "external_physical_unit_source", "satisfied": False, "detail": "no length/action/clock unit is exported beyond P3036-blocked representatives"},
        {"obligation": "scale_orbit_quotient_or_absolute_calibration", "satisfied": False, "detail": "U rescales the score; normalization fixes representatives only"},
        {"obligation": "readout_coupling_theorem", "satisfied": False, "detail": "no theorem maps the integrated density to a named observer-independent physical readout row"},
        {"obligation": "nonproxy_variational_or_action_installation", "satisfied": False, "detail": "no unit-bearing density/action/EOM/Hamiltonian insertion is exported"},
    ]
    return {
        "object": "IntegratedDensityUnitReadoutCoupling_ObstructionMatrix",
        "target_formula": "S_b = sum_i b * U_readout * Delta_i * chi_i from P3038",
        "density_norms": {k: round(v, 15) for k, v in norms.items()},
        "unit_grid": UNIT_GRID,
        "branch_scores_at_U1": {"plus": round(plus_base, 15), "minus": round(minus_base, 15)},
        "receiver_rows": receivers,
        "proof_obligations": obligations,
        "finite_certificate": {
            "receiver_rows": len(receivers),
            "finite_nonzero_rows": sum(1 for row in receivers if row["finite_nonzero"]),
            "accepted_physical_unit_readout_coupling_rows": sum(1 for row in receivers if row["accepted_as_physical_unit_readout_coupling"]),
            "unit_grid_rows": len(scaled),
            "unit_scaling_exact_rows": sum(1 for row in scaled if abs(row["score_over_unit"] - round(plus_base, 15)) <= 1e-12),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "physical_unit_readout_coupling_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_physical_unit_readout_coupling"] for row in receivers),
        },
    }


def build_payload() -> dict[str, Any]:
    for path in (P3036, P3038, P3039, P3040):
        read_json(path)
    matrix = build_matrix()
    return {
        "status": "P3041_INTEGRATED_DENSITY_UNIT_READOUT_COUPLING_OBSTRUCTION_NO_SOURCE_EXPORT",
        "input_hashes": {
            "P3036": hashlib.sha256(P3036.read_bytes()).hexdigest() if P3036.exists() else None,
            "P3038": hashlib.sha256(P3038.read_bytes()).hexdigest() if P3038.exists() else None,
            "P3039": hashlib.sha256(P3039.read_bytes()).hexdigest() if P3039.exists() else None,
            "P3040": hashlib.sha256(P3040.read_bytes()).hexdigest() if P3040.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The integrated density has a nonzero finite readout and a fully explicit unit-scale orbit: multiplying U rescales the score exactly.  However U=1, norm-to-one, score-to-one, and imported P3036 unit candidates are all internal gauges or blocked representatives.  No external physical unit, absolute scale quotient, readout-coupling theorem, or nonproxy variational/action installation is exported.",
            "negative_export_flags": {k: False for k in ["physical_unit_readout_coupling_exported", "p3038_promoted_to_selector_source", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay U=1, norm-to-one, score-to-one, or P3036-blocked unit imports as physical readout coupling.  Since P3039-P3041 have now attacked the direct P3038 source premises, the next honest move is a P3038-P3041 integrated selector-candidate reconciliation/no-selector-export certificate unless a genuinely new strict source law is supplied outside these exhausted receiver classes.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3041/S1991 integrated-density unit/readout coupling obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- receiver rows: `{c['receiver_rows']}`",
        f"- finite nonzero rows: `{c['finite_nonzero_rows']}`",
        f"- accepted physical unit/readout coupling rows: `{c['accepted_physical_unit_readout_coupling_rows']}`",
        f"- unit grid rows: `{c['unit_grid_rows']}`",
        f"- unit scaling exact rows: `{c['unit_scaling_exact_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- physical unit/readout coupling exported: `{c['physical_unit_readout_coupling_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3041/S1991 integrated-density unit/readout coupling obstruction", "## P3041/S1991 integrated-density unit/readout coupling obstruction\n\n`P3041/S1991` attacks the remaining direct P3038 premise: a physical unit/readout coupling theorem for `S_b = sum_i b * U_readout * Delta_i * chi_i`.  The integrated density has a nonzero finite readout and the unit-scale orbit is exact: multiplying `U` rescales the score.  But `U=1`, norm-to-one/score-to-one couplings, and imported P3036 unit candidates are internal gauges or blocked representatives.  No external physical unit, absolute scale quotient, readout-coupling theorem, nonproxy action/EOM installation, selector closure, `QW-2191` discharge, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3041/S1991 integrated density unit/readout `L_total` guard", "## P3041/S1991 integrated density unit/readout `L_total` guard\n\n`P3041/S1991` adds no physical `L_total` term.  It verifies exact unit-scale action on the P3038 readout score, but all tested couplings are dimensionless placeholders, internal normalizations, target gauges, or P3036-blocked unit imports rather than unit-bearing variational/action/EOM insertions.\n")
    append_once(AGENTS, "Current integrated-density unit/readout coupling guardrail (P3041/S1991, 2026-06-23)", "## Current integrated-density unit/readout coupling guardrail (P3041/S1991, 2026-06-23)\n\n- P3041 attacks the remaining direct P3038 source premise: a physical unit/readout coupling theorem for the integrated density.\n- The finite unit-scale orbit is exact, but `U=1`, norm-to-one/score-to-one couplings, and P3036 unit imports are internal gauges or blocked representatives rather than external physical unit/readout source laws.\n- Do not promote unit-scale orbits, placeholder `U_readout`, norm/score normalizations, or P3036-blocked unit imports to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move should reconcile P3038-P3041 as an integrated selector-candidate no-export certificate unless a genuinely new strict source law outside these receiver classes is supplied.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
