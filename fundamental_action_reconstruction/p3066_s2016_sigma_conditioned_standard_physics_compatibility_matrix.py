#!/usr/bin/env python3
"""P3066/S2016: sigma-conditioned standard-physics compatibility matrix.

P3065 allowed `sigma_selector` as an explicit boundary axiom so selector search
no longer blocks informational-nadsoliton work.  P3066 takes the next honest
step: do not ask again whether sigma is derived; instead build a finite matrix
for what can be studied under T_sigma and what would be needed for compatibility
with standard theoretical physics.

The constructed object is a SigmaConditionedStandardPhysicsCompatibilityMatrix.
It separates sigma-invariant rows from sigma-odd rows and checks six physics
obligations without claiming Lorentz/gauge/unitarity/renormalization/EOM/
empirical closure.
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3065_s2015_fundamental_selector_axiom_boundary_certificate import OUT as P3065

OUT = GEN / "p3066_s2016_sigma_conditioned_standard_physics_compatibility_matrix.json"
MD = GEN / "p3066_s2016_sigma_conditioned_standard_physics_compatibility_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "sigma_conditioned_standard_physics_matrix": r"sigma[-_ ]conditioned|conditioned-on-`?sigma|T_sigma.*standard|standard.*T_sigma|orientation-conditioned.*physics",
    "observable_invariant_under_selector_axiom": r"observable[-_ ]invariant|sigma[-_ ]invariant|orientation.*invariant|state-transition.*observable",
    "standard_physics_compatibility_obligations": r"Lorentz|gauge invariance|unitarity|renormalization|variational.*EOM|standard theoretical physics|observed physics",
}

SIGMA_BRANCHES = ("sigma_plus", "sigma_minus")
OBSERVABLE_ROWS = (
    {"id": "kernel_shape_scalar", "description": "K_strict_gate scalar-shape bookkeeping under fixed sigma", "sigma_parity": "even", "exported_currently": True},
    {"id": "entropy_count_scalar", "description": "selector-independent entropy/counting scalar", "sigma_parity": "even", "exported_currently": True},
    {"id": "orientation_boundary_datum", "description": "recorded sigma orientation boundary parameter", "sigma_parity": "odd", "exported_currently": True},
    {"id": "chiral_bispectrum_sign", "description": "signed/chiral marker candidate inherited only as orientation-sensitive evidence", "sigma_parity": "odd", "exported_currently": False},
    {"id": "light_emergence_interface", "description": "nadsoliton -> light transition interface", "sigma_parity": "unknown", "exported_currently": False},
    {"id": "matter_coupling_interface", "description": "light -> matter coupling interface", "sigma_parity": "unknown", "exported_currently": False},
)
PHYSICS_OBLIGATIONS = (
    "Lorentz_covariance",
    "gauge_invariance",
    "unitarity_or_positive_probability",
    "renormalization_or_scale_control",
    "variational_EOM_or_dynamics",
    "empirical_constant_map",
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def evaluate_observable(row: dict[str, Any]) -> dict[str, Any]:
    parity = row["sigma_parity"]
    invariant = parity == "even"
    axiom_safe = parity in {"even", "odd"}
    strict_physics_ready = False
    blockers = []
    if parity == "odd":
        blockers.append("requires_explicit_pseudoscalar_or_orientation_coupling_for_standard_physics_use")
    if parity == "unknown":
        blockers.append("requires_exported_transition_law_before_physics_obligation_check")
    if not row["exported_currently"]:
        blockers.append("object_not_exported_currently")
    blockers.extend(PHYSICS_OBLIGATIONS)
    return {
        **row,
        "sigma_invariant": invariant,
        "axiom_conditioned_bookkeeping_safe": axiom_safe,
        "strict_standard_physics_ready": strict_physics_ready,
        "remaining_blockers": blockers,
    }


def evaluate_branch_matrix() -> list[dict[str, Any]]:
    obs = [evaluate_observable(row) for row in OBSERVABLE_ROWS]
    matrix = []
    for sigma in SIGMA_BRANCHES:
        sign = 1 if sigma == "sigma_plus" else -1
        for row in obs:
            sigma_value = 1 if row["sigma_parity"] == "even" else sign if row["sigma_parity"] == "odd" else None
            matrix.append({
                "sigma_branch": sigma,
                "observable_id": row["id"],
                "sigma_value_proxy": sigma_value,
                "sigma_invariant": row["sigma_invariant"],
                "axiom_conditioned_bookkeeping_safe": row["axiom_conditioned_bookkeeping_safe"],
                "strict_standard_physics_ready": row["strict_standard_physics_ready"],
                "remaining_blockers": row["remaining_blockers"],
            })
    return matrix


def evaluate_physics_obligation_matrix() -> list[dict[str, Any]]:
    rows = []
    for obs in [evaluate_observable(row) for row in OBSERVABLE_ROWS]:
        for obligation in PHYSICS_OBLIGATIONS:
            rows.append({
                "observable_id": obs["id"],
                "obligation": obligation,
                "sigma_invariant_input": obs["sigma_invariant"],
                "accepted_now": False,
                "reason": "P3066 is an acceptance matrix, not a closure proof; the required standard-physics theorem is not exported for this observable",
            })
    return rows


def build_payload() -> dict[str, Any]:
    p3065 = read_json(P3065)
    grep_rows = content_grep()
    observable_rows = [evaluate_observable(row) for row in OBSERVABLE_ROWS]
    branch_matrix = evaluate_branch_matrix()
    physics_matrix = evaluate_physics_obligation_matrix()
    proof_obligations = [
        {"obligation": "content_first_grep_before_sigma_conditioned_physics_matrix", "satisfied": True, "detail": "searched by sigma-conditioned, observable-invariant, and standard-physics-obligation content"},
        {"obligation": "construct_sigma_conditioned_observable_matrix", "satisfied": True, "detail": "built 2 sigma branches times 6 observables = 12 branch rows"},
        {"obligation": "separate_sigma_invariant_from_sigma_odd_or_unknown_rows", "satisfied": True, "detail": "even rows are invariant, odd rows require pseudoscalar/orientation coupling, unknown rows require transition laws"},
        {"obligation": "audit_standard_physics_compatibility_obligations", "satisfied": True, "detail": "checked 6 obligations for each observable without promoting closure"},
        {"obligation": "export_standard_physics_closure", "satisfied": False, "detail": "no Lorentz/gauge/unitarity/renormalization/EOM/empirical closure theorem is exported"},
    ]
    return {
        "status": "P3066_SIGMA_CONDITIONED_STANDARD_PHYSICS_COMPATIBILITY_MATRIX_NO_CLOSURE",
        "input_hashes": {"P3065": hashlib.sha256(P3065.read_bytes()).hexdigest() if P3065.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "sigma_conditioned_standard_physics_compatibility_matrix": {
                "object": "SigmaConditionedStandardPhysicsCompatibilityMatrix",
                "definition": "Under T_sigma, classify each informational-nadsoliton observable by sigma parity and by six standard-physics compatibility obligations.",
                "sigma_branches": list(SIGMA_BRANCHES),
                "physics_obligations": list(PHYSICS_OBLIGATIONS),
                "observable_rows": observable_rows,
                "branch_matrix": branch_matrix,
                "physics_obligation_matrix": physics_matrix,
            },
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "sigma_branches": len(SIGMA_BRANCHES),
            "observable_rows": len(OBSERVABLE_ROWS),
            "branch_matrix_rows": len(branch_matrix),
            "sigma_invariant_observables": sum(1 for r in observable_rows if r["sigma_invariant"]),
            "sigma_odd_observables": sum(1 for r in observable_rows if r["sigma_parity"] == "odd"),
            "unknown_parity_observables": sum(1 for r in observable_rows if r["sigma_parity"] == "unknown"),
            "physics_obligations": len(PHYSICS_OBLIGATIONS),
            "physics_obligation_rows": len(physics_matrix),
            "accepted_physics_obligation_rows": sum(1 for r in physics_matrix if r["accepted_now"]),
            "p3065_axiom_augmented_admitted_rows": p3065.get("finite_certificate", {}).get("axiom_augmented_admitted_rows"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3066 advances beyond selector search by accepting T_sigma as a boundary condition and constructing a finite sigma-conditioned compatibility matrix for informational-nadsoliton observables.  The matrix separates 2 sigma-invariant scalar rows, 2 sigma-odd rows, and 2 unknown transition-interface rows, then checks 36 standard-physics obligation rows.  Current artifacts export 0 standard-physics closures, so this is a proof-grade road map rather than a claim of compatibility with observed physics.",
            "negative_export_flags": {k: False for k in ["standard_physics_compatibility_exported", "Lorentz_closure_exported", "gauge_closure_exported", "unitarity_closure_exported", "renormalization_closure_exported", "variational_EOM_closure_exported", "empirical_constant_map_exported", "observed_physics_exported", "ltotal_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"sigma_conditioned_physics_matrix_constructed": True, "sigma_invariant_rows_identified": True, "sigma_odd_rows_blocked_pending_coupling": True, "selector_search_no_longer_blocks_conditioned_analysis": True},
            "next_honest_step": "Choose exactly one row from the P3066 matrix and prove or obstruct it.  The best next move is the light_emergence_interface row: construct a sigma-conditioned nadsoliton-to-light transition law and test Lorentz covariance on the finite branch matrix.  If that transition law cannot be exported, pivot to a sigma-invariant scalar row and build a finite conservation/scale-control theorem before discussing empirical compatibility.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3066/S2016 sigma-conditioned standard-physics compatibility matrix", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- sigma branches: `{c['sigma_branches']}`",
        f"- observable rows: `{c['observable_rows']}`",
        f"- branch matrix rows: `{c['branch_matrix_rows']}`",
        f"- sigma-invariant observables: `{c['sigma_invariant_observables']}`",
        f"- sigma-odd observables: `{c['sigma_odd_observables']}`",
        f"- unknown-parity observables: `{c['unknown_parity_observables']}`",
        f"- physics obligation rows: `{c['physics_obligation_rows']}`",
        f"- accepted physics obligation rows: `{c['accepted_physics_obligation_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3066/S2016 sigma-conditioned standard-physics compatibility matrix", "## P3066/S2016 sigma-conditioned standard-physics compatibility matrix\n\n`P3066/S2016` stops replaying selector derivation and works under the explicit `T_sigma` axiom boundary from `P3065`.  It constructs a `SigmaConditionedStandardPhysicsCompatibilityMatrix` with `2` sigma branches, `6` informational-nadsoliton observable rows, `12` branch rows, and `36` standard-physics obligation rows covering Lorentz covariance, gauge invariance, unitarity/positive probability, renormalization/scale control, variational EOM/dynamics, and empirical constant mapping.  The current matrix identifies `2` sigma-invariant scalar rows, `2` sigma-odd rows needing explicit orientation/pseudoscalar coupling, and `2` unknown transition-interface rows needing exported transition laws; it exports `0` standard-physics closure rows and no observed-physics, `L_total`, bridge/role-transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3066/S2016 standard-physics matrix is not an `L_total` closure", "## P3066/S2016 standard-physics matrix is not an `L_total` closure\n\n`P3066/S2016` is a sigma-conditioned acceptance/obstruction matrix for compatibility obligations.  It does not add a unit-bearing action, gauge sector, propagator, counterterm flow, EOM, Hamiltonian, or empirical constant map; every standard-physics row remains a future theorem obligation.\n")
    append_once(AGENTS, "Current sigma-conditioned standard-physics compatibility guardrail (P3066/S2016, 2026-06-24)", "## Current sigma-conditioned standard-physics compatibility guardrail (P3066/S2016, 2026-06-24)\n\n- P3066 proceeds under the P3065 `T_sigma` boundary instead of replaying strict selector derivation: selector search no longer blocks conditioned informational-nadsoliton analysis.\n- The constructed matrix has `2` sigma branches, `6` observable rows, `12` branch rows, and `36` standard-physics obligation rows; it separates sigma-invariant scalar rows from sigma-odd rows and unknown transition-interface rows.\n- Current artifacts export `0` Lorentz/gauge/unitarity/renormalization/EOM/empirical compatibility closures.  Do not promote this matrix to observed physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next proof-grade move should pick exactly one matrix row, preferably `light_emergence_interface`, construct a sigma-conditioned nadsoliton-to-light transition law, and test Lorentz covariance; otherwise pivot to a sigma-invariant scalar conservation/scale-control theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
