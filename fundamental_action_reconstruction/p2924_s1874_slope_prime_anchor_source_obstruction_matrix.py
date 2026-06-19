#!/usr/bin/env python3
"""P2924/S1874: slope/prime anchor source obstruction matrix.

P2923 exported formal prime-log character readiness but left two residual atoms:
strict prime-log value source and strict slope/prime anchor source.  P2924 attacks
exactly one of them: the slope/prime anchor.  It verifies that the audited
unital multiplicative log-character family y_d = delta * log(d) is homogeneous
in the free slope delta, so finite monoid/prime-log readiness does not select
`delta=4/5` or `eta=9/5` unless a new strict nadsoliton source law supplies the
anchor.
"""
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2923 = GEN / "p2923_s1873_prime_log_proportionality_source_matrix.json"
OUT = GEN / "p2924_s1874_slope_prime_anchor_source_obstruction_matrix.json"
MD = GEN / "p2924_s1874_slope_prime_anchor_source_obstruction_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NODES = list(range(1, 12))
SLOPE_CANDIDATES = [Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(3, 2)]
STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)


def frac_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def y_value(d: int, delta: Fraction) -> float:
    return float(delta) * math.log(d)


def slope_family_rows() -> list[dict[str, Any]]:
    rows = []
    for delta in SLOPE_CANDIDATES:
        defects = []
        ratios = []
        for d in NODES[1:]:
            ratios.append(y_value(d, delta) / math.log(d))
        for d in NODES:
            for e in NODES:
                if d * e <= NODES[-1]:
                    defects.append(y_value(d * e, delta) - y_value(d, delta) - y_value(e, delta))
        rows.append({
            "delta": frac_text(delta),
            "eta_if_eta_equals_1_plus_delta": frac_text(1 + delta),
            "is_strict_numeric_target": delta == STRICT_DELTA,
            "unital_y1": y_value(1, delta),
            "max_abs_multiplicative_defect": max(abs(x) for x in defects),
            "ratio_y_d_over_log_d_constant": max(ratios) - min(ratios) < 1e-12,
            "ratio_value": float(delta),
            "passes_prime_log_character_shape": max(abs(x) for x in defects) < 1e-12,
        })
    return rows


def anchor_obstruction_rows() -> list[dict[str, Any]]:
    return [
        {
            "row": "unital_multiplicative_character",
            "equation": "y_1=0 and y_de=y_d+y_e",
            "selects_delta_4_5": False,
            "free_parameter_left": "delta",
            "witness": "all audited delta candidates satisfy the shape up to floating tolerance",
        },
        {
            "row": "prime_log_proportionality_shape",
            "equation": "y_d/log(d)=delta for d=2..11 once log atoms are admitted",
            "selects_delta_4_5": False,
            "free_parameter_left": "delta",
            "witness": "delta=1/2, 4/5, 1, and 3/2 all remain proportional families",
        },
        {
            "row": "eta_translation",
            "equation": "eta=1+delta",
            "selects_delta_4_5": False,
            "free_parameter_left": "eta or delta scale",
            "witness": "the relation maps any delta to an eta and does not source delta=4/5",
        },
        {
            "row": "finite_node_support",
            "equation": "nodes 1..11 with prime basis {2,3,5,7,11}",
            "selects_delta_4_5": False,
            "free_parameter_left": "global slope multiplier",
            "witness": "support cardinality/factorization fixes basis shape, not the slope value",
        },
    ]


def candidate_anchor_sources() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "literal_delta_4_5_token",
            "computes_delta_4_5": True,
            "strict_nadsoliton_provenance": False,
            "nonconventional_anchor": False,
            "couples_to_prime_log_character": True,
            "accepted_as_slope_prime_anchor_source": False,
            "failure": "naming the target value is an imported numeric token, not a source theorem",
        },
        {
            "candidate": "eta_9_5_minus_one_relation",
            "computes_delta_4_5": True,
            "strict_nadsoliton_provenance": False,
            "nonconventional_anchor": False,
            "couples_to_prime_log_character": True,
            "accepted_as_slope_prime_anchor_source": False,
            "failure": "eta=9/5 is equivalent to delta=4/5 and cannot source itself",
        },
        {
            "candidate": "two_node_ratio_fit",
            "computes_delta_4_5": True,
            "strict_nadsoliton_provenance": False,
            "nonconventional_anchor": False,
            "couples_to_prime_log_character": True,
            "accepted_as_slope_prime_anchor_source": False,
            "failure": "a fitted two-node anchor can encode the value but remains a boundary condition unless strict provenance is added",
        },
        {
            "candidate": "P2923_factorization_readiness",
            "computes_delta_4_5": False,
            "strict_nadsoliton_provenance": True,
            "nonconventional_anchor": False,
            "couples_to_prime_log_character": True,
            "accepted_as_slope_prime_anchor_source": False,
            "failure": "factorization readiness supplies the carrier, not the slope multiplier",
        },
        {
            "candidate": "external_log_scale_calibration",
            "computes_delta_4_5": True,
            "strict_nadsoliton_provenance": False,
            "nonconventional_anchor": False,
            "couples_to_prime_log_character": True,
            "accepted_as_slope_prime_anchor_source": False,
            "failure": "external calibration is empirical/imported rather than a strict nadsoliton source law",
        },
        {
            "candidate": "Strict_Damping_Slope_Prime_Anchor_Source_Theorem",
            "computes_delta_4_5": False,
            "strict_nadsoliton_provenance": False,
            "nonconventional_anchor": False,
            "couples_to_prime_log_character": True,
            "accepted_as_slope_prime_anchor_source": False,
            "failure": "this is the missing theorem schema, not an exported theorem instance",
        },
    ]


def build_payload(p2923: dict[str, Any]) -> dict[str, Any]:
    family = slope_family_rows()
    candidates = candidate_anchor_sources()
    accepted = [c for c in candidates if c["accepted_as_slope_prime_anchor_source"]]
    strict_row = next(row for row in family if row["delta"] == "4/5")
    return {
        "status": "P2924_SLOPE_PRIME_ANCHOR_SOURCE_OBSTRUCTION_MATRIX_NO_ACCEPTED_ANCHOR",
        "input_hashes": {"P2923": hashlib.sha256(P2923.read_bytes()).hexdigest() if P2923.exists() else None},
        "constructed_theoretical_objects": {
            "missing_theorem_name": "Strict_Damping_Slope_Prime_Anchor_Source_Theorem",
            "target_statement": "strict_nadsoliton_data -> delta=4/5, eta=9/5, coupled to y_d=delta*sum_p v_p(d)*L_p",
            "slope_family_rows": family,
            "anchor_obstruction_rows": anchor_obstruction_rows(),
            "candidate_anchor_sources": candidates,
            "layman_closure_potential_note": "The theory has many working gears: finite carriers, multiplicative shape, and obstruction tests.  The missing piece here is not arithmetic but a principled reason why nature should choose the 4/5 damping slope rather than another allowed slope.",
        },
        "acceptance_matrix": {
            "p2923_formal_log_character_readiness_inherited": p2923.get("acceptance_matrix", {}).get("formal_log_character_readiness_exported") is True,
            "audited_slope_family_count": len(family),
            "slope_family_rows_passing_shape": sum(1 for row in family if row["passes_prime_log_character_shape"]),
            "strict_delta_row_passes_shape": strict_row["passes_prime_log_character_shape"],
            "strict_delta_selected_by_finite_shape": False,
            "anchor_obstruction_row_count": len(anchor_obstruction_rows()),
            "candidate_anchor_source_count": len(candidates),
            "candidate_anchor_sources_computing_target_value": sum(1 for c in candidates if c["computes_delta_4_5"]),
            "accepted_anchor_source_count": len(accepted),
            "strict_slope_prime_anchor_source_exported": False,
            "strict_damping_beta_eta_source_exported": False,
        },
        "decision": {
            "positive_witnesses": {
                "strict_delta_4_5_is_compatible_with_prime_log_character_shape": strict_row["passes_prime_log_character_shape"],
                "homogeneous_family_obstruction_computed": True,
                "anchor_acceptance_schema_exported": True,
            },
            "negative_export_flags": {
                "strict_slope_prime_anchor_source_exported": False,
                "strict_damping_beta_eta_source_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2924 audits the slope/prime anchor atom left by P2923.  The strict numeric target delta=4/5 is compatible with the finite unital multiplicative prime-log character, but the same finite equations also admit other slope values.  Thus the finite shape is homogeneous in delta and cannot select eta=9/5 without a new strict nadsoliton source theorem.",
            "next_honest_step": "Either construct an explicit strict source law deriving delta=4/5/eta=9/5 from nadsoliton data and coupling it to the P2923 prime-log character, or stop the damping lane with a post-P2923/P2924 no-new-live-frontier certificate.  Do not promote fitted/named slope tokens, factorization readiness, bridge prose, role transfer, L_total, or ToE closure.",
            "layman_toe_potential": "For a layperson: this is promising as a disciplined research program because it keeps finding precise gears and precise missing bolts.  But it is not yet a Theory of Everything closure.  A ToE-level claim would need the missing source laws to be derived internally rather than chosen, fitted, or imported.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2924/S1874 slope/prime anchor source obstruction matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite anchor gate",
        f"- inherited P2923 log-character readiness: `{acc['p2923_formal_log_character_readiness_inherited']}`",
        f"- audited slope-family rows: `{acc['audited_slope_family_count']}`",
        f"- rows passing prime-log character shape: `{acc['slope_family_rows_passing_shape']}`",
        f"- strict delta=4/5 row passes shape: `{acc['strict_delta_row_passes_shape']}`",
        f"- strict delta selected by finite shape: `{acc['strict_delta_selected_by_finite_shape']}`",
        f"- candidate anchor sources: `{acc['candidate_anchor_source_count']}`",
        f"- accepted anchor sources: `{acc['accepted_anchor_source_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
        "",
        "## Layman ToE potential note",
        payload["decision"]["layman_toe_potential"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2923))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2924/S1874 slope/prime anchor source obstruction matrix", "## P2924/S1874 slope/prime anchor source obstruction matrix\n\n`P2924/S1874` attacks exactly one residual atom after P2923: the strict slope/prime anchor.  The audited family `y_d = delta*log(d)` satisfies the unital multiplicative prime-log character shape for every tested `delta`, including `delta=4/5` and non-target slopes.  Therefore the finite shape is homogeneous in the free slope and does not select `delta=4/5` or `eta=9/5`.  Six candidate anchors are rejected as target naming, equivalent eta naming, boundary fitting, factorization readiness, external calibration, or a theorem schema without an instance.  No strict damping `beta/eta` source, nonproxy `L_total`, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2924/S1874 slope/prime anchor `L_total` guard", "## P2924/S1874 slope/prime anchor `L_total` guard\n\n`P2924/S1874` shows that the strict numeric damping target `delta=4/5` is compatible with the P2923 prime-log character but not selected by it.  Because other slopes satisfy the same finite character equations, the damping term remains non-role-bearing until a strict nadsoliton source law derives the slope/prime anchor and couples it to the prime-log carrier.  No EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(AGENTS, "Current slope/prime anchor obstruction guardrail (P2924/S1874, 2026-06-19)", "## Current slope/prime anchor obstruction guardrail (P2924/S1874, 2026-06-19)\n\n- P2924 attacks exactly one post-P2923 residual atom: the strict slope/prime anchor for `delta=4/5` / `eta=9/5`.\n- The finite unital multiplicative prime-log character shape admits the target slope but also admits other slope values; it is homogeneous in the free `delta` parameter and does not source the target value.\n- Rejected anchors include named `4/5`, named `eta=9/5`, two-node fitting, P2923 factorization readiness, external calibration, and the missing theorem name itself.\n- Do not promote P2923/P2924 damping readiness to strict damping `beta/eta`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n- A next admissible move must either provide an explicit strict source law deriving `delta=4/5`/`eta=9/5` and coupling it to the prime-log character, or emit a post-P2923/P2924 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
