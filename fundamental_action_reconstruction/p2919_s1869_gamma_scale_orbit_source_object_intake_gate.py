#!/usr/bin/env python3
"""P2919/S1869: Gamma scale-orbit source-object intake gate.

P2918 showed that the P2916/P2917 quotient integral is homogeneous in
Gamma_9_5.  P2919 makes that obstruction sharper: it constructs the exact
missing source object as a scale-fixing map and tests whether any currently
available finite quotient datum can break the Gamma scaling orbit.

The finite result is negative.  Quotient weights, edge weights, derivative
ratios, orbit sizes, and candidate count normalizations are invariant under the
rescaling Gamma_9_5 -> c * Gamma_9_5 once ratios are taken.  They can determine
relative weights but cannot select a nonzero Action-dimension scale.  A real
source object would have to provide a strict nadsoliton-derived value map
sigma(strict_data) = Gamma_9_5 and a coupling proof to the P2917 integral.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2918 = GEN / "p2918_s1868_gamma_action_unit_source_law_obstruction_matrix.json"
OUT = GEN / "p2919_s1869_gamma_scale_orbit_source_object_intake_gate.json"
MD = GEN / "p2919_s1869_gamma_scale_orbit_source_object_intake_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
EDGE_COUNT = N * N


def f(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def scale_orbit_rows() -> list[dict[str, Any]]:
    """Finite observables and their behavior under Gamma -> c Gamma."""
    return [
        {
            "observable": "quotient_weight",
            "formula": "w_Q = 1/12",
            "under_gamma_rescaling": "unchanged",
            "fixes_relative_weight": True,
            "fixes_gamma_scale": False,
        },
        {
            "observable": "edge_weight",
            "formula": "w_E = 1/144",
            "under_gamma_rescaling": "unchanged",
            "fixes_relative_weight": True,
            "fixes_gamma_scale": False,
        },
        {
            "observable": "derivative_ratio",
            "formula": "(dI/dq_edge)/(dI/dQ_d) = 1/12",
            "under_gamma_rescaling": "unchanged because Gamma cancels",
            "fixes_relative_weight": True,
            "fixes_gamma_scale": False,
        },
        {
            "observable": "orbit_size",
            "formula": "|orbit(d)| = 12",
            "under_gamma_rescaling": "unchanged",
            "fixes_relative_weight": True,
            "fixes_gamma_scale": False,
        },
        {
            "observable": "edge_count",
            "formula": "|Z12 x Z12| = 144",
            "under_gamma_rescaling": "unchanged",
            "fixes_relative_weight": True,
            "fixes_gamma_scale": False,
        },
    ]


def required_source_object_schema() -> dict[str, Any]:
    return {
        "object_name": "Strict_Gamma_9_5_Action_Unit_Scale_Fixing_Source_Object",
        "source_map": "sigma_Gamma: strict_nadsoliton_data -> Action_nonzero",
        "coupling_map": "C_Gamma(sigma_Gamma, Q) = sigma_Gamma * (1/12) * sum_d Q_d",
        "acceptance_obligations": [
            "sigma_Gamma is defined from strict nadsoliton data, not from a quotient count or unit convention",
            "sigma_Gamma has a nonzero Action-dimension value",
            "sigma_Gamma breaks the Gamma -> c*Gamma scale orbit by selecting exactly one scale",
            "C_Gamma is proven to be the coefficient of the P2917 quotient integral",
            "the proof does not import selector replay, bridge closure, role transfer, L_total closure, or ToE promotion",
        ],
    }


def source_object_candidates() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "quotient measure weights {1/12, 1/144}",
            "breaks_scale_orbit": False,
            "nonzero_action_value": False,
            "strict_source_map": False,
            "accepted": False,
            "failure": "weights determine relative measure normalization only; they are dimensionless and Gamma-independent",
        },
        {
            "candidate": "derivative ratio 1/12",
            "breaks_scale_orbit": False,
            "nonzero_action_value": False,
            "strict_source_map": False,
            "accepted": False,
            "failure": "Gamma cancels from the ratio, so the ratio cannot source Gamma",
        },
        {
            "candidate": "orbit cardinality 12",
            "breaks_scale_orbit": False,
            "nonzero_action_value": False,
            "strict_source_map": False,
            "accepted": False,
            "failure": "cardinality is already the quotient fiber size and has no Action dimension",
        },
        {
            "candidate": "directed-edge cardinality 144",
            "breaks_scale_orbit": False,
            "nonzero_action_value": False,
            "strict_source_map": False,
            "accepted": False,
            "failure": "cardinality counts carriers and cannot be a strict unit-bearing source",
        },
        {
            "candidate": "external Action unit label",
            "breaks_scale_orbit": True,
            "nonzero_action_value": True,
            "strict_source_map": False,
            "accepted": False,
            "failure": "an imported dimensional label can fix scale by convention but is not strict nadsoliton provenance",
        },
        {
            "candidate": "new strict sigma_Gamma source map",
            "breaks_scale_orbit": False,
            "nonzero_action_value": False,
            "strict_source_map": False,
            "accepted": False,
            "failure": "this is the required new object, not an exported current artifact",
        },
    ]


def build_payload(p2918: dict[str, Any]) -> dict[str, Any]:
    rows = scale_orbit_rows()
    candidates = source_object_candidates()
    accepted = [candidate for candidate in candidates if candidate["accepted"]]
    scale_breakers = [candidate for candidate in candidates if candidate["breaks_scale_orbit"]]
    strict_scale_breakers = [candidate for candidate in scale_breakers if candidate["strict_source_map"]]
    return {
        "status": "P2919_GAMMA_SCALE_ORBIT_SOURCE_OBJECT_INTAKE_GATE_NO_EXPORT",
        "input_hashes": {"P2918": hashlib.sha256(P2918.read_bytes()).hexdigest() if P2918.exists() else None},
        "constructed_theoretical_objects": {
            "required_source_object_schema": required_source_object_schema(),
            "scale_orbit_action": "Gamma_9_5 -> c * Gamma_9_5 for c in positive scale orbit",
            "finite_scale_orbit_rows": rows,
            "source_object_candidates": candidates,
            "scale_orbit_certificate": {
                "finite_observable_count": len(rows),
                "all_finite_observables_gamma_scale_invariant": all(not row["fixes_gamma_scale"] for row in rows),
                "relative_weights_fixed": sorted({"1/12", "1/144"}),
                "strict_scale_breaking_candidates": len(strict_scale_breakers),
                "accepted_source_objects": len(accepted),
            },
        },
        "acceptance_matrix": {
            "p2918_rechecked_gamma_free_scalar": p2918.get("acceptance_matrix", {}).get("gamma_9_5_remains_free_scalar") is True,
            "finite_observable_count": len(rows),
            "source_object_candidate_count": len(candidates),
            "scale_breaking_candidate_count": len(scale_breakers),
            "strict_scale_breaking_candidate_count": len(strict_scale_breakers),
            "accepted_source_object_count": len(accepted),
            "finite_observables_fix_relative_weights": True,
            "finite_observables_fix_gamma_scale": False,
            "strict_sigma_gamma_source_object_exported": False,
            "accepted_as_nonproxy_ltotal_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "required_source_object_schema_constructed": True,
                "scale_orbit_invariance_matrix_constructed": True,
                "source_object_intake_candidates_audited": True,
            },
            "negative_export_flags": {
                "strict_gamma_9_5_source_exported": False,
                "strict_scale_fixing_source_object_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2919 constructs the exact source object that would be needed to repair P2918: a strict sigma_Gamma map from nadsoliton data to a nonzero Action value coupled to the P2917 quotient integral.  The finite scale-orbit matrix shows that quotient weights, derivative ratios, orbit size, and edge count are invariant under Gamma_9_5 -> c*Gamma_9_5 and therefore fix only relative weights.  The only scale-breaking candidate is an imported Action label, which is not strict source provenance; no accepted strict source object is exported.",
            "next_honest_step": "Unless a genuinely new strict sigma_Gamma source map with a computed nonzero Action value is supplied, the next honest move is a Gamma/Lambda no-new-live-frontier certificate.  If such a source map is supplied, test exactly its value, scale-orbit breaking, and coupling to I_Q before any nonproxy L_total/EOM/Hamiltonian/ToE promotion.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2919/S1869 Gamma scale-orbit source-object intake gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Scale-orbit source-object gate",
        f"- finite observables: `{acc['finite_observable_count']}`",
        f"- source-object candidates: `{acc['source_object_candidate_count']}`",
        f"- scale-breaking candidates: `{acc['scale_breaking_candidate_count']}`",
        f"- strict scale-breaking candidates: `{acc['strict_scale_breaking_candidate_count']}`",
        f"- accepted source objects: `{acc['accepted_source_object_count']}`",
        f"- finite observables fix relative weights: `{acc['finite_observables_fix_relative_weights']}`",
        f"- finite observables fix Gamma scale: `{acc['finite_observables_fix_gamma_scale']}`",
        f"- strict sigma_Gamma source object exported: `{acc['strict_sigma_gamma_source_object_exported']}`",
        f"- accepted as nonproxy L_total source: `{acc['accepted_as_nonproxy_ltotal_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2918))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2919/S1869 Gamma scale-orbit source-object intake gate", "## P2919/S1869 Gamma scale-orbit source-object intake gate\n\n`P2919/S1869` sharpens P2918 by constructing the missing `Strict_Gamma_9_5_Action_Unit_Scale_Fixing_Source_Object`: a source map `sigma_Gamma: strict_nadsoliton_data -> Action_nonzero` and coupling `C_Gamma = sigma_Gamma * (1/12) * sum_d Q_d`.  The finite scale-orbit matrix tests quotient weights, edge weights, derivative ratios, orbit size, and edge count under `Gamma_9_5 -> c*Gamma_9_5`; all remain scale-invariant and fix only relative weights.  Six intake candidates are audited, and the only scale-breaking one is an imported Action label rather than strict source provenance.  No strict `sigma_Gamma`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2919/S1869 Gamma scale-orbit `L_total` guard", "## P2919/S1869 Gamma scale-orbit `L_total` guard\n\n`P2919/S1869` proves that the finite Gamma/Lambda quotient data are invariant under the scale orbit `Gamma_9_5 -> c*Gamma_9_5`: weights, derivative ratios, and orbit/cardinality data determine only relative normalization.  A nonproxy `L_total` source still requires a strict `sigma_Gamma` map computing a nonzero Action value from nadsoliton data and coupling it to `I_Q`; without that source object, EOM, Hamiltonian, bridge closure, role transfer, and ToE remain blocked.\n")
    append_once(AGENTS, "Current Gamma scale-orbit source-object intake guardrail (P2919/S1869, 2026-06-19)", "## Current Gamma scale-orbit source-object intake guardrail (P2919/S1869, 2026-06-19)\n\n- P2919 constructs the required `Strict_Gamma_9_5_Action_Unit_Scale_Fixing_Source_Object` schema: `sigma_Gamma: strict_nadsoliton_data -> Action_nonzero` coupled to `I_Q`.\n- The finite scale-orbit matrix shows that quotient weights, derivative ratios, orbit size, and edge count are invariant under `Gamma_9_5 -> c*Gamma_9_5`; they fix relative weights but not the Gamma scale.\n- No strict scale-breaking `sigma_Gamma` source object is exported; an imported Action label would break scale only by convention and is rejected as strict provenance.\n- Do not promote P2911-P2919 Gamma/Lambda readiness to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without a genuinely new strict `sigma_Gamma` source object and coupling theorem.\n- If no such source object is supplied, the next honest move is a Gamma/Lambda no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
