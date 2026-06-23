#!/usr/bin/env python3
"""P3047/S1997: memory-lag polarity source-law obstruction.

P3046 found exactly two conditional couplings from the signed memory-lag torsor
to the orientation/selector torsor.  The remaining honest premise is a strict
polarity-selection/source law for the memory-lag sign.  P3047 attacks only that
premise.

The finite theorem is a source-domain obstruction.  A nonpremise strict source
with trivial Aut(Z12) action cannot select a point of the memory-lag sign torsor,
because inversion units 7 and 11 flip the sign.  An inversion-odd source domain
would have the correct representation type and gives two equivariant maps, but
current artifacts do not export a nonzero strict odd source value; therefore the
polarity remains unsourced.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3046_s1996_memory_lag_torsor_coupling_polarity_audit import OUT as P3046, SIGN_TORSOR, UNITS, INVERSION_UNITS

OUT = GEN / "p3047_s1997_memory_lag_polarity_source_law_obstruction.json"
MD = GEN / "p3047_s1997_memory_lag_polarity_source_law_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def target_action(unit: int, sign: int) -> int:
    return -sign if unit in INVERSION_UNITS else sign


def domain_action(domain: str, unit: int, value: int) -> int:
    if domain == "inversion_odd_sign_domain" and unit in INVERSION_UNITS:
        return -value
    return value


def equivariant_maps(domain: str, values: list[int]) -> list[dict[int, int]]:
    maps = []
    if len(values) == 1:
        possible = [{values[0]: target} for target in SIGN_TORSOR]
    else:
        possible = [
            {values[0]: SIGN_TORSOR[a], values[1]: SIGN_TORSOR[b]}
            for a in range(len(SIGN_TORSOR))
            for b in range(len(SIGN_TORSOR))
        ]
    for mapping in possible:
        ok = True
        for unit in UNITS:
            for value in values:
                lhs = mapping[domain_action(domain, unit, value)]
                rhs = target_action(unit, mapping[value])
                if lhs != rhs:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            maps.append(mapping)
    return maps


def source_domain_rows() -> list[dict[str, Any]]:
    specs = [
        ("trivial_singleton_source", [0], False, "Aut-trivial source datum cannot map equivariantly to an inversion-odd sign torsor"),
        ("trivial_binary_even_source", [-1, 1], False, "two labels with trivial Aut action are still sign-blind; equivariance would require fixed target signs"),
        ("p3044_even_receiver_magnitude_source", [0], False, "commutator magnitude/norm data are even receiver signatures, not signed source values"),
        ("p3045_lag2_receiver_winner_source", [0], False, "lag-2 winner is a chart receiver label and is Aut-trivial after forgetting the inserted polarity"),
        ("inversion_odd_sign_domain", [-1, 1], False, "correct representation type exists abstractly, but no nonzero strict odd source value is exported"),
    ]
    rows = []
    for domain, values, exported_input, failure in specs:
        maps = equivariant_maps(domain, values)
        rows.append({
            "source_domain": domain,
            "domain_size": len(values),
            "equivariant_maps_to_memory_lag_sign_torsor": len(maps),
            "has_correct_representation_type": len(maps) > 0,
            "nonzero_strict_source_value_exported": exported_input,
            "accepted_polarity_source_law": exported_input and len(maps) == 1,
            "failure": failure,
        })
    return rows


def candidate_law_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate_law": "choose_positive_memory_lag_sign",
            "computable": True,
            "aut_compatible": False,
            "nonpremise": False,
            "accepted": False,
            "reason": "chooses a target sign by convention; inversion sends it to the opposite sign",
        },
        {
            "candidate_law": "choose_lag2_winner_coupling_polarity",
            "computable": True,
            "aut_compatible": False,
            "nonpremise": False,
            "accepted": False,
            "reason": "imports the P3045 receiver winner as polarity source rather than deriving a strict sign",
        },
        {
            "candidate_law": "abstract_inversion_odd_strict_source_value",
            "computable": False,
            "aut_compatible": True,
            "nonpremise": False,
            "accepted": False,
            "reason": "records the missing object type but no concrete strict signed value is present",
        },
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3046)
    source_rows = source_domain_rows()
    laws = candidate_law_rows()
    obligations = [
        {"obligation": "p3046_coupling_polarity_gap_used", "satisfied": True, "detail": "P3047 attacks only the source law needed to select one P3046 coupling polarity"},
        {"obligation": "trivial_source_no_section_theorem", "satisfied": all(row["equivariant_maps_to_memory_lag_sign_torsor"] == 0 for row in source_rows[:4]), "detail": "Aut-trivial source domains have zero equivariant maps to the sign torsor"},
        {"obligation": "odd_source_representation_identified", "satisfied": source_rows[-1]["equivariant_maps_to_memory_lag_sign_torsor"] == 2, "detail": "an inversion-odd source domain would have two equivariant maps"},
        {"obligation": "nonzero_strict_odd_source_value", "satisfied": False, "detail": "no concrete nonzero strict odd source value for memory-lag polarity is exported"},
        {"obligation": "unique_polarity_selection_law", "satisfied": False, "detail": "no current row selects exactly one polarity nonconventionally"},
        {"obligation": "selector_readout_installation", "satisfied": False, "detail": "no selected polarity is installed as QW-2191 discharge or physical readout"},
    ]
    return {
        "object": "MemoryLagPolarity_SourceLawObstructionMatrix",
        "tested_premise": "strict polarity-selection/source law for the memory-lag sign",
        "sign_torsor": SIGN_TORSOR,
        "source_domain_rows": source_rows,
        "candidate_law_rows": laws,
        "proof_obligations": obligations,
        "finite_certificate": {
            "source_domain_rows": len(source_rows),
            "trivial_source_rows": 4,
            "trivial_source_equivariant_maps": sum(row["equivariant_maps_to_memory_lag_sign_torsor"] for row in source_rows[:4]),
            "odd_source_equivariant_maps": source_rows[-1]["equivariant_maps_to_memory_lag_sign_torsor"],
            "candidate_law_rows": len(laws),
            "accepted_polarity_source_law_rows": sum(1 for row in source_rows if row["accepted_polarity_source_law"]) + sum(1 for row in laws if row["accepted"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_polarity_source_law_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3047_MEMORY_LAG_POLARITY_SOURCE_LAW_OBSTRUCTION_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3046": hashlib.sha256(P3046.read_bytes()).hexdigest() if P3046.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3047 proves the finite source-domain obstruction for memory-lag polarity.  Aut-trivial source data have zero equivariant maps to the inversion-odd memory-lag sign torsor.  An inversion-odd source domain has the right representation type and two equivariant maps, but current artifacts export no nonzero strict odd source value, so no unique polarity-selection law is available.",
            "negative_export_flags": {k: False for k in ["strict_polarity_source_law_exported", "coupling_polarity_selected", "selector_readout_coupling_exported", "strict_memory_lag_source_law_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay positive sign, lag-2 winner, or even commutator magnitudes as a polarity source.  A next admissible move must supply one concrete nonzero strict inversion-odd source value coupled to the memory-lag sign, or pivot to a different genuinely new typed object; otherwise preserve the P3044-P3047 bounded no-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3047/S1997 memory-lag polarity source-law obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- source domain rows: `{c['source_domain_rows']}`",
        f"- trivial source rows: `{c['trivial_source_rows']}`",
        f"- trivial-source equivariant maps: `{c['trivial_source_equivariant_maps']}`",
        f"- odd-source equivariant maps: `{c['odd_source_equivariant_maps']}`",
        f"- candidate law rows: `{c['candidate_law_rows']}`",
        f"- accepted polarity source-law rows: `{c['accepted_polarity_source_law_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict polarity source law exported: `{c['strict_polarity_source_law_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3047/S1997 memory-lag polarity source-law obstruction", "## P3047/S1997 memory-lag polarity source-law obstruction\n\n`P3047/S1997` attacks exactly one P3046 remaining premise: a strict polarity-selection/source law for the memory-lag sign.  The finite source-domain theorem shows that Aut-trivial source data have zero equivariant maps to the inversion-odd memory-lag sign torsor.  An abstract inversion-odd source domain has the right representation type and two equivariant maps, but no concrete nonzero strict odd source value is exported.  Therefore no unique polarity-selection law, selected P3046 coupling polarity, selector/readout installation, `QW-2191` discharge, observed physics, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3047/S1997 memory-lag polarity source `L_total` guard", "## P3047/S1997 memory-lag polarity source `L_total` guard\n\n`P3047/S1997` adds no physical `L_total` term.  It identifies the required inversion-odd source representation for memory-lag polarity, but current artifacts export no nonzero strict odd source value, no unique polarity law, and no unit-bearing variational/action/EOM installation.\n")
    append_once(AGENTS, "Current memory-lag polarity source-law guardrail (P3047/S1997, 2026-06-23)", "## Current memory-lag polarity source-law guardrail (P3047/S1997, 2026-06-23)\n\n- P3047 attacks exactly one P3046 remaining premise: strict polarity-selection/source law for the memory-lag sign.\n- Aut-trivial source data have zero equivariant maps to the inversion-odd memory-lag sign torsor; an inversion-odd domain has the right representation type but no concrete nonzero strict odd source value is exported.\n- Do not promote positive sign conventions, lag-2 receiver winners, even commutator magnitudes, abstract odd source slots, or conditional P3046 coupling maps to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move requires one concrete nonzero strict inversion-odd source value coupled to memory-lag polarity, or a different genuinely new typed object; otherwise preserve the P3044-P3047 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
