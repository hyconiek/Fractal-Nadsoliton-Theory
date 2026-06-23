#!/usr/bin/env python3
"""P3046/S1996: memory-lag torsor coupling-polarity audit.

P3045 left an explicit next premise: a coupling theorem from the signed
memory-lag torsor to a selector/readout torsor.  P3046 attacks only that
coupling premise.  It does not replay lag-localizer scoring and it does not
claim a strict source law for the commutator.

The finite object is the sign torsor of P3045 integrated lag scores.  Under
Aut(Z12) inversion, the torsor sign flips.  This representation type can couple
to the orientation/selector torsor: exactly two equivariant polarity maps exist.
The obstruction is that they are opposite polarity choices, and no strict law in
current artifacts selects one coupling polarity or installs it as a readout row.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3045_s1995_memory_lag_localizer_equivariance_obstruction import OUT as P3045

OUT = GEN / "p3046_s1996_memory_lag_torsor_coupling_polarity_audit.json"
MD = GEN / "p3046_s1996_memory_lag_torsor_coupling_polarity_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SIGN_TORSOR = [-1, 1]
ORIENTATION_TORSOR = ["-omega", "+omega"]
UNITS = [1, 5, 7, 11]
INVERSION_UNITS = [7, 11]
PRESERVING_UNITS = [1, 5]


def sign_action(unit: int, sign: int) -> int:
    return -sign if unit in INVERSION_UNITS else sign


def orientation_action(unit: int, orientation: str) -> str:
    if unit not in INVERSION_UNITS:
        return orientation
    return "+omega" if orientation == "-omega" else "-omega"


def candidate_couplings() -> list[dict[str, Any]]:
    return [
        {"polarity": "+", "map": {-1: "-omega", 1: "+omega"}},
        {"polarity": "-", "map": {-1: "+omega", 1: "-omega"}},
    ]


def is_equivariant(mapping: dict[int, str]) -> bool:
    for unit in UNITS:
        for sign in SIGN_TORSOR:
            lhs = mapping[sign_action(unit, sign)]
            rhs = orientation_action(unit, mapping[sign])
            if lhs != rhs:
                return False
    return True


def coupling_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate in candidate_couplings():
        mapping = candidate["map"]
        rows.append({
            "polarity": candidate["polarity"],
            "maps_negative_sign_to": mapping[-1],
            "maps_positive_sign_to": mapping[1],
            "aut_equivariant": is_equivariant(mapping),
            "polarity_selected_by_current_artifacts": False,
            "accepted_selector_readout_coupling": False,
            "failure": "equivariant coupling exists only as one of two opposite polarity choices; no strict source law selects this row",
        })
    return rows


def source_readout_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "p3045_signed_lag_torsor_to_orientation_torsor",
            "finite_positive": True,
            "strict_polarity_selection_law": False,
            "nonpremise_selector_readout_installation": False,
            "accepted_as_coupling_theorem": False,
            "reason": "representation-level coupling exists, but coupling polarity and readout installation remain unsourced",
        },
        {
            "candidate": "lag2_receiver_winner_to_selector_row",
            "finite_positive": True,
            "strict_polarity_selection_law": False,
            "nonpremise_selector_readout_installation": False,
            "accepted_as_coupling_theorem": False,
            "reason": "lag-2 is a receiver winner from P3045, not a theorem selecting a torsor polarity",
        },
        {
            "candidate": "abstract_signed_source_counterfactual",
            "finite_positive": False,
            "strict_polarity_selection_law": False,
            "nonpremise_selector_readout_installation": False,
            "accepted_as_coupling_theorem": False,
            "reason": "records the missing object type; no concrete strict signed source artifact is supplied",
        },
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3045)
    couplings = coupling_rows()
    readouts = source_readout_rows()
    obligations = [
        {"obligation": "p3045_signed_lag_torsor_used", "satisfied": True, "detail": "P3046 attacks only the selector/readout coupling premise left by P3045"},
        {"obligation": "finite_equivariant_coupling_pair_enumerated", "satisfied": all(row["aut_equivariant"] for row in couplings), "detail": "two Aut-equivariant maps from sign torsor to orientation torsor are explicit"},
        {"obligation": "coupling_polarity_pair_exposed", "satisfied": len(couplings) == 2, "detail": "the two maps are opposite polarity choices"},
        {"obligation": "strict_polarity_selection_law", "satisfied": False, "detail": "current artifacts do not select one of the two equivariant polarity maps"},
        {"obligation": "nonpremise_selector_readout_installation", "satisfied": False, "detail": "no row installs the coupling as QW-2191 discharge or physical readout"},
        {"obligation": "strict_memory_lag_source_law", "satisfied": False, "detail": "P3044/P3045 commutator and lag polarity still lack strict nadsoliton provenance"},
    ]
    return {
        "object": "MemoryLagTorsor_CouplingPolarityAudit",
        "tested_premise": "explicit coupling theorem from signed memory-lag torsor to selector/readout torsor",
        "sign_torsor": SIGN_TORSOR,
        "orientation_torsor": ORIENTATION_TORSOR,
        "aut_units": UNITS,
        "coupling_rows": couplings,
        "source_readout_rows": readouts,
        "proof_obligations": obligations,
        "finite_certificate": {
            "sign_torsor_size": len(SIGN_TORSOR),
            "orientation_torsor_size": len(ORIENTATION_TORSOR),
            "candidate_coupling_rows": len(couplings),
            "aut_equivariant_coupling_rows": sum(1 for row in couplings if row["aut_equivariant"]),
            "polarity_selected_rows": sum(1 for row in couplings if row["polarity_selected_by_current_artifacts"]),
            "accepted_selector_readout_coupling_rows": sum(1 for row in couplings if row["accepted_selector_readout_coupling"]),
            "source_readout_rows": len(readouts),
            "accepted_source_readout_rows": sum(1 for row in readouts if row["accepted_as_coupling_theorem"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "selector_readout_coupling_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3046_MEMORY_LAG_TORSOR_COUPLING_POLARITY_AUDIT_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3045": hashlib.sha256(P3045.read_bytes()).hexdigest() if P3045.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3046 finds real representation-level progress: the signed memory-lag torsor has exactly two Aut-equivariant maps to the orientation/selector torsor.  But the maps are opposite coupling-polarity choices, and current artifacts export no strict law selecting one polarity or installing it as a nonpremise selector/readout row.",
            "negative_export_flags": {k: False for k in ["selector_readout_coupling_exported", "coupling_polarity_selected", "strict_memory_lag_source_law_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote the two conditional coupling maps to selector closure.  A next proof-grade move must supply a strict polarity-selection/source law for the memory-lag sign, or pivot to a different new typed object; otherwise preserve the P3044-P3046 bounded no-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3046/S1996 memory-lag torsor coupling-polarity audit", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- sign torsor size: `{c['sign_torsor_size']}`",
        f"- orientation torsor size: `{c['orientation_torsor_size']}`",
        f"- candidate coupling rows: `{c['candidate_coupling_rows']}`",
        f"- Aut-equivariant coupling rows: `{c['aut_equivariant_coupling_rows']}`",
        f"- polarity-selected rows: `{c['polarity_selected_rows']}`",
        f"- accepted selector/readout coupling rows: `{c['accepted_selector_readout_coupling_rows']}`",
        f"- source/readout rows: `{c['source_readout_rows']}`",
        f"- accepted source/readout rows: `{c['accepted_source_readout_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- selector/readout coupling exported: `{c['selector_readout_coupling_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3046/S1996 memory-lag torsor coupling-polarity audit", "## P3046/S1996 memory-lag torsor coupling-polarity audit\n\n`P3046/S1996` attacks exactly one P3045 remaining premise: an explicit coupling theorem from the signed memory-lag torsor to a selector/readout torsor.  The finite audit finds real representation-level progress: exactly two Aut-equivariant maps exist from the lag-sign torsor to the orientation torsor.  They are opposite polarity choices, and current artifacts export no strict polarity-selection law, no nonpremise selector/readout installation, and no strict memory-lag source law.  Therefore the conditional coupling pair does not discharge `QW-2191`, export selector closure, install observed physics, add unit-bearing action/EOM or `L_total`, close the bridge, transfer roles, or close ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3046/S1996 memory-lag torsor coupling `L_total` guard", "## P3046/S1996 memory-lag torsor coupling `L_total` guard\n\n`P3046/S1996` adds no physical `L_total` term.  The signed lag torsor has two Aut-equivariant orientation-torsor couplings, but neither polarity is strictly selected and no selector/readout row or unit-bearing variational/action/EOM installation is exported.\n")
    append_once(AGENTS, "Current memory-lag torsor coupling-polarity guardrail (P3046/S1996, 2026-06-23)", "## Current memory-lag torsor coupling-polarity guardrail (P3046/S1996, 2026-06-23)\n\n- P3046 attacks exactly one P3045 remaining premise: explicit coupling from the signed memory-lag torsor to a selector/readout torsor.\n- Exactly two Aut-equivariant coupling maps exist, but they are opposite polarity choices; no strict source law selects one polarity or installs it as a nonpremise selector/readout row.\n- Do not promote the conditional coupling pair, lag-sign torsor, lag-2 receiver winner, or memory-lag commutator sign to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move requires a strict polarity-selection/source law for the memory-lag sign, or a different genuinely new typed object; otherwise preserve the P3044-P3046 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
