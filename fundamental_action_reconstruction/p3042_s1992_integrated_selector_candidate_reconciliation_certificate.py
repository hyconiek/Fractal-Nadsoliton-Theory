#!/usr/bin/env python3
"""P3042/S1992: integrated selector-candidate reconciliation certificate.

P3038 constructed a branch-separating integrated selector candidate.  P3039,
P3040, and P3041 then audited its direct source premises.  P3042 is the finite
reconciliation step: build the Boolean source-premise lattice and decide whether
the candidate can be promoted to a strict selector mechanism on current
artifacts.

The result is proof-grade but bounded.  The operator/score atom is real, while
all source-export atoms remain open/false, so only the all-source-closed
counterfactual profile would accept.  The current profile is a no-selector-export
certificate, not a new closure claim.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import OUT as P3038
from p3039_s1989_chiral_projector_sign_localizer_obstruction import OUT as P3039
from p3040_s1990_retarded_path_anisotropy_source_obstruction import OUT as P3040
from p3041_s1991_integrated_density_unit_readout_coupling_obstruction import OUT as P3041

OUT = GEN / "p3042_s1992_integrated_selector_candidate_reconciliation_certificate.json"
MD = GEN / "p3042_s1992_integrated_selector_candidate_reconciliation_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ATOM_ORDER = [
    "branch_separating_integrated_operator",
    "nonpremise_chiral_sign_localizer",
    "sourced_retardation_path_anisotropy",
    "physical_unit_readout_coupling",
    "nonproxy_variational_or_readout_installation",
]


def load_inputs() -> dict[str, Any]:
    return {"P3038": read_json(P3038), "P3039": read_json(P3039), "P3040": read_json(P3040), "P3041": read_json(P3041)}


def atom_ledger(inputs: dict[str, Any]) -> list[dict[str, Any]]:
    p3038 = inputs["P3038"]
    p3039 = inputs["P3039"]
    p3040 = inputs["P3040"]
    p3041 = inputs["P3041"]
    return [
        {
            "atom": "branch_separating_integrated_operator",
            "source": "P3038",
            "closed": bool(p3038["finite_certificate"]["branch_score_nonzero"]),
            "evidence": "P3038 finite branch score is nonzero and opposite across branches",
            "closure_kind": "computational_positive_not_source_export",
        },
        {
            "atom": "nonpremise_chiral_sign_localizer",
            "source": "P3039",
            "closed": bool(p3039["finite_certificate"]["nonpremise_chiral_sign_localizer_exported"]),
            "evidence": "P3039 leaves phase origin and polarity unsourced",
            "closure_kind": "required_source_atom_open",
        },
        {
            "atom": "sourced_retardation_path_anisotropy",
            "source": "P3040",
            "closed": bool(p3040["finite_certificate"]["sourced_retardation_path_anisotropy_exported"]),
            "evidence": "P3040 leaves rho and retarded path sectors unsourced",
            "closure_kind": "required_source_atom_open",
        },
        {
            "atom": "physical_unit_readout_coupling",
            "source": "P3041",
            "closed": bool(p3041["finite_certificate"]["physical_unit_readout_coupling_exported"]),
            "evidence": "P3041 leaves external unit and readout coupling unsourced",
            "closure_kind": "required_source_atom_open",
        },
        {
            "atom": "nonproxy_variational_or_readout_installation",
            "source": "P3038/P3041",
            "closed": False,
            "evidence": "P3038 keeps nonproxy variational/readout theorem false and P3041 keeps action/readout installation false",
            "closure_kind": "required_installation_atom_open",
        },
    ]


def profile_lattice() -> list[dict[str, Any]]:
    profiles = []
    for bits in itertools.product([False, True], repeat=len(ATOM_ORDER)):
        profile = dict(zip(ATOM_ORDER, bits))
        accepted = all(profile.values())
        profiles.append({"profile": profile, "accepted_as_strict_selector_mechanism": accepted})
    return profiles


def build_matrix() -> dict[str, Any]:
    inputs = load_inputs()
    ledger = atom_ledger(inputs)
    current_profile = {row["atom"]: row["closed"] for row in ledger}
    profiles = profile_lattice()
    current_accepted = all(current_profile.values())
    obligations = [
        {"obligation": "p3038_integrated_operator_accounted", "satisfied": current_profile["branch_separating_integrated_operator"], "detail": "finite branch separation is retained as a positive computational atom"},
        {"obligation": "p3039_chiral_source_audited", "satisfied": True, "detail": "chiral sign/localizer source premise has a dedicated obstruction audit"},
        {"obligation": "p3040_path_source_audited", "satisfied": True, "detail": "retarded path-anisotropy source premise has a dedicated obstruction audit"},
        {"obligation": "p3041_unit_readout_audited", "satisfied": True, "detail": "physical unit/readout coupling premise has a dedicated obstruction audit"},
        {"obligation": "all_required_source_atoms_closed", "satisfied": all(row["closed"] for row in ledger[1:]), "detail": "chiral/path/unit/nonproxy source atoms must all close"},
        {"obligation": "current_profile_accepts_selector_export", "satisfied": current_accepted, "detail": "only all five atoms closed would accept strict selector export"},
    ]
    return {
        "object": "IntegratedSelectorCandidate_ReconciliationNoExportCertificate",
        "atom_order": ATOM_ORDER,
        "atom_ledger": ledger,
        "current_profile": current_profile,
        "profile_lattice": profiles,
        "proof_obligations": obligations,
        "finite_certificate": {
            "atom_count": len(ledger),
            "closed_atoms": sum(1 for row in ledger if row["closed"]),
            "open_required_source_atoms": sum(1 for row in ledger[1:] if not row["closed"]),
            "profile_count": len(profiles),
            "accepted_counterfactual_profiles": sum(1 for row in profiles if row["accepted_as_strict_selector_mechanism"]),
            "current_profile_accepted": current_accepted,
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_selector_mechanism_exported": current_accepted,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3042_INTEGRATED_SELECTOR_CANDIDATE_RECONCILIATION_NO_SELECTOR_EXPORT",
        "input_hashes": {
            "P3038": hashlib.sha256(P3038.read_bytes()).hexdigest() if P3038.exists() else None,
            "P3039": hashlib.sha256(P3039.read_bytes()).hexdigest() if P3039.exists() else None,
            "P3040": hashlib.sha256(P3040.read_bytes()).hexdigest() if P3040.exists() else None,
            "P3041": hashlib.sha256(P3041.read_bytes()).hexdigest() if P3041.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3042 reconciles the P3038 integrated selector candidate after its direct source-premise audits.  The finite branch-separating operator is real, but the chiral sign/localizer, retarded path anisotropy, physical unit/readout coupling, and nonproxy variational/readout installation atoms remain open.  The 32-profile lattice has exactly one accepting counterfactual profile, the all-five-atoms-closed profile; the current profile is not accepted.",
            "negative_export_flags": {k: False for k in ["strict_selector_mechanism_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not continue replaying P3038 receiver classes.  A next admissible move requires a genuinely new strict source law outside sine/chiral phase receivers, rho/path tuning, unit normalization/imports, and nonproxy placeholders; otherwise preserve this P3038-P3042 no-selector-export certificate and pivot to a broad state-map/new-object intake.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3042/S1992 integrated selector-candidate reconciliation certificate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- atom count: `{c['atom_count']}`",
        f"- closed atoms: `{c['closed_atoms']}`",
        f"- open required source atoms: `{c['open_required_source_atoms']}`",
        f"- profile count: `{c['profile_count']}`",
        f"- accepted counterfactual profiles: `{c['accepted_counterfactual_profiles']}`",
        f"- current profile accepted: `{c['current_profile_accepted']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict selector mechanism exported: `{c['strict_selector_mechanism_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3042/S1992 integrated selector-candidate reconciliation certificate", "## P3042/S1992 integrated selector-candidate reconciliation certificate\n\n`P3042/S1992` reconciles P3038-P3041 by constructing the five-atom selector-candidate source lattice: branch-separating integrated operator, nonpremise chiral sign/localizer, sourced retarded path anisotropy, physical unit/readout coupling, and nonproxy variational/readout installation.  The branch-separating operator atom is closed as a computational positive, but the four required source/installation atoms remain open.  The finite Boolean lattice has `32` profiles and only the all-five-atoms-closed counterfactual profile accepts selector export; the current profile does not.  No `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3042/S1992 integrated selector candidate reconciliation `L_total` guard", "## P3042/S1992 integrated selector candidate reconciliation `L_total` guard\n\n`P3042/S1992` adds no physical `L_total` term.  Its lattice records that P3038's finite branch-separating score remains real, but chiral source, path source, unit/readout source, and nonproxy variational/readout installation are all still open.\n")
    append_once(AGENTS, "Current integrated selector-candidate reconciliation guardrail (P3042/S1992, 2026-06-23)", "## Current integrated selector-candidate reconciliation guardrail (P3042/S1992, 2026-06-23)\n\n- P3042 reconciles P3038-P3041 in a five-atom source lattice: branch-separating integrated operator, chiral sign/localizer, retarded path anisotropy, physical unit/readout coupling, and nonproxy variational/readout installation.\n- Only the computational branch-separation atom is closed; the four required source/installation atoms remain open, and only the all-five-atoms-closed counterfactual profile would export a strict selector mechanism.\n- Do not replay P3038 receiver classes, sine/chiral phase receivers, rho/path tuning, unit normalization/imports, or nonproxy placeholders as `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move requires a genuinely new strict source law outside these exhausted receiver classes, or a broad state-map/new-object intake rather than more replay.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
