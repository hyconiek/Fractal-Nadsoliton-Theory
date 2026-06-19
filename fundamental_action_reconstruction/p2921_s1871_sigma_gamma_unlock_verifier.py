#!/usr/bin/env python3
"""P2921/S1871: executable sigma_Gamma unlock verifier.

P2920 closed the current Gamma/Lambda finite chain unless exactly one new typed
object is supplied: a strict sigma_Gamma Action source map with a computed
nonzero value and coupling to the P2917 quotient integral.  P2921 does not
replay quotient/Jacobian arithmetic.  It constructs the missing verifier that
any proposed source map must pass before it can unlock nonproxy L_total work.

The verifier is deliberately finite and proof-facing: it formalizes the minimal
unlock packet into five Boolean obligations and audits a small candidate family.
Current candidates still fail, so no L_total/EOM/Hamiltonian/ToE closure is
exported.  The useful output is an executable acceptance object for the next
honest source proposal.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2920 = GEN / "p2920_s1870_gamma_lambda_no_new_live_frontier_certificate.json"
OUT = GEN / "p2921_s1871_sigma_gamma_unlock_verifier.json"
MD = GEN / "p2921_s1871_sigma_gamma_unlock_verifier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

OBLIGATIONS = [
    "computed_nonzero_value",
    "strict_nadsoliton_provenance",
    "action_dimension_exported",
    "scale_orbit_broken_nonconventionally",
    "explicit_IQ_coupling",
]


def candidate_sigma_maps() -> list[dict[str, Any]]:
    """Candidate source maps to run through the P2920 minimal unlock verifier."""
    return [
        {
            "name": "sigma_Gamma_from_quotient_weight",
            "formula": "sigma_Gamma = 1/12",
            "computed_nonzero_value": True,
            "strict_nadsoliton_provenance": False,
            "action_dimension_exported": False,
            "scale_orbit_broken_nonconventionally": False,
            "explicit_IQ_coupling": True,
            "failure": "dimensionless quotient weight is not strict Action provenance and fixes scale only by convention",
        },
        {
            "name": "sigma_Gamma_from_edge_weight",
            "formula": "sigma_Gamma = 1/144",
            "computed_nonzero_value": True,
            "strict_nadsoliton_provenance": False,
            "action_dimension_exported": False,
            "scale_orbit_broken_nonconventionally": False,
            "explicit_IQ_coupling": True,
            "failure": "edge-normalization weight is dimensionless quotient arithmetic, not an Action source",
        },
        {
            "name": "sigma_Gamma_from_orbit_size",
            "formula": "sigma_Gamma = 12",
            "computed_nonzero_value": True,
            "strict_nadsoliton_provenance": False,
            "action_dimension_exported": False,
            "scale_orbit_broken_nonconventionally": False,
            "explicit_IQ_coupling": True,
            "failure": "orbit size is a counting invariant already consumed by quotienting",
        },
        {
            "name": "sigma_Gamma_imported_action_unit",
            "formula": "sigma_Gamma = ActionUnit_external",
            "computed_nonzero_value": True,
            "strict_nadsoliton_provenance": False,
            "action_dimension_exported": True,
            "scale_orbit_broken_nonconventionally": False,
            "explicit_IQ_coupling": True,
            "failure": "external unit has dimension but lacks strict nadsoliton provenance and nonconventional scale breaking",
        },
        {
            "name": "sigma_Gamma_zero",
            "formula": "sigma_Gamma = 0",
            "computed_nonzero_value": False,
            "strict_nadsoliton_provenance": False,
            "action_dimension_exported": True,
            "scale_orbit_broken_nonconventionally": False,
            "explicit_IQ_coupling": True,
            "failure": "zero value cannot source a nonzero action coefficient",
        },
        {
            "name": "sigma_Gamma_strict_source_placeholder",
            "formula": "sigma_Gamma(strict_nadsoliton_data) = gamma_*",
            "computed_nonzero_value": False,
            "strict_nadsoliton_provenance": False,
            "action_dimension_exported": False,
            "scale_orbit_broken_nonconventionally": False,
            "explicit_IQ_coupling": True,
            "failure": "placeholder names the required theorem but supplies no computed value or provenance",
        },
    ]


def evaluate_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    passed = {obligation: bool(candidate[obligation]) for obligation in OBLIGATIONS}
    accepted = all(passed.values())
    missing = [name for name, ok in passed.items() if not ok]
    return {
        "name": candidate["name"],
        "formula": candidate["formula"],
        "obligation_results": passed,
        "accepted": accepted,
        "missing_obligations": missing,
        "failure": None if accepted else candidate["failure"],
    }


def verifier_schema() -> dict[str, Any]:
    return {
        "verifier_name": "Strict_sigma_Gamma_Action_Source_Map_Unlock_Verifier",
        "input_type": "candidate sigma_Gamma source map plus coupling to I_Q",
        "output_type": "accepted unlock witness or rejected candidate with missing obligations",
        "acceptance_obligations": OBLIGATIONS,
        "accepted_unlock_consequence": "only after acceptance may a later step attempt nonproxy L_total coupling; this verifier itself does not export L_total",
    }


def build_payload(p2920: dict[str, Any]) -> dict[str, Any]:
    candidates = [evaluate_candidate(candidate) for candidate in candidate_sigma_maps()]
    accepted = [candidate for candidate in candidates if candidate["accepted"]]
    obligation_failure_counts = {
        obligation: sum(1 for candidate in candidates if not candidate["obligation_results"][obligation])
        for obligation in OBLIGATIONS
    }
    return {
        "status": "P2921_SIGMA_GAMMA_UNLOCK_VERIFIER_EXECUTED_NO_ACCEPTED_SOURCE",
        "input_hashes": {"P2920": hashlib.sha256(P2920.read_bytes()).hexdigest() if P2920.exists() else None},
        "constructed_theoretical_objects": {
            "verifier_schema": verifier_schema(),
            "candidate_evaluations": candidates,
            "obligation_failure_counts": obligation_failure_counts,
            "minimal_success_shape": {
                "required_formula": "sigma_Gamma(strict_nadsoliton_data) = gamma_* in Action_nonzero",
                "required_coupling": "I_Q = sigma_Gamma * (1/12) * sum_d Q_d",
                "required_status": "all five obligations true",
            },
        },
        "acceptance_matrix": {
            "p2920_rechecked_certificate": p2920.get("acceptance_matrix", {}).get("no_new_live_frontier_certificate_exported") is True,
            "obligation_count": len(OBLIGATIONS),
            "candidate_count": len(candidates),
            "accepted_candidate_count": len(accepted),
            "computed_nonzero_value_pass_count": sum(1 for c in candidates if c["obligation_results"]["computed_nonzero_value"]),
            "strict_nadsoliton_provenance_pass_count": sum(1 for c in candidates if c["obligation_results"]["strict_nadsoliton_provenance"]),
            "action_dimension_pass_count": sum(1 for c in candidates if c["obligation_results"]["action_dimension_exported"]),
            "scale_orbit_breaking_pass_count": sum(1 for c in candidates if c["obligation_results"]["scale_orbit_broken_nonconventionally"]),
            "explicit_IQ_coupling_pass_count": sum(1 for c in candidates if c["obligation_results"]["explicit_IQ_coupling"]),
            "strict_sigma_gamma_unlock_verifier_exported": True,
            "strict_sigma_gamma_source_accepted": False,
            "accepted_as_nonproxy_ltotal_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "unlock_verifier_constructed": True,
                "candidate_family_evaluated": True,
                "missing_obligations_reported": True,
            },
            "negative_export_flags": {
                "strict_sigma_gamma_source_accepted": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2921 turns the P2920 minimal unlock packet into an executable verifier.  Six candidate sigma_Gamma maps are evaluated against five obligations.  Existing quotient/count/imported/zero/placeholder candidates fail strict provenance, Action sourcehood, nonconventional scale breaking, or nonzero value obligations; accepted candidates remain zero.  The verifier is exported, not a source theorem or L_total closure.",
            "next_honest_step": "A future move should supply exactly one concrete sigma_Gamma formula and run this verifier.  If it passes all five obligations, only then audit nonproxy L_total coupling; if it fails or no formula is supplied, pivot outside the Gamma/Lambda lane.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2921/S1871 sigma_Gamma unlock verifier",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Verifier gate",
        f"- obligations: `{acc['obligation_count']}`",
        f"- candidates: `{acc['candidate_count']}`",
        f"- accepted candidates: `{acc['accepted_candidate_count']}`",
        f"- computed nonzero value pass count: `{acc['computed_nonzero_value_pass_count']}`",
        f"- strict nadsoliton provenance pass count: `{acc['strict_nadsoliton_provenance_pass_count']}`",
        f"- Action dimension pass count: `{acc['action_dimension_pass_count']}`",
        f"- scale-orbit breaking pass count: `{acc['scale_orbit_breaking_pass_count']}`",
        f"- explicit I_Q coupling pass count: `{acc['explicit_IQ_coupling_pass_count']}`",
        f"- verifier exported: `{acc['strict_sigma_gamma_unlock_verifier_exported']}`",
        f"- strict sigma_Gamma source accepted: `{acc['strict_sigma_gamma_source_accepted']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2920))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2921/S1871 sigma_Gamma unlock verifier", "## P2921/S1871 sigma_Gamma unlock verifier\n\n`P2921/S1871` constructs the executable `Strict_sigma_Gamma_Action_Source_Map_Unlock_Verifier` required by the P2920 minimal unlock packet.  It tests six candidate `sigma_Gamma` maps against five obligations: computed nonzero value, strict nadsoliton provenance, Action dimension, nonconventional scale-orbit breaking, and explicit coupling to `I_Q`.  No candidate passes all obligations; quotient/count/imported/zero/placeholder candidates remain non-strict or incomplete.  The verifier is exported as a future acceptance object, but no strict `sigma_Gamma` source, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2921/S1871 sigma_Gamma verifier `L_total` guard", "## P2921/S1871 sigma_Gamma verifier `L_total` guard\n\n`P2921/S1871` provides an executable verifier for any proposed `sigma_Gamma` source map.  The current candidate family has `0` accepted sources, so the Gamma/Lambda lane remains blocked for nonproxy `L_total`.  A later source proposal must pass all five verifier obligations before any EOM, Hamiltonian, bridge, role-transfer, or ToE promotion is considered.\n")
    append_once(AGENTS, "Current sigma_Gamma unlock verifier guardrail (P2921/S1871, 2026-06-19)", "## Current sigma_Gamma unlock verifier guardrail (P2921/S1871, 2026-06-19)\n\n- P2921 exports an executable `Strict_sigma_Gamma_Action_Source_Map_Unlock_Verifier` for the P2920 minimal unlock packet.\n- A candidate `sigma_Gamma` source must pass five obligations: computed nonzero value, strict nadsoliton provenance, Action dimension, nonconventional scale-orbit breaking, and explicit `I_Q` coupling.\n- The current six-candidate family has `0` accepted sources; quotient/count/imported/zero/placeholder candidates do not unlock nonproxy `L_total`.\n- Do not promote Gamma/Lambda readiness to `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE unless a concrete `sigma_Gamma` formula passes the verifier.\n- If no such formula is supplied, pivot outside the Gamma/Lambda lane rather than replaying quotient/Jacobian/source-name analyses.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
