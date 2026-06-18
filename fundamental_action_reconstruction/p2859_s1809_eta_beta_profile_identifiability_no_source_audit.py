#!/usr/bin/env python3
"""P2859/S1809: eta/beta compression-profile identifiability no-source audit.

P2858 closed the current phase-bit/source replay lane by proving an open phase
cell.  This packet pivots to the other recommended atom: strict damping /
compression beta-eta sourcing.  It asks a narrow question: if the strict
compression profile C(d)=1+beta*d^eta is already given on the finite positive
Z12 distances, does it identify beta=1 and eta=9/5?

Result: yes, the profile is locally and two-point identifiable.  But inverse
identifiability is not a source theorem: it recovers parameters from an already
supplied strict profile and does not export a pre-profile law producing beta,
eta, units, localization, L_total, EOM, Hamiltonian, bridge closure, or ToE.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN, fraction_payload, prime_factors

P2858 = GEN / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.json"
OUT = GEN / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.json"
MD = GEN / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

STRICT_BETA = Fraction(1, 1)
STRICT_ETA = Fraction(9, 5)
POSITIVE_DOMAIN = tuple(range(1, 12))
IDENTIFY_PAIR = (2, 3)


def compression_value(d: int, beta: Fraction = STRICT_BETA, eta: Fraction = STRICT_ETA) -> float:
    return 1.0 + float(beta) * (float(d) ** float(eta))


def compression_rows() -> list[dict[str, Any]]:
    return [
        {
            "d": d,
            "compression_value": compression_value(d),
            "tail_value": compression_value(d) - 1.0,
            "log_tail": math.log(compression_value(d) - 1.0),
            "log_d": math.log(d),
        }
        for d in POSITIVE_DOMAIN
    ]


def two_point_inverse(d1: int, d2: int) -> dict[str, Any]:
    y1 = compression_value(d1) - 1.0
    y2 = compression_value(d2) - 1.0
    eta_hat = math.log(y2 / y1) / math.log(d2 / d1)
    beta_hat = y1 / (d1**eta_hat)
    return {
        "d_pair": [d1, d2],
        "eta_hat": eta_hat,
        "beta_hat": beta_hat,
        "eta_abs_error": abs(eta_hat - float(STRICT_ETA)),
        "beta_abs_error": abs(beta_hat - float(STRICT_BETA)),
    }


def jacobian_certificate(d1: int, d2: int) -> dict[str, Any]:
    beta = float(STRICT_BETA)
    eta = float(STRICT_ETA)
    row1 = [d1**eta, beta * d1**eta * math.log(d1)]
    row2 = [d2**eta, beta * d2**eta * math.log(d2)]
    determinant = row1[0] * row2[1] - row1[1] * row2[0]
    closed_form = beta * (d1**eta) * (d2**eta) * math.log(d2 / d1)
    return {
        "d_pair": [d1, d2],
        "jacobian_rows_beta_eta": [row1, row2],
        "determinant": determinant,
        "closed_form_determinant": closed_form,
        "nonzero_determinant": abs(determinant) > 1e-12,
        "local_identifiability": abs(determinant) > 1e-12,
    }


def eta_denominator_payload() -> dict[str, Any]:
    factors = prime_factors(STRICT_ETA.denominator)
    return {
        "eta": fraction_payload(STRICT_ETA),
        "denominator_prime_factors": factors,
        "requires_prime5": 5 in factors,
        "pure_z12_denominator_source": set(factors).issubset({2, 3}),
    }


def candidate_matrix(inverse: dict[str, Any], jacobian: dict[str, Any], eta_payload: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "two_point_profile_inverse",
            "finite_witness_passes": inverse["eta_abs_error"] < 1e-12 and inverse["beta_abs_error"] < 1e-12,
            "exports_pre_profile_source_law": False,
            "verdict": "identifies beta/eta from an already supplied strict compression profile; inverse fit is not a source law.",
        },
        {
            "candidate": "jacobian_local_identifiability",
            "finite_witness_passes": jacobian["local_identifiability"],
            "exports_pre_profile_source_law": False,
            "verdict": "nonzero Jacobian proves local identifiability, not origin/source of the profile.",
        },
        {
            "candidate": "pure_z12_eta_denominator_source",
            "finite_witness_passes": not eta_payload["pure_z12_denominator_source"] and eta_payload["requires_prime5"],
            "exports_pre_profile_source_law": False,
            "verdict": "eta=9/5 again imports prime 5; denominator representation is not a source theorem.",
        },
        {
            "candidate": "beta_one_normalization",
            "finite_witness_passes": STRICT_BETA == 1,
            "exports_pre_profile_source_law": False,
            "verdict": "beta=1 is a normalization in the supplied strict profile unless a target-independent unit/coupling law exports it.",
        },
    ]


def build_payload(p2858: dict[str, Any]) -> dict[str, Any]:
    rows = compression_rows()
    inverse = two_point_inverse(*IDENTIFY_PAIR)
    jacobian = jacobian_certificate(*IDENTIFY_PAIR)
    eta_payload = eta_denominator_payload()
    matrix = candidate_matrix(inverse, jacobian, eta_payload)
    accepted_count = sum(1 for row in matrix if row["exports_pre_profile_source_law"])
    facts = {
        "p2858_rechecked": p2858.get("status") == "P2858_PHASE_BIT_CELL_CONTINUUM_NO_SOURCE_AUDIT_NO_CLOSURE",
        "two_point_inverse_recovers_eta": inverse["eta_abs_error"] < 1e-12,
        "two_point_inverse_recovers_beta": inverse["beta_abs_error"] < 1e-12,
        "jacobian_is_nonzero": jacobian["nonzero_determinant"],
        "eta_requires_prime5": eta_payload["requires_prime5"],
        "accepted_count_zero": accepted_count == 0,
    }
    return {
        "status": "P2859_ETA_BETA_PROFILE_IDENTIFIABILITY_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2858": sha(P2858)},
        "eta_beta_profile_identifiability_no_source_audit": {
            "input_status_rechecked": p2858.get("status"),
            "strict_beta": fraction_payload(STRICT_BETA),
            "strict_eta": eta_payload,
            "compression_profile": rows,
            "two_point_inverse": inverse,
            "jacobian_certificate": jacobian,
            "candidate_matrix": matrix,
            "accepted_candidate_count": accepted_count,
            "proof_certificate": {
                "inverse_step": "The positive strict compression profile C(d)=1+beta*d^eta determines beta and eta from d=2,3 when the profile is already supplied.",
                "jacobian_step": "The beta/eta Jacobian on d=2,3 has determinant beta*2^eta*3^eta*log(3/2), hence is nonzero at beta=1, eta=9/5.",
                "source_boundary_step": "Identifiability of parameters from a supplied profile is weaker than a pre-profile source law for beta, eta, units, coupling, or dynamics.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_eta_beta_identifiability_no_source_audit": all(facts.values()),
            "exports_eta_beta_source_law": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "eta_beta_source_law_exported": False,
                "target_independent_beta_source_exported": False,
                "eta_source_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2859 shows that beta=1 and eta=9/5 are identifiable from an already supplied strict compression profile: a two-point inverse recovers them and the beta/eta Jacobian is nonzero.  This is not a source law; it presupposes the strict profile and does not export a target-independent beta source, eta source, unit/coupling law, L_total term, EOM, Hamiltonian, bridge, or ToE closure.",
            "next_honest_step": "Do not replay profile fitting or local identifiability as eta/beta sourcehood.  The next proof-grade move must supply a pre-profile strict source law for the compression profile or a unit-bearing coupling/localization theorem; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["eta_beta_profile_identifiability_no_source_audit"]
    inv = audit["two_point_inverse"]
    jac = audit["jacobian_certificate"]
    lines = [
        "# P2859/S1809 eta/beta profile identifiability no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Identifiability witnesses",
        f"- two-point pair: `{inv['d_pair']}`",
        f"- eta_hat: `{inv['eta_hat']}`; eta_abs_error: `{inv['eta_abs_error']}`",
        f"- beta_hat: `{inv['beta_hat']}`; beta_abs_error: `{inv['beta_abs_error']}`",
        f"- Jacobian determinant: `{jac['determinant']}`; nonzero: `{jac['nonzero_determinant']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2858))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2859/S1809 eta/beta profile identifiability no-source audit",
        "## P2859/S1809 eta/beta profile identifiability no-source audit\n\n"
        "`P2859/S1809` pivots after the P2858 phase-bit open-cell no-source result to strict damping/compression.  The supplied strict profile `C(d)=1+beta*d^eta` is locally identifiable: the `d=2,3` inverse recovers `beta=1`, `eta=9/5`, and the beta/eta Jacobian determinant is nonzero.  This proves inverse identifiability from an already supplied profile, not a pre-profile source law for `eta`, target-independent `beta`, units, coupling, localization, `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2859/S1809 eta/beta identifiability `L_total` guard",
        "## P2859/S1809 eta/beta identifiability `L_total` guard\n\n"
        "`P2859/S1809` adds no action term.  Recovering `beta=1` and `eta=9/5` from an already supplied compression profile does not provide a pre-profile unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current eta/beta profile-identifiability no-source guardrail (P2859/S1809, 2026-06-18)",
        "## Current eta/beta profile-identifiability no-source guardrail (P2859/S1809, 2026-06-18)\n\n"
        "- P2859 shows that the supplied strict compression profile locally identifies `beta=1` and `eta=9/5`; a two-point inverse recovers them and the beta/eta Jacobian is nonzero.\n"
        "- This is inverse identifiability, not a pre-profile source law; `eta=9/5` still imports prime `5`, and `beta=1` remains a supplied normalization absent a target-independent unit/coupling theorem.\n"
        "- Do not promote profile fitting, local identifiability, beta normalization, or eta denominator representability to strict damping/compression bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a pre-profile strict compression source law or a unit-bearing coupling/localization theorem; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
