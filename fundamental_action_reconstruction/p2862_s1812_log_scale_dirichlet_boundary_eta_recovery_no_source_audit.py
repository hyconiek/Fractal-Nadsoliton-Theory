#!/usr/bin/env python3
"""P2862/S1812: log-scale Dirichlet boundary eta-recovery no-source audit.

P2861 showed that the log-scale affine/harmonic variational equation is
eta-blind unless an extra eta-selecting premise is supplied.  This packet tests
the narrow next candidate: add Dirichlet boundary data at d=1 and d=11 and ask
whether the exponent is selected.

Result: yes, a log-scale affine solution with endpoint data recovers a unique
slope eta.  But when the endpoint value is the strict tail value T(11)=11^(9/5),
the boundary datum already contains eta=9/5.  Varying the endpoint value yields
other eta values.  Therefore Dirichlet recovery is a conditional inverse fit,
not a pre-profile source law, coupling/localization theorem, or closure.
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

P2861 = GEN / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.json"
OUT = GEN / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.json"
MD = GEN / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

STRICT_BETA = Fraction(1, 1)
STRICT_ETA = Fraction(9, 5)
LEFT_D = 1
RIGHT_D = 11
DOMAIN = list(range(1, 12))
ALTERNATIVE_ETAS = [Fraction(1, 1), Fraction(6, 5), Fraction(7, 5), Fraction(2, 1), Fraction(12, 5), Fraction(3, 1)]
TOL = 1e-12


def strict_tail(d: int, beta: Fraction = STRICT_BETA, eta: Fraction = STRICT_ETA) -> float:
    return float(beta) * float(d) ** float(eta)


def endpoint_log_value(beta: Fraction, eta: Fraction, d: int = RIGHT_D) -> float:
    return math.log(float(beta)) + float(eta) * math.log(d)


def recover_eta_from_dirichlet(left_log_tail: float, right_log_tail: float, left_d: int = LEFT_D, right_d: int = RIGHT_D) -> float:
    return (right_log_tail - left_log_tail) / (math.log(right_d) - math.log(left_d))


def recover_beta_from_left(left_log_tail: float) -> float:
    # left_d=1, so T(1)=beta and log(T(1))=log(beta).
    return math.exp(left_log_tail)


def dirichlet_solution_rows(left_log_tail: float, right_log_tail: float) -> list[dict[str, Any]]:
    eta = recover_eta_from_dirichlet(left_log_tail, right_log_tail)
    beta = recover_beta_from_left(left_log_tail)
    rows = []
    for d in DOMAIN:
        x = math.log(d)
        y = left_log_tail + eta * (x - math.log(LEFT_D))
        tail_value = math.exp(y)
        rows.append(
            {
                "d": d,
                "log_tail": y,
                "tail_value": tail_value,
                "recovered_beta_d_power_eta": beta * (d**eta),
                "abs_reconstruction_residual": abs(tail_value - beta * (d**eta)),
            }
        )
    return rows


def alternative_boundary_rows() -> list[dict[str, Any]]:
    rows = []
    left_log_tail = endpoint_log_value(STRICT_BETA, Fraction(0, 1), LEFT_D)
    for eta in ALTERNATIVE_ETAS:
        right_log_tail = endpoint_log_value(STRICT_BETA, eta, RIGHT_D)
        recovered_eta = recover_eta_from_dirichlet(left_log_tail, right_log_tail)
        rows.append(
            {
                "imposed_endpoint_eta": fraction_payload(eta),
                "right_endpoint_log_tail": right_log_tail,
                "right_endpoint_tail_value": math.exp(right_log_tail),
                "recovered_eta": recovered_eta,
                "abs_eta_error": abs(recovered_eta - float(eta)),
                "is_strict_eta": eta == STRICT_ETA,
                "requires_prime5": 5 in prime_factors(eta.denominator),
            }
        )
    return rows


def build_payload(p2861: dict[str, Any]) -> dict[str, Any]:
    left_log_tail = endpoint_log_value(STRICT_BETA, Fraction(0, 1), LEFT_D)
    strict_right_log_tail = endpoint_log_value(STRICT_BETA, STRICT_ETA, RIGHT_D)
    recovered_eta = recover_eta_from_dirichlet(left_log_tail, strict_right_log_tail)
    recovered_beta = recover_beta_from_left(left_log_tail)
    rows = dirichlet_solution_rows(left_log_tail, strict_right_log_tail)
    alternative_rows = alternative_boundary_rows()
    non_strict_boundary_count = sum(1 for row in alternative_rows if not row["is_strict_eta"] and row["abs_eta_error"] < TOL)
    max_reconstruction_residual = max(row["abs_reconstruction_residual"] for row in rows)
    facts = {
        "p2861_rechecked": p2861.get("status") == "P2861_LOG_SCALE_AFFINE_VARIATIONAL_ETA_BLIND_NO_SOURCE_AUDIT_NO_CLOSURE",
        "strict_dirichlet_data_recovers_eta": abs(recovered_eta - float(STRICT_ETA)) < TOL,
        "strict_dirichlet_data_recovers_beta": abs(recovered_beta - float(STRICT_BETA)) < TOL,
        "dirichlet_solution_reconstructs_supplied_profile": max_reconstruction_residual < TOL,
        "alternative_endpoint_data_select_other_eta_values": non_strict_boundary_count > 0,
        "accepted_count_zero": True,
    }
    return {
        "status": "P2862_LOG_SCALE_DIRICHLET_BOUNDARY_ETA_RECOVERY_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2861": sha(P2861)},
        "log_scale_dirichlet_boundary_eta_recovery_no_source_audit": {
            "input_status_rechecked": p2861.get("status"),
            "candidate_premise": "log-scale affine variational equation plus Dirichlet boundary data at d=1 and d=11",
            "strict_beta": fraction_payload(STRICT_BETA),
            "strict_eta": fraction_payload(STRICT_ETA),
            "left_boundary": {"d": LEFT_D, "log_tail": left_log_tail, "tail_value": math.exp(left_log_tail)},
            "right_boundary": {"d": RIGHT_D, "log_tail": strict_right_log_tail, "tail_value": math.exp(strict_right_log_tail)},
            "recovered_eta_from_strict_boundaries": recovered_eta,
            "recovered_beta_from_left_boundary": recovered_beta,
            "max_reconstruction_residual": max_reconstruction_residual,
            "dirichlet_solution_rows": rows,
            "alternative_endpoint_rows": alternative_rows,
            "non_strict_boundary_eta_count": non_strict_boundary_count,
            "candidate_matrix": [
                {
                    "candidate": "strict_endpoint_dirichlet_eta_recovery",
                    "finite_witness_passes": abs(recovered_eta - float(STRICT_ETA)) < TOL,
                    "exports_eta_source_law": False,
                    "exports_boundary_source_law": False,
                    "verdict": "strict endpoint data recovers eta=9/5, but the right endpoint value already contains the strict exponent.",
                },
                {
                    "candidate": "alternative_endpoint_counterfamily",
                    "finite_witness_passes": non_strict_boundary_count > 0,
                    "exports_eta_source_law": False,
                    "exports_boundary_source_law": False,
                    "verdict": "changing only the endpoint tail value selects different eta values with the same log-affine equation.",
                },
                {
                    "candidate": "unit_boundary_beta_recovery",
                    "finite_witness_passes": abs(recovered_beta - float(STRICT_BETA)) < TOL,
                    "exports_eta_source_law": False,
                    "exports_boundary_source_law": False,
                    "verdict": "beta=1 is recovered from the imposed T(1)=1 boundary, not sourced as a target-independent unit law.",
                },
            ],
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "recovery_step": "For log-affine y=A+eta*x, two Dirichlet endpoint values determine A and eta uniquely.",
                "circularity_step": "Choosing the right endpoint as log(11^(9/5)) imports eta=9/5 in the boundary datum.",
                "counterfamily_step": "Alternative endpoint values log(11^eta') recover eta' just as exactly, so the boundary condition must itself be sourced before it can select eta.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_log_scale_dirichlet_boundary_eta_recovery_no_source_audit": all(facts.values()),
            "exports_eta_source_law": False,
            "exports_boundary_source_law": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "eta_source_exported": False,
                "boundary_source_law_exported": False,
                "target_independent_beta_unit_source_exported": False,
                "unit_bearing_coupling_localization_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2862 shows that Dirichlet boundary data can recover eta only conditionally.  The strict right endpoint log(11^(9/5)) already imports eta=9/5, and alternative endpoint data recover alternative eta values.  Thus boundary recovery is inverse fitting, not an eta source law or unit-bearing coupling/localization theorem.",
            "next_honest_step": "Do not replay endpoint fitting, strict tail boundary values, or Dirichlet recovery as eta sourcehood.  The next proof-grade move must source the boundary datum itself via a unit-bearing coupling/localization theorem, introduce a genuine exponent-fixing source term, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["log_scale_dirichlet_boundary_eta_recovery_no_source_audit"]
    lines = [
        "# P2862/S1812 log-scale Dirichlet boundary eta-recovery no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Dirichlet recovery witness",
        f"- recovered eta from strict endpoints: `{audit['recovered_eta_from_strict_boundaries']}`",
        f"- recovered beta from left endpoint: `{audit['recovered_beta_from_left_boundary']}`",
        f"- max reconstruction residual: `{audit['max_reconstruction_residual']}`",
        f"- non-strict boundary eta count: `{audit['non_strict_boundary_eta_count']}`",
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
    payload = build_payload(read_json(P2861))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2862/S1812 log-scale Dirichlet boundary eta-recovery no-source audit",
        "## P2862/S1812 log-scale Dirichlet boundary eta-recovery no-source audit\n\n"
        "`P2862/S1812` adds Dirichlet endpoint data to the log-scale affine variational equation after `P2861`.  The strict endpoints recover `eta=9/5` and `beta=1`, but this is conditional inverse recovery: the right endpoint value `log(11^(9/5))` already contains the strict exponent, and alternative endpoint values recover alternative exponents.  Endpoint fitting does not source the boundary datum, an eta law, a target-independent beta unit, a strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2862/S1812 Dirichlet boundary recovery `L_total` guard",
        "## P2862/S1812 Dirichlet boundary recovery `L_total` guard\n\n"
        "`P2862/S1812` adds no strict action term.  Dirichlet endpoint recovery fixes the log-affine fit only after importing endpoint tail values; it is not a unit-bearing boundary/source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current log-scale Dirichlet boundary eta-recovery guardrail (P2862/S1812, 2026-06-18)",
        "## Current log-scale Dirichlet boundary eta-recovery guardrail (P2862/S1812, 2026-06-18)\n\n"
        "- P2862 tests Dirichlet endpoint data as the next eta-selecting candidate after P2861.\n"
        "- The strict endpoints recover `eta=9/5` and `beta=1`, but the right endpoint already imports `log(11^(9/5))`; alternative endpoint data recover alternative eta values.\n"
        "- Do not promote endpoint fitting, strict tail boundary values, Dirichlet recovery, log-scale harmonicity, or multiplicative scale covariance to strict damping/compression bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must source the boundary datum itself through a unit-bearing coupling/localization theorem, introduce a genuine exponent-fixing source term, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
