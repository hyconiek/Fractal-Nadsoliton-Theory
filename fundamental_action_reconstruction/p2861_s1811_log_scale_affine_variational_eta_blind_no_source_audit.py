#!/usr/bin/env python3
"""P2861/S1811: log-scale affine variational eta-blind no-source audit.

P2860 found that multiplicative tail covariance is real but eta-blind.  This
packet tests the natural next eta-selecting candidate: promote the tail to a
log-scale variational object and ask whether a scale-affine/harmonic action
selects the strict exponent eta=9/5.

Let T(d)=beta*d^eta and y(d)=log(T(d)) on positive Z12 distances.  In the
scale coordinate x=log(d), y=log(beta)+eta*x is exactly affine.  Therefore the
finite scale-affine residuals and the associated squared residual action vanish
for every eta and every positive beta.  The variational equation is a real
structural witness, but it is eta-blind and beta-scale-blind; it cannot source
eta=9/5 or a target-independent beta unit.
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

P2860 = GEN / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.json"
OUT = GEN / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.json"
MD = GEN / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

STRICT_BETA = Fraction(1, 1)
STRICT_ETA = Fraction(9, 5)
DOMAIN = list(range(1, 12))
SAMPLE_ETAS = [Fraction(n, 5) for n in range(5, 16)]
SAMPLE_BETAS = [Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(3, 1)]
TOL = 1e-12


def log_tail(d: int, beta: Fraction, eta: Fraction) -> float:
    return math.log(float(beta)) + float(eta) * math.log(d)


def scale_affine_residual_rows(beta: Fraction = STRICT_BETA, eta: Fraction = STRICT_ETA) -> list[dict[str, Any]]:
    """Residuals of affine interpolation in x=log(d) across all triples a<b<c."""
    rows = []
    for i, a in enumerate(DOMAIN):
        for j, b in enumerate(DOMAIN[i + 1 :], start=i + 1):
            for c in DOMAIN[j + 1 :]:
                xa, xb, xc = math.log(a), math.log(b), math.log(c)
                ya, yb, yc = log_tail(a, beta, eta), log_tail(b, beta, eta), log_tail(c, beta, eta)
                # y_b should equal affine interpolation between (x_a,y_a) and (x_c,y_c).
                weight = (xb - xa) / (xc - xa)
                interpolated = ya + weight * (yc - ya)
                residual = yb - interpolated
                rows.append(
                    {
                        "a": a,
                        "b": b,
                        "c": c,
                        "residual": residual,
                        "abs_residual": abs(residual),
                    }
                )
    return rows


def action_value(beta: Fraction, eta: Fraction) -> float:
    return sum(row["residual"] ** 2 for row in scale_affine_residual_rows(beta, eta))


def eta_action_samples() -> list[dict[str, Any]]:
    rows = []
    for eta in SAMPLE_ETAS:
        action = action_value(STRICT_BETA, eta)
        rows.append(
            {
                "eta": fraction_payload(eta),
                "requires_prime5": 5 in prime_factors(eta.denominator),
                "scale_affine_action": action,
                "passes_variational_equation": action < TOL,
                "is_strict_eta": eta == STRICT_ETA,
            }
        )
    return rows


def beta_action_samples() -> list[dict[str, Any]]:
    rows = []
    for beta in SAMPLE_BETAS:
        action = action_value(beta, STRICT_ETA)
        rows.append(
            {
                "beta": fraction_payload(beta),
                "scale_affine_action": action,
                "passes_variational_equation": action < TOL,
                "is_strict_beta": beta == STRICT_BETA,
            }
        )
    return rows


def build_payload(p2860: dict[str, Any]) -> dict[str, Any]:
    strict_rows = scale_affine_residual_rows()
    strict_action = sum(row["residual"] ** 2 for row in strict_rows)
    eta_rows = eta_action_samples()
    beta_rows = beta_action_samples()
    non_strict_eta_pass_count = sum(1 for row in eta_rows if row["passes_variational_equation"] and not row["is_strict_eta"])
    non_strict_beta_pass_count = sum(1 for row in beta_rows if row["passes_variational_equation"] and not row["is_strict_beta"])
    facts = {
        "p2860_rechecked": p2860.get("status") == "P2860_COMPRESSION_TAIL_MULTIPLICATIVE_SCALE_LAW_NO_ETA_SOURCE_AUDIT_NO_CLOSURE",
        "strict_tail_satisfies_log_scale_affine_variational_equation": strict_action < TOL,
        "non_strict_eta_samples_also_satisfy": non_strict_eta_pass_count > 0,
        "non_strict_beta_samples_also_satisfy": non_strict_beta_pass_count > 0,
        "accepted_count_zero": True,
    }
    return {
        "status": "P2861_LOG_SCALE_AFFINE_VARIATIONAL_ETA_BLIND_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2860": sha(P2860)},
        "log_scale_affine_variational_eta_blind_no_source_audit": {
            "input_status_rechecked": p2860.get("status"),
            "candidate_variational_principle": "minimize sum of squared affine-interpolation residuals of y=log(T(d)) in x=log(d)",
            "strict_beta": fraction_payload(STRICT_BETA),
            "strict_eta": fraction_payload(STRICT_ETA),
            "triple_count": len(strict_rows),
            "strict_scale_affine_action": strict_action,
            "strict_residual_rows_first16": strict_rows[:16],
            "eta_action_samples": eta_rows,
            "beta_action_samples": beta_rows,
            "non_strict_eta_pass_count": non_strict_eta_pass_count,
            "non_strict_beta_pass_count": non_strict_beta_pass_count,
            "candidate_matrix": [
                {
                    "candidate": "log_scale_affine_variational_equation",
                    "finite_witness_passes": strict_action < TOL,
                    "exports_eta_source_law": False,
                    "exports_beta_unit_source_law": False,
                    "verdict": "the strict tail is stationary/zero-action because log(T)=log(beta)+eta*log(d) is scale-affine.",
                },
                {
                    "candidate": "eta_selection_by_scale_affine_action",
                    "finite_witness_passes": non_strict_eta_pass_count > 0,
                    "exports_eta_source_law": False,
                    "exports_beta_unit_source_law": False,
                    "verdict": "non-strict eta samples have the same zero action, so eta=9/5 is not selected.",
                },
                {
                    "candidate": "beta_unit_selection_by_log_scale_affine_action",
                    "finite_witness_passes": non_strict_beta_pass_count > 0,
                    "exports_eta_source_law": False,
                    "exports_beta_unit_source_law": False,
                    "verdict": "positive beta only shifts log(T) by a constant, so the action does not fix beta=1 as a unit-bearing source.",
                },
            ],
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "affine_step": "For T(d)=beta*d^eta, y=log(T(d))=log(beta)+eta*log(d) is affine in x=log(d).",
                "variational_step": "Every affine-interpolation residual across triples a<b<c vanishes, so the squared residual action is zero.",
                "no_source_step": "The zero-action condition holds for all sampled eta and all sampled positive beta; it is a structural consistency equation, not an eta- or unit-selecting source law.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_log_scale_affine_variational_eta_blind_no_source_audit": all(facts.values()),
            "exports_eta_source_law": False,
            "exports_target_independent_beta_unit_source": False,
            "exports_unit_bearing_variational_source_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "eta_source_exported": False,
                "target_independent_beta_unit_source_exported": False,
                "unit_bearing_variational_source_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2861 tests the natural log-scale variational upgrade after P2860.  The strict tail has zero scale-affine action, but so do many non-strict eta values and many positive beta normalizations.  Therefore the variational equation is eta- and unit-blind, not a source law for eta=9/5 or beta=1.",
            "next_honest_step": "Do not replay log-scale affine/harmonic variational language as eta sourcehood.  A next proof-grade move must introduce a genuinely eta-selecting source term, boundary condition, unit-bearing coupling/localization theorem, or a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["log_scale_affine_variational_eta_blind_no_source_audit"]
    lines = [
        "# P2861/S1811 log-scale affine variational eta-blind no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Variational witness",
        f"- triple count: `{audit['triple_count']}`",
        f"- strict scale-affine action: `{audit['strict_scale_affine_action']}`",
        f"- non-strict eta pass count: `{audit['non_strict_eta_pass_count']}`",
        f"- non-strict beta pass count: `{audit['non_strict_beta_pass_count']}`",
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
    payload = build_payload(read_json(P2860))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2861/S1811 log-scale affine variational eta-blind no-source audit",
        "## P2861/S1811 log-scale affine variational eta-blind no-source audit\n\n"
        "`P2861/S1811` tests the natural variational upgrade after `P2860`: require `y=log(T(d))` to be affine/harmonic in the scale coordinate `x=log(d)`, equivalently minimize squared affine-interpolation residuals across positive `Z12` distance triples.  The strict tail has zero action, but so do sampled non-strict `eta` values and sampled positive `beta` normalizations.  Thus log-scale variational consistency is eta- and unit-blind; it does not source `eta=9/5`, target-independent `beta=1`, a strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2861/S1811 log-scale variational `L_total` guard",
        "## P2861/S1811 log-scale variational `L_total` guard\n\n"
        "`P2861/S1811` adds no action term to the strict theory.  A zero auxiliary log-scale affine residual action is a fitted consistency condition on the supplied tail, not a unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current log-scale affine variational eta-blind guardrail (P2861/S1811, 2026-06-18)",
        "## Current log-scale affine variational eta-blind guardrail (P2861/S1811, 2026-06-18)\n\n"
        "- P2861 tests a log-scale affine/harmonic variational candidate for the compression tail after P2860.\n"
        "- The strict tail has zero scale-affine action, but non-strict eta samples and positive beta rescalings also have zero action; the condition is eta- and unit-blind.\n"
        "- Do not promote log-scale harmonicity, affine residual minimization, beta rescaling invariance, multiplicative scale covariance, or profile identifiability to strict damping/compression bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must add a genuinely eta-selecting source term, boundary condition, unit-bearing coupling/localization theorem, or a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
