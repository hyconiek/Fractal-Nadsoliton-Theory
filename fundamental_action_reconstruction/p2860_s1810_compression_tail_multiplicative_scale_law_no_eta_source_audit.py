#!/usr/bin/env python3
"""P2860/S1810: compression-tail multiplicative scale-law no-eta-source audit.

P2859 showed that the supplied strict compression profile identifies beta and
eta, but inverse identifiability is not sourcehood.  This packet tests a more
structural pre-profile candidate: the strict compression tail
T(d)=C(d)-1=beta*d^eta might be sourced by multiplicative scale covariance
T(ab)=T(a)T(b) on the finite positive Z12 distance semigroup where ab <= 11.

Result: the strict tail passes this scale-law exactly, and the law constrains
nonzero beta to beta=1.  But it leaves eta free: every eta in a tested rational
continuum also satisfies the same law with beta=1.  Therefore multiplicative
scale covariance is a real structural witness, not an eta source law.
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

P2859 = GEN / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.json"
OUT = GEN / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.json"
MD = GEN / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

STRICT_BETA = Fraction(1, 1)
STRICT_ETA = Fraction(9, 5)
DOMAIN_MAX = 11
SAMPLE_ETAS = [Fraction(n, 5) for n in range(5, 16)]  # 1, 6/5, ..., 3
TOL = 1e-12


def product_pairs(limit: int = DOMAIN_MAX) -> list[tuple[int, int, int]]:
    return [(a, b, a * b) for a in range(1, limit + 1) for b in range(1, limit + 1) if a * b <= limit]


def tail(d: int, beta: Fraction = STRICT_BETA, eta: Fraction = STRICT_ETA) -> float:
    return float(beta) * (float(d) ** float(eta))


def multiplicative_residual_rows(beta: Fraction = STRICT_BETA, eta: Fraction = STRICT_ETA) -> list[dict[str, Any]]:
    rows = []
    for a, b, ab in product_pairs():
        left = tail(ab, beta, eta)
        right = tail(a, beta, eta) * tail(b, beta, eta)
        residual = left - right
        rows.append(
            {
                "a": a,
                "b": b,
                "ab": ab,
                "left_T_ab": left,
                "right_T_a_T_b": right,
                "residual": residual,
                "abs_residual": abs(residual),
            }
        )
    return rows


def sample_eta_family() -> list[dict[str, Any]]:
    rows = []
    for eta in SAMPLE_ETAS:
        residuals = multiplicative_residual_rows(STRICT_BETA, eta)
        max_abs = max(row["abs_residual"] for row in residuals)
        rows.append(
            {
                "eta": fraction_payload(eta),
                "requires_prime5": 5 in prime_factors(eta.denominator),
                "max_abs_multiplicative_residual": max_abs,
                "passes_scale_law": max_abs < TOL,
                "is_strict_eta": eta == STRICT_ETA,
            }
        )
    return rows


def beta_constraint_rows() -> list[dict[str, Any]]:
    rows = []
    for beta in [Fraction(0, 1), Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(3, 1)]:
        residuals = multiplicative_residual_rows(beta, STRICT_ETA)
        nontrivial = beta != 0
        rows.append(
            {
                "beta": fraction_payload(beta),
                "nontrivial_tail": nontrivial,
                "max_abs_multiplicative_residual": max(row["abs_residual"] for row in residuals),
                "passes_scale_law": max(row["abs_residual"] for row in residuals) < TOL,
            }
        )
    return rows


def build_payload(p2859: dict[str, Any]) -> dict[str, Any]:
    strict_rows = multiplicative_residual_rows()
    max_abs = max(row["abs_residual"] for row in strict_rows)
    eta_family = sample_eta_family()
    beta_rows = beta_constraint_rows()
    non_strict_eta_pass_count = sum(1 for row in eta_family if row["passes_scale_law"] and not row["is_strict_eta"])
    beta_one_unique_nonzero_sample = [row for row in beta_rows if row["nontrivial_tail"] and row["passes_scale_law"]]
    facts = {
        "p2859_rechecked": p2859.get("status") == "P2859_ETA_BETA_PROFILE_IDENTIFIABILITY_NO_SOURCE_AUDIT_NO_CLOSURE",
        "strict_tail_passes_multiplicative_scale_law": max_abs < TOL,
        "beta_one_is_unique_nonzero_sample_solution": len(beta_one_unique_nonzero_sample) == 1 and beta_one_unique_nonzero_sample[0]["beta"]["fraction"] == "1/1",
        "non_strict_eta_samples_also_pass": non_strict_eta_pass_count > 0,
        "accepted_count_zero": True,
    }
    return {
        "status": "P2860_COMPRESSION_TAIL_MULTIPLICATIVE_SCALE_LAW_NO_ETA_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2859": sha(P2859)},
        "compression_tail_multiplicative_scale_law_no_eta_source_audit": {
            "input_status_rechecked": p2859.get("status"),
            "scale_law_candidate": "T(ab)=T(a)T(b) for T(d)=C(d)-1 on positive Z12 distances with ab<=11",
            "strict_beta": fraction_payload(STRICT_BETA),
            "strict_eta": fraction_payload(STRICT_ETA),
            "product_pair_count": len(product_pairs()),
            "strict_scale_law_max_abs_residual": max_abs,
            "strict_scale_law_rows_first16": strict_rows[:16],
            "eta_family_samples": eta_family,
            "non_strict_eta_pass_count": non_strict_eta_pass_count,
            "beta_constraint_samples": beta_rows,
            "candidate_matrix": [
                {
                    "candidate": "multiplicative_tail_scale_law",
                    "finite_witness_passes": max_abs < TOL,
                    "exports_eta_source_law": False,
                    "exports_beta_source_law": False,
                    "verdict": "exact scale covariance is real for the supplied strict tail, but it is shared by a continuum of eta values.",
                },
                {
                    "candidate": "nonzero_beta_multiplicativity_constraint",
                    "finite_witness_passes": facts["beta_one_is_unique_nonzero_sample_solution"],
                    "exports_eta_source_law": False,
                    "exports_beta_source_law": False,
                    "verdict": "the law constrains nonzero multiplicative tails toward beta=1, but does not export target-independent units or coupling sourcehood.",
                },
                {
                    "candidate": "sampled_eta_continuum_counterexample",
                    "finite_witness_passes": non_strict_eta_pass_count > 0,
                    "exports_eta_source_law": False,
                    "exports_beta_source_law": False,
                    "verdict": "many non-strict eta samples satisfy the same scale law, so eta=9/5 is not selected.",
                },
            ],
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "scale_step": "For beta=1, T(d)=d^eta obeys T(ab)=T(a)T(b) for every eta and every audited product ab<=11.",
                "beta_step": "For T(d)=beta*d^eta, nonzero multiplicativity forces beta=1 algebraically because beta*(ab)^eta = beta^2*(ab)^eta.",
                "eta_boundary_step": "The same law is eta-blind; it cannot select eta=9/5 without an additional source principle.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_multiplicative_scale_law_no_eta_source_audit": all(facts.values()),
            "exports_eta_source_law": False,
            "exports_target_independent_beta_unit_source": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "eta_source_exported": False,
                "target_independent_beta_unit_source_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2860 finds a real structural witness: the compression tail is multiplicative and this algebraically singles out beta=1 among nonzero tails.  However, the same law is eta-blind: every sampled eta in the rational family also satisfies it.  Therefore scale covariance does not source eta=9/5 or export units, coupling, L_total, EOM, Hamiltonian, bridge, or ToE closure.",
            "next_honest_step": "Do not replay multiplicative scale covariance as eta sourcehood.  The next proof-grade move must add an eta-selecting strict principle, such as a unit-bearing variational/coupling/localization theorem with a real exponent-fixing equation; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["compression_tail_multiplicative_scale_law_no_eta_source_audit"]
    lines = [
        "# P2860/S1810 compression-tail multiplicative scale-law no-eta-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Scale-law witness",
        f"- product pair count: `{audit['product_pair_count']}`",
        f"- strict max residual: `{audit['strict_scale_law_max_abs_residual']}`",
        f"- non-strict eta pass count: `{audit['non_strict_eta_pass_count']}`",
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
    payload = build_payload(read_json(P2859))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2860/S1810 compression-tail multiplicative scale-law no-eta-source audit",
        "## P2860/S1810 compression-tail multiplicative scale-law no-eta-source audit\n\n"
        "`P2860/S1810` tests a structural pre-profile candidate for damping/compression: multiplicative tail scale covariance `T(ab)=T(a)T(b)` for `T(d)=C(d)-1`.  The strict tail passes exactly on all audited positive `Z12` product pairs with `ab<=11`, and nonzero multiplicativity algebraically constrains `beta=1`.  But the law is eta-blind: sampled non-strict rational exponents satisfy the same scale law, so `eta=9/5` is not selected.  No eta source, target-independent beta/unit source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2860/S1810 multiplicative scale-law `L_total` guard",
        "## P2860/S1810 multiplicative scale-law `L_total` guard\n\n"
        "`P2860/S1810` adds no action term.  Multiplicative tail scale covariance constrains a supplied tail shape but does not provide an eta-selecting unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current compression-tail multiplicative scale-law guardrail (P2860/S1810, 2026-06-18)",
        "## Current compression-tail multiplicative scale-law guardrail (P2860/S1810, 2026-06-18)\n\n"
        "- P2860 verifies that the strict compression tail satisfies multiplicative scale covariance and that nonzero multiplicativity constrains `beta=1`.\n"
        "- The same scale law is eta-blind: non-strict rational eta samples satisfy it, so `eta=9/5` is not selected or sourced.\n"
        "- Do not promote multiplicative scale covariance, beta normalization, profile identifiability, or eta denominator representability to strict damping/compression bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply an eta-selecting strict principle, such as a unit-bearing variational/coupling/localization theorem with an exponent-fixing equation; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
