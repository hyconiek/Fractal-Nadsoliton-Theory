#!/usr/bin/env python3
"""Scratch analytic-domain audit for the legacy/strict kernel comparison.

Motivation:
- The strict kernel is numerically/gate useful, but its d^(9/5) damping changes
  the analytic domain relative to the legacy meromorphic kernel.
- The D_f-1 bridge-prep route is numerically strong, but D_f-1=4 ln 2 - 1 is
  not the same analytic class as eta=9/5.

This packet records a sharper theorem-prep obstruction/choice:
positive-real operational bridge vs complex analytic/branched-cover bridge.
It does not export a bridge theorem or role transfer.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLACEMENT_AUDIT = HERE / "bridge_df_transport_placement_discriminator_report.json"
F151 = ROOT / "fundamental_action_reconstruction" / "F151_FIRST_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_PACKET.md"
K1 = ROOT / "fundamental_action_reconstruction" / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
K2 = ROOT / "fundamental_action_reconstruction" / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md"
OUT_JSON = HERE / "bridge_kernel_analytic_domain_audit_report.json"
OUT_MD = HERE / "bridge_kernel_analytic_domain_audit_report.md"

ALPHA_GEO = 4.0 * math.log(2.0)
GAMMA_F = ALPHA_GEO - 1.0
STRICT_ETA = Fraction(9, 5)
LEGACY_BETA_TORS = Fraction(1, 100)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rational_approx_rows(x: float, denominators: list[int]) -> list[dict[str, Any]]:
    rows = []
    for max_den in denominators:
        best = Fraction(x).limit_denominator(max_den)
        strict_gap = abs(float(STRICT_ETA) - x)
        best_gap = abs(float(best) - x)
        rows.append(
            {
                "max_denominator": max_den,
                "best_fraction": f"{best.numerator}/{best.denominator}",
                "best_value": float(best),
                "best_abs_gap_to_D_f_minus_1": best_gap,
                "strict_9_over_5_abs_gap_to_D_f_minus_1": strict_gap,
                "strict_9_over_5_is_best_at_this_bound": best == STRICT_ETA,
            }
        )
    return rows


def continued_fraction(x: float, n_terms: int = 10) -> list[int]:
    terms = []
    value = x
    for _ in range(n_terms):
        term = math.floor(value)
        terms.append(term)
        frac = value - term
        if abs(frac) < 1e-15:
            break
        value = 1.0 / frac
    return terms


def symbolic_domain_facts() -> dict[str, Any]:
    z, t = sp.symbols("z t")
    legacy_den = 1 + sp.Rational(1, 100) * z
    strict_den_z = 1 + z ** sp.Rational(9, 5)
    strict_cover_den = sp.factor(1 + t**9)  # z=t^5
    strict_cover_roots = [sp.sstr(root) for root in sp.solve(sp.Eq(strict_cover_den, 0), t)]
    return {
        "legacy_denominator": sp.sstr(legacy_den),
        "legacy_pole": "z = -100",
        "strict_denominator_principal_symbolic": sp.sstr(strict_den_z),
        "strict_eta_exact": "9/5",
        "strict_branch_points_on_z_plane": ["z=0", "z=infinity"],
        "strict_finite_cover_substitution": "z = t^5",
        "strict_denominator_on_cover": sp.sstr(strict_cover_den),
        "strict_cover_pole_equation": "t^9 = -1",
        "strict_cover_roots_sample": strict_cover_roots,
        "D_f_minus_1_exact_expression": "4*log(2) - 1",
        "D_f_minus_1_domain_note": "By Hermite-Lindemann, log(2) is transcendental; hence 4*log(2)-1 is not a rational finite-sheet exponent.  A d^(D_f-1) bridge is naturally positive-real or infinite-branched, not a finite algebraic cover like eta=9/5.",
    }


def text_contains(path: Path, needles: list[str]) -> dict[str, bool]:
    text = path.read_text(encoding="utf-8") if path.exists() else ""
    return {needle: needle in text for needle in needles}


def main() -> None:
    placement = load_json(PLACEMENT_AUDIT)
    approximation_rows = rational_approx_rows(GAMMA_F, [5, 9, 10, 12, 20, 30, 50, 100])
    low_den_counterexamples = [row for row in approximation_rows if not row["strict_9_over_5_is_best_at_this_bound"]]
    symbolic = symbolic_domain_facts()
    provenance = {
        "F151": text_contains(F151, ["B_legacy_strict_bridge_target_v1", "R_damp_renorm_target_v1", "no false pass"]),
        "K1": text_contains(K1, ["K_legacy_ont", "K_strict_gate", "not yet rigorously identified"]),
        "K2": text_contains(K2, ["omega = 0.18575", "eta   = 1.80000", "not directly derived"]),
    }

    report = {
        "status": "OPEN_ANALYTIC_DOMAIN_TRADEOFF_TRACE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_KERNEL_ANALYTIC_DOMAIN_AUDIT__NOT_A_THEOREM",
        "source_reports": {
            "placement_discriminator": str(PLACEMENT_AUDIT.relative_to(ROOT)),
            "F151_bridge_target_packet": str(F151.relative_to(ROOT)),
            "K1_kernel_split_note": str(K1.relative_to(ROOT)),
            "K2_strict_derivation_chain_note": str(K2.relative_to(ROOT)),
        },
        "constants": {
            "legacy_beta_tors": float(LEGACY_BETA_TORS),
            "strict_eta_fraction": "9/5",
            "strict_eta_float": float(STRICT_ETA),
            "D_f": ALPHA_GEO,
            "D_f_minus_1": GAMMA_F,
            "strict_eta_minus_D_f_minus_1": float(STRICT_ETA) - GAMMA_F,
        },
        "provenance_checks": provenance,
        "symbolic_domain_facts": symbolic,
        "rationalization_audit": {
            "continued_fraction_D_f_minus_1_first_terms": continued_fraction(GAMMA_F),
            "approximation_rows": approximation_rows,
            "strict_9_over_5_best_for_all_tested_denominator_bounds": all(
                row["strict_9_over_5_is_best_at_this_bound"] for row in approximation_rows
            ),
            "counterexample_best_approximants": low_den_counterexamples,
            "verdict": "eta=9/5 cannot be justified merely as the best small-denominator rationalization of D_f-1; e.g. 7/4 is closer for max denominator 5 and 16/9 is closer for max denominator 9.",
        },
        "bridge_domain_choice": {
            "positive_real_operational_route": "Compatible with prior D_f-1 denominator-placement evidence and avoids complex branch-cover obligations, but does not preserve legacy meromorphic C-domain semantics.",
            "complex_finite_cover_route": "Compatible with strict eta=9/5 via z=t^5, but D_f-1 itself is not a rational finite-cover exponent; requires an extra rationalization/renormalization theorem.",
            "legacy_meromorphic_route": "Preserves legacy one-pole meromorphic structure but cannot silently produce d^(9/5) or d^(D_f-1) damping.",
        },
        "replayed_numeric_context": {
            "placement_denominator_wins_all_windows": placement["aggregate_summary"]["all_windows_denominator_beats_numerator_and_inverse_by_bic"],
            "placement_min_bic_numerator_minus_denominator": placement["aggregate_summary"]["min_bic_numerator_minus_denominator"],
            "placement_min_bic_inverse_minus_denominator": placement["aggregate_summary"]["min_bic_inverse_numerator_minus_denominator"],
        },
        "honest_interpretation": [
            "The strict kernel is computationally/gate selected and has a controlled finite algebraic cover because eta=9/5, but it is not meromorphic on the original complex d-plane.",
            "The D_f-1 route is the strongest current real-positive damping clue, but it is analytically more singular than strict eta=9/5 if promoted as an exact complex exponent.",
            "Therefore the next bridge theorem must explicitly choose a domain: positive real operational kernel, finite branched-cover strict kernel, or a rationalization map from D_f-1 to 9/5.  It cannot claim all three silently.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No proof that D_f-1 rationalizes to eta=9/5 is exported.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Build a real-positive transport bridge lemma or a separate finite-cover rationalization lemma; do not mix real-positive D_f-1 evidence with complex meromorphic claims without an explicit domain theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch analytic-domain audit: legacy vs strict kernel\n\n"
        "Status: analytic-domain tradeoff recorded; no bridge theorem.\n\n"
        f"- Legacy denominator is meromorphic with pole `z=-100`; strict eta is exact `9/5` with branch points `z=0, infinity` and finite cover `z=t^5`.\n"
        f"- `D_f-1={GAMMA_F:.15f}` differs from strict `9/5` by `{float(STRICT_ETA) - GAMMA_F:.15f}` and is not a finite-sheet rational exponent.\n"
        f"- `9/5` is not the best tested rational approximation to `D_f-1`: all-tested-best flag `{report['rationalization_audit']['strict_9_over_5_best_for_all_tested_denominator_bounds']}`.\n"
        f"- Prior denominator-placement win replayed: `{report['replayed_numeric_context']['placement_denominator_wins_all_windows']}`.\n"
        "- No false pass: no kernel identity, no D_f-1→9/5 rationalization theorem, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
