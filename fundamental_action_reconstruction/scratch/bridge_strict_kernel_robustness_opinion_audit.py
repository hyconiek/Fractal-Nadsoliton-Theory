#!/usr/bin/env python3
"""Scratch audit of the opinion that strict-kernel robustness proves physical law.

The opinion cites P2086/P2260/P2270/P1620-style evidence and concludes that
eta=1.8 is not a numerical artifact but a physical fractal law, with legacy as
near-field and strict as whole-ring generalization.  This audit checks the cited
repo artifacts and separates supported operational robustness from unproven
bridge/ontology claims.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GEN = ROOT / "fundamental_action_reconstruction" / "generated"
OUT_JSON = HERE / "bridge_strict_kernel_robustness_opinion_audit_report.json"
OUT_MD = HERE / "bridge_strict_kernel_robustness_opinion_audit_report.md"

SOURCES = {
    "P2086_EOM": GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json",
    "P2260_STOCHASTIC": GEN / "p2260_s1210_strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe.json",
    "P2270_LIPSCHITZ": GEN / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.json",
    "P1620_POSTERIOR": GEN / "p1620_s570_strict_measured_covariance_posterior_summary.json",
    "DF_MARGINAL_RG": HERE / "bridge_df_marginal_rg_flow_report.json",
}

ALPHA_GEO = 4.0 * math.log(2.0)
GAMMA_F = ALPHA_GEO - 1.0
STRICT_ETA = 1.8
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def symbolic_limit_checks() -> dict[str, Any]:
    d, alpha, omega_l, phi_l, beta_t, omega_s, phi_s, beta, eta = sp.symbols(
        "d alpha omega_l phi_l beta_t omega_s phi_s beta eta", positive=True, real=True
    )
    k_legacy = alpha * sp.cos(omega_l * d + phi_l) / (1 + beta_t * d)
    k_strict = sp.cos(omega_s * d + phi_s) / (1 + beta * d**eta)
    ratio_limit = sp.limit(k_legacy / k_strict, d, 0, dir="+")
    damping_ratio_limit = sp.limit((1 / (1 + beta_t * d)) / (1 / (1 + beta * d**eta)), d, 0, dir="+")
    numeric_ratio = float(
        ALPHA_GEO * math.cos(LEGACY["phi"]) / math.cos(STRICT["phi"])
    )
    return {
        "symbolic_full_kernel_limit_Klegacy_over_Kstrict_d_to_0_plus": sp.sstr(ratio_limit),
        "numeric_full_kernel_limit_with_repo_parameters": numeric_ratio,
        "alpha_geo": ALPHA_GEO,
        "numeric_limit_equals_alpha_geo": abs(numeric_ratio - ALPHA_GEO) < 1e-12,
        "phase_factor_ratio_cos_phi_legacy_over_cos_phi_strict": math.cos(LEGACY["phi"]) / math.cos(STRICT["phi"]),
        "damping_only_ratio_limit": sp.sstr(damping_ratio_limit),
        "verdict": "The claimed alpha_geo near-field limit is true only for amplitude/damping-normalized or phase-aligned comparison; the full kernels give alpha_geo*cos(phi_legacy)/cos(phi_strict), not alpha_geo.",
    }


def main() -> None:
    p2086 = load_json(SOURCES["P2086_EOM"])
    p2260 = load_json(SOURCES["P2260_STOCHASTIC"])
    p2270 = load_json(SOURCES["P2270_LIPSCHITZ"])
    p1620 = load_json(SOURCES["P1620_POSTERIOR"])
    rg_flow = load_json(SOURCES["DF_MARGINAL_RG"])

    p2260_probe = p2260["strict_nu_branch_group_policy_stochastic_boundary_hit_validation_probe"]
    p2270_probe = p2270["strict_nu_branch_group_policy_symbolic_derivative_bound_probe"]
    posterior_eta = p1620["posterior_update"]["eta"]
    eta_ci_contains_strict = posterior_eta["q05"] <= STRICT_ETA <= posterior_eta["q95"]
    eta_ci_contains_df_minus_one = posterior_eta["q05"] <= GAMMA_F <= posterior_eta["q95"]
    stochastic_draws_total = int(p2260_probe["inputs"]["draws_per_horizon"] * len(p2260_probe["inputs"]["horizon_scales"]))

    evidence_bundle = {
        "P2086_EOM_residual_evidence": {
            "result_kind": p2086["result_kind"],
            "all_symbolic_residual_zero": p2086["gatekeeper_checks"]["all_symbolic_residual_zero"],
            "all_numeric_probe_residual_zero": p2086["gatekeeper_checks"]["all_numeric_probe_residual_zero"],
            "global_status": p2086["global_status"],
            "c3_theorem_proven": p2086["gatekeeper_checks"]["c3_theorem_proven"],
        },
        "P2260_stochastic_boundary_evidence": {
            "result_kind": p2260["result_kind"],
            "draws_per_horizon": p2260_probe["inputs"]["draws_per_horizon"],
            "horizon_count": len(p2260_probe["inputs"]["horizon_scales"]),
            "total_draws": stochastic_draws_total,
            "global_first_hit_rate": p2260_probe["stochastic_hit_rate_comparison"]["global_first_hit_rate"],
            "global_second_hit_rate": p2260_probe["stochastic_hit_rate_comparison"]["global_second_hit_rate"],
        },
        "P2270_lipschitz_evidence": {
            "result_kind": p2270["result_kind"],
            "admissible_box": p2270_probe["admissible_box"],
            "l1_lipschitz_over_box": p2270_probe["symbolic_box_certificates"]["l1_lipschitz_over_box"],
            "monotone_in_rho_over_box": p2270_probe["symbolic_box_certificates"]["monotone_in_rho_over_box"],
            "monotone_in_kappa_over_box": p2270_probe["symbolic_box_certificates"]["monotone_in_kappa_over_box"],
        },
        "P1620_posterior_evidence": {
            "status": p1620["status"],
            "eta_posterior": posterior_eta,
            "eta_ci_contains_strict_eta_1p8": eta_ci_contains_strict,
            "eta_ci_contains_D_f_minus_1": eta_ci_contains_df_minus_one,
            "legacy_bridge_used": p1620["legacy_bridge_used"],
            "strict_only": p1620["strict_only"],
        },
        "DF_marginal_rg_bridge_prep_replay": {
            "result_kind": rg_flow["result_kind"],
            "denominator_residual": rg_flow["symbolic_rg_flow"]["denominator_residual"],
            "max_exact_rg_abs_residual_vs_strict": rg_flow["aggregate_summary"]["max_exact_rg_abs_residual_vs_strict"],
            "no_bridge_theorem_status": rg_flow["status"],
        },
    }

    claim_audit = {
        "strict_is_not_single_point_overfit": {
            "verdict": "supported_operationally_not_theorem_level",
            "reason": "P2086 residual zeros, P2260 finite-sample zero hit rates, P2270 Lipschitz/monotonic box bounds, and P1620 posterior concentration jointly argue against a single calibration-point artifact.",
        },
        "eta_1p8_is_physical_law": {
            "verdict": "not_discharged",
            "reason": "The evidence supports an operational strict parameter, but no theorem derives eta=1.8 from strict-side nadsoliton dynamics; P1620 is posterior evidence and the D_f RG flow remains a target.",
        },
        "legacy_near_field_strict_whole_ring_generalization": {
            "verdict": "plausible_bridge_story_but_not_exported",
            "reason": "Scratch D_f probes give a denominator/RG candidate, but K1/F2 still forbid silent identification or role transfer.",
        },
        "alpha_geo_limit_claim": symbolic_limit_checks(),
    }

    support_score = {
        "operational_robustness_supported": all(
            [
                evidence_bundle["P2086_EOM_residual_evidence"]["all_symbolic_residual_zero"],
                evidence_bundle["P2086_EOM_residual_evidence"]["all_numeric_probe_residual_zero"],
                evidence_bundle["P2260_stochastic_boundary_evidence"]["global_first_hit_rate"] == 0.0,
                evidence_bundle["P2260_stochastic_boundary_evidence"]["global_second_hit_rate"] == 0.0,
                evidence_bundle["P2270_lipschitz_evidence"]["l1_lipschitz_over_box"] < 0.5,
                evidence_bundle["P1620_posterior_evidence"]["eta_ci_contains_strict_eta_1p8"],
            ]
        ),
        "physical_law_or_bridge_discharged": False,
        "near_field_alpha_geo_full_kernel_limit_disconfirmed": not claim_audit["alpha_geo_limit_claim"]["numeric_limit_equals_alpha_geo"],
    }

    report = {
        "status": "OPEN_ROBUSTNESS_OPINION_AUDIT_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_STRICT_KERNEL_ROBUSTNESS_OPINION_AUDIT__NOT_A_THEOREM",
        "source_reports": {key: str(path.relative_to(ROOT)) for key, path in SOURCES.items()},
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "D_f_minus_1": GAMMA_F,
            "strict_eta": STRICT_ETA,
        },
        "evidence_bundle": evidence_bundle,
        "claim_audit": claim_audit,
        "support_score": support_score,
        "honest_interpretation": [
            "The cited artifacts support strict-kernel operational robustness beyond a single fitting point.",
            "They do not by themselves prove that eta=1.8 is a physical fractal law or that legacy and strict kernels are bridged.",
            "The alpha_geo near-field limit must be stated with normalization/phase caveats: the full-kernel d->0 ratio is not alpha_geo for the repo parameter values.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No theorem deriving eta=1.8 from nadsoliton/fractal dynamics is exported.",
            "No theorem deriving the D_f marginal RG law is exported.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Use the robustness evidence only as operational support; the next proof task remains deriving the D_f marginal RG beta law and a phase-normalized near-field map, not declaring eta=1.8 a physical law by posterior robustness alone.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-kernel robustness opinion audit\n\n"
        "Status: robustness evidence checked; no bridge theorem.\n\n"
        f"- Operational robustness supported: `{support_score['operational_robustness_supported']}` from P2086/P2260/P2270/P1620.\n"
        f"- P2270 Lipschitz `L1={evidence_bundle['P2270_lipschitz_evidence']['l1_lipschitz_over_box']:.12f}`; P2260 total stochastic draws `{stochastic_draws_total}` with global hit rates `0/0`.\n"
        f"- P1620 eta posterior mean `{posterior_eta['mean']}`, q05 `{posterior_eta['q05']}`, q95 `{posterior_eta['q95']}`; contains strict eta `{eta_ci_contains_strict}` and D_f-1 `{eta_ci_contains_df_minus_one}`.\n"
        f"- Full-kernel d->0 ratio `Klegacy/Kstrict={claim_audit['alpha_geo_limit_claim']['numeric_full_kernel_limit_with_repo_parameters']:.12f}` vs `alpha_geo={ALPHA_GEO:.12f}`; alpha limit claim full-kernel-valid `{claim_audit['alpha_geo_limit_claim']['numeric_limit_equals_alpha_geo']}`.\n"
        "- No false pass: no kernel identity, no eta physical-law theorem, no D_f RG theorem, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
