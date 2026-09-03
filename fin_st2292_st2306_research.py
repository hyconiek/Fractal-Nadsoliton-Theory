#!/usr/bin/env python3
"""FIN ST2292--ST2306: multiplicative valuations and conditional p=1 theorem."""
import numpy as np

from fin_st2262_st2351_common import (
    mean_normalize, orbit_cycle_vector, power_map, write_packet, write_round,
)


LO, HI = 2292, 2306
NAMES = [
    "PositiveMultiplicativeValuation", "PowerLawClassification", "NormalizedPowerMap",
    "PowerMapComposition", "IdempotenceEquation", "NonconstantProfileTheorem",
    "TrivialErasureBranch", "NontrivialIdentityBranch", "ScaleStationaryPathPrinciple",
    "PrincipleNecessityAudit", "ContinuousSemigroupFreedom", "CoassociativityBoundary",
    "ConditionalPSelectionGate", "RoundThreeVerdict", "RoundThreeRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    x0 = mean_normalize(orbit_cycle_vector())
    p, q = 0.7, 1.8
    composition_residual = float(np.linalg.norm(power_map(power_map(x0, q), p) - power_map(x0, p*q), np.inf))
    idem = {str(r): float(np.linalg.norm(power_map(power_map(x0, r), r) - power_map(x0, r), np.inf))
            for r in (0.0, 0.5, 1.0, 1.8, 2.0)}

    out = {}
    out["ST2292"] = packet(2292, "Defined", "Assumes a positive continuous valuation of multiplicative edge/loop composition.", {
        "law": "g(xy)=g(x)g(y)", "domain": "positive reals"})
    out["ST2293"] = packet(2293, "Proven classification", "Continuity excludes discontinuous Cauchy characters.", {
        "theorem": "Every continuous positive multiplicative g is g(x)=x^p for a real p; monotonicity gives p>=0."})
    out["ST2294"] = packet(2294, "Constructed", "Mean normalization removes common amplitude but not exponent.", {
        "map": "T_p(x)=x^p/mean(x^p)", "strict_profile_nonconstant": bool(np.ptp(x0) > 0)})
    out["ST2295"] = packet(2295, "Proven", "Exact algebra; finite residual shown for the strict orbit vector.", {
        "composition": "T_p composed with T_q equals T_(pq)", "residual_inf": composition_residual})
    out["ST2296"] = packet(2296, "Proven", "Scale-stationary two-step equals one-step requires T_(p^2)(x)=T_p(x).", {
        "idempotence_residuals": idem})
    out["ST2297"] = packet(2297, "Proven theorem", "Use any unequal pair xi/xj and compare its transformed ratio.", {
        "theorem": "For nonconstant positive x, T_p(T_p(x))=T_p(x) implies p^2=p, hence p in {0,1}."})
    out["ST2298"] = packet(2298, "Proven", "p=0 sends every positive profile to the uniform vector and erases all orbit information.", {
        "p": 0, "information_preserving": False, "output_distinct_values": 1})
    out["ST2299"] = packet(2299, "Proven", "p=1 is the identity on mean-normalized profiles.", {
        "p": 1, "information_preserving": True, "nontrivial_compression": False})
    out["ST2300"] = packet(2300, "Conditional axiom", "Stronger than ordinary coassociativity: the same direct map must represent one or any number of refinement steps.", {
        "axiom": "scale-stationary path independence/idempotence", "strict_exported": False})
    out["ST2301"] = packet(2301, "Axiom-removal audit", "Every premise is needed for the p=1 conclusion.", {
        "necessary": ["positive multiplicative power class", "nonconstant profile", "idempotence", "information preservation"],
        "removals": {"idempotence": "all p survive", "information preservation": "p=0 survives", "nonconstant": "p unidentifiable", "power class": "other maps survive"}})
    out["ST2302"] = packet(2302, "Proven counterfamily", "Continuous refinement semigroups remain possible without idempotence.", {
        "family": "p(t)=exp(-c t), c>=0", "composition": "p(t+s)=p(t)p(s)", "free_rate": "c"})
    out["ST2303"] = packet(2303, "Guardrail confirmed", "ST289/ST619 already prove that plain coassociativity leaves a rate modulus/per-level sequence.", {
        "ordinary_coassociativity_selects_p": False})
    out["ST2304"] = packet(2304, "Conditional pass", "p=1 follows only after adding scale-stationary idempotence and excluding erasure.", {
        "p_selected": 1, "strict_selection": False})
    out["ST2305"] = packet(2305, "Round verdict", "Fractal path independence can select p=1, but only in a strong idempotent form absent from current FIN.", {
        "breakthrough_conditional": True, "strict_theorem": False})
    out["ST2306"] = packet(2306, "Recommendation", "Test the strong principle on the actual nonconstant orbit vector, including relabellings, entropy flow and level-dependent rates.", {
        "next": ["permutation self-similarity", "entropy monotonicity", "fixed-profile no-go", "vertical-rate escape", "legacy sign boundary"]})
    write_round(LO, HI, out)


if __name__ == "__main__":
    main()
