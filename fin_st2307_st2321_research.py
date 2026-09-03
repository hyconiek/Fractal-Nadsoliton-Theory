#!/usr/bin/env python3
"""FIN ST2307--ST2321: exact self-similar refinement falsification."""
import math
import numpy as np

from fin_st2262_st2351_common import (
    entropy_of_mean_one_weights, mean_normalize, orbit_cycle_vector, power_map,
    write_packet, write_round,
)


LO, HI = 2307, 2321
NAMES = [
    "StrictOrbitProfile", "CenteredLogCoordinates", "PermutationSelfSimilarityTheorem",
    "IdempotenceNumericalWitness", "CompressionEntropyFlow", "SharpeningEntropyFlow",
    "FixedProfileNoGo", "IdentityBranchBoundary", "ErasureBranchBoundary",
    "SemigroupRateFreedom", "LevelDependentEscape", "VerticalRateCompatibility",
    "LegacySignedBoundary", "RoundFourVerdict", "RoundFourRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    profile = mean_normalize(orbit_cycle_vector())
    y = np.log(profile) - np.mean(np.log(profile))
    ps = (0.0, 0.5, 1.0, 1.8, 2.0)
    idem = {str(p): float(np.linalg.norm(power_map(power_map(profile, p), p) - power_map(profile, p), np.inf)) for p in ps}
    ent = {str(p): entropy_of_mean_one_weights(power_map(profile, p)) for p in ps}
    h0 = entropy_of_mean_one_weights(profile)

    x = {}
    x["ST2307"] = packet(2307, "Computed", "Twelve strictly positive and distinct orbit values.", {
        "dimension": len(profile), "distinct": len(set(np.round(profile, 14))), "minimum": float(profile.min()), "maximum": float(profile.max())})
    x["ST2308"] = packet(2308, "Constructed", "Mean-zero log-ratio coordinates turn T_p into scalar multiplication y->p y.", {
        "centered_log_norm": float(np.linalg.norm(y)), "nonzero": bool(np.linalg.norm(y) > 0)})
    x["ST2309"] = packet(2309, "Proven theorem", "If T_p(x)=P x for a permutation P, then p||y||=||Py||=||y||; for p>0 and nonconstant x, p=1.", {
        "self_similarity_up_to_orbit_permutation_requires_p": 1})
    x["ST2310"] = packet(2310, "Computed witness", "Large residuals reject nontrivial idempotence on the strict profile.", {
        "idempotence_residuals_inf": idem})
    x["ST2311"] = packet(2311, "Proven/Computed", "For 0<p<1 the normalized power map moves the finite profile toward uniformity.", {
        "base_entropy_nats": h0, "p_half_entropy_nats": ent["0.5"], "entropy_increases": ent["0.5"] > h0})
    x["ST2312"] = packet(2312, "Proven/Computed", "For p>1 the map sharpens the largest orbit and lowers entropy.", {
        "p_1_8_entropy_nats": ent["1.8"], "p_2_entropy_nats": ent["2.0"], "entropy_decreases": ent["2.0"] < h0})
    x["ST2313"] = packet(2313, "Proven no-go", "A nonconstant finite positive profile cannot be an exact fixed shape of T_p for p not equal to 1, even up to relabelling.", {
        "nontrivial_exact_power_self_similarity": False})
    x["ST2314"] = packet(2314, "Boundary", "p=1 preserves the profile but performs no compression or renormalization evolution.", {
        "p": 1, "map": "identity"})
    x["ST2315"] = packet(2315, "Boundary", "p=0 is idempotent but collapses all orbit information to uniformity.", {
        "p": 0, "map": "erasure"})
    x["ST2316"] = packet(2316, "Proven", "Without fixed-profile/idempotence, p_n may vary and cumulative exponent is product_n p_n.", {
        "free_data": "one exponent/rate per level"})
    x["ST2317"] = packet(2317, "Proven escape", "A direct level-n map may use the cumulative exponent; ordinary path composition remains exact for arbitrary sequences.", {
        "path_consistency_without_stationarity": True, "rate_selected": False})
    x["ST2318"] = packet(2318, "Guardrail alignment", "Independent vertical q_N modes from ST552--ST641 supply further unsourced refinement freedom.", {
        "vertical_rates_per_level": True, "selected": False})
    x["ST2319"] = packet(2319, "Kernel-split boundary", "Legacy loop products have both signs; real noninteger powers require abs/sign choices.", {
        "noninteger_power_on_signed_legacy_real_without_extra_rule": False, "legacy_bridge_completed": False})
    x["ST2320"] = packet(2320, "Round verdict", "Exact finite self-similarity selects only identity or erasure; nontrivial fractal evolution requires new level data.", {
        "strict_nontrivial_self_similar_refinement": False})
    x["ST2321"] = packet(2321, "Recommendation", "Construct the minimal instrument and non-Gaussian state law that would make hidden ternary information operational, then mark every added premise.", {
        "next": ["parity observable", "Fisher information", "sample bound", "three-body exponential family", "adaptive source parameter"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
