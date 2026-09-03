#!/usr/bin/env python3
"""FIN ST2217--ST2231: exact hidden-synergy/refinement obstruction."""
import itertools
import math

from fin_st2172_st2261_common import hidden_synergy_distribution, moments, write_packet, write_round


LO, HI = 2217, 2231
NAMES = [
    "BinaryFiberSynergyFamily", "ProbabilityPositivity", "OnePointMomentEquality",
    "TwoPointMomentEquality", "ThirdMomentSeparation", "ContinuumPairIndistinguishability",
    "PairReconstructionNoGo", "EntropySynergy", "MarginalizationErasure",
    "CoassociativeErasure", "RefinementEmbedding", "RefinementLiftNonuniqueness",
    "InternalLayerObserverBoundary", "RoundFourVerdict", "RoundFourRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def entropy(p):
    return -sum(v * math.log2(v) for v in p.values() if v > 0)


def marginal(p, indices):
    out = {}
    for x, mass in p.items():
        key = tuple(x[i] for i in indices)
        out[key] = out.get(key, 0.0) + mass
    return out


def main():
    samples = [-1.0, -0.75, -0.25, 0.0, 0.25, 0.75, 1.0]
    distributions = {e: hidden_synergy_distribution(e) for e in samples}
    stats = {e: moments(p) for e, p in distributions.items()}
    pair_marginals = {e: [marginal(p, pair) for pair in ((0, 1), (0, 2), (1, 2))]
                      for e, p in distributions.items()}
    reference = pair_marginals[0.0]
    max_pair_difference = max(
        abs(row[key] - reference[j][key])
        for e in samples for j, row in enumerate(pair_marginals[e]) for key in row
    )
    pplus, pminus = distributions[0.75], distributions[-0.75]
    tv = 0.5 * sum(abs(pplus[x] - pminus[x]) for x in pplus)
    # Marginalizing in either order is ordinary finite summation.
    order_a = marginal(marginal(pplus, (0, 1)), (0,))
    direct = marginal(pplus, (0,))
    coassoc_residual = max(abs(order_a[k] - direct[k]) for k in direct)

    x = {}
    x["ST2217"] = packet(2217, "Constructed", "Abstract binary refinement fiber; not asserted as an already exported FIN state law.", {
        "formula": "p_epsilon(x)=1/8(1+epsilon x1 x2 x3)", "support_size": 8, "epsilon_interval": [-1, 1]})
    x["ST2218"] = packet(2218, "Proven", "Exact for every |epsilon|<=1.", {
        "sample_minimum_probabilities": {str(e): stats[e]["minimum_probability"] for e in samples},
        "sample_masses": {str(e): stats[e]["mass"] for e in samples}})
    x["ST2219"] = packet(2219, "Proven", "Character orthogonality on the cube.", {
        "all_sample_means": {str(e): stats[e]["means"] for e in samples}})
    x["ST2220"] = packet(2220, "Proven", "Every two-variable marginal is uniform and every pair moment is zero.", {
        "maximum_pair_marginal_difference": max_pair_difference,
        "all_sample_pair_moments": {str(e): stats[e]["pairs"] for e in samples}})
    x["ST2221"] = packet(2221, "Proven", "Direct exact summation gives E[x1 x2 x3]=epsilon.", {
        "sample_triples": {str(e): stats[e]["triple"] for e in samples}})
    x["ST2222"] = packet(2222, "Proven", "A continuum of probability laws has identical complete one/two-point records.", {
        "parameter_cardinality": "continuum", "pair_record_depends_on_epsilon": False})
    x["ST2223"] = packet(2223, "Proven impossibility", "No deterministic function of the complete one/two-point record can equal epsilon for all family members.", {
        "proof": "identical input record with distinct required outputs", "pair_to_triple_reconstruction_exists": False})
    x["ST2224"] = packet(2224, "Proven", "All pair entropies are two bits while full entropy detects |epsilon|, not its sign.", {
        "pair_entropy_bits": 2.0, "full_entropy_bits": {str(e): entropy(distributions[e]) for e in samples},
        "sign_selected_by_entropy": False})
    x["ST2225"] = packet(2225, "Proven", "Tracing out any one of the three fine variables returns the same uniform pair state.", {
        "fine_total_variation_plus_vs_minus": tv, "coarse_total_variation": 0.0})
    x["ST2226"] = packet(2226, "Proven", "Finite marginalization is associative; changing the order does not restore erased synergy.", {
        "coassociative_marginal_residual": coassoc_residual})
    x["ST2227"] = packet(2227, "Conditional construction", "Three binary fiber labels on a refined triangle realize the exact family.", {
        "compatible_with_12_to24_binary_fibers": True, "strict_FIN_probability_law_exported": False})
    x["ST2228"] = packet(2228, "Proven nonuniqueness", "A fixed coarse/pair record admits every epsilon in [-1,1]; refinement intertwining alone cannot choose a lift.", {
        "lift_modulus_dimension_at_least": 1, "coassociativity_selects_epsilon": False})
    x["ST2229"] = packet(2229, "Conditional implication", "An observer restricted to pair/coarse instruments is exactly blind to fine third-order synergy.", {
        "requires": ["binary fiber model", "pair-limited instrument"],
        "does_not_derive": ["physical scale", "cosmology", "actual observer"]})
    x["ST2230"] = packet(2230, "Round verdict", "Irreducible ternary information is mathematically independent of all one/two-point data and can be erased by refinement coarse-graining.", {
        "strict_pair_kernel_can_reconstruct_general_ternary_state": False})
    x["ST2231"] = packet(2231, "Recommendation", "Test the strongest remaining positive construction: use strict loop products as a target-free flag-complex Hodge measure, then attack uniqueness.", {
        "next": ["support clique complex", "tau Hodge", "power-family counterexample", "spectral divergence", "continuum/locality boundary"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
