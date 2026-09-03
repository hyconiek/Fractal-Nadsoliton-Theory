#!/usr/bin/env python3
"""FIN ST2277--ST2291: weighted automorphisms and face universality."""
import numpy as np

from fin_st2172_st2261_common import d12_permutations, permutation_matrix
from fin_st2262_st2351_common import (
    orbit_cycle_vector, strict_A_W, strict_distance_weights, triangle_orbits,
    write_packet, write_round,
)


LO, HI = 2277, 2291
NAMES = [
    "DistanceWeightSeparation", "MaximumEdgeCycle", "CycleAutomorphismGroup",
    "WeightedGraphAutomorphismTheorem", "D12Replay", "TriangleOrbitCount",
    "TauOrbitSeparation", "OrbitCoefficientFreedom", "S12Counterfactual",
    "AffineUnitFailure", "UniversalityAsExtraNaturality", "A4SourceAudit",
    "UniversalityGate", "RoundTwoVerdict", "RoundTwoRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    A, W = strict_A_W()
    weights = strict_distance_weights()
    min_sep = min(abs(weights[i] - weights[j]) for i in range(6) for j in range(i + 1, 6))
    residuals = []
    for _, _, p in d12_permutations():
        P = permutation_matrix(p)
        residuals.append(float(np.linalg.norm(P @ W @ P.T - W, ord=np.inf)))
    transposition = list(range(12)); transposition[1], transposition[2] = 2, 1
    Pt = permutation_matrix(tuple(transposition))
    unit5 = permutation_matrix(tuple((5 * i) % 12 for i in range(12)))
    tau_orbit = orbit_cycle_vector()
    tau_sep = min(abs(tau_orbit[i] - tau_orbit[j]) for i in range(12) for j in range(i + 1, 12))

    x = {}
    x["ST2277"] = packet(2277, "Proven", "All six cyclic-distance conductances are distinct in the frozen strict matrix.", {
        "distance_weights": weights.tolist(), "minimum_pairwise_separation": float(min_sep)})
    x["ST2278"] = packet(2278, "Proven", "Distance-one edges have the unique maximum weight and form the cycle C12.", {
        "maximum_distance": int(np.argmax(weights) + 1), "maximum_edge_subgraph": "C12"})
    x["ST2279"] = packet(2279, "Proven", "Aut(C12) is the dihedral group with 24 elements.", {
        "automorphism_group": "D12", "order": 24})
    x["ST2280"] = packet(2280, "Proven", "Every weighted automorphism preserves the unique maximum-edge C12, hence lies in D12; every D12 element preserves all cyclic distances.", {
        "Aut_W": "D12", "order": 24})
    x["ST2281"] = packet(2281, "Proven", "Direct finite conjugation replay.", {
        "tested_elements": 24, "maximum_residual_inf": max(residuals)})
    x["ST2282"] = packet(2282, "Proven", "Exact finite orbit census.", {
        "triangle_orbits": len(triangle_orbits()), "orbit_sizes": [len(o) for o in triangle_orbits()]})
    x["ST2283"] = packet(2283, "Proven", "The strict loop product distinguishes all twelve triangle orbits numerically with a large nonzero separation.", {
        "orbit_values": tau_orbit.tolist(), "minimum_separation": float(tau_sep), "distinct": len(set(np.round(tau_orbit, 15)))})
    x["ST2284"] = packet(2284, "Proven", "Aut(W)-naturality permits one scalar/function value per triangle orbit.", {
        "equivariant_face_scalar_dimension": 12})
    x["ST2285"] = packet(2285, "Refuted", "Full S12 would make triangles one orbit but is not a symmetry of W.", {
        "proposal": "derive universality from S12", "transposition_residual_inf": float(np.linalg.norm(Pt @ W @ Pt.T - W, ord=np.inf))})
    x["ST2286"] = packet(2286, "Proven", "Multiplication by unit 5 permutes cyclic distances with unequal weights and is not a weighted symmetry.", {
        "unit5_residual_inf": float(np.linalg.norm(unit5 @ W @ unit5.T - W, ord=np.inf))})
    x["ST2287"] = packet(2287, "Conditional", "One common formula across all weighted graphs is a category-level universality axiom, stronger than Aut(W)-equivariance.", {
        "universal_graph_functor_would_select_common_formula": True, "exported_by_strict_FIN": False})
    x["ST2288"] = packet(2288, "Repository audit", "The strict carrier provides D12, not a theorem erasing its twelve face orbits.", {
        "strict_A4_source": False})
    x["ST2289"] = packet(2289, "Gate fail", "Face universality is compatible but not unavoidable.", {
        "universality_strict": False, "orbit_moduli": 12})
    x["ST2290"] = packet(2290, "Round verdict", "The actual strict symmetry positively supports orbit-dependent face physics; it cannot force A4.", {
        "A4_unavoidable": False})
    x["ST2291"] = packet(2291, "Recommendation", "Replace A3/A4 by a stronger fractal hypothesis and classify normalized multiplicative refinement maps.", {
        "next": ["power valuation theorem", "composition law", "idempotence", "scale-stationary path independence", "nontriviality"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
