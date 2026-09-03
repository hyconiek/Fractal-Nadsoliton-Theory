#!/usr/bin/env python3
"""FIN ST2172--ST2186: strict third-order object classification."""
import numpy as np

from fin_st2172_st2261_common import (
    TRIANGLES, cycle_field, d12_permutations, permutation_matrix,
    strict_A_W, triangle_orbits, write_packet, write_round,
)


LO, HI = 2172, 2186
NAMES = [
    "StrictD12SymmetryAudit", "OddScalarFixedPointNoGo", "QuadraticGaussianCore",
    "GaussianThirdCumulantZero", "LinearSemigroupGaussianClosure",
    "FunctionalCalculusCommutativity", "FunctionalAssociatorZero",
    "SingleGeneratorCurvatureNoGo", "StrictTriangleCycleField", "TriangleOrbitCensus",
    "CycleFieldCovariance", "ReducibleVersusConnectedThirdOrder",
    "CompetingNaturalCycleFunctions", "RoundOneVerdict", "RoundOneRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    A, W = strict_A_W()
    residuals = []
    for _, _, p in d12_permutations():
        P = permutation_matrix(p)
        residuals.append(float(np.linalg.norm(P @ A @ P.T - A, ord=np.inf)))
    powers = [np.linalg.matrix_power(A, k) for k in range(6)]
    comm = max(float(np.linalg.norm(X @ Y - Y @ X, ord=np.inf)) for X in powers for Y in powers)
    assoc = max(float(np.linalg.norm((X @ Y) @ Z - X @ (Y @ Z), ord=np.inf))
                for X in powers[:3] for Y in powers[:3] for Z in powers[:3])
    tau = cycle_field(W)
    orbits = triangle_orbits()
    orbit_spread = max(max(float(tau[TRIANGLES.index(t)]) for t in o) -
                       min(float(tau[TRIANGLES.index(t)]) for t in o) for o in orbits)
    tau1 = tau / np.mean(tau)
    tau2 = tau**2 / np.mean(tau**2)
    proportional_residual = float(np.linalg.norm(tau1 - tau2, ord=np.inf))

    x = {}
    x["ST2172"] = packet(2172, "Proven", "Numerical replay of the frozen strict matrix.", {
        "group_elements": 24, "max_D12_conjugation_residual_inf": max(residuals)})
    x["ST2173"] = packet(2173, "Proven", "Applies to a global reflection-odd scalar derived equivariantly from D12-fixed input.", {
        "theorem": "If r.x=x and F(r.x)=-F(x), then F(x)=0.", "strict_input_fixed": True})
    x["ST2174"] = packet(2174, "Constructed", "Dimensionless centered Gaussian state with covariance (A+I)^-1.", {
        "covariance_lambda_min": float(np.linalg.eigvalsh(np.linalg.inv(A + np.eye(12))).min())})
    x["ST2175"] = packet(2175, "Proven", "Wick theorem for the quadratic/zero-mean Gaussian core.", {
        "connected_third_cumulant": 0, "all_component_third_moments": 0})
    x["ST2176"] = packet(2176, "Proven", "Linear unitary, wave, or heat evolution of Gaussian data remains Gaussian.", {
        "new_connected_third_order_generated": False})
    x["ST2177"] = packet(2177, "Proven", "Audited powers 0 through 5; theorem holds for all Borel functions of one A.", {
        "max_polynomial_commutator_residual_inf": comm})
    x["ST2178"] = packet(2178, "Proven", "Ordinary matrix multiplication is associative; C*(A) cannot source a nonzero associator.", {
        "max_polynomial_associator_residual_inf": assoc})
    x["ST2179"] = packet(2179, "Refuted", "Only the single-generator functional-calculus proposal is refuted.", {
        "proposal": "derive curvature from [f(A),g(A)] or an associator", "commutator_and_associator_nonzero": False})
    x["ST2180"] = packet(2180, "Proven", "A reducible graph 3-cycle statistic, not a probabilistic connected cumulant.", {
        "formula": "tau_ijk=W_ij W_jk W_ki", "faces": len(tau), "minimum": float(tau.min()),
        "maximum": float(tau.max()), "nonzero": int(np.sum(tau != 0)), "positive": int(np.sum(tau > 0))})
    x["ST2181"] = packet(2181, "Proven", "Exact finite D12 orbit census on 220 unoriented triangle supports.", {
        "orbit_count": len(orbits), "orbit_sizes": [len(o) for o in orbits]})
    x["ST2182"] = packet(2182, "Proven", "The field is D12-natural and reflection-even as a support weight.", {
        "maximum_within_orbit_spread": orbit_spread, "orientation_selector": False})
    x["ST2183"] = packet(2183, "Proven", "Corrects the overly coarse earlier gate: a cycle observable exists, but no independent third cumulant follows.", {
        "reducible_cycle_exists": True, "irreducible_connected_third_order_from_pair_kernel": False})
    x["ST2184"] = packet(2184, "Refuted", "Two target-free D12-natural positive rules already disagree.", {
        "rules": ["tau/mean(tau)", "tau^2/mean(tau^2)"],
        "nonproportional_residual_inf": proportional_residual})
    x["ST2185"] = packet(2185, "Round verdict", "Strict pair data produce a positive even 3-cycle field, but not chirality, an associator, or a connected third cumulant.", {
        "strict_reducible_cycle": True, "strict_irreducible_ternary_source": False})
    x["ST2186"] = packet(2186, "Recommendation", "Test state-dependent signed candidates and whether degeneracy can be broken without imported preparation.", {
        "next": ["Fourier branch pair", "current sign", "invariant mixture", "equivariant-channel theorem", "spontaneous branch test"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
