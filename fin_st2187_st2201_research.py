#!/usr/bin/env python3
"""FIN ST2187--ST2201: state-dependent chirality and selector audit."""
import hashlib
import json
from pathlib import Path

import numpy as np

from fin_st2172_st2261_common import (
    ROOT, d12_permutations, edge_current, fourier_state, permutation_matrix,
    strict_A_W, write_packet, write_round,
)


LO, HI = 2187, 2201
NAMES = [
    "FourierBranchPair", "ReflectionExchange", "DegenerateEnergy", "SignedStateCurrent",
    "InvariantMixtureCancellation", "D12TwirlCancellation", "EquivariantChannelInvariantState",
    "ThermalStateNoChirality", "PitchforkBranchPair", "SymmetricNoiseNoUniqueBranch",
    "BiasImportsSelector", "ChiralBispectrumBoundary", "MemoryAreaBoundary",
    "RoundTwoVerdict", "RoundTwoRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    A, W = strict_A_W()
    psi_p, psi_m = fourier_state(1), fourier_state(-1)
    rho_p, rho_m = np.outer(psi_p, psi_p.conj()), np.outer(psi_m, psi_m.conj())
    reflection = permutation_matrix(tuple((-i) % 12 for i in range(12)))
    exchange_residual = float(np.linalg.norm(reflection @ rho_p @ reflection.T - rho_m, ord=np.inf))
    energies = [float(np.real(v.conj() @ A @ v)) for v in (psi_p, psi_m)]
    currents = [edge_current(r, W) for r in (rho_p, rho_m)]
    rho_mix = (rho_p + rho_m) / 2
    twirl = np.zeros_like(rho_p)
    for _, _, p in d12_permutations():
        P = permutation_matrix(p)
        twirl += P @ rho_p @ P.T
    twirl /= 24
    beta = 1.0
    ev, U = np.linalg.eigh(A)
    gibbs = (U * np.exp(-beta * ev)) @ U.conj().T
    gibbs /= np.trace(gibbs)
    source = ROOT / "fundamental_action_reconstruction/generated/p3048_s1998_memory_phase_area_odd_source_candidate.json"
    p3048 = json.loads(source.read_text())
    cert = p3048["finite_certificate"]

    x = {}
    x["ST2187"] = packet(2187, "Constructed", "The branch label is a state/preparation datum, not selected by A.", {
        "branches": [1, -1], "density_rank": 1})
    x["ST2188"] = packet(2188, "Proven", "Exact group action up to floating Fourier evaluation.", {
        "reflection_exchange_residual_inf": exchange_residual})
    x["ST2189"] = packet(2189, "Proven", "Real circulant A has lambda_k=lambda_-k.", {
        "energy_plus": energies[0], "energy_minus": energies[1], "difference": abs(energies[0] - energies[1])})
    x["ST2190"] = packet(2190, "Proven", "A valid receiver once a chiral state is prepared.", {
        "formula": "Lambda=sum_i 2 Im(rho_i,i+1 W_i+1,i)",
        "current_plus": currents[0], "current_minus": currents[1],
        "odd_sum_residual": abs(currents[0] + currents[1])})
    x["ST2191"] = packet(2191, "Proven", "The reflection-invariant state has no signed current.", {
        "mixture_current": edge_current(rho_mix, W)})
    x["ST2192"] = packet(2192, "Proven", "Group averaging removes the branch sign.", {
        "twirled_current": edge_current(twirl, W), "twirl_reflection_residual": float(np.linalg.norm(reflection @ twirl @ reflection.T - twirl))})
    x["ST2193"] = packet(2193, "Proven", "For a D12-covariant channel Phi, g rho=rho implies g Phi(rho)=Phi(rho).", {
        "invariant_input_can_generate_odd_expectation": False})
    x["ST2194"] = packet(2194, "Proven", "The Gibbs state is a function of A and is reflection invariant.", {
        "gibbs_current": edge_current(gibbs, W),
        "reflection_residual_inf": float(np.linalg.norm(reflection @ gibbs @ reflection.T - gibbs, ord=np.inf))})
    x["ST2195"] = packet(2195, "Constructed", "For dq/dt=(1-q^2)q, q=+-1 are stable and q=0 is unstable.", {
        "fixed_points": [-1, 0, 1], "stable": [-1, 1], "selected_unique_branch": False})
    x["ST2196"] = packet(2196, "Proven", "Reflection-symmetric initial/noise laws give equal branch probabilities by path-law symmetry.", {
        "ensemble_mean_branch": 0, "probability_plus": 0.5, "probability_minus": 0.5})
    x["ST2197"] = packet(2197, "Refuted", "A term h is an extra inversion-odd premise even if it makes one branch stable/preferred.", {
        "proposal": "use a biased pitchfork as strict selector", "bias_required": True})
    x["ST2198"] = packet(2198, "Strong evidence", "Reuses the exact P2718 finite certificate; it is a receiver, not a source law.", {
        "marker_values": [-2, 2], "translation_or_source_localizer_exported": False})
    x["ST2199"] = packet(2199, "Bounded no-go", "P3048 has a real three-point odd area but lacks strict provenance and polarity selection.", {
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "area_rows": cert["area_rows"], "finite_nonzero_rows": cert["finite_nonzero_area_rows"],
        "accepted_strict_rows": cert["accepted_strict_odd_source_value_rows"]})
    x["ST2200"] = packet(2200, "Round verdict", "State-dependent chirality exists, but every audited nonzero sign is carried by a prepared branch, chart orientation, or unselected polarity pair.", {
        "state_receiver_exists": True, "strict_state_source_exists": False, "QW2191_discharged": False})
    x["ST2201"] = packet(2201, "Recommendation", "Audit the actual adaptive law: an equivariant update may amplify a state but cannot create an odd source from invariant data.", {
        "next": ["gradient operator flow", "self-generated covariance", "commuting invariant manifold", "Hebbian fixed point", "non-Gaussian escape condition"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
