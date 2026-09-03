#!/usr/bin/env python3
"""FIN ST2202--ST2216: adaptive-law ternary-source audit."""
import numpy as np

from fin_st2172_st2261_common import fourier_state, strict_A_W, write_packet, write_round


LO, HI = 2202, 2216
NAMES = [
    "DeclaredOperatorGradientFlow", "SpectralSourceInvariantManifold",
    "CommutingFlowEigenbasisBoundary", "SelfGeneratedDephasedCovariance",
    "DephasedCovarianceCommutator", "InitialStateOriginTransfer",
    "LinearHebbFixedPoint", "StrictWIndefinitenessNoGo", "OjaCompressionBoundary",
    "SymmetricPreparationOutcome", "PotentialReconstructionFreedom",
    "GaussianAdaptiveClosureBoundary", "EquivariantAdaptiveSymmetryNoGo",
    "RoundThreeVerdict", "RoundThreeRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def dephased_covariance(A, psi):
    values, vectors = np.linalg.eigh(A)
    C = np.zeros_like(A, dtype=complex)
    groups = []
    for value in sorted(set(np.round(values, 10))):
        idx = np.where(abs(values - value) < 1e-8)[0]
        P = vectors[:, idx] @ vectors[:, idx].T
        z = P @ psi
        C += np.outer(z, z.conj())
        groups.append(len(idx))
    return C, groups


def main():
    A, W = strict_A_W()
    delta = np.zeros(12, dtype=complex)
    delta[0] = 1.0
    Cdelta, multiplicities = dephased_covariance(A, delta)
    uniform = np.ones(12, dtype=complex) / np.sqrt(12)
    Cuniform, _ = dephased_covariance(A, uniform)
    comm = float(np.linalg.norm(A @ Cdelta - Cdelta @ A, ord=np.inf))
    eigW = np.linalg.eigvalsh(W)
    fixed_residual = float(np.linalg.norm(Cdelta.real - W, ord="fro"))
    origin_shift = np.roll(Cdelta, 1, axis=0)
    origin_shift = np.roll(origin_shift, 1, axis=1)

    x = {}
    x["ST2202"] = packet(2202, "Classified", "The repository formula is a family until V, C, projection and state law are sourced.", {
        "law": "dot K=Pi[C-V'(K)]", "free_objects": ["V", "C", "Pi", "initial state", "time scale"]})
    x["ST2203"] = packet(2203, "Proven", "If C and V'(K) lie in C*(K), the vector field remains in the same commutative spectral algebra.", {
        "can_generate_noncommuting_associator": False})
    x["ST2204"] = packet(2204, "Proven", "For a simple spectrum, commuting flow changes eigenvalues but not eigenvectors; degeneracies require state data to split.", {
        "new_orientation_without_degenerate_block_source": False})
    x["ST2205"] = packet(2205, "Constructed", "Infinite-time unitary covariance of a localized state.", {
        "spectral_multiplicities": multiplicities, "trace": float(np.trace(Cdelta).real),
        "rank": int(np.linalg.matrix_rank(Cdelta, tol=1e-10))})
    x["ST2206"] = packet(2206, "Proven", "Time dephasing is block diagonal in A's spectral projections.", {
        "commutator_residual_inf": comm})
    x["ST2207"] = packet(2207, "Proven", "The localized covariance changes when the prepared delta origin is translated.", {
        "translation_difference_fro": float(np.linalg.norm(Cdelta - origin_shift)),
        "origin_is_sourced_by_preparation": True})
    x["ST2208"] = packet(2208, "Proven", "For dot W=eta(C-gamma W), the only fixed point is W=C/gamma.", {
        "gamma_positive": True, "strict_fixed_point_residual_at_gamma1_fro": fixed_residual})
    x["ST2209"] = packet(2209, "Refuted", "This targets positive-covariance linear Hebb with positive decay, not all nonlinear learning laws.", {
        "proposal": "strict W is a fixed point of linear covariance Hebb", "W_lambda_min": float(eigW.min()),
        "W_lambda_max": float(eigW.max()), "PSD_covariance_can_equal_positive_gamma_W": False})
    x["ST2210"] = packet(2210, "Proven", "Uniform preparation gives a rank-one learned covariance; localized dephasing gives rank seven, neither derives W.", {
        "uniform_covariance_rank": int(np.linalg.matrix_rank(Cuniform, tol=1e-10)),
        "localized_dephased_rank": int(np.linalg.matrix_rank(Cdelta, tol=1e-10))})
    x["ST2211"] = packet(2211, "Proven", "The uniform state is the zero mode and carries no chiral/origin datum.", {
        "uniform_energy": float(np.real(uniform.conj() @ A @ uniform)),
        "overlap_with_k1": float(abs(np.vdot(fourier_state(1), uniform)))})
    x["ST2212"] = packet(2212, "Proven underdetermination", "Infinitely many V can be made stationary at a chosen target K; stationarity is not derivation.", {
        "target_fitting_needed": True, "unique_variational_law_from_K": False})
    x["ST2213"] = packet(2213, "Proven", "Linear Gaussian state/operator evolution preserves zero connected third cumulants.", {
        "escape_requires": ["non-Gaussian initial law", "nonlinear state interaction", "non-Gaussian noise"]})
    x["ST2214"] = packet(2214, "Proven", "An equivariant deterministic recursion started in the fixed set stays in the fixed set by induction.", {
        "roundoff_or_noise_branch_is_strict_selector": False, "asymmetric_input_needed": True})
    x["ST2215"] = packet(2215, "Round verdict", "Current adaptive laws amplify or learn supplied covariance/state structure; they do not export an intrinsic odd or irreducible ternary source.", {
        "adaptive_source_gate_pass": False, "linear_Hebb_strictW_refuted": True})
    x["ST2216"] = packet(2216, "Recommendation", "Use an exact refinement counterexample to decide whether pair/coarse records can ever reconstruct hidden third-order information.", {
        "next": ["binary-fiber synergy family", "equal one/two-point moments", "different triple moment", "coarse marginal", "refinement no-lift theorem"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
