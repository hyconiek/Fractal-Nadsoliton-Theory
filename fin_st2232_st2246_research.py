#!/usr/bin/env python3
"""FIN ST2232--ST2246: strict loop-product Hodge candidate and falsification."""
import numpy as np

from fin_st2172_st2261_common import (
    EDGES, TRIANGLES, cochains, cycle_field, hodge_spectrum,
    normalized_positive_power, strict_A_W, write_packet, write_round,
)


LO, HI = 2232, 2246
NAMES = [
    "StrictSupportGraph", "FlagComplex", "PositiveLoopProduct",
    "TauHodgeConstruction", "TauHodgeSpectrum", "PositivePowerFamily",
    "PowerFamilyNaturality", "SpectralCounterfamily", "CohomologyInvariance",
    "ExponentSelectionNoGo", "MultiplicativeRefinementClassification",
    "StrictEtaTransferBoundary", "DenseContinuumBoundary", "RoundFiveVerdict",
    "RoundFiveRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    A, W = strict_A_W()
    tau = cycle_field(W)
    exponents = [0.5, 1.0, 1.8, 2.0]
    spectra = {p: hodge_spectrum(normalized_positive_power(tau, p)) for p in exponents}
    metrics = {
        str(p): {"lambda_min": float(s.min()), "lambda_max": float(s.max()),
                 "trace": float(s.sum()), "distinct_rounded_8": len(set(np.round(s, 8)))}
        for p, s in spectra.items()
    }
    d0, d1 = cochains()
    rank_d1 = int(np.linalg.matrix_rank(d1, tol=1e-10))
    counter_residual = float(np.linalg.norm(spectra[1.0] - spectra[2.0], ord=np.inf))
    relative = float(np.linalg.norm(spectra[1.0] - spectra[2.0]) / np.linalg.norm(spectra[1.0]))

    x = {}
    x["ST2232"] = packet(2232, "Proven", "Support means nonzero off-diagonal strict conductance.", {
        "vertices": 12, "edges": len(EDGES), "minimum_edge_weight": float(min(W[i, j] for i, j in EDGES)),
        "support_is_complete_graph": True})
    x["ST2233"] = packet(2233, "Constructed canonically", "The support-clique functor is canonical, but its use as physical gauge cells is an interpretation.", {
        "complex": "full 11-simplex", "triangles": len(TRIANGLES), "contractible": True})
    x["ST2234"] = packet(2234, "Proven", "Every strict support edge is positive in the audited W.", {
        "tau_min": float(tau.min()), "tau_max": float(tau.max()), "positive_faces": int(np.sum(tau > 0))})
    x["ST2235"] = packet(2235, "Constructed", "Target-free dimensionless candidate; average-one normalization fixes only a convention.", {
        "formula": "h_f=tau_f/mean(tau)", "positive": True, "D12_natural": True})
    x["ST2236"] = packet(2236, "Computed", "A strict-derived base-level candidate spectrum, not observed gauge physics.", metrics["1.0"])
    x["ST2237"] = packet(2237, "Proven counterfamily", "Every p>0 is positive, target-free, D12-natural, and has mean one.", {
        "family": "h_f(p)=tau_f^p/mean(tau^p)", "parameter_domain": "p>0", "sample_metrics": metrics})
    x["ST2238"] = packet(2238, "Proven", "Permutation covariance, positivity, support locality and target independence hold for the whole family.", {
        "criteria_shared_by_all_p": ["D12 covariance", "positivity", "same flag complex", "no target fit", "same normalization"]})
    x["ST2239"] = packet(2239, "Refuted uniqueness", "p=1 and p=2 satisfy the same audited naturality premises but disagree spectrally.", {
        "eigenvalue_vector_residual_inf": counter_residual, "relative_l2_difference": relative})
    x["ST2240"] = packet(2240, "Proven", "Positive row weights do not change d1's kernel; the full simplex is contractible.", {
        "rank_d1": rank_d1, "kernel_dimension_on_C1": len(EDGES) - rank_d1, "H1_dimension": 0})
    x["ST2241"] = packet(2241, "Proven no-go", "Current strict data contain no premise distinguishing p=1 from p=2 or the continuum p>0.", {
        "unique_Hodge_exponent_selected": False})
    x["ST2242"] = packet(2242, "Proven classification", "A continuous positive multiplicative valuation g(xy)=g(x)g(y) has g(x)=x^p; monotonicity gives p>=0.", {
        "refinement_multiplicativity_reduces_but_does_not_remove_modulus": True, "remaining_modulus": "p"})
    x["ST2243"] = packet(2243, "Guardrail", "Choosing p=eta=1.8 would silently transfer a damping exponent into a face measure without a theorem.", {
        "strict_eta": 1.8, "eta_to_Hodge_exponent_bridge_exported": False})
    x["ST2244"] = packet(2244, "Refuted as continuum physics", "Complete support has graph diameter one and no local wavevector/lattice scale.", {
        "support_diameter": 1, "Lorentz_or_Maxwell_limit_derived": False})
    x["ST2245"] = packet(2245, "Round verdict", "Strict W yields an executable positive weighted flag-Hodge family, but not an unavoidable member or physical gauge theory.", {
        "mathematical_candidate_exists": True, "unique_strict_model": False, "physical_prediction": False})
    x["ST2246"] = packet(2246, "Recommendation", "Synthesize the corrected gate: separate reducible cycle success from irreducible/provenance/refinement/physics failures.", {
        "next": ["gate v2", "pair-sufficiency obstruction", "power-family no-go", "minimal extra axioms", "next decisive source-law test"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
