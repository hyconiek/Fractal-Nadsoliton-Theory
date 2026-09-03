#!/usr/bin/env python3
"""FIN ST2322--ST2336: minimal hidden-synergy instrument and extension."""
import math

from fin_st2262_st2351_common import write_packet, write_round


LO, HI = 2322, 2336
NAMES = [
    "ParityObservable", "ThreeBodyMinimality", "PairInstrumentFisherZero",
    "ParityLikelihood", "ParityFisherInformation", "EfficientEstimator",
    "FiniteSampleBound", "ThreeBodyExponentialFamily", "ThreeBodyGibbsHamiltonian",
    "AdaptiveLikelihoodFlow", "ZeroSourceInvariant", "PolarityBoundary",
    "OperationalRecordSpec", "RoundFiveVerdict", "RoundFiveRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    epsilon0, delta = 0.75, 0.01
    n_bound = math.ceil(2 * math.log(2 / delta) / epsilon0**2)
    fisher = 1 / (1 - epsilon0**2)
    variance = (1 - epsilon0**2) / n_bound

    x = {}
    x["ST2322"] = packet(2322, "Constructed", "Dimensionless three-body observable on the exact hidden-synergy family.", {
        "observable": "Y=x1*x2*x3", "expectation": "epsilon", "outcomes": [-1, 1]})
    x["ST2323"] = packet(2323, "Proven", "Fourier characters on the binary cube are orthogonal; parity is not a sum of one/two-body functions.", {
        "minimum_interaction_order": 3})
    x["ST2324"] = packet(2324, "Proven", "All pair marginals are epsilon-independent, so their score and Fisher information vanish.", {
        "pair_only_Fisher_information": 0})
    x["ST2325"] = packet(2325, "Proven", "Y is Bernoulli on +-1 with P(Y=y)=(1+epsilon*y)/2.", {
        "likelihood": "(1+epsilon*y)/2"})
    x["ST2326"] = packet(2326, "Proven", "Direct score expectation.", {
        "single_event_Fisher": "1/(1-epsilon^2)", "at_epsilon_0": fisher})
    x["ST2327"] = packet(2327, "Proven", "The sample mean of parity is unbiased and reaches the Cramer-Rao variance.", {
        "estimator": "epsilon_hat=mean(Y)", "variance": "(1-epsilon^2)/n", "at_n_bound": variance})
    x["ST2328"] = packet(2328, "Proven bound", "Hoeffding two-sided bound for deciding the sign when |epsilon|>=epsilon0.", {
        "epsilon0": epsilon0, "delta": delta, "sufficient_events": n_bound})
    x["ST2329"] = packet(2329, "Constructed", "Minimal non-Gaussian extension carrying exactly the missing sufficient statistic.", {
        "family": "p_theta(x)=exp(theta*x1*x2*x3)/(8*cosh(theta))", "epsilon": "tanh(theta)"})
    x["ST2330"] = packet(2330, "Conditional physical form", "A three-body Gibbs interaction realizes theta=beta*J, but beta and J are added dimensional/thermodynamic data.", {
        "Hamiltonian": "H=-J*x1*x2*x3", "strict_J_source": False})
    x["ST2331"] = packet(2331, "Constructed", "Expected likelihood gradient has a unique stable stationary point for a supplied target epsilon_star.", {
        "flow": "dot(theta)=epsilon_star-tanh(theta)", "stationary": "theta=artanh(epsilon_star)"})
    x["ST2332"] = packet(2332, "Proven", "With epsilon_star=0 and theta(0)=0 the adaptive flow remains at the pair/Gaussian baseline.", {
        "spontaneous_nonzero_source": False})
    x["ST2333"] = packet(2333, "Boundary", "The sign of theta/epsilon is the same unsourced polarity problem in a minimal exact model.", {
        "QW2191_discharged": False, "odd_source_needed": True})
    x["ST2334"] = packet(2334, "Executable specification", "Raw records must retain all three outcomes from one jointly prepared trial.", {
        "record_fields": ["run_id", "x1", "x2", "x3", "parity", "timestamp", "configuration"],
        "pairwise_histograms_sufficient": False, "physical_platform_supplied": False})
    x["ST2335"] = packet(2335, "Round verdict", "A minimal ternary instrument and model exist, but both explicitly add the three-body statistic/source absent from strict FIN.", {
        "mathematical_bridge": True, "strict_source_bridge": False, "laboratory_evidence": False})
    x["ST2336"] = packet(2336, "Recommendation", "Synthesize the two conditional routes and issue a no-new-live-frontier result for the current pair/refinement core.", {
        "next": ["gate v3", "minimal principle lattice", "strict/conditional split", "stop rules", "next source class"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
