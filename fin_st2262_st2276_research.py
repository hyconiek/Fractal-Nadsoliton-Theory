#!/usr/bin/env python3
"""FIN ST2262--ST2276: does strict composition force separate linearity?"""
import numpy as np

from fin_st2262_st2351_common import strict_A_W, strict_distance_weights, write_packet, write_round


LO, HI = 2262, 2276
NAMES = [
    "CompositionLawInventory", "DirichletParallelAdditivity", "LiftNonImplication",
    "SeparateAdditivityTheorem", "ProductRuleTest", "GeometricMeanCounterrule",
    "HarmonicMeanCounterrule", "MinimumCounterrule", "HomogeneityConflict",
    "SupportZeroClassification", "ParallelAdditivityPremise", "StrictSourceAudit",
    "SeparateLinearityGate", "RoundOneVerdict", "RoundOneRecommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def defects(F, a, ap, b, c):
    return float(F(a + ap, b, c) - F(a, b, c) - F(ap, b, c))


def main():
    A, W = strict_A_W()
    w = strict_distance_weights()
    a, ap, b, c = map(float, (w[0], w[1], w[2], w[3]))
    product = lambda x, y, z: x * y * z
    geometric = lambda x, y, z: (x * y * z) ** (1 / 3)
    harmonic = lambda x, y, z: 3 / (1 / x + 1 / y + 1 / z)
    minimum = lambda x, y, z: min(x, y, z)
    rules = {"product": product, "geometric_mean": geometric,
             "harmonic_mean": harmonic, "minimum": minimum}
    add_defects = {name: defects(F, a, ap, b, c) for name, F in rules.items()}
    rng = np.random.default_rng(2262)
    f = rng.normal(size=12)
    W2 = 0.37 * W
    def energy(Q):
        return 0.5 * sum(Q[i, j] * (f[i] - f[j])**2 for i in range(12) for j in range(12))
    parallel_residual = abs(energy(W + W2) - energy(W) - energy(W2))

    x = {}
    x["ST2262"] = packet(2262, "Classified", "Current strict data expose edge energy, spectral calculus, adaptive potential and refinement; no face-composition law is exported.", {
        "candidate_compositions": ["parallel conductance addition", "series/harmonic composition", "path multiplication", "mean aggregation"]})
    x["ST2263"] = packet(2263, "Proven", "Dirichlet energy is linear in W for a fixed field f.", {
        "parallel_energy_additivity_residual": parallel_residual})
    x["ST2264"] = packet(2264, "Proven nonimplication", "Linearity of a scalar edge-energy sum does not determine a nonlinear map from three edge weights to one face weight.", {
        "edge_energy_linear": True, "face_rule_logically_determined": False})
    x["ST2265"] = packet(2265, "Conditional theorem", "Continuity rules out pathological Cauchy solutions.", {
        "theorem": "A continuous support-zero rule additive in each positive edge variable is F(a,b,c)=kappa*a*b*c."})
    x["ST2266"] = packet(2266, "Proven", "The product satisfies separate parallel additivity.", {
        "additivity_defect": add_defects["product"]})
    x["ST2267"] = packet(2267, "Counterexample", "Local, symmetric, positive, support-zero and degree-one, but not separately additive.", {
        "rule": "(abc)^(1/3)", "additivity_defect": add_defects["geometric_mean"]})
    x["ST2268"] = packet(2268, "Counterexample", "Natural series-conductance aggregation violates separate additivity.", {
        "rule": "3/(1/a+1/b+1/c)", "additivity_defect": add_defects["harmonic_mean"]})
    x["ST2269"] = packet(2269, "Counterexample", "Local, symmetric, monotone and support-zero, but nonsmooth/nonadditive.", {
        "rule": "min(a,b,c)", "additivity_defect": add_defects["minimum"]})
    x["ST2270"] = packet(2270, "Proven ambiguity", "Product has total degree three; geometric/harmonic/minimum have degree one.", {
        "homogeneity_degrees": {"product": 3, "geometric_mean": 1, "harmonic_mean": 1, "minimum": 1},
        "strict_face_dimension_exported": False})
    x["ST2271"] = packet(2271, "Proven", "Support-zero alone permits infinitely many symmetric positive rules.", {
        "examples": ["(abc)^p", "harmonic mean", "min", "abc/(a+b+c)^2"]})
    x["ST2272"] = packet(2272, "Identified missing premise", "A3 is equivalent to requiring the face lift to preserve parallel addition independently in every edge slot.", {
        "premise": "face-lift parallel additivity", "present_in_edge_energy": True, "present_in_face_lift": False})
    x["ST2273"] = packet(2273, "Repository audit", "No current strict packet exports face-lift parallel additivity as a source theorem.", {
        "strict_A3_source": False})
    x["ST2274"] = packet(2274, "Gate fail", "A3 selects the loop product conditionally but is not derived by the strict pair operator.", {
        "separate_linearity_strict": False, "conditional_selection": True})
    x["ST2275"] = packet(2275, "Round verdict", "Dirichlet linearity does not propagate automatically one categorical level upward to triangle weights.", {
        "A3_unavoidable": False})
    x["ST2276"] = packet(2276, "Recommendation", "Audit A4 using the exact weighted-graph automorphism group and triangle orbit separation.", {
        "next": ["Aut(W)", "D12 completeness", "triangle orbit count", "tau orbit separation", "S12 counterfactual"]})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
