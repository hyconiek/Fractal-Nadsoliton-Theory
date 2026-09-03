#!/usr/bin/env python3
"""FIN ST2247--ST2261: synthesis, axiom removal, gate v2, report trigger."""
import hashlib
import json
import math

import matplotlib.pyplot as plt
import numpy as np

from fin_st2172_st2261_common import (
    ROOT, TRIANGLES, cycle_field, hodge_spectrum, normalized_positive_power,
    strict_A_W, triangle_orbits, write_packet, write_round,
)


LO, HI = 2247, 2261
NAMES = [
    "LoopProductFlagLift", "ConditionalLoopProductUniqueness", "LocalityRemoval",
    "SupportZeroRemoval", "SeparateLinearityRemoval", "UniversalityRemoval",
    "SymmetryAxiomRedundancy", "PairSufficiencyObstruction", "TernarySourceGateV2",
    "GateScore", "MinimalAxiomLedger", "LegacyStrictBoundary", "PhysicsBoundary",
    "HighestValueNextTheorem", "FinalCycleReportTrigger",
]
FIGDIR = ROOT / "FIN_ST2172_ST2261_Figures"
FIG = FIGDIR / "ternary_source_and_flag_hodge_audit.png"
GATE = ROOT / "FIN_ST2247_ST2261_TernarySourceGateV2.json"
AXIOMS = ROOT / "FIN_ST2247_ST2261_LoopProductAxiomLedger.json"


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def legacy_cycle_stats():
    n = 12
    W = np.zeros((n, n))
    alpha, omega, phi, beta = 4 * math.log(2), math.pi / 4, math.pi / 6, 0.01
    for i in range(n):
        for j in range(i + 1, n):
            d = min((i - j) % n, (j - i) % n)
            W[i, j] = W[j, i] = alpha * math.cos(omega * d + phi) / (1 + beta * d)
    tau = cycle_field(W)
    return {"minimum": float(tau.min()), "maximum": float(tau.max()),
            "positive": int(np.sum(tau > 0)), "negative": int(np.sum(tau < 0))}


def make_figure(tau, spectra, gate):
    FIGDIR.mkdir(exist_ok=True)
    orbits = triangle_orbits()
    orbit_values = [np.mean([tau[TRIANGLES.index(t)] for t in o]) for o in orbits]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.2))
    ax[0, 0].bar(range(1, 13), orbit_values)
    ax[0, 0].set_yscale("log")
    ax[0, 0].set(title="Strict loop-product field", xlabel="$D_{12}$ triangle orbit", ylabel="$\\tau$ (log scale)")
    for p, ev in spectra.items():
        ax[0, 1].plot(range(1, len(ev) + 1), ev, label=f"p={p}")
    ax[0, 1].set(title="Hodge spectral counterfamily", xlabel="ordered mode", ylabel="eigenvalue")
    ax[0, 1].legend()
    eps = np.linspace(-1, 1, 101)
    ax[1, 0].plot(eps, np.zeros_like(eps), label="all pair moments")
    ax[1, 0].plot(eps, eps, label="third moment")
    ax[1, 0].set(title="Hidden-synergy obstruction", xlabel="$\\varepsilon$", ylabel="observable")
    ax[1, 0].legend()
    strict_math = sum(bool(r["strict_mathematical_pass"]) for r in gate["rows"])
    physics = sum(bool(r["physical_source_pass"]) for r in gate["rows"])
    ax[1, 1].bar(["strict math", "physical source"], [strict_math, physics], color=["#4c78a8", "#e45756"])
    ax[1, 1].set_ylim(0, len(gate["rows"]))
    ax[1, 1].set(title="Updated gate", ylabel=f"passes / {len(gate['rows'])}")
    for a in ax.flat:
        a.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG, dpi=180)
    plt.close(fig)


def main():
    _, W = strict_A_W()
    tau = cycle_field(W)
    spectra = {p: hodge_spectrum(normalized_positive_power(tau, p)) for p in (0.5, 1.0, 1.8, 2.0)}
    gate = {"schema": "FIN-TERNARY-SOURCE-GATE-v2", "rows": [
        {"id": "G1", "name": "strict reducible three-cycle observable", "strict_mathematical_pass": True, "physical_source_pass": False, "evidence": "ST2180"},
        {"id": "G2", "name": "D12 covariance and nonzero orbit witness", "strict_mathematical_pass": True, "physical_source_pass": False, "evidence": "ST2181-ST2182"},
        {"id": "G3", "name": "support-clique face complex", "strict_mathematical_pass": True, "physical_source_pass": False, "evidence": "ST2232-ST2233"},
        {"id": "G4", "name": "irreducible connected ternary source", "strict_mathematical_pass": False, "physical_source_pass": False, "evidence": "ST2175/ST2223"},
        {"id": "G5", "name": "nonpremise reflection-odd polarity", "strict_mathematical_pass": False, "physical_source_pass": False, "evidence": "ST2187-ST2200"},
        {"id": "G6", "name": "unique positive Hodge law", "strict_mathematical_pass": False, "physical_source_pass": False, "evidence": "ST2237-ST2242"},
        {"id": "G7", "name": "unique 12-24-48 ternary lift", "strict_mathematical_pass": False, "physical_source_pass": False, "evidence": "ST2225-ST2228"},
        {"id": "G8", "name": "dimensional local continuum scaling", "strict_mathematical_pass": False, "physical_source_pass": False, "evidence": "ST2244"},
        {"id": "G9", "name": "OA platform and independent record", "strict_mathematical_pass": False, "physical_source_pass": False, "evidence": "not exported"},
    ]}
    GATE.write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n")
    axioms = {"schema": "FIN-LOOP-PRODUCT-AXIOMS-v1", "axioms": [
        {"id": "L", "name": "triangle-edge locality", "necessary_for_uniqueness": True,
         "removal": "global invariant multipliers and nonlocal face rules enter"},
        {"id": "Z", "name": "vanish when any support edge vanishes", "necessary_for_uniqueness": True,
         "removal": "constant, linear, and pair terms enter"},
        {"id": "M", "name": "separate linearity in each edge coupling", "necessary_for_uniqueness": True,
         "removal": "the positive power family tau^p enters"},
        {"id": "U", "name": "one universal face rule rather than one coefficient per D12 orbit", "necessary_for_uniqueness": True,
         "removal": "twelve independent orbit coefficients enter"},
        {"id": "S", "name": "explicit edge-permutation symmetry", "necessary_for_uniqueness": False,
         "removal": "redundant once L,Z,M,U imply kappa*a*b*c"},
    ], "conditional_theorem": "L+Z+M+U imply F(a,b,c)=kappa*a*b*c; mean normalization fixes descriptive kappa."}
    AXIOMS.write_text(json.dumps(axioms, indent=2, sort_keys=True) + "\n")
    make_figure(tau, spectra, gate)
    legacy = legacy_cycle_stats()
    math_passes = sum(r["strict_mathematical_pass"] for r in gate["rows"])
    physical_passes = sum(r["physical_source_pass"] for r in gate["rows"])

    x = {}
    x["ST2247"] = packet(2247, "Constructed", "Named mathematical lift, not a physical-law claim.", {
        "object": "strict loop-product flag lift", "map": "W -> (Clique(supp W), tau_ijk=W_ij W_jk W_ki)"})
    x["ST2248"] = packet(2248, "Conditional theorem", "Requires L,Z,M,U; overall scale is descriptive until physical calibration.", {
        "theorem": "Every universal local separately-linear support-zero triangle rule is kappa*a*b*c."})
    x["ST2249"] = packet(2249, "Axiom necessary", "Removing locality permits multiplication by arbitrary global invariants such as Tr(W^2).", {
        "axiom": "L", "counterfamily_exists": True})
    x["ST2250"] = packet(2250, "Axiom necessary", "Removing support-zero permits symmetric lower-degree terms.", {
        "axiom": "Z", "counterfamily": "c0+c1(a+b+c)+c2(ab+bc+ca)+c3abc"})
    x["ST2251"] = packet(2251, "Axiom necessary", "Removing separate linearity admits tau^p and general positive functions.", {
        "axiom": "M", "counterfamily": "(abc)^p, p>0"})
    x["ST2252"] = packet(2252, "Axiom necessary", "D12 equivariance alone permits one multiplier for each of twelve triangle orbits.", {
        "axiom": "U", "counterfamily_dimension": 12})
    x["ST2253"] = packet(2253, "Axiom redundant", "The surviving monomial abc is already symmetric.", {
        "axiom": "S", "independent_necessity": False})
    x["ST2254"] = packet(2254, "Proven obstruction", "Identical complete pair records can have arbitrary epsilon third moment.", {
        "general_pair_sufficiency_for_ternary_information": False, "counterexample": "ST2217-ST2228"})
    x["ST2255"] = packet(2255, "Gate constructed", "Separates mathematical lifts from physical source claims.", {
        "file": GATE.name, "rows": len(gate["rows"])})
    x["ST2256"] = packet(2256, "Gate result", "Three mathematical rows pass; no physical-source row passes.", {
        "strict_mathematical_passes": math_passes, "physical_source_passes": physical_passes})
    x["ST2257"] = packet(2257, "Axiom ledger", "Only L,Z,M,U give conditional loop-product uniqueness; S is redundant.", {
        "file": AXIOMS.name, "necessary": ["L", "Z", "M", "U"], "redundant": ["S"]})
    x["ST2258"] = packet(2258, "Kernel-split boundary", "Legacy raw cycle products have both signs, so strict positivity cannot be silently transferred.", {
        "legacy_cycle_stats": legacy, "legacy_role_transfer": False})
    x["ST2259"] = packet(2259, "Physics no-go", "No unique Hodge rule, local continuum, units, polarity, state law, apparatus, or record follows.", {
        "fundamental_physics_closed": False, "ToE_closed": False})
    x["ST2260"] = packet(2260, "Recommendation", "Test whether an already strict adaptive/refinement law forces separate linearity M and universality U; otherwise classify them as minimal Sector/Hodge axioms.", {
        "highest_value_next_theorem": "derive or refute M+U from strict composition/refinement rather than assume p=1",
        "stop_rule": "do not enumerate more cycle functions unless they attack M or U"})
    x["ST2261"] = packet(2261, "Cycle complete", "Release 11.08 report required.", {
        "programs": 90, "rounds": 6, "release": "11.08",
        "main_result": "strict reducible loop-product flag lift exists; irreducible ternary reconstruction and unique Hodge physics fail",
        "gate_sha256": hashlib.sha256(GATE.read_bytes()).hexdigest(),
        "axiom_ledger_sha256": hashlib.sha256(AXIOMS.read_bytes()).hexdigest(),
        "figure": str(FIG.relative_to(ROOT)), "strict_ToE_closure": False})
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
