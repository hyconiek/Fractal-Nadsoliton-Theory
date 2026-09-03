#!/usr/bin/env python3
"""FIN ST2337--ST2351: synthesis, gate v3 and Release 11.09 trigger."""
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np

from fin_st2262_st2351_common import (
    ROOT, entropy_of_mean_one_weights, mean_normalize, orbit_cycle_vector,
    power_map, write_packet, write_round,
)


LO, HI = 2337, 2351
NAMES = [
    "A3A4JointFailure", "ConditionalRouteOne", "ConditionalRouteTwo",
    "PowerSelfSimilarityNoGo", "TernaryInstrumentBoundary", "MinimalPrincipleLattice",
    "UpdatedGateV3", "GateV3Score", "KernelSplitBoundary", "PhysicsDecision",
    "StopRules", "HighestValueNextObject", "RecommendedNextPrograms",
    "FinalInterpretation", "CycleSixReportTrigger",
]
FIGDIR = ROOT / "FIN_ST2262_ST2351_Figures"
FIG = FIGDIR / "composition_refinement_and_synergy_audit.png"
GATE = ROOT / "FIN_ST2337_ST2351_HigherOrderRefinementGateV3.json"
LATTICE = ROOT / "FIN_ST2337_ST2351_MinimalPrincipleLattice.json"


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def make_figure(profile, gate):
    FIGDIR.mkdir(exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.2))
    names = ["product", "geometric", "harmonic", "minimum"]
    defects = [0.0, -0.0785221321, -0.0785966273, -0.0470291687]
    ax[0, 0].bar(names, np.abs(defects)); ax[0, 0].tick_params(axis='x', rotation=18)
    ax[0, 0].set(title="Separate-additivity test", ylabel="absolute defect")
    for p in (0.5, 1.0, 1.8, 2.0):
        ax[0, 1].plot(power_map(profile, p), marker='o', label=f"p={p}")
    ax[0, 1].set(title="Normalized orbit profiles", xlabel="orbit", ylabel="weight"); ax[0, 1].legend()
    ps = np.linspace(0, 2, 101)
    resid = [np.linalg.norm(power_map(power_map(profile, p), p)-power_map(profile, p), np.inf) for p in ps]
    ax[1, 0].plot(ps, resid); ax[1, 0].set(title="Path-independence residual", xlabel="p", ylabel="infinity norm")
    vals = [sum(r[key] for r in gate["rows"]) for key in ("mathematical_pass", "strict_source_pass", "physical_pass")]
    ax[1, 1].bar(["math", "strict source", "physical"], vals, color=["#4c78a8", "#f2cf5b", "#e45756"])
    ax[1, 1].set_ylim(0, len(gate["rows"])); ax[1, 1].set(title="Gate v3", ylabel=f"passes / {len(gate['rows'])}")
    for a in ax.flat: a.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(FIG, dpi=180); plt.close(fig)


def main():
    profile = mean_normalize(orbit_cycle_vector())
    gate = {"schema": "FIN-HIGHER-ORDER-REFINEMENT-GATE-v3", "rows": [
        {"id":"R1","name":"strict reducible loop-product field","mathematical_pass":True,"strict_source_pass":True,"physical_pass":False},
        {"id":"R2","name":"exact Aut(W)=D12 and twelve face orbits","mathematical_pass":True,"strict_source_pass":True,"physical_pass":False},
        {"id":"R3","name":"support-clique complex","mathematical_pass":True,"strict_source_pass":True,"physical_pass":False},
        {"id":"R4","name":"strict separate edge linearity","mathematical_pass":True,"strict_source_pass":False,"physical_pass":False},
        {"id":"R5","name":"strict face universality","mathematical_pass":True,"strict_source_pass":False,"physical_pass":False},
        {"id":"R6","name":"idempotent p=1 selection theorem","mathematical_pass":True,"strict_source_pass":False,"physical_pass":False},
        {"id":"R7","name":"nontrivial exact self-similar refinement","mathematical_pass":False,"strict_source_pass":False,"physical_pass":False},
        {"id":"R8","name":"irreducible ternary state source","mathematical_pass":False,"strict_source_pass":False,"physical_pass":False},
        {"id":"R9","name":"minimal parity instrument","mathematical_pass":True,"strict_source_pass":False,"physical_pass":False},
        {"id":"R10","name":"nonpremise odd polarity","mathematical_pass":False,"strict_source_pass":False,"physical_pass":False},
        {"id":"R11","name":"dimensional continuum OA realization","mathematical_pass":False,"strict_source_pass":False,"physical_pass":False}]}
    GATE.write_text(json.dumps(gate, indent=2, sort_keys=True)+"\n")
    lattice = {"schema":"FIN-MINIMAL-PRINCIPLE-LATTICE-v1","routes":[
        {"route":"unique loop Hodge by multilinearity","premises":["A3 separate edge additivity","A4 face universality"],"strict":False},
        {"route":"unique p=1 by refinement","premises":["positive multiplicative valuation","scale-stationary idempotence","nonconstant profile","information preservation"],"strict":False,"result":"identity, not nontrivial compression"},
        {"route":"irreducible ternary state","premises":["three-body sufficient statistic/source","joint parity instrument","polarity law"],"strict":False},
        {"route":"physical theory","premises":["one prior route","CA units","local continuum law","OA apparatus and independent record"],"strict":False}]}
    LATTICE.write_text(json.dumps(lattice, indent=2, sort_keys=True)+"\n")
    make_figure(profile, gate)
    scores={key:sum(r[key] for r in gate["rows"]) for key in ("mathematical_pass","strict_source_pass","physical_pass")}

    x={}
    x["ST2337"]=packet(2337,"Proven synthesis","A3 and A4 are compatible additions but neither follows from strict edge energy or Aut(W).",{"A3_strict":False,"A4_strict":False})
    x["ST2338"]=packet(2338,"Conditional route","A3+A4 select tau up to scale.",{"premises":["separate edge additivity","face universality"],"strict":False})
    x["ST2339"]=packet(2339,"Conditional route","Multiplicativity plus scale-stationary idempotence and information preservation select p=1.",{"selected_p":1,"strict":False,"nontrivial_compression":False})
    x["ST2340"]=packet(2340,"Proven no-go","On the nonconstant finite orbit profile, exact power self-similarity up to relabelling permits only identity p=1; erasure p=0 is the other idempotent boundary.",{"nontrivial_exact_power_fractal":False})
    x["ST2341"]=packet(2341,"Constructed boundary","Parity is the minimal sufficient instrument, but the three-body state law and physical platform are added.",{"sufficient_events_epsilon075_delta001":19,"strict_source":False})
    x["ST2342"]=packet(2342,"Constructed","Separates mathematical lift, unique-model premises, ternary state and physical realization.",{"file":LATTICE.name,"routes":4})
    x["ST2343"]=packet(2343,"Gate constructed","Version 3 distinguishes mathematical constructions, strict sources and physical passes.",{"file":GATE.name,"rows":len(gate['rows'])})
    x["ST2344"]=packet(2344,"Gate result","Only the original graph lift is strict-sourced; no physical row passes.",scores)
    x["ST2345"]=packet(2345,"Kernel-split boundary","Signed legacy loops do not admit real noninteger power refinement without an extra abs/sign/completion rule.",{"legacy_transfer":False})
    x["ST2346"]=packet(2346,"Physics verdict","Neither conditional route generates units, locality, Lorentz/Maxwell structure, apparatus or evidence.",{"fundamental_physics":False,"ToE":False})
    x["ST2347"]=packet(2347,"Stop rules","Freeze cycle-function enumeration and plain coassociativity replay.",{"stop":["more positive face functions","p=eta transfer","S12 universality import","plain coassociativity as scale selector","pair records as ternary evidence"]})
    x["ST2348"]=packet(2348,"Highest-value next object","A genuinely strict non-Gaussian three-body sufficient statistic/state law with fixed polarity and refinement transport is the only object that changes the source verdict.",{"required":["formula","strict provenance","nonzero connected third cumulant","polarity","12-24-48 law"]})
    x["ST2349"]=packet(2349,"Recommended programmes","ST2352--ST2366 should audit one strict nonlinear state law, not another receiver family.",{"range":"ST2352-ST2366","count":15})
    x["ST2350"]=packet(2350,"Final interpretation","Fractal self-similarity of the finite loop profile collapses to identity or erasure; nontrivial hierarchy needs unsourced level data.",{"strict_pair_refinement_frontier":"bounded no-go"})
    x["ST2351"]=packet(2351,"Cycle complete","Release 11.09 report required.",{"programs":90,"rounds":6,"release":"11.09","main_result":"A3/A4 not strict; exact power self-similarity gives identity or erasure; minimal parity bridge is conditional","strict_ToE_closure":False,"gate_sha256":hashlib.sha256(GATE.read_bytes()).hexdigest(),"lattice_sha256":hashlib.sha256(LATTICE.read_bytes()).hexdigest(),"figure":str(FIG.relative_to(ROOT))})
    write_round(LO,HI,x)


if __name__=="__main__":main()
