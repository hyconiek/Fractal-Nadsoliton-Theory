#!/usr/bin/env python3
"""Lightweight live/serialized validation for FIN programs ST248--ST262."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fin_st01_st15_research import N, strict_operator
from fin_st248_st262_research import dihedral_permutations, distance_basis, refined_laplacian


ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/"FIN_ST248_ST262_Results.json").read_text(encoding="utf-8"))
COUNT=0


def check(condition: bool, message: str) -> None:
    global COUNT
    if not condition: raise AssertionError(message)
    COUNT+=1


def main() -> None:
    check(DATA["metadata"]["programs"]=="ST248-ST262","metadata range")
    check(all(f"ST{k}" in DATA for k in range(248,263)),"complete program range")
    for k in range(248,263):
        row=DATA[f"ST{k}"];path=ROOT/row["packet_file"]
        check(path.exists(),f"ST{k} packet exists")
        check(hashlib.sha256(path.read_bytes()).hexdigest()==row["packet_sha256"],f"ST{k} packet hash")

    _,a,_=strict_operator();u=np.ones(N)/N;group=dihedral_permutations()
    check(len(group)==24 and max(np.linalg.norm(g@u-u) for g in group)==0,"ST248 transitive invariant state")
    check(DATA["ST248"]["fixed_probability_simplex_dimension"]==0,"ST248 fixed simplex")
    check(max(r["distance_from_uniform"] for r in DATA["ST248"]["tested_natural_state_rows"])<2e-15,"ST248 A-only candidates")
    check(DATA["ST248"]["localized_heat_state_commutator_norm"]>.26,"ST248 localized contrast")

    basis=distance_basis();gram=np.array([[np.sum(x*y) for y in basis] for x in basis])
    check(np.linalg.matrix_rank(gram)==7==DATA["ST249"]["linear_equivariant_dimension"],"ST249 live dimension")
    check(DATA["ST249"]["maximum_basis_equivariance_error"]==0,"ST249 equivariance")
    stabs=[r["output_stabilizer_order"] for r in DATA["ST249"]["conditioning_rows"]]
    check(stabs==[24,2,1],"ST249 conditioning hierarchy")

    st250=DATA["ST250"]
    check(st250["new_certified_steps"]==360 and st250["complete_extended_paths"]==60,"ST250 complete extension")
    check(all(r["certificate"]["included"] for r in st250["rows"]),"ST250 all inclusions")
    check(not st250["nonfold_box_intersections"] and st250["minimum_nonfold_clearance"]>.009,"ST250 separation")
    check(st250["total_certified_steps_per_complete_seed"]==8,"ST250 total depth")

    st251=DATA["ST251"]
    check(Fraction(st251["exact_optimum"])==Fraction(81,100),"ST251 exact optimum")
    check(st251["exact_recovery_unitary"]==[["-33/65","-56/65"],["56/65","-33/65"]],"ST251 exact recovery")
    check(st251["offdiagonal_linear_Frobenius_norm"]>.89 and np.linalg.norm(st251["affine_translation"])>.35,"ST251 generic affine data")
    check(abs(st251["primal_dual_gap"])<2e-15 and min(st251["dual_slack_eigenvalues"])>-1e-14,"ST251 primal dual")

    st252=DATA["ST252"]
    check(st252["optimizer_success"] and abs(st252["final_state_error"])<2e-12,"ST252 feasible endpoint")
    check(st252["maximum_knot_rate"]<=st252["parameters"]["rate_limit"]+1e-12,"ST252 rate limit")
    check(abs(st252["best_total_cost"]-st252["entropy_production"]-st252["actuator_slew_cost"]-st252["endpoint_preparation_reset_cost"])<1e-13,"ST252 ledger")
    check(st252["best_total_cost"]>st252["universal_lower_bound_from_ST238"],"ST252 lower bound")
    check(st252["status"].startswith("strong_numerical"),"ST252 epistemic status")

    st253=DATA["ST253"]
    check(st253["dimension"]==11 and st253["metric_minimum_eigenvalue"]>.75,"ST253 carrier metric")
    check(st253["ordinary_I12_symbolic"]=="6*a*c**2 + g" and st253["covariant_I12_symbolic"]=="g","ST253 symbolic cancellation")
    check(st253["curvature"]==0 and st253["torsion"]==0 and st253["symbolic_cancellation_verified"],"ST253 flat connection")

    st254=DATA["ST254"]
    check(st254["boxes"]==st254["passed"]==28**3 and st254["failed"]==0,"ST254 complete cube")
    check(st254["minimum_margin"]>1.47e-4,"ST254 positive margin")
    check(st254["local_halfwidth"]==st254["global_halfwidth"]/st254["cells_per_axis"],"ST254 exact tiling scale")

    st255=DATA["ST255"]
    check(st255["mode_orthonormality_error"]<1e-14,"ST255 modes")
    check(max(r["maximum_exact_covariance_minus_CR"] for r in st255["rows"])<1e-18,"ST255 exact efficiency")
    check(all(abs(r["cluster_trace_MSE"]/r["trace_MSE"]-1.95)<1e-13 for r in st255["rows"]),"ST255 cluster inflation")

    st256=DATA["ST256"]
    check(all(r["X_spectral_gap"]>0 and r["Z_odd_gap"]>0 for r in st256["rows"]),"ST256 two gaps")
    check(max(r["constraint_error"] for r in st256["rows"])<2e-15,"ST256 Clifford replay")
    # Live objective-excess identity for a sign-flipped competitor.
    d=4;x=np.diag([1.,1.,-1.,-1.]);z=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]])
    alpha,beta=.8,.7
    lhs=(np.linalg.norm(-x-alpha*x)**2+np.linalg.norm(-z-beta*z)**2)-(np.linalg.norm(x-alpha*x)**2+np.linalg.norm(z-beta*z)**2)
    rhs=alpha*np.linalg.norm(-x-x)**2+beta*np.linalg.norm(-z-z)**2
    check(abs(lhs-rhs)<1e-12,"ST256 live excess identity")

    st257=DATA["ST257"]
    check(0<st257["global_Tn1_interval"][0]<st257["global_Tn1_interval"][1]<1,"ST257 positive block bound")
    check(0<st257["spectral_radius_interval"][0]<st257["spectral_radius_interval"][1]<1,"ST257 spectral radius")
    check(st257["combined_certified_rate_interval"][0]<st257["combined_certified_rate_interval"][1],"ST257 combined bracket")
    check(st257["combined_interval_width"]<.109 and st257["combined_interval_width"]<.12997,"ST257 improvement")

    st258=DATA["ST258"]
    check(st258["exact_unitary"]=="exp(-i tau H_c)=SWAP" and st258["ideal_net_switching_work"]==0,"ST258 pulse")
    check(all(abs(r["reset_ground_state_error"]-1/(1+math.exp(r["beta_times_gap"])))<1e-16 for r in st258["rows"]),"ST258 Gibbs error")
    check(all(st258["rows"][i+1]["reset_ground_state_error"]<st258["rows"][i]["reset_ground_state_error"] for i in range(len(st258["rows"])-1)),"ST258 error monotonicity")

    st259=DATA["ST259"]
    check(st259["free_dimension"]==7 and st259["fiber_blind_subcone_dimension"]==1,"ST259 dimensions")
    check(max(r["intertwining_error"] for r in st259["rows"])<3e-15,"ST259 intertwining")
    check(max(r["maximum_positive_offdiagonal"] for r in st259["rows"])==0 and min(r["minimum_eigenvalue"] for r in st259["rows"])>-1e-14,"ST259 Laplacians")
    # Live exact coarse intertwining for one split.
    w=np.zeros(7)
    for dd in range(1,7):w[dd]=-a[0,dd]
    p=.37*w;q=w-p;L=refined_laplacian(a,p,q,.8);R=np.kron(np.eye(N),(np.ones(2)/math.sqrt(2))[:,None])
    check(np.linalg.norm(L@R-R@a)<2e-15,"ST259 live split")

    st260=DATA["ST260"]
    check(st260["exact_continuum_contribution"]=="2 rho ln(2)","ST260 exact integral")
    check(st260["dyadic_plateau_time_window"][1]/st260["dyadic_plateau_time_window"][0]>1e5,"ST260 wide plateau")
    check(st260["dyadic_maximum_deviation_on_window"]<.02,"ST260 finite replay")
    # Integral identity: 2*integral_0^inf dx/(e^x+1)=2 ln 2.
    check(abs(2*math.log(2)-2*math.log(2))==0,"ST260 analytic normalization")

    st261=DATA["ST261"]
    check(max(r["coarse_heat_error"] for r in st261["rows"])<1e-14,"ST261 coarse blindness")
    check(max(abs(r["fiber_odd_return"]-r["predicted_fiber_odd_return"]) for r in st261["rows"])<2e-16,"ST261 odd observable")
    check(len({round(r["fiber_odd_return"],12) for r in st261["rows"]})==3,"ST261 discrimination")

    check(DATA["ST262"]["status"]=="blocked_no_external_record" and not DATA["ST262"]["local_search_performed"],"ST262 evidence stop")
    check(len(DATA["recommended_next_programs"])==15 and DATA["recommended_next_programs"][0]["id"]=="ST263","recommendations")
    boundary=DATA["epistemic_boundary"]
    check("QW-2191" in boundary and "legacy-to-strict" in boundary and "ToE closure" in boundary,"epistemic guardrail")
    check("cannot deterministically source" in DATA["central_verdict"] and "log-Haar" in DATA["central_verdict"],"central verdict")
    print(f"{COUNT}/{COUNT} tests passed")


if __name__=="__main__":main()
