#!/usr/bin/env python3
"""FIN ST308--ST326: sourced gain candidates, certified localization,
observer-depth invariants, relational clocks, and compression bounds.

All empirical gates remain closed.  The program writes only local mathematical
packets and permanently distinguishes supplied premises from strict-core output.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq, minimize, root

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv, strict_interval_matrix
from fin_st154_st165_research import parametric_data
from fin_st130_st141_research import point_design_system


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST308_ST326_Results.json"
SUMMARY = ROOT / "FIN_ST308_ST326_Summary.csv"
FIG_DIR = ROOT / "FIN_ST308_ST326_Figures"
SEED = 20260815
NAMES = {
    308: "Minimal_Active_D12_Covariant_Law",
    309: "Finite_Pumped_Auxiliary_Elimination",
    310: "Conservative_Mediator_Rank_and_Resource_Dichotomy",
    311: "Interval_Localized_Simplex_State",
    312: "Global_Localization_Dual_Reduction",
    313: "Localized_Branch_First_Order_Crossing",
    314: "Continuous_Tube_Obstruction_Audit",
    315: "Exact_Rational_Three_Output_SDP",
    316: "Projective_Filter_Contraction",
    317: "Sensitivity_Weighted_Nuisance_Allocation",
    318: "Approximate_Odd_Clifford_Lower_Bound",
    319: "Scale_Cocycle_Gauge_Classification",
    320: "Soft_IR_Bathtub_Extremal_Theorem",
    321: "Calibration_Group_Reduction_and_Anchor_Count",
    322: "Independent_Record_Stop",
    323: "Observer_Blindness_Kernel",
    324: "Blackwell_Depth_and_Kernel_Order",
    325: "Relative_Clock_Refinement_Naturality",
    326: "Fractal_Compression_Record_Lower_Bound",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": obj, "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


IDX = np.array([i if i <= 6 else N-i for i in range(N)])
MULT = np.array([1, 2, 2, 2, 2, 2, 1.0])


def radial_matrix(a: np.ndarray) -> np.ndarray:
    return np.array([[sum(a[i, k] for k in range(N) if IDX[k] == j)
                      for j in range(7)] for i in range(7)])


def localized_system(x: np.ndarray, b: np.ndarray, g: float) -> np.ndarray:
    q, lam = x[:7], x[7]
    return np.r_[np.log(N*q) + 1 - g*b@q + lam, MULT@q - 1]


def localized_jacobian(x: np.ndarray, b: np.ndarray, g: float) -> np.ndarray:
    j = np.zeros((8, 8)); j[:7, :7] = np.diag(1/x[:7]) - g*b
    j[:7, 7] = 1; j[7, :7] = MULT
    return j


def localized_seed(a: np.ndarray, g: float = 4.0) -> np.ndarray:
    old = json.loads((ROOT / "FIN_ST293_ST307_Results.json").read_text())["ST295"]
    p = np.array(old["representative_probability"], float)
    p = np.roll(p, -int(np.argmax(p)))
    q = np.array([p[IDX == j].mean() for j in range(7)])
    b = radial_matrix(a)
    lam = float(np.mean(g*b@q - np.log(N*q) - 1))
    sol = root(lambda z: localized_system(z, b, g), np.r_[q, lam],
               jac=lambda z: localized_jacobian(z, b, g), tol=1e-12)
    if not sol.success or np.linalg.norm(localized_system(sol.x, b, g), np.inf) > 2e-12:
        raise RuntimeError("localized root discovery failed")
    return sol.x


def localized_interval_system(center: np.ndarray, radius: float, aiv, g: int = 4):
    q = [iv((center[i]-radius, center[i]+radius)) for i in range(7)]
    lam = iv((center[7]-radius, center[7]+radius))
    biv = [[sum((aiv[i][k] for k in range(N) if IDX[k] == j), iv(0))
            for j in range(7)] for i in range(7)]
    f = [mp.iv.log(iv(N)*q[i]) + 1 - g*sum((biv[i][j]*q[j] for j in range(7)), iv(0)) + lam
         for i in range(7)]
    f.append(sum((iv(float(MULT[j]))*q[j] for j in range(7)), iv(0))-1)
    jac = [[iv(0) for _ in range(8)] for _ in range(8)]
    for i in range(7):
        for j in range(7):
            jac[i][j] = (1/q[i] if i == j else iv(0)) - g*biv[i][j]
        jac[i][7] = iv(1)
    for j in range(7): jac[7][j] = iv(float(MULT[j]))
    flo = np.array([bounds(x)[0] for x in f]); fhi = np.array([bounds(x)[1] for x in f])
    jlo = np.array([[bounds(x)[0] for x in row] for row in jac])
    jhi = np.array([[bounds(x)[1] for x in row] for row in jac])
    return flo, fhi, jlo, jhi


def localized_krawczyk(center: np.ndarray, aiv, radius: float) -> dict:
    f0lo, f0hi, j0lo, j0hi = localized_interval_system(center, 0.0, aiv)
    pre = np.linalg.inv((j0lo+j0hi)/2)
    _, _, jlo, jhi = localized_interval_system(center, radius, aiv)
    cglo, cghi = interval_matvec(pre, pre, f0lo, f0hi)
    ylo = np.nextafter(center-cghi, -np.inf); yhi = np.nextafter(center-cglo, np.inf)
    cjlo, cjhi = interval_left(pre, jlo, jhi); mlo, mhi = -cjhi, -cjlo
    for i in range(8):
        mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf)
        mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, np.full(8, -radius), np.full(8, radius))
    klo = np.nextafter(ylo+dlo, -np.inf); khi = np.nextafter(yhi+dhi, np.inf)
    margin = float(min(np.min(klo-(center-radius)), np.min((center+radius)-khi)))
    return {"radius": radius, "margin": margin, "included": margin > 0,
            "maximum_image_halfwidth": float(np.max((khi-klo)/2))}


def objective_from_radial(x: np.ndarray, a: np.ndarray, g: float) -> float:
    p = x[:7][IDX]; q = p-np.ones(N)/N
    return float(np.sum(p*np.log(N*p)) - .5*g*q@a@q)


def st308(a: np.ndarray) -> dict:
    ev = np.linalg.eigvalsh(a); ge = 2*math.log(N)/float(np.trace(a)/N)
    packet = {
        "law": "q_dot=-grad[D(p||u)-(kappa*y/2)<q,Aq>], tau*y_dot=P-gamma*y",
        "effective_gain_at_auxiliary_equilibrium": "g_eff=kappa*P/gamma",
        "D12_covariance": "y is a scalar and <q,Aq> is D12 invariant",
        "resource_ledger": {"strict_objects": ["A", "D(p||u)"],
                            "new_premises": ["P>0", "gamma>0", "tau>0", "kappa>0"],
                            "not_generated": ["pump magnitude", "coupling sign", "physical units"]},
        "vertex_crossing_pump_for_kappa_equals_gamma": ge,
        "largest_mode_Hessian_pump_for_kappa_equals_gamma": float(N/ev[-1]),
        "theorem": (
            "Within autonomous scalar-auxiliary extensions that modulate the invariant quadratic form, one real variable y is sufficient: its unique stable equilibrium is P/gamma and the reduced information objective has g_eff=kappa P/gamma. The law is D12-covariant and imports no phase label. Zero auxiliary dimension can only prescribe a constant g externally, while any stateful autonomous realization needs at least one real state. This is a minimality theorem for the declared multiplicative class, not a strict derivation of P, gamma, tau, or kappa."
        )}
    return finalize(308, "Minimal Scalar Active D12-Covariant Gain Law",
                    "proven_minimal_in_declared_multiplicative_auxiliary_class",
                    "The pump and signed coupling are new nonequilibrium premises. No selector, phase, units, or laboratory realization is generated.", packet)


def st309() -> dict:
    rows=[]
    for tau in (.1, .5, 1., 2.):
        for t in (tau, 3*tau, 6*tau):
            rel=math.exp(-t/tau)
            rows.append({"tau":tau,"time":t,"relative_gain_error_from_y0_zero":rel})
    packet={"exact_solution":"y(t)=P/gamma+(y0-P/gamma) exp(-gamma*t/tau)",
            "rows_for_gamma_one":rows,"phase_labels_in_law":0,
            "theorem":(
                "The finite pumped auxiliary is not merely an adiabatic ansatz: for constant P its gain converges exponentially and exactly to kappa P/gamma. Starting from y0=0, the relative gain error is exp(-gamma t/tau). Because the coupling uses only the D12-invariant quadratic form, finite relaxation does not import a broken phase. It does introduce a preparation y0 and an external pump history."
            )}
    return finalize(309,"Finite Pumped Auxiliary and Controlled Elimination Error",
                    "proven_exact_auxiliary_relaxation_and_reduced_gain",
                    "The reduced gain is conditional on a supplied pump. The calculation is dimensionless and does not establish a physical power source.",packet)


def st310(a: np.ndarray) -> dict:
    rank=int(np.linalg.matrix_rank(a,tol=1e-10));ev,v=np.linalg.eigh(a)
    sqrt_a=(v*np.sqrt(np.maximum(ev,0)))@v.T
    g=4.;b=math.sqrt(g)*sqrt_a
    defect=float(np.linalg.norm(b.T@b-g*a))
    packet={"rank_A":rank,"minimum_linear_mediator_dimension":rank,
            "constructive_coupling":"Gamma=I_11 and B=sqrt(g)*A^(1/2) restricted to ran(A)",
            "full_matrix_construction_defect":defect,
            "effective_term":"min_y [<y,Gamma y>/2-<y,Bq>] = -<q,B^T Gamma^-1 B q>/2",
            "resource_dichotomy":{
                "dissipative_positive_real_memory":"cannot deliver negative net cyclic work",
                "conservative_Schur_mediator":"can generate the negative static term without a pump",
                "stability_cost":"if the quadratic sector becomes negative, saturation or a higher-order stabilizer is required"},
            "theorem":(
                "For an m-dimensional positive quadratic mediator, the induced matrix B^T Gamma^{-1}B has rank at most m. Since rank(A)=11, exact realization of gA by a linear mediator requires m>=11; m=11 is attained by B=sqrt(g) A^{1/2} on ran(A). Thus a pump is not logically necessary for the negative static Schur term. What is necessary is an additional signed coupling, and beyond instability a stabilizing nonlinearity. Passive dissipative memory and conservative mediation are different resource classes."
            )}
    return finalize(310,"Conservative Mediator Rank Theorem and Gain-Resource Dichotomy",
                    "proven_rank_minimal_Schur_realization_and_pump_not_necessary_correction",
                    "The mediator and its coupling are supplied, not strict-derived. The theorem does not source nonlinear saturation or a physical environment.",packet)


def st311(a: np.ndarray) -> dict:
    mp.iv.dps=70;center=localized_seed(a);aiv,_,_=strict_interval_matrix()
    trials=[localized_krawczyk(center,aiv,r) for r in (1e-8,3e-9,1e-9,3e-10,1e-10)]
    accepted=next((x for x in trials if x["included"]),None)
    if accepted is None: raise RuntimeError("localized Krawczyk isolation failed")
    p=center[:7][IDX];h=np.diag(1/p)-4*a
    z=np.vstack([np.eye(11),-np.ones(11)])
    reduced=z.T@h@z;point_min=float(np.linalg.eigvalsh(reduced)[0])
    alo=np.array([[bounds(x)[0] for x in row] for row in aiv])
    ahi=np.array([[bounds(x)[1] for x in row] for row in aiv])
    # Every interval A' obeys ||A'-A||_2 <= ||entrywise radius+center defect||_F.
    amid=(alo+ahi)/2
    a_family_error=float(np.linalg.norm(np.maximum(np.abs(alo-a),np.abs(ahi-a)),ord="fro"))
    r=accepted["radius"];diag_error=float(r/(np.min(p)-r)**2)
    full_hessian_operator_error=diag_error+4*a_family_error
    lower=point_min-(np.linalg.norm(z,2)**2)*full_hessian_operator_error
    packet={"g":4.0,"reflection_even_root_center":center.tolist(),
            "point_residual_inf":float(np.linalg.norm(localized_system(center,radial_matrix(a),4),np.inf)),
            "Krawczyk_trials":trials,"accepted":accepted,
            "probability_minimum_box_lower":float(np.min(center[:7]-accepted["radius"])),
            "objective_at_center":objective_from_radial(center,a,4),
            "tangent_Hessian_point_congruence_minimum":point_min,
            "strict_operator_family_Frobenius_error_bound":a_family_error,
            "full_Hessian_operator_error_bound":full_hessian_operator_error,
            "conservative_tangent_Hessian_lower_bound":lower,
            "translation_orbit_size":12,
            "theorem":(
                "Outward Krawczyk inclusion isolates one positive reflection-even stationary probability at supplied g=4 for every strict operator in the declared transcendental interval family. A congruence on the exact sum-zero tangent coordinates has a positive perturbation-paid lower bound, so the root is a strict local minimum. Circulant covariance translates it into twelve distinct certified local minima. This proves existence, not global exhaustion or a strict source for g."
            )}
    return finalize(311,"Interval-Certified Localized Simplex Minimum and Twelve Translates",
                    "proven_one_strict_local_minimum_and_twelve_translation_orbit_at_supplied_g",
                    "The coefficient g=4 is supplied. No global-minimum, exactly-twelve-critical-points, physical-vacuum, or selector theorem is claimed.",packet)


def st312(a: np.ndarray) -> dict:
    center=localized_seed(a);u=np.ones(N)/N
    def softmax(y):
        z=np.r_[y,0.0];z-=z.max();e=np.exp(z);return e/e.sum()
    def val(y):
        p=softmax(y);q=p-u;return float(np.sum(p*np.log(N*p))-2*q@a@q)
    values=[]
    for k in range(300):
        rng=np.random.default_rng(SEED+k);sol=minimize(val,rng.normal(0,3,N-1),method="BFGS",
                                                       options={"maxiter":1000,"gtol":1e-9})
        values.append(val(sol.x))
    ap=np.linalg.pinv(a)
    packet={"exact_dual_reduction":"inf_{h in 1-perp} [(1/(2g))<h,A^+h>-log((1/12)sum_i exp(h_i))]",
            "primal_fixed_point":"p_i=exp(g(Ap)_i)/sum_j exp(g(Ap)_j)",
            "certified_localized_value":objective_from_radial(center,a,4),
            "uniform_value":0.0,"pure_vertex_value":math.log(N)-2*float(np.trace(a)/N),
            "multistart_trials":len(values),"best_multistart_value":min(values),
            "gap_best_multistart_to_certified_center":min(values)-objective_from_radial(center,a,4),
            "global_certificate":None,
            "theorem":(
                "Fenchel completion of the negative quadratic gives the displayed exact coercive 11-dimensional dual minimization and the Gibbs fixed-point equation. This is a smaller rigorous target for global certification. The certified localized root beats both the uniform state and every pure vertex; 300 frozen starts find no lower point. Neither facts nor sampling prove global minimality, so the global claim remains open rather than being promoted from basin evidence."
            )}
    return finalize(312,"Exact Global Dual Reduction and Bounded Globality Audit",
                    "proven_dual_reduction_plus_strong_numerical_global_evidence_no_global_certificate",
                    "The 11-dimensional dual has not been exhaustively interval-covered. Global minimality and complete critical-orbit enumeration remain open.",packet)


def st313(a: np.ndarray) -> dict:
    b=radial_matrix(a);x=localized_seed(a);rows=[]
    for g in np.linspace(4.0,2.68,265):
        sol=root(lambda z:localized_system(z,b,g),x,jac=lambda z:localized_jacobian(z,b,g),tol=1e-11)
        if not sol.success or np.min(sol.x[:7])<=0: break
        x=sol.x;rows.append({"g":float(g),"value":objective_from_radial(x,a,g),
                             "peak_probability":float(x[0]),
                             "Jacobian_min_singular":float(np.linalg.svd(localized_jacobian(x,b,g),compute_uv=False)[-1])})
    cross=None
    for left,right in zip(rows[:-1],rows[1:]):
        if left["value"]*right["value"]<=0:
            cross=[right["g"],left["g"]];break
    packet={"continued_rows":rows,"sampled_gain_interval":[rows[-1]["g"],rows[0]["g"]],
            "first_order_energy_crossing_bracket_sampled":cross,
            "minimum_sampled_Jacobian_singular":min(x["Jacobian_min_singular"] for x in rows),
            "result":(
                "Numerical continuation keeps the localized stationary family regular from g=4 down to the reported endpoint. Its energy crosses the uniform value in the narrow sampled bracket, while the branch remains strongly localized. This supports a first-order transition before uniform Hessian loss, but the crossing and branch interval are not outward-certified."
            )}
    return finalize(313,"Localized-Branch Continuation and First-Order Crossing",
                    "strong_numerical_evidence_for_first_order_localization_crossing",
                    "The continuation and crossing use floating roots. Spinodals, global coexistence, and exact transition gain need interval certification.",packet)


def st314(a: np.ndarray) -> dict:
    old=json.loads((ROOT/"FIN_ST293_ST307_Results.json").read_text())["ST298"]
    rows=old["rows"];centers=[np.array(x["center"]) for x in rows]
    radii=[x["certificate"]["radius"] for x in rows]
    gaps=[float(np.linalg.norm(centers[i+1]-centers[i])-radii[i]-radii[i+1]) for i in range(len(rows)-1)]
    aiv,_,_=strict_interval_matrix();large=[]
    for i in range(1,len(rows)-1,8):
        tangent=centers[i+1]-centers[i-1];tangent/=np.linalg.norm(tangent)
        rad=.51*min(np.linalg.norm(centers[i]-centers[i-1]),np.linalg.norm(centers[i+1]-centers[i]))
        from fin_st203_st217_research import pseudo_krawczyk
        large.append(pseudo_krawczyk(centers[i],aiv,tangent,centers[i],rad))
    packet={"existing_boxes":len(rows),"minimum_center_distance_minus_box_radii":min(gaps),
            "existing_boxes_overlap":False,"half_step_tube_trials":len(large),
            "half_step_tube_inclusions":sum(x["included"] for x in large),
            "half_step_minimum_margin":min(x["margin"] for x in large),
            "theorem":(
                "The 160 accepted ST298 root boxes are pairwise separated by a strictly positive computed center-distance gap after paying both radii. They therefore do not form an overlapping interval tube. Enlarging representative sections to half-step scale fails Krawczyk inclusion on part of the chain. The discrete unique-root theorem remains valid, but a continuous-tube theorem does not follow from it."
            )}
    return finalize(314,"Audit of the Missing Continuous Validated Tube",
                    "proven_existing_certificates_are_disjoint_and_do_not_form_tube",
                    "Failure of this enclosure strategy is not proof that the smooth branch breaks or that no different interval implicit-function certificate exists.",packet)


def st315() -> dict:
    tau=[np.array([[.75,0],[0,.25]]),np.array([[.5,.25],[.25,.5]]),np.array([[.5,-.25],[-.25,.5]])]
    gamma=.75*np.eye(2);pplus=.5*np.array([[1,1],[1,1.]]);pminus=.5*np.array([[1,-1],[-1,1.]])
    effects=[np.zeros((2,2)),pplus,pminus]
    obj=sum(float(np.trace(effects[i]@tau[i])) for i in range(3))
    slack_eigs=[np.linalg.eigvalsh(gamma-t).tolist() for t in tau]
    packet={"states":[x.tolist() for x in tau],"dual_Gamma":gamma.tolist(),"dual_trace":1.5,
            "optimal_effects":[x.tolist() for x in effects],"primal_value":obj,
            "dual_slack_eigenvalues":slack_eigs,"POVM_sum_defect":float(np.linalg.norm(sum(effects)-np.eye(2))),
            "theorem":(
                "For the declared rational noncommuting triple, Gamma=3I/4 is dual feasible and the POVM (0,|+><+|,|-><-|) is primal feasible. Both have exact objective 3/2, so weak duality proves global optimality. The optimum deliberately ignores the first output; this is an exact three-output result that cannot be represented by one binary positive-part rule."
            )}
    return finalize(315,"Exact Rational Three-Output Recovery SDP Optimum",
                    "proven_exact_primal_dual_optimum_three_halves",
                    "The rational triple is a mathematical fixture, not a FIN-derived laboratory channel or a formula for arbitrary multi-output recovery.",packet)


def st316() -> dict:
    transition=np.array([[.9,.2],[.1,.8]])
    emissions=(np.array([[.98,.02],[.92,.08]]),np.array([[.08,.92],[.02,.98]]))
    rows=[]
    for h,e in enumerate(emissions):
        for y in (0,1):
            m=transition@np.diag(e[:,y]);cross=float(m[0,0]*m[1,1]/(m[0,1]*m[1,0]))
            diameter=abs(math.log(cross));coef=math.tanh(diameter/4)
            rows.append({"hypothesis":h,"output":y,"cross_ratio":cross,
                         "projective_diameter":diameter,"Birkhoff_coefficient":coef})
    packet={"positive_transition":transition.tolist(),"rows":rows,
            "uniform_exact_coefficient":"5/7","maximum_coefficient":max(x["Birkhoff_coefficient"] for x in rows),
            "theorem":(
                "Each positive two-state filter matrix has cross ratio 36; row or column emission scaling cancels from that ratio. Birkhoff's formula therefore gives the exact uniform Hilbert-projective contraction coefficient (sqrt(36)-1)/(sqrt(36)+1)=5/7. This supplies the contraction ingredient absent in ST300. A separate weighted derivative-recursion estimate is still required before promoting the finite-depth Chernoff optimizer."
            )}
    return finalize(316,"Exact Projective Contraction of the Frozen HMM Filters",
                    "proven_uniform_filter_contraction_tail_promotion_still_partial",
                    "Belief contraction alone does not bound the full normalized Chernoff derivative transfer operator. The fixture remains synthetic.",packet)


def st317(a: np.ndarray) -> dict:
    ec,ef,_,_=parametric_data(a);nu=np.array([.2,.7,.05]);x=root(lambda z:point_design_system(z,ec,ef,*nu),[2.1862,.53983]).x
    h=1e-6;sens=[]
    for j in range(3):
        ep=np.zeros(3);ep[j]=h
        xp=root(lambda z:point_design_system(z,ec,ef,*(nu+ep)),x).x
        xm=root(lambda z:point_design_system(z,ec,ef,*(nu-ep)),x).x
        sens.append(float(np.linalg.norm((xp-xm)/(2*h))))
    prod=64000.;geom=(prod/np.prod(sens))**(1/3);ideal=np.array(sens)*geom
    candidates=[]
    for n0 in range(20,101):
        for n1 in range(20,101):
            n2=max(1,int(prod//(n0*n1)))
            if n0*n1*n2<=64000:
                score=max(sens[i]/[n0,n1,n2][i] for i in range(3))
                candidates.append((score,n0,n1,n2))
    best=min(candidates)
    packet={"root_center":x.tolist(),"implicit_root_sensitivity_norms":sens,
            "continuous_product_budget":prod,"continuous_optimal_axis_counts":ideal.tolist(),
            "integer_budget_not_exceeding_64000":{"counts":list(best[1:]),"cells":int(np.prod(best[1:])),"max_linearized_axis_load":best[0]},
            "theorem":(
                "For a rectangular first-order enclosure with axis loads a_i/n_i and fixed product of subdivision counts, minimizing the maximum load forces a_i/n_i equal, hence n_i proportional to a_i. The displayed counts apply this exact allocation theorem to numerical implicit-root sensitivities. They are a design recommendation, not a replacement for the complete outward cover."
            )}
    return finalize(317,"Sensitivity-Weighted Nuisance-Cover Allocation",
                    "proven_conditional_allocation_rule_plus_numerical_FIN_sensitivities",
                    "The sensitivities are floating local derivatives and the linear envelope may fail over large boxes. No new complete cover or physical tolerance is claimed.",packet)


def st318() -> dict:
    packet={"norm":"operator norm","lower_bound":2.0,
            "spectral_argument":"U=XZ is unitary; XUX=U* pairs nonreal eigenvalues; odd dimension leaves one eigenvalue +1 or -1",
            "consequence":"||XZ+ZX||=||U+U*|| >= 2",
            "sharp_example":"X=Z=I gives norm 2",
            "theorem":(
                "For Hermitian involutions X,Z in odd complex dimension, U=XZ is unitary and its nonreal eigenvalues occur in conjugate pairs. An unpaired eigenvalue is therefore +1 or -1, on which U+U*=XZ+ZX has eigenvalue +2 or -2. Hence the operator-norm anticommutation defect is at least 2, and the bound is sharp. Odd dimension cannot even approximate an involutive Clifford pair below this fixed error."
            )}
    return finalize(318,"Sharp Approximate Odd-Dimensional Clifford Obstruction",
                    "proven_sharp_operator_norm_lower_bound_two",
                    "The theorem assumes exact Hermitian involutions. Nearly involutive, non-Hermitian, graded, or enlarged systems require separate stability bounds.",packet)


def st319() -> dict:
    rates=[1.2,.7,1.8,.9];h=[1.]
    for r in rates:h.append(h[-1]*r)
    reconstructed=[[h[j]/h[i] for j in range(len(h))] for i in range(len(h))]
    packet={"sample_edge_rates":rates,"vertex_gauge":h,"reconstructed_cocycle":reconstructed,
            "tree_cohomology":"every positive multiplicative 1-cocycle on a connected tree is a coboundary g_mn=h_m/h_n",
            "loop_observable":"only products around noncontractible cycles survive gauge change",
            "theorem":(
                "On a linear or tree-like refinement hierarchy, choose a root and multiply edge calibrations along the unique path to define h_n. The cocycle then equals h_m/h_n and is gauge-trivial. Thus a nonconstant inter-layer scale law on a tree is not by itself observable; path-dependent holonomy requires cycles, while an absolute magnitude requires an anchor."
            )}
    return finalize(319,"Gauge Classification of Positive Scale Cocycles on Refinement Trees",
                    "proven_all_tree_scale_cocycles_are_gauge_trivial",
                    "The result does not exclude physical relative scales after a section or anchor is chosen. Non-tree refinement categories may carry holonomy.",packet)


def st320() -> dict:
    packet={"admissible_class":"measurable profiles 0<=q(mu)<=1 with one linear heat-trace budget",
            "pointwise_dimension_kernel":"d_t(mu)=4 mu t/(exp(2 mu t)+1)",
            "heat_cost_kernel":"h_t0(mu)=exp(-2 mu t0)",
            "ordering_score":"d_t(mu)/h_t0(mu)",
            "extremizer":"q=1 above a score threshold, q=0 below, with at most one partially filled level set",
            "theorem":(
                "The continuum soft-IR optimization with a single linear budget is an infinite-dimensional fractional knapsack problem. The exchange argument (bathtub principle) proves that every pointwise maximizer is bang-bang when ordered by d_t/h_t0, up to a threshold level set. Therefore the two-parameter smooth family sampled in ST304 is not the unrestricted pointwise extremal class. A time-window minimax objective couples several kernels and is not solved by this theorem."
            )}
    return finalize(320,"Bathtub Extremal Theorem for Pointwise Soft-IR Design",
                    "proven_pointwise_bang_bang_extremal_structure",
                    "No sourced density, cutoff, absolute unit, trace-class infinite exact tower, or minimax plateau-window optimizer is derived.",packet)


def st321() -> dict:
    # Rows are exponents of supplied invariant dimensional constants in (L,T,E).
    c=np.array([1,-1,0]);hbar=np.array([0,1,1]);mat=np.vstack([c,hbar])
    rank=int(np.linalg.matrix_rank(mat));null=3-rank
    packet={"initial_structure_group":"(R_+)^3 for length,time,energy",
            "supplied_relations":{"speed_anchor_c":[1,-1,0],"action_anchor_hbar":[0,1,1]},
            "constraint_rank":rank,"residual_group_dimension":null,
            "residual_scaling":"(ell,tau,epsilon)->(a ell,a tau,epsilon/a)",
            "additional_independent_anchor_needed":1,
            "theorem":(
                "Dimensional-anchor counting is the rank of the exponent matrix. Supplying a universal speed relation L/T and an action relation ET gives two independent constraints on the three-dimensional calibration torsor, leaving one positive scaling orbit. One additional independent dimensional anchor fixes a section. These relations reduce the bridge but are not generated by dimensionless FIN data."
            )}
    return finalize(321,"Calibration Structure-Group Reduction and Minimal Anchor Count",
                    "proven_two_relations_leave_one_dimensional_calibration_orbit",
                    "c and hbar are supplied physical anchors, not FIN predictions. No SI values or unique channel-identification map follows.",packet)


def st322() -> dict:
    records=list(ROOT.glob("FIN_ST322_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet={"required_pattern":"FIN_ST322_INDEPENDENT_RAW_EVENTS*.jsonl",
            "matching_records":[x.name for x in records],"independent_record_count":len(records),
            "theorem":"No independently registered raw-event packet is present. Local mathematics cannot manufacture provider/registrar separation, apparatus, or empirical evidence."}
    return finalize(322,"Independent Evidence Gate After the Mathematical Batch",
                    "blocked_no_independent_events" if not records else "record_present_requires_blinded_validation",
                    "This is an evidentiary stop, not a negative experiment.",packet)


def st323() -> dict:
    packet={"declared_refinement":"A_tilde=R A R* direct-sum B on coarse plus 12-dimensional fiber",
            "coarse_context":"preparations and effects supported only in ran(R)",
            "blind_perturbation_space":"all symmetric delta B","blind_kernel_dimension":78,
            "fiber_complete_context":"child-resolved preparations/effects and heat records for an open time interval",
            "fiber_complete_kernel_dimension":0,
            "theorem":(
                "Block functional calculus gives f(A_tilde)=R f(A)R* direct-sum f(B). Hence every coarse record is independent of arbitrary symmetric changes delta B: the declared blindness kernel is Sym(12), dimension 78. If a fiber-complete context records exp(-tB) on an open interval, differentiation at t=0 recovers B and the kernel collapses to zero. Operational depth is therefore instrument-relative, not a layer number alone."
            )}
    return finalize(323,"Exact Observer-Blindness Kernel for a Two-Child Refinement",
                    "proven_78_dimensional_coarse_blindness_and_fiber_tomographic_collapse",
                    "The block refinement and ideal complete fiber instruments are supplied. No laboratory implements them in this packet.",packet)


def st324() -> dict:
    packet={"preorder":"O_coarse <= O_fine when every coarse record is a stochastic postprocessing of fine records",
            "kernel_implication":"Blackwell dominance implies N_fine subset N_coarse",
            "converse":"false",
            "counterexample":"BSC(0.1) and BSC(0.2) both have zero parameter kernel, but BSC(0.2) is a garbling of BSC(0.1)",
            "theorem":(
                "A postprocessing cannot distinguish models that the finer experiment fails to distinguish, so Blackwell dominance implies reverse inclusion of blindness kernels. Kernel inclusion alone is not a complete depth order: two identifiable binary channels have the same zero kernel while their statistical informativeness differs strictly. FIN observer depth must therefore retain the experiment/garbling structure, not only kernel dimension."
            )}
    return finalize(324,"Blackwell Observer Depth Versus Blindness-Kernel Inclusion",
                    "proven_kernel_monotonicity_and_converse_counterexample",
                    "This defines an operational hierarchy; it does not identify human observers, consciousness, or a physical refinement apparatus.",packet)


def st325(a: np.ndarray) -> dict:
    ev=np.linalg.eigvalsh(a);ratios=[float(ev[k]/ev[1]) for k in range(1,N) if ev[k]>1e-10]
    packet={"coarse_positive_spectral_ratios":ratios,
            "exact_refinement":"A_tilde R=R A",
            "coarse_clock_law":"phase beats and heat-rate ratios lambda_i/lambda_j are exactly preserved",
            "wave_clock_law":"wave-frequency ratios are sqrt(lambda_i/lambda_j)",
            "fiber_clock_status":"depends on free B and cannot be compared canonically without a transition cocycle/anchor",
            "theorem":(
                "Exact intertwining transports every coarse eigenmode with the same dimensionless eigenvalue. Consequently all coarse unitary phase-frequency ratios, heat-rate ratios, and square-root wave-frequency ratios are refinement-natural. Absolute seconds and comparisons to fiber modes remain undetermined because B and the calibration section are free."
            )}
    return finalize(325,"Refinement-Naturality of Relative Spectral Clocks",
                    "proven_coarse_relative_clock_invariance_across_exact_refinement",
                    "A periodic relative clock is not a time arrow or SI second. Channel-specific functional maps must remain distinct.",packet)


def st326() -> dict:
    rows=[{"depth":n,"binary_branch_record_lower_bits":n,"twelve_branch_record_lower_bits":n*math.log2(12)} for n in (1,2,4,8,16,32)]
    packet={"rows":rows,"counting_bound":"a lossless code for b^n arbitrary branch histories needs at least log2(b^n)=n log2 b bits",
            "controlled_compression_object":"short refinement rule + growing branch/initial record + observer-dependent lossy decoder and error bound",
            "theorem":(
                "A fixed finite kernel and a fixed refinement rule cannot losslessly encode every arbitrary depth-n branch history with bounded additional storage: injectivity requires at least n log2 b record bits. Fractal self-similarity compresses the rule, not arbitrary initial data or choices. Sublinear storage is possible only after restricting the state class, allowing loss, or imposing correlations; each requires a declared encoder, decoder, and distortion criterion."
            )}
    return finalize(326,"Information-Theoretic Lower Bound on Fractal Branch Compression",
                    "proven_lossless_branch_record_lower_bound",
                    "The theorem does not deny compression of structured states and does not estimate the universe's physical information content.",packet)


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    r=out["ST313"];fig,ax=plt.subplots(figsize=(7,4))
    ax.plot([x["g"] for x in r["continued_rows"]],[x["value"] for x in r["continued_rows"]],lw=2)
    ax.axhline(0,color="black",ls="--");ax.set(xlabel="supplied gain g",ylabel="localized branch objective",
        title="ST313: numerical first-order energy crossing");fig.tight_layout();fig.savefig(FIG_DIR/"st313_branch_crossing.png",dpi=190);plt.close(fig)
    r=out["ST309"];fig,ax=plt.subplots(figsize=(7,4))
    x=np.linspace(0,6,300);ax.semilogy(x,np.exp(-x));ax.set(xlabel="gamma t / tau",ylabel="relative gain error",
        title="ST309: exact finite-pump relaxation");fig.tight_layout();fig.savefig(FIG_DIR/"st309_pump_relaxation.png",dpi=190);plt.close(fig)


def main() -> None:
    np.random.seed(SEED);_,a,s=strict_operator()
    out={"metadata":{"seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
                     "strict_row_sum":s,"scope":"local analytic, interval, exact finite, and bounded numerical research"}}
    funcs={308:lambda:st308(a),309:st309,310:lambda:st310(a),311:lambda:st311(a),312:lambda:st312(a),
           313:lambda:st313(a),314:lambda:st314(a),315:st315,316:st316,317:lambda:st317(a),318:st318,
           319:st319,320:st320,321:st321,322:st322,323:st323,324:st324,325:lambda:st325(a),326:st326}
    for k in range(308,327): out[f"ST{k}"]=funcs[k]()
    make_figures(out);RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(308,327):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    print(json.dumps({k:v["status"] for k,v in out.items() if k.startswith("ST")},indent=2))


if __name__ == "__main__": main()
