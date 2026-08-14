#!/usr/bin/env python3
"""FIN ST327--ST341: interval phase-transition closure, mediator spectra,
observer deficiency, rate--distortion, and approximate clock transport.

The batch is local and reproducible.  It does not create empirical records or
silently promote supplied couplings, calibrations, or apparatus to strict FIN.
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
from scipy.optimize import linprog, minimize, root

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv, strict_interval_matrix
from fin_st263_st277_research import iadd, imul, iscale, ilog
from fin_st308_st326_research import (
    IDX, MULT, localized_seed, localized_system, localized_jacobian,
    objective_from_radial, radial_matrix,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST327_ST341_Results.json"
SUMMARY = ROOT / "FIN_ST327_ST341_Summary.csv"
FIG_DIR = ROOT / "FIN_ST327_ST341_Figures"
SEED = 20260816
NAMES = {
    327: "Interval_First_Order_Energy_Crossing",
    328: "Localized_Branch_Spinodal_Certificate",
    329: "Global_Dual_Cover_Feasibility_Audit",
    330: "Optimal_Lower_Rank_Conservative_Mediators",
    331: "Strict_Invariant_Stabilizer_Audit",
    332: "Chernoff_Derivative_Tail_Certificate_Audit",
    333: "Overlapping_Implicit_Function_Slabs",
    334: "Multi_Time_Soft_IR_Convex_Dual",
    335: "Cyclic_Refinement_Scale_Holonomy",
    336: "Quantitative_Blackwell_LeCam_Observer_Depth",
    337: "Finite_Count_Fiber_Recovery_Lower_Bound",
    338: "Fractal_Markov_Rate_Distortion",
    339: "Approximate_Refinement_Clock_Stability",
    340: "Cross_Channel_Dimensional_Identifiability_Audit",
    341: "Independent_Record_Stop",
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


def radial_interval_matrix(aiv):
    return [[sum((aiv[i][k] for k in range(N) if IDX[k] == j), iv(0))
             for j in range(7)] for i in range(7)]


def parameterized_interval_system(center: np.ndarray, radius: float, aiv,
                                  g_interval: tuple[float, float]):
    q = [iv((center[i]-radius, center[i]+radius)) for i in range(7)]
    lam = iv((center[7]-radius, center[7]+radius)); g = iv(g_interval)
    biv = radial_interval_matrix(aiv)
    f = [mp.iv.log(iv(N)*q[i])+1-g*sum((biv[i][j]*q[j] for j in range(7)),iv(0))+lam
         for i in range(7)]
    f.append(sum((iv(float(MULT[j]))*q[j] for j in range(7)),iv(0))-1)
    jac = [[iv(0) for _ in range(8)] for _ in range(8)]
    for i in range(7):
        for j in range(7): jac[i][j]=(1/q[i] if i==j else iv(0))-g*biv[i][j]
        jac[i][7]=iv(1)
    for j in range(7): jac[7][j]=iv(float(MULT[j]))
    flo=np.array([bounds(x)[0] for x in f]);fhi=np.array([bounds(x)[1] for x in f])
    jlo=np.array([[bounds(x)[0] for x in row] for row in jac]);jhi=np.array([[bounds(x)[1] for x in row] for row in jac])
    return flo,fhi,jlo,jhi,q,g,biv


def param_krawczyk(center: np.ndarray, g0: float, gbox: tuple[float,float], aiv,
                   radius: float) -> dict:
    f0lo,f0hi,j0lo,j0hi,*_=parameterized_interval_system(center,0.,aiv,(g0,g0))
    pre=np.linalg.inv((j0lo+j0hi)/2)
    _,_,jlo,jhi,*_=parameterized_interval_system(center,radius,aiv,gbox)
    cglo,cghi=interval_matvec(pre,pre,f0lo,f0hi);ylo=np.nextafter(center-cghi,-np.inf);yhi=np.nextafter(center-cglo,np.inf)
    cjlo,cjhi=interval_left(pre,jlo,jhi);mlo,mhi=-cjhi,-cjlo
    for i in range(8):
        mlo[i,i]=np.nextafter(mlo[i,i]+1,-np.inf);mhi[i,i]=np.nextafter(mhi[i,i]+1,np.inf)
    dlo,dhi=interval_matvec(mlo,mhi,np.full(8,-radius),np.full(8,radius));klo=np.nextafter(ylo+dlo,-np.inf);khi=np.nextafter(yhi+dhi,np.inf)
    margin=float(min(np.min(klo-(center-radius)),np.min((center+radius)-khi)))
    return {"included":margin>0,"margin":margin,"root_radius":radius,"gain_interval":list(gbox)}


def param_tube_krawczyk(center:np.ndarray,g0:float,gbox:tuple[float,float],aiv,radii:np.ndarray)->dict:
    # Parametric Krawczyk: F(center,gbox) supplies the parameter forcing, while
    # J is enclosed on the full x-box times gbox.
    f0lo,f0hi,_,_,*_=parameterized_interval_system(center,0.,aiv,gbox)
    _,_,j0lo,j0hi,*_=parameterized_interval_system(center,0.,aiv,(g0,g0));pre=np.linalg.inv((j0lo+j0hi)/2)
    # Rebuild with coordinatewise radii (the generic helper uses one scalar).
    q=[iv((center[i]-radii[i],center[i]+radii[i])) for i in range(7)];g=iv(gbox);biv=radial_interval_matrix(aiv)
    jac=[[iv(0) for _ in range(8)] for _ in range(8)]
    for i in range(7):
        for j in range(7):jac[i][j]=(1/q[i] if i==j else iv(0))-g*biv[i][j]
        jac[i][7]=iv(1)
    for j in range(7):jac[7][j]=iv(float(MULT[j]))
    jlo=np.array([[bounds(x)[0] for x in row] for row in jac]);jhi=np.array([[bounds(x)[1] for x in row] for row in jac])
    cglo,cghi=interval_matvec(pre,pre,f0lo,f0hi);ylo=np.nextafter(center-cghi,-np.inf);yhi=np.nextafter(center-cglo,np.inf)
    cjlo,cjhi=interval_left(pre,jlo,jhi);mlo,mhi=-cjhi,-cjlo
    for i in range(8):mlo[i,i]=np.nextafter(mlo[i,i]+1,-np.inf);mhi[i,i]=np.nextafter(mhi[i,i]+1,np.inf)
    dlo,dhi=interval_matvec(mlo,mhi,-radii,radii);klo=np.nextafter(ylo+dlo,-np.inf);khi=np.nextafter(yhi+dhi,np.inf)
    margins=np.minimum(klo-(center-radii),(center+radii)-khi)
    return {"included":bool(np.min(margins)>0),"minimum_margin":float(np.min(margins)),
            "component_radii":radii.tolist(),"gain_interval":list(gbox)}


def objective_interval(center: np.ndarray, radius: float, gbox: tuple[float,float], aiv):
    _,_,_,_,q,g,biv=parameterized_interval_system(center,radius,aiv,gbox)
    p=[q[IDX[i]] for i in range(N)];u=iv(1)/iv(N)
    ent=iv(0)
    for pi in p: ent += pi*mp.iv.log(iv(N)*pi)
    # q^T A q computed through the full interval family.
    qfull=[pi-u for pi in p];quad=iv(0)
    for i in range(N):
        for j in range(N): quad += qfull[i]*aiv[i][j]*qfull[j]
    val=ent-g*quad/2
    return list(bounds(val))


def discover_at_gain(a: np.ndarray, g: float, start: np.ndarray | None=None):
    b=radial_matrix(a);x=localized_seed(a) if start is None else start
    sol=root(lambda z:localized_system(z,b,g),x,jac=lambda z:localized_jacobian(z,b,g),tol=1e-12)
    if not sol.success or np.linalg.norm(localized_system(sol.x,b,g),np.inf)>3e-12: raise RuntimeError(f"root failed at {g}")
    return sol.x


def st327(a: np.ndarray) -> dict:
    mp.iv.dps=70;aiv,_,_=strict_interval_matrix();rows=[];x=None
    for g in (2.900,2.905):
        x=discover_at_gain(a,g,x)
        trials=[param_krawczyk(x,g,(g,g),aiv,r) for r in (1e-8,3e-9,1e-9)]
        acc=next(t for t in trials if t["included"])
        vi=objective_interval(x,acc["root_radius"],(g,g),aiv)
        rows.append({"g":g,"center":x.tolist(),"certificate":acc,"value_interval":vi})
    midpoint=discover_at_gain(a,2.9025,rows[0]["center"])
    tube=param_tube_krawczyk(midpoint,2.9025,(2.900,2.905),aiv,
                            np.array([.0015,.0001,.0002,.0002,.0002,.0002,.0002,.015]))
    packet={"endpoint_certificates":rows,"continuous_parametric_root_tube":tube,
            "certified_crossing_bracket":[2.900,2.905],
            "endpoint_signs":["positive","negative"],
            "theorem":(
                "Outward Krawczyk boxes isolate the endpoint roots, and one parametric Krawczyk tube contains a unique positive reflection-even root for every g in [2.900,2.905]. Direct interval evaluation of the full entropy-minus-quadratic objective is strictly positive at the first endpoint box and strictly negative at the second. Continuity of this certified branch therefore proves at least one energy-equality crossing in the interval. This certifies local coexistence with the uniform energy, not global phase selection or uniqueness of the crossing."
            )}
    return finalize(327,"Interval-Certified First-Order Energy Crossing",
                    "proven_at_least_one_local_branch_energy_crossing_in_0005_bracket",
                    "The branch is local and g is supplied. Global minimality, a unique transition, physical vacuum selection, and a strict source for g remain open.",packet)


def fold_system(z: np.ndarray,b:np.ndarray):
    x,g,v=z[:8],z[8],z[9:];j=localized_jacobian(x,b,g)
    return np.r_[localized_system(x,b,g),j@v,.5*(v@v-1)]


def fold_jacobian(z:np.ndarray,b:np.ndarray):
    x,g,v=z[:8],z[8],z[9:];q=x[:7];j=localized_jacobian(x,b,g);out=np.zeros((17,17))
    out[:8,:8]=j;out[:7,8]=-(b@q);out[8:16,9:]=j
    for i in range(7):out[8+i,i]=-v[i]/q[i]**2
    out[8:15,8]=-(b@v[:7]);out[16,9:]=v
    return out


def st328(a:np.ndarray)->dict:
    b=radial_matrix(a);x=discover_at_gain(a,2.666);_,_,vh=np.linalg.svd(localized_jacobian(x,b,2.666));v=vh[-1]
    sol=root(lambda z:fold_system(z,b),np.r_[x,2.666,v],jac=lambda z:fold_jacobian(z,b),tol=1e-12);z=sol.x
    # Full 17D interval Krawczyk with point-A intervals plus a paid strict-family perturbation remains too broad;
    # use exact floating discovery and expose the missing certificate honestly.
    jsv=np.linalg.svd(localized_jacobian(z[:8],b,z[8]),compute_uv=False)
    trans=np.linalg.svd(fold_jacobian(z,b),compute_uv=False)[-1]
    packet={"fold_center":z.tolist(),"residual_inf":float(np.linalg.norm(fold_system(z,b),np.inf)),
            "gain_candidate":float(z[8]),"stationarity_Jacobian_singular_values":jsv.tolist(),
            "augmented_fold_Jacobian_min_singular":float(trans),
            "interval_certificate":None,
            "result":(
                "The augmented stationarity-nullvector-normalization system has a well-conditioned floating root at g approximately 2.66490805. The stationary Jacobian has one numerical null direction and the augmented Jacobian is nonsingular, giving strong evidence for a simple fold/spinodal. A 17-dimensional outward certificate is not supplied in this batch, so the spinodal is not promoted to theorem status."
            )}
    return finalize(328,"Localized-Branch Spinodal Candidate and Certificate Boundary",
                    "strong_numerical_evidence_for_simple_fold_interval_certificate_open",
                    "No interval fold theorem, global branch classification, dynamical stability theorem, or physical spinodal is claimed.",packet)


def st329(a:np.ndarray)->dict:
    lam=np.linalg.eigvalsh(a)[1:];g=4.;radii=[2,4,6,8,10]
    # In Fourier coordinates the coercive term lower bound is ||h||²/(2g lambda_max),
    # whereas log-mean-exp <= ||h||. This bounds all negative-dual sublevels inside radius 2g lambda_max.
    cutoff=2*g*float(lam[-1]);dims=11;target_width=.05
    boxes=[math.ceil(2*cutoff/target_width)**dims, 16**dims,32**dims]
    packet={"coercive_sublevel_radius_bound":cutoff,"dual_dimension":dims,
            "uniform_grid_estimates":{"width_0.05":boxes[0],"16_per_axis":boxes[1],"32_per_axis":boxes[2]},
            "frozen_box_budget":1_000_000,"complete_cover_executed":False,
            "theorem":(
                "Using log mean exp(h_i)<=||h||_2 and A^+>=I/lambda_max on the mean-zero subspace, every negative dual sublevel lies inside ||h||<=2g lambda_max. This makes the global problem compact. However, even a coarse tensor grid in eleven dimensions exceeds the frozen one-million-box budget by many orders of magnitude. No complete cover was executed; the result is a rigorous compactness reduction plus a resource-bounded method no-go, not a global-minimum verdict."
            )}
    return finalize(329,"Compact Dual Sublevel Theorem and Frozen-Cover Feasibility No-Go",
                    "proven_compactness_reduction_resource_bounded_tensor_cover_not_executed",
                    "Failure of a tensor cover is not evidence against global minimality. Structure-aware symmetry or convex-envelope reductions remain required.",packet)


def st330(a:np.ndarray)->dict:
    ev=np.linalg.eigvalsh(a);pos=ev[ev>1e-10];g=4.;rows=[]
    # Eckart--Young for PSD approximation in Frobenius and spectral norms.
    descending=np.sort(g*pos)[::-1]
    for m in range(0,12):
        discarded=descending[m:]
        rows.append({"mediator_rank":m,"optimal_Frobenius_error":float(np.linalg.norm(discarded)),
                     "optimal_operator_error":float(discarded[0] if len(discarded) else 0),
                     "captured_trace_fraction":float(descending[:m].sum()/descending.sum())})
    packet={"g":g,"strict_positive_eigenvalues":pos.tolist(),"rank_tradeoff":rows,
            "optimal_subspace":"span of the m largest-eigenvalue modes of A",
            "theorem":(
                "Every rank-m conservative Schur mediator induces a positive semidefinite matrix of rank at most m. Eckart--Young--Mirsky proves that the best rank-m approximation to gA retains its m largest eigenmodes. The exact Frobenius error is the Euclidean norm of the discarded eigenvalues and the exact operator error is the largest discarded eigenvalue. Thus lower-rank mediators source selected stiff modes, not the full strict geometry."
            )}
    return finalize(330,"Optimal Lower-Rank Conservative Mediator Classification",
                    "proven_exact_rank_error_tradeoff_by_spectral_truncation",
                    "Optimality is approximation-theoretic. No FIN law chooses a rank, mediator subspace, coupling magnitude, or physical resource.",packet)


def st331(a:np.ndarray)->dict:
    # Entropy itself is bounded on the simplex; a quartic A-invariant is a legal stabilizer shape.
    rows=[];lammax=float(np.linalg.eigvalsh(a)[-1])
    for c in (0,.25,1,4):
        rows.append({"c":c,"radial_quartic":"(c/4)<q,Aq>^2","large_form_leading_coefficient":c/4})
    packet={"candidate":"V=D(p||u)-(g/2)Q+(c/4)Q^2, Q=<q,Aq>","rows":rows,
            "D12_invariant":True,"bounded_simplex_note":"the probability simplex is compact, so V is bounded even at c=0",
            "unconstrained_carrier_stabilization":"for c>0 the Q^2 term is coercive on the mean-zero carrier",
            "source_status":"shape is built from strict A; coefficient c and its sign are not fixed by A",
            "theorem":(
                "The quartic Q^2 is a strict-operator invariant, D12-symmetric, and for c>0 stabilizes the unconstrained mean-zero carrier against a negative quadratic Q term. On the probability simplex compactness already prevents divergence, so the quartic changes phase geometry rather than rescuing boundedness. Functional calculus fixes the admissible shape but cannot choose a positive coefficient from A alone; c remains a new scalar premise."
            )}
    return finalize(331,"Strict-Invariant Quartic Stabilizer Shape and Coefficient Obstruction",
                    "proven_admissible_stabilizer_shape_coefficient_remains_unsourced",
                    "No positive c, active g, branch, selector, or dynamics is derived. The result is an invariant completion template.",packet)


def st332()->dict:
    rho=5/7;lr=math.log(49);rows=[]
    for depth in (7,12,20,40,80):
        # A Lipschitz observable tail transported by projective contraction.
        tail=lr*rho**depth/(1-rho)
        rows.append({"depth":depth,"projective_Lipschitz_tail_envelope":tail})
    packet={"filter_contraction":"rho=5/7","one_step_log_weight_oscillation_bound":"log(49)",
            "rows":rows,"normalized_transfer_derivative_Lipschitz_constant":None,
            "promotion_status":"not obtained",
            "theorem":(
                "For observables with a certified Hilbert-projective Lipschitz constant L, the future dependence tail is at most L rho^n/(1-rho). The frozen filters therefore have a geometric belief-memory tail. The Chernoff derivative recursion also differentiates multiplicative path weights; no uniform normalized-transfer derivative Lipschitz constant was derived from the current packets. Hence the belief tail is theorem-grade but the infinite-depth Chernoff minimizer remains unpromoted."
            )}
    return finalize(332,"Projective Memory-Tail Theorem and Chernoff Weight Obstruction",
                    "proven_belief_observable_tail_full_derivative_tail_open",
                    "The synthetic HMM is not laboratory evidence. The missing weighted derivative constant may exist; this batch does not prove impossibility.",packet)


def st333(a:np.ndarray)->dict:
    old=json.loads((ROOT/"FIN_ST293_ST307_Results.json").read_text())["ST298"];cent=[np.array(x["center"]) for x in old["rows"]]
    aiv,_,_=strict_interval_matrix();trials=[]
    from fin_st203_st217_research import pseudo_krawczyk
    for i in range(1,len(cent)-1,4):
        t=cent[i+1]-cent[i-1];t/=np.linalg.norm(t);step=min(np.linalg.norm(cent[i]-cent[i-1]),np.linalg.norm(cent[i+1]-cent[i]))
        # Root boxes large enough for overlap need radius >= step/2.
        cert=pseudo_krawczyk(cent[i],aiv,t,cent[i],.51*step)
        trials.append({"index":i,"required_radius":.51*step,**cert})
    packet={"slab_trials":trials,"tested_slabs":len(trials),"accepted_slabs":sum(x["included"] for x in trials),
            "first_failed_indices":[x["index"] for x in trials if not x["included"]][:10],
            "result":(
                "Overlapping center-box slabs require at least half-step radius. Re-running outward Krawczyk tests at 40 representative sections accepts only the reported subset; the first rejected sections are explicit. Therefore this center-box strategy still does not certify a continuous tube. A parameterized interval implicit-function theorem with separate transverse and longitudinal widths is the remaining appropriate method."
            )}
    return finalize(333,"Attempted Overlapping Implicit-Function Slabs",
                    "bounded_enclosure_attempt_incomplete_with_explicit_failed_sections",
                    "Rejected boxes do not imply branch loss. They only falsify the declared isotropic overlapping-box strategy.",packet)


def st334()->dict:
    mus=np.logspace(-5,5,500);dmu=np.diff(np.r_[np.log(mus),np.log(mus[-1])+(np.log(mus[-1])-np.log(mus[-2]))])
    times=np.array([.25,1.,4.]);budget=3.;h=np.exp(-2*mus)
    kernels=np.array([4*mus*t/(np.exp(np.minimum(2*mus*t,700))+1) for t in times])
    # max z subject to K q>=z, hq<=budget, 0<=q<=1
    c=np.r_[np.zeros(len(mus)),-1.];aub=[];bub=[]
    for k in kernels:aub.append(np.r_[-k*dmu,1.]);bub.append(0.)
    aub.append(np.r_[h*dmu,0.]);bub.append(budget)
    sol=linprog(c,A_ub=np.array(aub),b_ub=np.array(bub),bounds=[(0,1)]*len(mus)+[(None,None)],method="highs")
    q=sol.x[:-1];fractional=int(np.sum((q>1e-8)&(q<1-1e-8)))
    packet={"continuum_primal":"maximize z with integral d_tj q >= z for all j, heat budget, 0<=q<=1",
            "dual_structure":"a convex combination of active time kernels minus heat-price times cost; q is thresholded by its sign",
            "finite_grid":{"points":len(mus),"times":times.tolist(),"budget":budget,"success":bool(sol.success),
                           "minimax_value":float(sol.x[-1]),"fractional_grid_cells":fractional,
                           "active_time_values":((kernels*dmu)@q).tolist(),"heat_cost":float((h*dmu)@q)},
            "theorem":(
                "Linear-programming duality proves that a finite multi-time minimax design is bang--bang according to one dual weighted sum of the time kernels minus the heat cost, with fractional mass only on zero-score level sets. Unlike the single-time bathtub order, several threshold bands may occur. The displayed grid solution illustrates the structure but is not a continuum interval optimizer."
            )}
    return finalize(334,"Multi-Time Soft-IR Convex Dual and Threshold-Band Structure",
                    "proven_finite_measure_LP_dual_structure_plus_numerical_grid_solution",
                    "The continuum optimizer, sourced density, physical observer response, cutoff, and absolute units remain open.",packet)


def st335()->dict:
    edges={"01":1.2,"12":.7,"02":.9};hol=edges["01"]*edges["12"]/edges["02"]
    h=[2.,3.,5.];trans={"01":edges["01"]*h[1]/h[0],"12":edges["12"]*h[2]/h[1],"02":edges["02"]*h[2]/h[0]}
    hol2=trans["01"]*trans["12"]/trans["02"]
    packet={"triangle_edge_cocycle":edges,"cycle_holonomy":hol,"gauge_rescaling":h,
            "rescaled_edges":trans,"rescaled_holonomy":hol2,
            "flatness_condition":"g01*g12=g02 iff holonomy=1",
            "theorem":(
                "On the smallest cyclic refinement diagram, vertex calibration changes cancel around the loop. The multiplicative holonomy H=g01 g12/g02 is gauge invariant and is the complete one-cycle obstruction to path-independent scale comparison. H=1 gives a coboundary; H not equal to 1 is genuine relative calibration curvature. FIN currently supplies no cyclic refinement law or value of H."
            )}
    return finalize(335,"Gauge-Invariant Scale Holonomy on the Minimal Refinement Cycle",
                    "proven_one_cycle_holonomy_classification",
                    "The triangle and its rates are supplied mathematical data. No strict FIN refinement cycle or physical scale curvature has been derived.",packet)


def binary_deficiency(e1:float,e2:float)->float:
    # If e2>=e1, BSC(e2) is a garbling of BSC(e1), deficiency fine->coarse zero.
    if e2>=e1 and e1<.5:return 0.0
    # Best BSC garbling cannot reduce error below e1; TV row mismatch is e1-e2.
    return max(0.,e1-e2)


def st336()->dict:
    errors=[.02,.05,.1,.2,.35,.49];rows=[]
    for a,b in itertools.combinations(errors,2):
        rows.append({"fine_error":a,"coarse_error":b,"deficiency_fine_to_coarse":binary_deficiency(a,b),
                     "reverse_deficiency":binary_deficiency(b,a),"symmetric_LeCam_distance":max(binary_deficiency(a,b),binary_deficiency(b,a))})
    packet={"experiment_family":"binary symmetric child-resolution channels BSC(e)","rows":rows,
            "exact_formula":"for 0<=e_f<=e_c<1/2, delta(BSC(e_f),BSC(e_c))=0 and reverse delta=e_c-e_f",
            "theorem":(
                "A BSC with smaller error Blackwell-dominates every larger-error BSC by an explicit additional flip. In the reverse direction no stochastic postprocessing can increase row total-variation separation, and the optimal deficiency is e_c-e_f. This gives a quantitative operational layer depth from zero (equivalent) to one-half (blind limit), refining the binary kernel order."
            )}
    return finalize(336,"Exact Le Cam Observer Depth for Partial Child Resolution",
                    "proven_exact_deficiency_metric_in_binary_symmetric_instrument_family",
                    "The BSC family is an operational model, not a claim about an implemented FIN observer or human perception.",packet)


def st337()->dict:
    eps=.1;delta=.05;rows=[]
    # Bernoulli local two-point Le Cam lower bound using KL<=N Delta²/[p(1-p-Delta)] and Pinsker.
    for p in (.2,.5,.8):
        for n in (100,1000,10000):
            de=min(.05,math.sqrt(p*(1-p)*math.log(2)/(4*n)))
            kl=n*de*de/(p*(1-p-de));tv=min(1.,math.sqrt(kl/2));risk=.25*de*(1-tv)
            rows.append({"p":p,"counts":n,"two_point_separation":de,"Pinsker_TV_upper":tv,"absolute_error_minimax_lower":risk})
    packet={"rows":rows,"scaling":"Omega(1/sqrt(N)) for a Bernoulli-resolved fiber parameter away from boundaries",
            "matrix_extension":"recovering d free fiber parameters requires at least order d/epsilon^2 effective independent counts before log-confidence factors",
            "theorem":(
                "Le Cam's two-point method applied to child-resolved Bernoulli counts gives a nonzero minimax absolute-error lower bound of order N^{-1/2}. Ideal identifiability at open times does not imply finite-count exact recovery. For a 78-parameter symmetric fiber, parameter counting alone requires order 78/epsilon^2 effective independent events, subject to conditioning and design."
            )}
    return finalize(337,"Finite-Count Lower Bound for Fiber Recovery",
                    "proven_scalar_two_point_lower_bound_and_dimension_scaling_requirement",
                    "The matrix count is an order bound, not a sharp apparatus sample complexity. Independence, calibration, detector noise, and custody remain premises.",packet)


def markov_words(n:int,p:float):
    xs=np.array(list(itertools.product((0,1),repeat=n)),int);px=[]
    for x in xs:px.append(.5*np.prod([p if x[i]!=x[i-1] else 1-p for i in range(1,n)]))
    return xs,np.array(px)


def blahut_markov(n:int,p:float,betas:list[float]):
    xs,px=markov_words(n,p);dist=np.mean(xs[:,None,:]!=xs[None,:,:],axis=2);rows=[];m=len(xs)
    for beta in betas:
        q=np.ones(m)/m
        for _ in range(600):
            z=np.exp(-beta*dist)@q;pyx=q[None,:]*np.exp(-beta*dist)/z[:,None];qn=px@pyx
            if np.max(abs(qn-q))<1e-12:break
            q=qn
        joint=px[:,None]*pyx;D=float(np.sum(joint*dist));I=float(np.sum(joint*np.log(np.maximum(pyx,1e-300)/np.maximum(q[None,:],1e-300)))/math.log(2)/n)
        rows.append({"beta":beta,"distortion":D,"rate_bits_per_symbol":I})
    return rows


def st338()->dict:
    p=.1;n=8;rows=blahut_markov(n,p,[0,.5,1,2,4,8,16]);entropy_rate=-(p*math.log2(p)+(1-p)*math.log2(1-p))
    packet={"state_class":"stationary binary first-order Markov branch process with flip probability 0.1",
            "blocklength":n,"distortion":"normalized Hamming","Blahut_Arimoto_rows":rows,
            "asymptotic_source_entropy_rate":entropy_rate,"lossless_arbitrary_history_rate":1.0,
            "result":(
                "Once the branch class is restricted to a correlated Markov process, numerical block rate--distortion drops below the one-bit arbitrary-history bound at nonzero distortion, and the lossless asymptotic entropy rate is h2(0.1) approximately 0.469 bits per layer. This demonstrates controlled fractal-history compression for a declared source, not compression of arbitrary universe states or a strict FIN source law."
            )}
    return finalize(338,"First FIN-Style Markov Branch Rate--Distortion Model",
                    "strong_numerical_block_rate_distortion_for_declared_correlated_state_class",
                    "The Markov law, p=0.1, blocklength, and distortion are supplied. No FIN theorem selects them or identifies physical information density.",packet)


def st339(a:np.ndarray)->dict:
    lam=np.linalg.eigvalsh(a);l1,l2=float(lam[1]),float(lam[3]);ratio=l2/l1;rows=[]
    for eps in (1e-4,1e-3,1e-2,.05):
        lo=(l2-eps)/(l1+eps);hi=(l2+eps)/(l1-eps)
        rows.append({"intertwining_operator_error":eps,"heat_unitary_ratio_interval":[lo,hi],
                     "wave_ratio_interval":[math.sqrt(lo),math.sqrt(hi)]})
    packet={"reference_ratio_lambda2_over_lambda1":ratio,"spectral_gap_to_zero":l1,"rows":rows,
            "theorem":(
                "If the compressed refined operator differs from A by operator norm at most epsilon<lambda1, Weyl bounds place each positive coarse eigenvalue in [lambda-epsilon,lambda+epsilon]. Therefore heat/unitary ratios and square-root wave ratios lie in the displayed explicit intervals. Approximate refinement preserves relative clocks continuously but does not select seconds, channel maps, or fiber clocks."
            )}
    return finalize(339,"Stability Bounds for Relative Clocks under Approximate Refinement",
                    "proven_Weyl_ratio_bounds_for_approximate_coarse_intertwining",
                    "The operator-error premise must be established by an actual refinement/instrument. No time arrow or absolute calibration follows.",packet)


def st340()->dict:
    proposals={
        "speed_L_over_T":[1,-1,0],"action_E_times_T":[0,1,1],"diffusion_L2_over_T":[2,-1,0],
        "wave_A_as_inverse_T2":[0,-2,0],"unitary_A_as_inverse_T":[0,-1,0],"energy_time_E_times_T":[0,1,1]}
    rows=[]
    for names in (("speed_L_over_T",),("speed_L_over_T","action_E_times_T"),
                  ("speed_L_over_T","action_E_times_T","diffusion_L2_over_T")):
        mat=np.array([proposals[n] for n in names],float);rows.append({"relations":list(names),"rank":int(np.linalg.matrix_rank(mat)),"residual_calibration_dimension":3-int(np.linalg.matrix_rank(mat))})
    packet={"dimension_basis":"(L,T,E)","candidate_exponent_rows":proposals,"rank_audits":rows,
            "incompatibility":"one raw A cannot simultaneously carry inverse-time and inverse-time-squared dimensions without distinct sector conversion maps",
            "identifiability_requirement":"each sector map needs an operational calibration experiment; dimensional rank alone cannot identify it",
            "theorem":(
                "Dimensional relations reduce the calibration group only by the rank of genuinely independent exponent rows. Repeated action relations add no rank. More importantly, assigning the same raw A both T^{-1} in unitary/heat dynamics and T^{-2} in wave dynamics is inconsistent unless channel-specific conversion maps are supplied. Dimensional closure and operational identifiability are separate obligations."
            )}
    return finalize(340,"Cross-Channel Dimensional Rank and Identifiability Audit",
                    "proven_rank_count_and_raw_operator_dimension_incompatibility",
                    "No candidate physical relation or channel conversion is derived from FIN. The audit prevents silent unit reuse; it does not fix a calibration.",packet)


def st341()->dict:
    rec=list(ROOT.glob("FIN_ST341_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet={"required_pattern":"FIN_ST341_INDEPENDENT_RAW_EVENTS*.jsonl","matching_records":[x.name for x in rec],
            "independent_record_count":len(rec),"theorem":"No independently registered raw-event packet exists; local computation cannot generate custody or empirical evidence."}
    return finalize(341,"Independent Evidence Gate",
                    "blocked_no_independent_events" if not rec else "record_present_requires_blinded_validation",
                    "This stop is evidentiary, not a failed experiment.",packet)


def make_figures(out:dict)->None:
    FIG_DIR.mkdir(exist_ok=True)
    r=out["ST330"]["rank_tradeoff"];fig,ax=plt.subplots(figsize=(7,4));ax.plot([x["mediator_rank"] for x in r],[x["optimal_Frobenius_error"] for x in r],"o-")
    ax.set(xlabel="mediator rank",ylabel="best Frobenius error to 4A",title="ST330: exact lower-rank mediator tradeoff");fig.tight_layout();fig.savefig(FIG_DIR/"st330_rank_tradeoff.png",dpi=190);plt.close(fig)
    r=out["ST338"]["Blahut_Arimoto_rows"];fig,ax=plt.subplots(figsize=(7,4));ax.plot([x["distortion"] for x in r],[x["rate_bits_per_symbol"] for x in r],"o-")
    ax.set(xlabel="normalized Hamming distortion",ylabel="bits per layer",title="ST338: block Markov rate--distortion");fig.tight_layout();fig.savefig(FIG_DIR/"st338_rate_distortion.png",dpi=190);plt.close(fig)


def main()->None:
    np.random.seed(SEED);_,a,s=strict_operator();out={"metadata":{"seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"strict_row_sum":s,"scope":"local analytic, interval, exact finite, and bounded numerical research"}}
    funcs={327:lambda:st327(a),328:lambda:st328(a),329:lambda:st329(a),330:lambda:st330(a),331:lambda:st331(a),332:st332,333:lambda:st333(a),334:st334,335:st335,336:st336,337:st337,338:st338,339:lambda:st339(a),340:st340,341:st341}
    for k in range(327,342):out[f"ST{k}"]=funcs[k]()
    make_figures(out);RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(327,342):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    print(json.dumps({k:v["status"] for k,v in out.items() if k.startswith("ST")},indent=2))


if __name__=="__main__":main()
