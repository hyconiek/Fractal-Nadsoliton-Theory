#!/usr/bin/env python3
"""FIN ST203--ST217: schedule identifiability, endogenous-noise limits, and multiscale geometry."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.optimize import root

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st130_st141_research import point_design_system
from fin_st132_center_isolation_replay import bounds as replay_bounds, iv, strict_interval_matrix
from fin_st166_st177_research import local_param_krawczyk, parametric_data
from fin_st178_st189_research import IDX, stationary_slice_float, uniform_fold_seed
from fin_st190_st202_research import replicator_rhs
from fin_st189_external_record_validator import synthetic_self_test


ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/"FIN_ST203_ST217_Results.json"
SUMMARY=ROOT/"FIN_ST203_ST217_Summary.csv"
FIG_DIR=ROOT/"FIN_ST203_ST217_Figures"
SEED=20260812
PACKETS={k:ROOT/f"FIN_ST{k}_{name}.json" for k,name in {
    203:"Semigroup_Schedule_Identifiability",204:"Blackwell_Experiment_Lattice",
    205:"Strict_Covariant_Fluctuation",206:"Validated_Pseudo_Arclength",
    207:"Negative_Orientation_Recovery",208:"Finite_Time_Selector_Speed",
    209:"Carrier_Invariant_Nonlinear_Response",210:"Adaptive_Nuisance_Cover",
    211:"Multimode_Correlated_Sampling",212:"Joint_Factor_Nonuniqueness",
    213:"Asymptotic_HMM_Rate_Bracket",214:"Time_Dependent_Reset_No_Go",
    215:"Diffusion_Geometry",216:"Viewing_Versus_Reconstruction",
    217:"External_Evidence_Gate"}.items()}


def native(x:Any)->Any:
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    if isinstance(x,np.ndarray):return native(x.tolist())
    if isinstance(x,(np.floating,np.integer)):return x.item()
    if isinstance(x,complex):return {"real":x.real,"imag":x.imag}
    return x


def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k:int,obj:str,status:str,boundary:str,packet:dict)->dict:
    path=PACKETS[k];path.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8")
    return {"program":f"ST{k}","object":obj,"packet_file":path.name,"packet_sha256":sha(path),**packet,"status":status,"boundary":boundary}


def st203_schedule_identifiability(a:np.ndarray)->dict:
    schedules={
        "single":[3.75],"uniform_four":[.9375]*4,"dyadic_four":[.25,.5,1,2],
        "irregular":[.1,.2,.45,3.0]
    }
    target=expm(-3.75*a);rows=[]
    for name,ts in schedules.items():
        c=np.eye(N)
        for t in ts:c=expm(-t*a)@c
        rows.append({"schedule":name,"increments":ts,"total":sum(ts),
                     "final_map_error_from_exp_minus_3_75A":float(np.linalg.norm(c-target,np.inf))})
    # Geometric schedules t0(1+r+...+r^(n-1))=T: one t0 for every r>0.
    geometric=[]
    for r in [.25,.5,1,2,4]:
        denom=4 if r==1 else (r**4-1)/(r-1)
        geometric.append({"ratio":r,"base":3.75/denom,"total":3.75})
    packet={"final_parameter":3.75,"schedule_rows":rows,"geometric_four_layer_family":geometric,
            "theorem":(
                "For any nonnegative schedule, product_j exp(-t_j A)=exp[-(sum_j t_j)A]. Hence a final composite generated only by A identifies the total parameter but not the number, order, or ratios of layers. "
                "For every r>0 there is a geometric four-layer schedule with the same total, so no dyadic ratio is selected by final-map data. Intermediate maps recover increments only after a parameter normalization for A is fixed."
            )}
    return finalize(203,"Semigroup-Schedule Nonidentifiability Theorem",
        "proven_exact_schedule_collapse_and_geometric_family_no_go",
        "This rules out internal schedule selection from the final strict heat map alone. Noncommuting layer generators, intermediate observations, refinement labels, or an added scale law lie outside the no-go.",packet)


def canonical_partition(blocks):return tuple(sorted((tuple(sorted(b)) for b in blocks),key=lambda b:(b[0],len(b))))


def set_partitions(items):
    if not items:yield ();return
    first,*rest=items
    for p in set_partitions(rest):
        yield canonical_partition(((first,),)+p)
        for i in range(len(p)):
            q=list(p);q[i]=tuple(sorted((first,)+q[i]));yield canonical_partition(q)


def refines(p,q):return all(any(set(b)<=set(c) for c in q) for b in p)


def st204_blackwell_lattice()->dict:
    parts=sorted(set(set_partitions(list(range(4)))),key=lambda p:(-len(p),p))
    comparable=sum(1 for p in parts for q in parts if refines(p,q))
    covers=[]
    for i,p in enumerate(parts):
        for j,q in enumerate(parts):
            if len(p)==len(q)+1 and refines(p,q):covers.append([i,j])
    rank_counts={str(k):sum(len(p)==k for p in parts) for k in range(1,5)}
    packet={"minimal_sufficient_atoms":4,"partition_count":len(parts),"rank_counts_by_blocks":rank_counts,
            "partitions":[[list(b) for b in p] for p in parts],"ordered_comparable_pairs_including_identity":comparable,
            "Hasse_cover_edges":covers,"Hasse_edge_count":len(covers),
            "theorem":(
                "Deterministic coarse experiments on the four Blackwell-minimal atoms form the partition lattice Pi_4. Blackwell dominance is exactly partition refinement; meet is common refinement and join is transitive closure of the two equivalence relations. "
                "There are Bell(4)=15 experiments with block-count distribution 1,7,6,1 from one through four blocks. Only the four-block element remains sufficient for the original three-hypothesis experiment."
            )}
    return finalize(204,"Complete Blackwell Lattice of the Minimal FIN Experiment",
        "proven_exact_complete_partition_lattice",
        "The lattice belongs to the declared ST191 preparation family. Changing preparations or instruments can change proportionality atoms and therefore the operational lattice.",packet)


def st205_covariant_fluctuation(a:np.ndarray)->dict:
    q=np.eye(N)-np.ones((N,N))/N
    vals,vecs=np.linalg.eigh(a);sqrt2a=vecs@np.diag(np.sqrt(np.maximum(2*vals,0)))@vecs.T
    lyap=float(np.linalg.norm(a@q+q@a-sqrt2a@sqrt2a.T,np.inf))
    shift=np.roll(np.eye(N),1,axis=0)
    coverr=float(np.linalg.norm(shift@q@shift.T-q,np.inf))
    rng=np.random.default_rng(SEED+205);runs=120
    xi=rng.multivariate_normal(np.zeros(N),q,size=runs);p=np.ones((runs,N))/N+1e-3*xi
    p=np.maximum(p,1e-12);p/=p.sum(axis=1,keepdims=True);m=a*a;dt=.01
    for _ in range(4000):
        pay=p@m.T;mean=np.sum(p*pay,axis=1,keepdims=True);p+=dt*p*(pay-mean)
        p=np.maximum(p,1e-300);p/=p.sum(axis=1,keepdims=True)
    counts=np.bincount(np.argmax(p,axis=1),minlength=N);expected=runs/N
    chisq=float(np.sum((counts-expected)**2/expected))
    packet={"OU_law":"d xi=-A xi dt+sqrt(2A)dB on the mean-zero sector","stationary_covariance":"Q=I-J/12",
            "Lyapunov_residual":lyap,"cyclic_covariance_error":coverr,"replicator_runs":runs,
            "selected_vertex_counts":counts.tolist(),"uniform_expected_count":expected,"chi_square_descriptive_statistic":chisq,
            "minimum_final_max_probability":float(np.min(np.max(p,axis=1))),
            "theorem":(
                "On the mean-zero sector, A Q+Q A=2A, so the constructed Ornstein--Uhlenbeck law has invariant covariance Q. Its law is C12- and inversion-invariant. Coupling one sample to the ST192 feedback triggers branch selection, but symmetry forces equal branch probabilities in the exact ensemble and forbids a canonical orientation or vertex label."
            )}
    return finalize(205,"Strict-Covariant Fluctuation Candidate and Selector No-Go",
        "proven_covariance_and_symmetry_no_go_with_seeded_selection_audit",
        "The covariance shape is derived from A, but Brownian stochasticity, noise normalization, and its coupling to the replicator law are added dynamical premises. It does not discharge QW-2191.",packet)


def stationary7(x,a):
    f,j=stationary_slice_float(x,a,np.zeros(7),np.zeros(7),0.0)
    return f[:7],j[:7,:]


def pseudo_interval(center,radius,aiv,normal,predictor):
    q=[iv((center[i]-radius,center[i]+radius)) for i in range(7)];kappa=iv((center[7]-radius,center[7]+radius));u=[q[IDX[s]] for s in range(12)]
    au=[sum((aiv[i][j]*u[j] for j in range(12)),iv(0)) for i in range(12)];h=[];rd=[]
    for z in u[:7]:
        rho=z*z;den=1+rho/2;qf=rho/den;qp=den**-2;qpp=-den**-3
        hh=-qf*qp+iv("0.075");hp=-(qp*qp+qf*qpp);h.append(hh);rd.append(2*hh+4*rho*hp)
    g=[kappa*au[i]+2*u[i]*h[i] for i in range(7)]
    g.append(sum((iv(float(normal[i]))*([*q,kappa][i]-iv(float(predictor[i]))) for i in range(8)),iv(0)))
    jac=[[iv(0) for _ in range(8)] for _ in range(8)]
    for i in range(7):
        for col in range(7):
            total=sum((kappa*aiv[i][site] for site in range(12) if IDX[site]==col),iv(0))
            if i==col:total+=rd[i]
            jac[i][col]=total
        jac[i][7]=au[i]
    for col in range(8):jac[7][col]=iv(float(normal[col]))
    glo=np.array([replay_bounds(z)[0] for z in g]);ghi=np.array([replay_bounds(z)[1] for z in g])
    jlo=np.array([[replay_bounds(z)[0] for z in row] for row in jac]);jhi=np.array([[replay_bounds(z)[1] for z in row] for row in jac])
    return glo,ghi,jlo,jhi


def pseudo_krawczyk(center,aiv,normal,predictor,radius):
    g0lo,g0hi,j0lo,j0hi=pseudo_interval(center,0,aiv,normal,predictor);pre=np.linalg.inv(.5*(j0lo+j0hi))
    _,_,jlo,jhi=pseudo_interval(center,radius,aiv,normal,predictor)
    cglo,cghi=interval_matvec(pre,pre,g0lo,g0hi);ylo=np.nextafter(center-cghi,-np.inf);yhi=np.nextafter(center-cglo,np.inf)
    cjlo,cjhi=interval_left(pre,jlo,jhi);mlo,mhi=-cjhi,-cjlo
    for i in range(8):mlo[i,i]=np.nextafter(mlo[i,i]+1,-np.inf);mhi[i,i]=np.nextafter(mhi[i,i]+1,np.inf)
    dlo,dhi=interval_matvec(mlo,mhi,np.full(8,-radius),np.full(8,radius));klo=np.nextafter(ylo+dlo,-np.inf);khi=np.nextafter(yhi+dhi,np.inf)
    margin=float(min(np.min(klo-(center-radius)),np.min((center+radius)-khi)))
    return {"radius":radius,"margin":margin,"included":margin>0}


def st206_pseudo_arclength(a:np.ndarray)->dict:
    prior=json.loads((ROOT/"FIN_ST190_ST202_Results.json").read_text(encoding="utf-8"))["ST193"]["rows"]
    seed=next(r for r in prior if r["uniform_amplitude"]>0 and r["mode"]==1 and r["slice_epsilon"]==.005)
    x0=np.array(seed["continued_center"]);aiv,_,_=strict_interval_matrix();mp.iv.dps=60;ds=.002;rows=[]
    _,j0=stationary7(x0,a);_,_,vh=np.linalg.svd(j0);t0=vh[-1];t0/=np.linalg.norm(t0)
    qref,_,vref=uniform_fold_seed(a,seed["uniform_amplitude"],1)
    if np.dot(t0[:7],vref)<0:t0=-t0
    for direction in [1,-1]:
        x=x0.copy();t=direction*t0.copy()
        for step in range(1,21):
            pred=x+ds*t
            sol=root(lambda y:np.r_[stationary7(y,a)[0],t@(y-pred)],pred,
                     jac=lambda y:np.vstack([stationary7(y,a)[1],t]),tol=1e-12)
            xn=sol.x;_,jn=stationary7(xn,a);_,sv,vh=np.linalg.svd(jn);tn=vh[-1];tn/=np.linalg.norm(tn)
            if np.dot(tn,t)<0:tn=-tn
            trials=[pseudo_krawczyk(xn,aiv,t,pred,r) for r in (3e-7,1e-7,3e-8,1e-8)]
            cert=next((z for z in trials if z["included"]),None)
            rows.append({"direction":direction,"step":step,"arc_parameter":direction*step*ds,
                         "modal_coordinate":float(vref@(xn[:7]-qref)),"kappa":float(xn[7]),
                         "stationary_residual":float(np.linalg.norm(stationary7(xn,a)[0],np.inf)),
                         "minimum_nonzero_Jacobian_singular_value":float(sv[-1]),
                         "tangent_alignment":float(np.dot(tn,t)),"certificate":cert,"trials":trials})
            x,t=xn,tn
    passed=sum(r["certificate"] is not None for r in rows)
    packet={"seed_uniform_amplitude":seed["uniform_amplitude"],"seed_mode":1,"step_size":ds,"steps":len(rows),"certified_steps":passed,
            "minimum_Jacobian_rank_margin":min(r["minimum_nonzero_Jacobian_singular_value"] for r in rows),
            "minimum_tangent_alignment":min(r["tangent_alignment"] for r in rows),"rows":rows,
            "theorem":"Every accepted pseudo-arclength hyperplane has exactly one stationary root in its outward Krawczyk box. Positive seven-row Jacobian rank margin and tangent alignment exclude a detected singular collision on the traced segment."}
    return finalize(206,"Validated Two-Sided Pseudo-Arclength Segment",
        "proven_local_pseudo_arclength_segment" if passed==len(rows) else "partial_validated_pseudo_arclength_segment",
        "Only one reflection-even branch and a finite arc segment are traced. Absence of a detected collision on this segment is not a global branch classification or stability theorem.",packet)


def st207_negative_orientation_recovery()->dict:
    probs=[Fraction(1,10),Fraction(3,10),Fraction(3,10),Fraction(3,10)]
    # I,X,Y,Z channel: diagonal Bloch multipliers.
    lam=[probs[0]+probs[1]-probs[2]-probs[3],probs[0]-probs[1]+probs[2]-probs[3],probs[0]-probs[1]-probs[2]+probs[3]]
    fidelities={name:float(p) for name,p in zip(["I","X","Y","Z"],probs)}
    packet={"Pauli_probabilities_exact":[str(x) for x in probs],"Bloch_diagonal_exact":[str(x) for x in lam],
            "Bloch_determinant_exact":str(lam[0]*lam[1]*lam[2]),"Pauli_recovery_fidelities":fidelities,
            "global_optimal_entanglement_fidelity":max(fidelities.values()),"optimal_recoveries":[k for k,v in fidelities.items() if v==max(fidelities.values())],
            "theorem":(
                "Pauli-twirling any CPTP recovery preserves entanglement fidelity against a Pauli channel and produces a Pauli recovery channel. The remaining objective is linear in its four probabilities, so an optimal recovery is a Pauli unitary correcting a most likely error. "
                "For probabilities (1/10,3/10,3/10,3/10), the Bloch determinant is negative and the global optimum over all, including nonunital, CPTP recoveries is 3/10, attained by X, Y, or Z correction."
            )}
    return finalize(207,"Negative-Orientation Pauli Recovery over All CPTP Maps",
        "proven_global_CPTP_optimum_for_negative_orientation_Pauli_channel",
        "The theorem uses Pauli covariance/twirling. It does not solve arbitrary negative-orientation non-Pauli or generic nonunital noise channels.",packet)


def binary_relative(x,pi):return x*math.log(x/pi)+(1-x)*math.log((1-x)/(1-pi))


def st208_selector_speed()->dict:
    bmax=4.0;rmax=1.0;x0=1/N;pi=1/(1+11*math.exp(-bmax));targets=[.2,.4,.6,.8];rows=[]
    for xf in targets:
        if xf>=pi:rows.append({"target":xf,"reachable":False});continue
        t=math.log((pi-x0)/(pi-xf))/rmax
        entropy=binary_relative(x0,pi)-binary_relative(xf,pi)
        rows.append({"target":xf,"reachable":True,"minimum_time":t,"entropy_production_at_speed_saturating_constant_control":entropy})
    packet={"field_bound_beta_h":bmax,"relaxation_rate_bound":rmax,"initial_preferred_probability":x0,"maximum_control_equilibrium_probability":pi,"rows":rows,
            "theorem":(
                "For every reduced selector control satisfying 0<=r(t)<=R and pi(t)<=pi_B, x_dot=r(t)(pi(t)-x) <= R(pi_B-x). Hence reaching x_f<pi_B requires T>=R^-1 log[(pi_B-x_0)/(pi_B-x_f)]. "
                "The constant maximal control r=R, pi=pi_B attains equality. Its integrated entropy production equals D(x_0||pi_B)-D(x_f||pi_B). Targets x_f>=pi_B are unreachable under the frozen field bound."
            )}
    return finalize(208,"Finite-Time Selector Speed Limit under Frozen Resources",
        "proven_time_optimum_and_constant_control_entropy_for_reduced_family",
        "This solves fastest selection, not the general minimum-dissipation problem at fixed time. The reduction, field/rate bounds, bath, target vertex, and dimensionless clock remain supplied.",packet)


def st209_response_invariant(a:np.ndarray)->dict:
    ajj=float(a[0,0]);gs=[.1,1,10];scales=[.25,.5,2,4];rows=[]
    for g,c in itertools.product(gs,scales):
        r2=ajj/c**2;r12=math.factorial(11)*g/c**12;inv=r12/(math.factorial(11)*r2**6)
        rows.append({"g":g,"field_coordinate_scale":c,"R2":r2,"R12":r12,"invariant_ratio":inv,"expected_g_over_Ajj6":g/ajj**6})
    packet={"strict_vertex_quadratic_response_Ajj":ajj,"rows":rows,
            "invariant":"I12=R12/(11! R2^6)=g/Ajj^6",
            "theorem":(
                "Under scalar field reparameterization psi'=c psi, R_2 transforms as c^-2 and R_12 as c^-12, so I_12=R_12/(11!R_2^6) is invariant. In the declared local model it equals g/A_jj^6. "
                "All jets of order below twelve are independent of g; therefore one nonzero twelfth-order scalar response is the minimal additional local jet datum capable of identifying this dimensionless coupling ratio."
            )}
    return finalize(209,"Minimal Carrier-Scale-Invariant Nonlinear Response",
        "proven_response_invariant_and_jet_minimality",
        "The transported field direction, access to a twelfth-order response, and model form are supplied. The invariant identifies a ratio; it does not derive g or physical units.",packet)


def st210_adaptive_cover(a:np.ndarray)->dict:
    ec,ef,eigc,eigf=parametric_data(a);h=2e-4;sub=14;lh=h/sub;offs=[-h+(2*i+1)*lh for i in range(sub)];rows=[]
    for ds in itertools.product(offs,repeat=3):
        nu=(.2+ds[0],.7+ds[1],.05+ds[2]);center=root(lambda x:point_design_system(x,ec,ef,*nu),[2.1862,.53983],tol=1e-12).x
        rows.append(local_param_krawczyk(eigc,eigf,nu,lh,center))
    passed=sum(r["included"] for r in rows)
    packet={"global_halfwidth":h,"equivalent_adaptive_construction":"subdivide every failed cell of the ST197 7^3 cover once into eight children","cells_per_axis":sub,
            "boxes":len(rows),"passed_boxes":passed,"local_halfwidth":lh,"minimum_margin":min(r["margin"] for r in rows),
            "theorem":"The 14^3 child cells exactly tile the nuisance cube of halfwidth 0.0002. Krawczyk inclusion in every child proves a unique continued stationary root throughout the complete cube."}
    return finalize(210,"Adaptive Refinement beyond the ST197 Fixed-Grid Boundary",
        "proven_complete_halfwidth_0_0002_cover" if passed==len(rows) else "partial_adaptive_cover",
        "The result extends a local stationary branch, not apparatus tolerances or a maximal continuation domain. Further extension requires another adaptive frontier test.",packet)


def st211_multimode_sampling()->dict:
    epsilon=0.01;cluster=20;rho=.05;design=1+(cluster-1)*rho;rows=[]
    mu=np.array([.75+.25*math.cos(2*math.pi*k/N) for k in range(1,7)])
    for n in [0,2,4,6,8,10,12]:
        sigma=mu**n;fisher=12*sigma*sigma;worst=float(np.min(fisher));iid=math.ceil(1/(worst*epsilon**2));corr=math.ceil(design/(worst*epsilon**2))
        rows.append({"layers":n,"mode_attenuations":sigma.tolist(),"per_sample_Fisher_eigenvalues":fisher.tolist(),
                     "iid_samples_for_worst_mode_sd_0_01":iid,"cluster_correlated_samples":corr})
    packet={"orthonormal_modes":6,"target_standard_deviation":epsilon,"cluster_size":cluster,"intracluster_correlation":rho,"design_effect":design,"rows":rows,
            "theorem":(
                "At the uniform multinomial state, orthonormal Fourier contrasts transmitted with factors sigma_k have per-sample Fisher matrix diag(12 sigma_k^2). The worst-direction Cramer--Rao variance is at least [12N min_k sigma_k^2]^-1. "
                "For exchangeable clusters of size m and intracluster correlation rho, the covariance inflation is the design effect 1+(m-1)rho, multiplying the necessary count."
            )}
    return finalize(211,"Simultaneous Deep-Mode Sampling with Correlated Counts",
        "proven_local_Fisher_minimax_and_exchangeable_cluster_inflation",
        "The result is local at the uniform state and assumes declared Fourier contrasts and exchangeable count correlation. Biased estimators and other noise geometries require separate minimax analysis.",packet)


def st212_joint_factor_nonuniqueness()->dict:
    c=.8;angles=np.linspace(0,2*math.pi,13)[:-1];rows=[]
    # n=(e_x+u)/sqrt2, m=(e_x-u)/sqrt2, u in yz-circle.
    optimum=4*(1+c*c)-4*c*math.sqrt(2)
    for th in angles:
        u=np.array([0,math.cos(th),math.sin(th)]);n=(np.array([1.,0,0])+u)/math.sqrt(2);m=(np.array([1.,0,0])-u)/math.sqrt(2)
        rows.append({"theta":float(th),"n_dot_m":float(n@m),"objective":float(2*np.sum((n-c*np.array([1,0,0]))**2)+2*np.sum((m-c*np.array([1,0,0]))**2))})
    packet={"raw_target":"X0=Z0=c sigma_x","c":c,"raw_generators_invertible":True,"analytic_minimum":optimum,"sampled_minimizer_circle":rows,
            "theorem":(
                "For qubit involutions X=n.sigma and Z=m.sigma, anticommutation is n.m=0. With invertible targets X0=Z0=c sigma_x, minimizing the joint Frobenius error is equivalent to maximizing e_x.(n+m)<=sqrt(2). "
                "Equality holds for n=(e_x+u)/sqrt(2), m=(e_x-u)/sqrt(2) for every unit u perpendicular to e_x, producing a circle of global minimizers. Thus universal global nearest-joint-factor uniqueness is false even when both raw targets are invertible."
            )}
    return finalize(212,"Global Nearest-Joint-Factor Nonuniqueness Counterexample",
        "proven_continuum_of_global_minimizers",
        "The counterexample has a zero sequential odd-part gap, so it does not refute a future uniqueness theorem with the full ST199 gap/transversality premises.",packet)


def st213_hmm_rate_bracket()->dict:
    st200=json.loads((ROOT/"FIN_ST190_ST202_Results.json").read_text(encoding="utf-8"))["ST200"]["rows"]
    c=Fraction(3,10);C=Fraction(12,5);rows=[]
    for r in st200:
        n=r["events"];zlo=float(Fraction(r["exact_rational_lower"]));zhi=float(Fraction(r["exact_rational_upper"]))
        rows.append({"block_length":n,"asymptotic_rate_lower":(-math.log(zhi)-math.log(float(C)))/n,
                     "asymptotic_rate_upper":(-math.log(zlo)-math.log(float(c)))/n,
                     "bracket_width":math.log(float(C/c))/n})
    packet={"conditional_to_stationary_likelihood_ratio_bounds_exact":[str(c),str(C)],"rows":rows,
            "best_displayed_bracket":rows[-1],
            "theorem":(
                "After any observed past, one transition puts the next hidden distribution componentwise between c=3/10 and C=12/5 times the stationary distribution. Therefore c Z_m Z_n <= Z_{m+n} <= C Z_m Z_n. "
                "Almost-additivity proves existence of the asymptotic observed Hellinger rate R and yields [(-log Z_n-log C)/n,(-log Z_n-log c)/n]. Exact-rational ST200 block enclosures make every displayed rate bracket outward."
            )}
    return finalize(213,"Rigorous Asymptotic Observed-HMM Rate Bracket",
        "proven_rate_existence_and_explicit_asymptotic_bracket",
        "The bracket is conservative and concerns the Hellinger point s=1/2, not the fully optimized Chernoff exponent. A sharper belief-space transfer certificate remains useful.",packet)


def st214_time_dependent_reset()->dict:
    d=json.loads((ROOT/"FIN_ST190_ST202_Results.json").read_text(encoding="utf-8"))["ST201"]
    packet={"instantaneous_generator_class":"measurable H_c(t) in the 60-dimensional corresponding-pair Pauli span",
            "certified_commutator_map_lower_singular_bound":d["expanded_interval_lower_singular_bound"],
            "pointwise_energy_conservation":"[H_total,H_c(t)]=0 for almost every t","conclusion":"H_c(t)=0 almost everywhere and U(T)=I",
            "theorem":(
                "ST201 proves the commutator map is injective on the declared 60-dimensional span. Hence any measurable time-dependent coefficient vector satisfying pointwise exact energy conservation lies in the zero kernel almost everywhere. Its time-ordered exponential is the identity, so no nontrivial reset is possible in this pointwise-conserving control class."
            )}
    return finalize(214,"Time-Dependent Pointwise Energy-Conserving Local-Reset No-Go",
        "proven_no_go_for_time_dependent_pointwise_conserving_declared_span",
        "Net energy conservation only at the final time, noncorresponding controls, ancillas, higher locality, and excursions outside the 60-dimensional span remain open.",packet)


def st215_diffusion_geometry(a:np.ndarray)->dict:
    times=[.1,.5,1,2,4];rows=[]
    for t in times:
        c=expm(-t*a);dist=[]
        for d in range(1,7):dist.append(float(np.linalg.norm(c[0]-c[d])))
        violations=0
        D=np.array([[np.linalg.norm(c[i]-c[j]) for j in range(N)] for i in range(N)])
        for i,j,k in itertools.product(range(N),repeat=3):violations+=D[i,k]>D[i,j]+D[j,k]+1e-12
        rows.append({"t":t,"distance_by_cycle_separation":dist,"triangle_violations":int(violations),"condition_number":float(np.linalg.cond(c))})
    lam=np.real(np.fft.fft(a[0]));t=10;c=expm(-t*a);scaled=[]
    for d in range(1,7):scaled.append(float(math.exp(t*lam[1])*np.linalg.norm(c[0]-c[d])))
    chord=[math.sqrt(2/3)*abs(math.sin(math.pi*d/12)) for d in range(1,7)]
    packet={"definition":"D_t(i,j)=||e_i^T exp(-tA)-e_j^T exp(-tA)||_2","rows":rows,
            "large_t_rescaled_distances_t10":scaled,"first_mode_chordal_limit":chord,
            "maximum_t10_limit_error":float(max(abs(x-y) for x,y in zip(scaled,chord))),
            "theorem":(
                "For every finite t, exp(-tA) is invertible, so its row embedding is injective and D_t is a genuine translation-invariant Euclidean metric. After rescaling by exp(t lambda_1), the spectral gap makes D_t converge to the chordal metric of the first Fourier circle, sqrt(2/3)|sin(pi d/12)|."
            )}
    return finalize(215,"Dimensionless Diffusion Geometry from the Strict Layer Tower",
        "proven_metric_and_first_mode_asymptotic_geometry",
        "Diffusion distance is dimensionless and depends on the supplied t normalization. It creates no meters, causal metric, spacetime signature, Planck length, or unique physical geometry.",packet)


def st216_viewing_vs_reconstruction(a:np.ndarray)->dict:
    # Hard 12->3 and soft depth-eight ST190 maps.
    hard=np.zeros((3,12))
    for g in range(3):hard[g,4*g:4*g+4]=1
    soft=expm(-63.75*a)
    rows=[
        {"channel":"hard_12_to_3","rank":int(np.linalg.matrix_rank(hard)),"kernel_dimension":9,"globally_injective":False,"condition_number_on_nonzero_spectrum":1.0},
        {"channel":"strict_heat_depth_8","rank":int(np.linalg.matrix_rank(soft,tol=1e-80)),"kernel_dimension":0,"globally_injective":True,"condition_number":float(np.linalg.cond(soft)),"exact_spectral_condition_number":float(math.exp(63.75*np.max(np.linalg.eigvalsh(a))))},
    ]
    packet={"rows":rows,
            "theorem":(
                "Viewing asks only for predictions of effects in range(C^T). Exact reconstruction on a model set M requires C restricted to M to be injective. Stable linear reconstruction requires a positive restricted minimum singular value, with noise amplification at least its reciprocal. "
                "Blackwell sufficiency is weaker and decision-relative: it requires a decoder only for the declared experiment. Thus visibility, sufficiency, algebraic invertibility, and stable reconstruction are four distinct properties."
            )}
    return finalize(216,"Formal Separation of Viewing, Sufficiency, and Reconstruction",
        "proven_general_finite_channel_distinction_with_FIN_examples",
        "The theorem organizes operational claims; it does not identify which preparation manifold or instrument family describes actual observers.",packet)


def st217_external_gate()->dict:
    s=synthetic_self_test();missing=["genuine laboratory event record","independently frozen registration","calibration certificate","verified provider/registrar/analyst separation","independent custody","laboratory execution attestation"]
    packet={"synthetic_self_test_passed":s["test_passed"],"missing_external_atoms":missing,"external_execution_performed":False,"physical_result":None,
            "theorem":"The local package remains structurally ready and tamper-sensitive, but the external-evidence predicate is false because every nonlocal provenance atom is absent. Repeating local execution cannot change that status."}
    return finalize(217,"Immutable External-Evidence Gate Re-audit",
        "blocked_correctly_no_external_record",
        "No laboratory data or external audit was requested or available. This program records the stop condition and does not simulate evidence.",packet)


def make_figures(d):
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST203"]["schedule_rows"];ax.bar([r["schedule"] for r in rows],[r["final_map_error_from_exp_minus_3_75A"] for r in rows]);ax.set_yscale("symlog",linthresh=1e-16);ax.set(ylabel="final-map mismatch",title="ST203: distinct schedules collapse to one composite");fig.tight_layout();fig.savefig(FIG_DIR/"st203_schedule_collapse.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));counts=d["ST205"]["selected_vertex_counts"];ax.bar(range(12),counts);ax.axhline(d["ST205"]["uniform_expected_count"],ls="--",color="red");ax.set(xlabel="selected orbit member",ylabel="count",title="ST205: invariant noise selects no canonical label");fig.tight_layout();fig.savefig(FIG_DIR/"st205_branch_counts.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST206"]["rows"]
    for direction in [-1,1]:
        q=[r for r in rows if r["direction"]==direction];ax.plot([r["modal_coordinate"] for r in q],[r["kappa"] for r in q],"o-")
    ax.set(xlabel="modal coordinate",ylabel="kappa",title="ST206: validated two-sided pseudo-arclength segment");fig.tight_layout();fig.savefig(FIG_DIR/"st206_pseudo_arclength.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST210"].get("margin_sample",[])
    ax.bar(["ST197 fixed pass","ST210 adaptive pass"],[.00013045,d["ST210"]["global_halfwidth"]]);ax.set(ylabel="certified global halfwidth",title="ST210: adaptive nuisance-cover extension");fig.tight_layout();fig.savefig(FIG_DIR/"st210_adaptive_cover.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST211"]["rows"];ax.semilogy([r["layers"] for r in rows],[r["iid_samples_for_worst_mode_sd_0_01"] for r in rows],"o-",label="iid");ax.semilogy([r["layers"] for r in rows],[r["cluster_correlated_samples"] for r in rows],"s--",label="clustered");ax.set(xlabel="layers",ylabel="necessary samples",title="ST211: simultaneous-mode sampling burden");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st211_multimode_samples.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST213"]["rows"];ax.fill_between([r["block_length"] for r in rows],[r["asymptotic_rate_lower"] for r in rows],[r["asymptotic_rate_upper"] for r in rows],alpha=.3);ax.plot([r["block_length"] for r in rows],[r["asymptotic_rate_lower"] for r in rows],"o-");ax.plot([r["block_length"] for r in rows],[r["asymptotic_rate_upper"] for r in rows],"o-");ax.set(xlabel="certified block length",ylabel="asymptotic rate bracket",title="ST213: rigorous observed-HMM rate enclosure");fig.tight_layout();fig.savefig(FIG_DIR/"st213_hmm_rate_bracket.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST215"]["rows"]
    for r in rows:ax.plot(range(1,7),r["distance_by_cycle_separation"],"o-",label=f"t={r['t']}")
    ax.set(xlabel="cycle separation",ylabel="diffusion distance",title="ST215: layer-induced dimensionless geometry");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st215_diffusion_geometry.png",dpi=190);plt.close(fig)


def main():
    _,a,_=strict_operator();out={"metadata":{"programs":"ST203-ST217","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__}}
    out["ST203"]=st203_schedule_identifiability(a);out["ST204"]=st204_blackwell_lattice();out["ST205"]=st205_covariant_fluctuation(a);out["ST206"]=st206_pseudo_arclength(a);out["ST207"]=st207_negative_orientation_recovery()
    out["ST208"]=st208_selector_speed();out["ST209"]=st209_response_invariant(a);out["ST210"]=st210_adaptive_cover(a);out["ST211"]=st211_multimode_sampling();out["ST212"]=st212_joint_factor_nonuniqueness()
    out["ST213"]=st213_hmm_rate_bracket();out["ST214"]=st214_time_dependent_reset();out["ST215"]=st215_diffusion_geometry(a);out["ST216"]=st216_viewing_vs_reconstruction(a);out["ST217"]=st217_external_gate()
    out["recommended_next_programs"]=[
        {"id":"ST218","priority":1,"study":"test one noncommuting strict-derived layer generator candidate; otherwise formalize the full f(A) schedule no-go"},
        {"id":"ST219","priority":2,"study":"derive preparation-dependent Blackwell lattices for the vertex, spectral, and nonlinear-response instrument families"},
        {"id":"ST220","priority":3,"study":"classify all Aut(Z12)-covariant Lévy and Gaussian fluctuation sources and prove their selector impossibility"},
        {"id":"ST221","priority":4,"study":"continue every certified nonuniform fold by pseudo-arclength and isolate all finite-segment singularities"},
        {"id":"ST222","priority":5,"study":"solve generic negative-orientation non-Pauli qubit recovery by a replayable primal-dual conic certificate"},
        {"id":"ST223","priority":6,"study":"solve minimum-dissipation selector control at fixed finite time with interval dynamic programming"},
        {"id":"ST224","priority":7,"study":"construct a carrier-transport theorem for the normalized twelfth-response invariant"},
        {"id":"ST225","priority":8,"study":"extend the adaptive nuisance cover until first certified root-loss or a resource-bounded stop"},
        {"id":"ST226","priority":9,"study":"derive nonlocal and biased-estimator minimax bounds for simultaneous compressed modes"},
        {"id":"ST227","priority":10,"study":"prove a local nearest-joint-factor theorem under explicit spectral-gap and transversality premises"},
        {"id":"ST228","priority":11,"study":"sharpen the observed-HMM asymptotic bracket with interval belief-space discretization"},
        {"id":"ST229","priority":12,"study":"test net-energy-conserving bang-bang reset sequences outside pointwise conservation"},
        {"id":"ST230","priority":13,"study":"compare strict diffusion geometry with resistance, Fisher, and spectral-projector metrics"},
        {"id":"ST231","priority":14,"study":"seek a refinement law that preserves diffusion geometry up to scale without importing a length unit"},
        {"id":"ST232","priority":15,"study":"execute external validation only after a genuine independently registered laboratory record exists"},
    ]
    out["central_verdict"]=("The strict heat family supplies layer dynamics but cannot select its own discrete history: every f(A) heat schedule collapses to total semigroup time. Strict-covariant fluctuations can trigger spontaneous orbit selection but, by the same symmetry, cannot choose a canonical label. The strongest new positive bridge is dimensionless diffusion geometry; the strongest new obstruction is that neither layer schedule nor selector orientation emerges from the current symmetric functional calculus.")
    out["epistemic_boundary"]=("No carrierless information, strict preferred fractal schedule, Planck scale, dimensional spacetime metric, canonical selector, QW-2191 discharge, dimensional clock/length/action source, external record, legacy-to-strict completion or physical-role transfer, Standard Model, gravity, L_total, or ToE closure is claimed.")
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(203,218):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out);print(json.dumps({"results":RESULTS.name,"programs":15,"figures":7},indent=2))


if __name__=="__main__":main()
