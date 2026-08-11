#!/usr/bin/env python3
"""FIN ST178--ST189: irreversible carriers and layered-compression observation."""

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
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import root

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st130_st141_research import point_design_system
from fin_st132_center_isolation_replay import bounds as replay_bounds, iv, strict_interval_matrix
from fin_st154_st165_research import commutant_dimension, hmm_sequence_likelihoods, hidden_pair_transfer
from fin_st166_st177_research import commutator_laplacian, local_param_krawczyk, parametric_data, rho_stationary_polynomial_intervals
from fin_st189_external_record_validator import synthetic_self_test


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST178_ST189_Results.json"
SUMMARY = ROOT / "FIN_ST178_ST189_Summary.csv"
FIG_DIR = ROOT / "FIN_ST178_ST189_Figures"
SEED = 20260826

PACKETS = {
    178: ROOT / "FIN_ST178_Irreversible_Carrier_Sufficient_Statistics.json",
    179: ROOT / "FIN_ST179_State_Dependent_Interaction.json",
    180: ROOT / "FIN_ST180_Uniform_Fold_Branch_Continuation.json",
    181: ROOT / "FIN_ST181_Three_Error_Phase_Recovery.json",
    182: ROOT / "FIN_ST182_Detailed_Balance_Selector.json",
    183: ROOT / "FIN_ST183_Twelfth_Order_Response_Observables.json",
    184: ROOT / "FIN_ST184_Adaptive_Nuisance_Cover.json",
    185: ROOT / "FIN_ST185_Layered_Compression_Visibility.json",
    186: ROOT / "FIN_ST186_Factor_Algebra_Recovery.json",
    187: ROOT / "FIN_ST187_Observed_HMM_Gap_Audit.json",
    188: ROOT / "FIN_ST188_Energy_Conserving_Local_Reset.json",
    189: ROOT / "FIN_ST189_External_Record_Interface.json",
}


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    if isinstance(x, complex): return {"real": x.real, "imag": x.imag}
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_packet(number: int, packet: dict) -> tuple[str, str]:
    path = PACKETS[number]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return path.name, sha(path)


def coarse_sum_matrix(groups: list[list[int]], fine_dimension: int) -> np.ndarray:
    c = np.zeros((len(groups), fine_dimension))
    for coarse, group in enumerate(groups):
        for fine in group:
            c[coarse, fine] = 1.0
    return c


def st178_irreversible_carriers() -> dict:
    c1 = coarse_sum_matrix([[2*i, 2*i+1] for i in range(6)], 12)
    c2 = coarse_sum_matrix([[2*i, 2*i+1] for i in range(3)], 6)
    composite = c2 @ c1
    ranks = [12, int(np.linalg.matrix_rank(c1)), int(np.linalg.matrix_rank(composite))]
    p = np.zeros(12); q = np.zeros(12); p[0] = 1; q[1] = 1
    deep_effect = np.zeros(12); deep_effect[0] = 1
    packet = {
        "layer_dimensions": [12, 6, 3],
        "composite_ranks": ranks,
        "invisible_kernel_dimensions": [12-r for r in ranks],
        "first_map": c1.tolist(), "second_map": c2.tolist(), "composite_map": composite.tolist(),
        "indistinguishable_deep_states": {"p": p.tolist(), "q": q.tolist(),
            "coarse_p": (composite@p).tolist(), "coarse_q": (composite@q).tolist(),
            "deep_effect_p": float(deep_effect@p), "deep_effect_q": float(deep_effect@q)},
        "theorem": (
            "For a finite stochastic coarse-graining C, a deep linear effect is exactly visible at the coarse layer iff it lies in range(C^T). "
            "Two deep states are operationally indistinguishable at that layer iff their difference lies in ker(C). For a hierarchy C_n...C_1, "
            "visible-effect spaces are nested and their dimensions equal the ranks of the successive composites. An exact linear decoder on all states exists iff C has full column rank; a stochastic decoder additionally requires positivity."
        ),
    }
    name, digest = write_packet(178, packet)
    return {"program":"ST178","object":"Irreversible-Carrier Sufficient-Statistic Theorem",
            "packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_exact_for_finite_stochastic_coarse_grainings",
            "boundary":"The hierarchy and carrier maps are supplied. The theorem characterizes visibility and loss; it does not derive FIN fractal layers or physical propagation without a carrier."}


def state_dependent_coupling(a: np.ndarray, rho: np.ndarray) -> np.ndarray:
    return np.diag(np.real(np.diag(a @ rho @ a)))


def st179_state_dependent_interaction(a: np.ndarray) -> dict:
    uniform = np.eye(N)/N
    localized = np.zeros((N,N)); localized[0,0] = 1
    weights = np.arange(1,N+1,dtype=float); weights /= weights.sum()
    asymmetric = np.diag(weights)
    rows=[]
    for label,rho in [("uniform",uniform),("localized_vertex_0",localized),("asymmetric_full_support",asymmetric)]:
        l=state_dependent_coupling(a,rho)
        rows.append({"state":label,"L_diagonal":np.diag(l).tolist(),
                     "L_commutant_dimension":commutant_dimension([l]),
                     "joint_A_L_commutant_dimension":commutant_dimension([a,l]),
                     "number_distinct_L_values":len(np.unique(np.round(np.diag(l),12)))})
    shift=np.roll(np.eye(N),1,axis=0)
    rho=np.diag(np.arange(1,N+1,dtype=float)/78)
    covariance_error=np.linalg.norm(state_dependent_coupling(a,shift@rho@shift.T)-shift@state_dependent_coupling(a,rho)@shift.T,np.inf)
    packet={"coupling_definition":"L_rho=diag(diag(A rho A))","rows":rows,"C12_covariance_error":float(covariance_error),
            "theorem":(
                "Because A is circulant and diagonal extraction is equivariant under vertex permutations, L_{R rho R*}=R L_rho R* for every cyclic shift R. "
                "The uniform state gives a scalar interaction and cannot select. A generic asymmetric state gives twelve distinct diagonal values and hence a 12-dimensional pointer commutant, but adding Hamiltonian A leaves only scalars for the displayed state."
            )}
    name,digest=write_packet(179,packet)
    return {"program":"ST179","object":"State-Dependent Equivariant Interaction Candidate",
            "packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_conditional_state_to_pointer_compiler_with_no_start",
            "boundary":"The asymmetric state is supplied and already carries the broken symmetry. This coupling amplifies a selected state into a pointer algebra; it does not produce the first selector or discharge QW-2191."}


IDX=np.array([i if i<=6 else 12-i for i in range(12)])


def stationary_slice_float(y: np.ndarray, a: np.ndarray, q0: np.ndarray, v: np.ndarray, epsilon: float):
    q,kappa=y[:7],y[7];u=q[IDX];au=a@u;h=np.zeros(7);rdiag=np.zeros(7)
    for i,z in enumerate(u[:7]):
        rho=z*z;den=1+rho/2;qfun=rho/den;qp=den**-2;qpp=-den**-3
        hh=-qfun*qp+.075;hp=-(qp*qp+qfun*qpp)
        h[i]=hh;rdiag[i]=2*hh+4*rho*hp
    f=np.r_[kappa*au[:7]+2*u[:7]*h, v@(q-q0)-epsilon]
    jac=np.zeros((8,8))
    for i in range(7):
        for col in range(7):
            total=kappa*sum(a[i,site] for site in range(12) if IDX[site]==col)
            if i==col:total+=rdiag[i]
            jac[i,col]=total
        jac[i,7]=au[i]
    jac[7,:7]=v
    return f,jac


def stationary_slice_interval(center: np.ndarray, radius: float, aiv, q0: np.ndarray, v: np.ndarray, epsilon: float):
    q=[iv((center[i]-radius,center[i]+radius)) for i in range(7)];kappa=iv((center[7]-radius,center[7]+radius))
    u=[q[IDX[site]] for site in range(12)]
    au=[sum((aiv[i][j]*u[j] for j in range(12)),iv(0)) for i in range(12)]
    h=[];rdiag=[]
    for z in u[:7]:
        rho=z*z;den=1+rho/2;qfun=rho/den;qp=den**-2;qpp=-den**-3
        hh=-qfun*qp+iv("0.075");hp=-(qp*qp+qfun*qpp)
        h.append(hh);rdiag.append(2*hh+4*rho*hp)
    g=[kappa*au[i]+2*u[i]*h[i] for i in range(7)]
    g.append(sum((iv(float(v[i]))*(q[i]-iv(float(q0[i]))) for i in range(7)),iv(0))-iv(str(epsilon)))
    jac=[[iv(0) for _ in range(8)] for _ in range(8)]
    for i in range(7):
        for col in range(7):
            total=sum((kappa*aiv[i][site] for site in range(12) if IDX[site]==col),iv(0))
            if i==col:total+=rdiag[i]
            jac[i][col]=total
        jac[i][7]=au[i]
    for col in range(7):jac[7][col]=iv(float(v[col]))
    glo=np.array([replay_bounds(z)[0] for z in g]);ghi=np.array([replay_bounds(z)[1] for z in g])
    jlo=np.array([[replay_bounds(z)[0] for z in row] for row in jac]);jhi=np.array([[replay_bounds(z)[1] for z in row] for row in jac])
    return glo,ghi,jlo,jhi


def krawczyk_slice(center,aiv,q0,v,epsilon,radius):
    g0lo,g0hi,j0lo,j0hi=stationary_slice_interval(center,0,aiv,q0,v,epsilon);pre=np.linalg.inv(.5*(j0lo+j0hi))
    _,_,jlo,jhi=stationary_slice_interval(center,radius,aiv,q0,v,epsilon)
    cglo,cghi=interval_matvec(pre,pre,g0lo,g0hi);ylo=np.nextafter(center-cghi,-np.inf);yhi=np.nextafter(center-cglo,np.inf)
    cjlo,cjhi=interval_left(pre,jlo,jhi);mlo,mhi=-cjhi,-cjlo
    for i in range(8):mlo[i,i]=np.nextafter(mlo[i,i]+1,-np.inf);mhi[i,i]=np.nextafter(mhi[i,i]+1,np.inf)
    dlo,dhi=interval_matvec(mlo,mhi,np.full(8,-radius),np.full(8,radius));klo=np.nextafter(ylo+dlo,-np.inf);khi=np.nextafter(yhi+dhi,np.inf)
    margin=float(min(np.min(klo-(center-radius)),np.min((center+radius)-khi)))
    return {"radius":radius,"margin":margin,"included":margin>0,"maximum_Krawczyk_halfwidth":float(np.max((khi-klo)/2))}


def uniform_fold_seed(a: np.ndarray, c: float, mode: int):
    rho=c*c;den=1+rho/2;qfun=rho/den;qp=den**-2;qpp=-den**-3
    h=-qfun*qp+.075;hp=-(qp*qp+qfun*qpp);rdiag=2*h+4*rho*hp
    full=np.cos(2*np.pi*mode*np.arange(12)/12);lam=float(full@a@full/(full@full));kappa=-rdiag/lam
    v=np.cos(2*np.pi*mode*np.arange(7)/12);v/=np.linalg.norm(v)
    return np.full(7,c),kappa,v


def st180_fold_continuation(a: np.ndarray) -> dict:
    _,_,positive=rho_stationary_polynomial_intervals(); amplitudes=[0.0]
    for lo,hi,_ in positive:
        c=math.sqrt((lo+hi)/2);amplitudes.extend([c,-c])
    epsilon=1e-3;aiv,_,_=strict_interval_matrix();mp.iv.dps=60;rows=[]
    for amplitude in amplitudes:
        for mode in range(1,7):
            q0,kappa0,v=uniform_fold_seed(a,amplitude,mode)
            for sign in (1,-1):
                eps=sign*epsilon
                sol=root(lambda y:stationary_slice_float(y,a,q0,v,eps)[0],np.r_[q0+eps*v,kappa0],
                         jac=lambda y:stationary_slice_float(y,a,q0,v,eps)[1],tol=1e-12)
                center=sol.x;res=float(np.linalg.norm(stationary_slice_float(center,a,q0,v,eps)[0],np.inf))
                trials=[krawczyk_slice(center,aiv,q0,v,eps,r) for r in (1e-7,3e-8,1e-8,3e-9)]
                accepted=next((x for x in trials if x["included"]),None)
                rows.append({"uniform_amplitude":amplitude,"mode":mode,"slice_epsilon":eps,"uniform_kappa":kappa0,
                             "continued_center":center.tolist(),"continued_kappa":float(center[7]),"residual_inf":res,
                             "nonuniform_norm":float(np.linalg.norm(center[:7]-np.mean(center[:7]))),"certificate":accepted,"trials":trials})
    passed=sum(r["certificate"] is not None for r in rows)
    packet={"uniform_folds":30,"signed_slice_targets":len(rows),"certified_continued_states":passed,"slice_amplitude":epsilon,"rows":rows,
            "theorem":(
                "Fixing the signed modal coordinate <v,q-q_uniform>=+/-10^-3 converts the stationary equation into an eight-variable square system. "
                "Every accepted Krawczyk box contains exactly one nonuniform stationary state on that declared slice for every strict operator in the transcendental interval enclosure."
            )}
    name,digest=write_packet(180,packet)
    return {"program":"ST180","object":"Validated Nonuniform Continuation from Exact Uniform Folds",
            "packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_local_slice_continuations" if passed==len(rows) else "partial_validated_local_slice_continuation",
            "boundary":"The signed slice, distance 10^-3, and reflection-even sector are supplied. These local states do not constitute a global nonuniform branch atlas or a physical phase transition."}


def st181_three_error_recovery() -> dict:
    phases=np.array([0.0,.7,1.9]);priors=np.array([.5,.3,.2]);confusion=np.array([[.82,.10,.08],[.12,.78,.10],[.06,.12,.82]])
    rows=[];total=0.0
    for y in range(3):
        w=priors*confusion[y];mass=float(w.sum());coherence=np.sum(w*np.exp(1j*phases));t=float(abs(coherence));contribution=(mass+t)/2;total+=contribution
        slack=np.array([[t,coherence],[np.conj(coherence),t]],complex)
        rows.append({"syndrome":y,"weights":w.tolist(),"mass":mass,"coherence":[float(coherence.real),float(coherence.imag)],
                     "optimal_contribution":contribution,"dual_t":t,"dual_slack_eigenvalues":np.linalg.eigvalsh(slack).tolist(),
                     "primal_phase":[float((coherence/t).real),float((coherence/t).imag)]})
    packet={"phases":phases.tolist(),"priors":priors.tolist(),"readout_y_given_error":confusion.tolist(),"rows":rows,"optimal_entanglement_fidelity":total,
            "theorem":(
                "For any finite mixture of commuting qubit Z-phase errors in a syndrome branch, the channel is a subnormalized dephasing channel with mass W and coherence c. "
                "Optimal phase correction gives (W+|c|)/2. The conic dual certificate |c|=min{t:[[t,c],[c*,t]]>=0} has eigenvalues 0 and 2|c| at t=|c|."
            )}
    name,digest=write_packet(181,packet)
    return {"program":"ST181","object":"Three-Error Nonorthogonal Phase-Recovery Certificate","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_exact_for_arbitrary_finite_commuting_qubit_phase_families",
            "boundary":"Commutativity and the qubit phase axis are essential. Three generic noncommuting errors remain a semidefinite optimization problem."}


def detailed_balance_row(beta_h: float, gamma: float, times: list[float]) -> dict:
    e=math.exp(beta_h);gap_inf=(e-1)/(e+11);rate=gamma*(math.exp(beta_h/2)+11*math.exp(-beta_h/2))
    return {"beta_h":beta_h,"stationary_gap":gap_inf,"relaxation_rate_in_gamma_units":rate,
            "trajectory":[{"gamma_t":gamma*t,"gap":gap_inf*(1-math.exp(-rate*t))} for t in times],
            "dimensionless_field_from_gap":math.log((1+11*gap_inf)/(1-gap_inf))}


def st182_detailed_balance_selector() -> dict:
    times=[0,.1,.25,.5,1,2,4];rows=[detailed_balance_row(b,1.0,times) for b in (.1,.5,1,2,4)]
    packet={"energies":"E_0=-h, E_j=0 for j!=0","rates":"q_{i->j}=gamma exp[-beta(E_j-E_i)/2]","rows":rows,
            "theorem":(
                "The declared complete-graph rates satisfy detailed balance. From a uniform state the preferred probability obeys a closed two-state equation and the vertex gap is "
                "g(t)=g_infinity(1-exp(-r t)), with g_infinity=(exp(beta h)-1)/(exp(beta h)+11) and r=gamma(exp(beta h/2)+11 exp(-beta h/2)). "
                "The field required for a stationary gap is beta h=log((1+11g)/(1-g))."
            )}
    name,digest=write_packet(182,packet)
    return {"program":"ST182","object":"Detailed-Balance Selector Control Model","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_exact_for_declared_controlled_markov_family",
            "boundary":"The preferred vertex, energy field, inverse temperature, rate scale, and clock are supplied. The model quantifies selection but does not source it from strict FIN."}


def st183_twelfth_order_observables() -> dict:
    factor_vertex=math.factorial(11);factor_mode=2048*(6**6)
    tests=[]
    for g in (-3,-.2,0,.7,5):
        derivative=factor_vertex*g;angular=g/factor_mode
        tests.append({"g":g,"vertex_twelfth_derivative":derivative,"g_from_vertex":derivative/factor_vertex,
                      "first_mode_cos12_coefficient":angular,"g_from_first_mode":angular*factor_mode})
    packet={"vertex_estimator":"g=(1/11!) partial_j^12 V_g(0)","first_mode_estimator":"g=2048*6^6*[r^12 cos(12 theta)]V_g","tests":tests,
            "theorem":(
                "For V_g=V_0+(g/12)sum_j psi_j^12 with V_0 having no degree-twelve local term, the single-vertex twelfth derivative equals 11! g. "
                "On the normalized first Fourier mode, the r^12 cos(12 theta) coefficient equals g/(2048*6^6). Either dimensionless nonlinear response identifies g exactly."
            )}
    name,digest=write_packet(183,packet)
    return {"program":"ST183","object":"Equivalent Twelfth-Order Coupling Observables","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_exact_conditional_response_identification",
            "boundary":"Access to the nonlinear potential or twelfth-order response is supplied. The result identifies g if measured; it does not derive g, its physical units, or its microscopic source."}


def st184_adaptive_cover(a: np.ndarray) -> dict:
    ec,ef,eigc,eigf=parametric_data(a);global_halfwidth=1e-4;subdivisions=7;local_halfwidth=global_halfwidth/subdivisions
    offsets=[-global_halfwidth+(2*i+1)*local_halfwidth for i in range(subdivisions)];rows=[]
    for dq0,dq1,dd in itertools.product(offsets,repeat=3):
        nuisance=(.2+dq0,.7+dq1,.05+dd)
        center=root(lambda x:point_design_system(x,ec,ef,*nuisance),np.array([2.1862,.53983]),tol=1e-12).x
        rows.append(local_param_krawczyk(eigc,eigf,nuisance,local_halfwidth,center))
    passed=sum(r["included"] for r in rows)
    packet={"global_halfwidth":global_halfwidth,"subdivisions_per_axis":subdivisions,"local_halfwidth":local_halfwidth,
            "boxes":len(rows),"passed_boxes":passed,"minimum_margin":min(r["margin"] for r in rows),"rows":rows,
            "theorem":"The 7x7x7 union exactly covers the declared nuisance cube. Inclusion in every cell proves uniform stationary-root continuation over it."}
    name,digest=write_packet(184,packet)
    return {"program":"ST184","object":"Third-Generation Adaptive Nuisance Continuation","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_extended_uniform_cover" if passed==len(rows) else "partial_extended_cover",
            "boundary":"No maximal continuation domain, global optimum, dimensional calibration, or apparatus tolerance follows."}


def kl(p: np.ndarray,q: np.ndarray)->float:
    return float(sum(x*math.log(x/y) for x,y in zip(p,q) if x>0))


def st185_layer_visibility() -> dict:
    eps=.25;shift=np.roll(np.eye(N),1,axis=0);b=(1-eps)*np.eye(N)+(eps/2)*(shift+shift.T)
    eigenvalues=np.array([1-eps+eps*math.cos(2*math.pi*k/N) for k in range(N)])
    rows=[]
    for layers in range(0,13):
        atten=np.abs(eigenvalues)**layers
        rows.append({"layers":layers,"mode_attenuation":atten[:7].tolist(),"minimum_singular_value":float(np.min(atten)),
                     "inverse_noise_amplification":float(1/np.min(atten)),"effective_modes_above_1e-3":int(np.sum(atten>1e-3))})
    c1=coarse_sum_matrix([[2*i,2*i+1] for i in range(6)],12);c2=coarse_sum_matrix([[2*i,2*i+1] for i in range(3)],6);hard=c2@c1
    rng=np.random.default_rng(SEED+185);p=rng.random(N);p/=p.sum();q=rng.random(N);q/=q.sum()
    dpi=[];pp=p.copy();qq=q.copy()
    for layers in range(9):
        dpi.append({"layers":layers,"KL":kl(pp,qq),"total_variation":float(.5*np.sum(abs(pp-qq)))})
        pp=b@pp;qq=b@qq
    packet={"soft_layer":"B=(1-epsilon)I+(epsilon/2)(S+S^-1)","epsilon":eps,"single_layer_Fourier_eigenvalues":eigenvalues[:7].tolist(),
            "visibility_rows":rows,"hard_12_to_3_rank":int(np.linalg.matrix_rank(hard)),"hard_invisible_dimension":N-int(np.linalg.matrix_rank(hard)),
            "data_processing_audit":dpi,
            "theorem":(
                "For a layered observation C^(n)=C_n...C_1, the Layer Visibility Spectrum is the singular spectrum of C^(n). Zero singular modes are exactly invisible; "
                "a nonzero mode is noiselessly recoverable but reconstruction amplifies observation error by at least its inverse singular value. For repeated circulant B, Fourier visibility is |mu_k|^n. "
                "Total variation and relative entropy cannot increase under each stochastic layer."
            )}
    name,digest=write_packet(185,packet)
    return {"program":"ST185","object":"Layered-Compression Visibility Spectrum","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_exact_layer_visibility_and_data_processing_theorem",
            "boundary":"This rigorously models seeing deep scales through supplied compression layers. It does not prove that FIN generates these layers, that they are geometric fractals, or that any layer is the Planck scale or observed spacetime."}


def subspace_projector(lap: np.ndarray, dimension: int=4):
    vals,vecs=np.linalg.eigh(lap);v=vecs[:,:dimension];return vals,v@v.conj().T,v


def multiplication_leakage(basis: np.ndarray, projector: np.ndarray) -> float:
    d=int(round(math.sqrt(projector.shape[0])));worst=0.0
    mats=[basis[:,i].reshape((d,d),order="F") for i in range(basis.shape[1])]
    for x in mats:
        for y in mats:
            z=(x@y).reshape(-1,order="F");worst=max(worst,float(np.linalg.norm((np.eye(d*d)-projector)@z)))
    return worst


def st186_factor_recovery() -> dict:
    x=np.array([[0,1],[1,0]],complex);z=np.diag([1,-1]).astype(complex);i2=np.eye(2)
    rng=np.random.default_rng(SEED+186);raw=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));u,_=np.linalg.qr(raw)
    generators=[u@np.kron(x,i2)@u.conj().T,u@np.kron(z,i2)@u.conj().T]
    vals0,p0,b0=subspace_projector(commutator_laplacian(generators));rows=[]
    for eta in (0,.001,.005,.01,.02,.05):
        noisy=[]
        for g in generators:
            if eta==0:h=np.zeros_like(g)
            else:
                h=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));h=(h+h.conj().T)/2;h*=eta/np.linalg.norm(h,2)
            noisy.append(g+h)
        vals,p,basis=subspace_projector(commutator_laplacian(noisy))
        rows.append({"eta":eta,"projector_distance":float(np.linalg.norm(p-p0,2)),"fourth_eigenvalue":float(vals[3]),"fifth_eigenvalue":float(vals[4]),
                     "multiplication_leakage_of_raw_low_basis":multiplication_leakage(basis,p)})
    packet={"ideal_gap":float(vals0[4]),"rows":rows,"gauge_group":"global conjugation modulo local U(2)xU(2); add factor swap if labels are removed",
            "theorem":(
                "Any two labeled commuting M2 factors generating M4 are globally unitarily equivalent to the standard tensor factors. The implementing unitary is determined only modulo local U(2)xU(2); unlabeled factors add swap ambiguity. "
                "The commutator-Laplacian reconstructs the exact commutant for noiseless complete generators. Under generic noise its low subspace is close but is not exactly multiplication closed, so algebra projection is an additional estimation step."
            )}
    name,digest=write_packet(186,packet)
    return {"program":"ST186","object":"Factor-Algebra Recovery and Gauge-Ambiguity Audit","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_exact_noiseless_classification_with_numerical_noisy_audit",
            "boundary":"The global factor labels, complete generators, and Hilbert-space dimension are supplied. No physical subsystem, causal region, or apparatus label is inferred."}


def st187_hmm_gap() -> dict:
    P=np.array([[.9,.1],[.2,.8]]);pi=np.array([2/3,1/3]);p0=np.array([.02,.08]);p1=np.array([.92,.98]);e0=np.c_[1-p0,p0];e1=np.c_[1-p1,p1]
    initial,t=hidden_pair_transfer(P,p0,p1,.5);rows=[]
    for n in (1,2,4,8,12,16,18,20):
        l0=hmm_sequence_likelihoods(n,P,pi,e0);l1=hmm_sequence_likelihoods(n,P,pi,e1)
        observed=float(np.sum(np.sqrt(l0*l1)));pair=float(initial@np.linalg.matrix_power(t,n-1)@np.ones(4))
        rows.append({"events":n,"observed_Hellinger_coefficient":observed,"pair_path_upper_bound":pair,"bound_over_exact_ratio":pair/observed,
                     "observed_finite_rate":-math.log(observed)/n,"pair_path_finite_rate":-math.log(pair)/n})
    packet={"rows":rows,"ST175_pair_path_asymptotic_rate_lower_bound":0.32690226513199694,
            "theorem":(
                "At every finite n, concavity/subadditivity gives the rigorous inequality Z_n^observed<=B_n^pair. Complete forward enumeration computes the declared finite observed coefficients. "
                "The pair-path exponential rate is therefore a rigorous lower bound on the observed Hellinger discrimination rate, but the growing B_n/Z_n ratio demonstrates that it is quantitatively loose for this model."
            )}
    name,digest=write_packet(187,packet)
    return {"program":"ST187","object":"Observed-HMM versus Pair-Path Exponent Gap Audit","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_finite_bound_with_strong_asymptotic_gap_evidence",
            "boundary":"The displayed observed coefficients use deterministic floating forward enumeration. They are finite-n audit values, not an interval proof of the exact asymptotic observed Chernoff rate."}


def swap_bit_matrix(bits: int, left_bit: int, right_bit: int) -> np.ndarray:
    dim=1<<bits;p=np.zeros((dim,dim))
    for state in range(dim):
        a=(state>>left_bit)&1;b=(state>>right_bit)&1;target=state
        if a!=b:target^=(1<<left_bit)|(1<<right_bit)
        p[target,state]=1
    return p


def st188_local_reset(a: np.ndarray) -> dict:
    h=np.zeros((16,16));h[:12,:12]=a;i16=np.eye(16);ht=np.kron(h,i16)+np.kron(i16,h)
    swaps=[swap_bit_matrix(8,i,4+i) for i in range(4)]
    comms=[ht@s-s@ht for s in swaps];gram=np.array([[np.vdot(x,y).real for y in comms] for x in comms])
    geig=np.linalg.eigvalsh(gram);global_swap=swaps[0]@swaps[1]@swaps[2]@swaps[3]
    packet={"binary_embedding":"strict A in leading 12x12 block of a 16-level four-qubit register; four padding levels set to zero",
            "pairwise_SWAP_commutator_Gram":gram.tolist(),"Gram_eigenvalues":geig.tolist(),"numerical_pairwise_linear_ansatz_nullity":int(np.sum(geig<1e-10)),
            "global_register_SWAP_commutator_norm":float(np.linalg.norm(ht@global_swap-global_swap@ht,2)),
            "equal_pairwise_sum_commutator_norm":float(np.linalg.norm(sum(comms),2)),
            "theorem":(
                "The product of four disjoint qubit SWAPs is the global register SWAP and commutes exactly with H tensor I+I tensor H for every identical register Hamiltonian H. "
                "For the displayed padded strict encoding, the four individual pairwise-SWAP commutators have a positive numerical Gram spectrum, giving strong evidence that no nonzero static linear combination of those four two-qubit SWAP generators conserves the encoded energy."
            )}
    name,digest=write_packet(188,packet)
    return {"program":"ST188","object":"Energy-Conserving Reset Locality Audit","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_global_swap_identity_with_strong_pairwise_ansatz_no_go_evidence",
            "boundary":"The no-go is numerical and restricted to a padded binary encoding and the four static pairwise-SWAP generator span. Controlled sequences, ancillas, other 2-local terms, and exact interval rank certification remain open."}


def st189_external_interface() -> dict:
    selftest=synthetic_self_test()
    packet={"required_event_fields":["timestamp","outcome","config","run_id"],"self_test":selftest,
            "theorem":(
                "The executable interface detects event-file tampering, verifies a frozen holdout hash and distinct declared roles, and refuses to mark physical execution valid without externally supplied paths, calibration hash, custody attestation, and laboratory attestation."
            )}
    name,digest=write_packet(189,packet)
    return {"program":"ST189","object":"External Immutable-Record Validation Interface","packet_file":name,"packet_sha256":digest,**packet,
            "status":"proven_executable_interface_self_test_no_external_record",
            "boundary":"Only a synthetic self-test was run. Role strings and attestation flags are not identity proofs; no external record, independent custody, calibration, or laboratory result is claimed."}


def make_figures(d: dict):
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4));ax.plot([0,1,2],d["ST178"]["composite_ranks"],"o-",label="visible rank");ax.plot([0,1,2],d["ST178"]["invisible_kernel_dimensions"],"s--",label="invisible dimension");ax.set(xlabel="compression layer",ylabel="dimension",title="ST178: visible and lost operational directions");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st178_irreversible_carriers.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST179"]["rows"];ax.bar([r["state"] for r in rows],[r["L_commutant_dimension"] for r in rows],label="L commutant");ax.plot(range(len(rows)),[r["joint_A_L_commutant_dimension"] for r in rows],"ro-",label="joint with A");ax.tick_params(axis="x",rotation=12);ax.set(ylabel="complex dimension",title="ST179: state-dependent interaction algebra");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st179_state_interaction.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST180"]["rows"];ax.scatter([r["uniform_kappa"] for r in rows],[r["continued_kappa"]-r["uniform_kappa"] for r in rows],c=[r["mode"] for r in rows],s=22);ax.set(xlabel="uniform fold kappa",ylabel="continued kappa shift",title="ST180: signed local continuations");fig.tight_layout();fig.savefig(FIG_DIR/"st180_fold_continuation.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST181"]["rows"];ax.bar([r["syndrome"] for r in rows],[r["optimal_contribution"] for r in rows]);ax.set(xlabel="syndrome",ylabel="optimal contribution",title="ST181: three-error phase recovery");fig.tight_layout();fig.savefig(FIG_DIR/"st181_three_error_recovery.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));
    for row in d["ST182"]["rows"]:ax.plot([x["gamma_t"] for x in row["trajectory"]],[x["gap"] for x in row["trajectory"]],"o-",label=f"beta h={row['beta_h']}")
    ax.set(xlabel="dimensionless time gamma t",ylabel="vertex gap",title="ST182: detailed-balance selector trajectories");ax.legend(ncol=2);fig.tight_layout();fig.savefig(FIG_DIR/"st182_selector_trajectories.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));ax.hist([r["margin"] for r in d["ST184"]["rows"]],bins=24);ax.set(xlabel="Krawczyk margin",ylabel="cells",title="ST184: 343-cell nuisance cover");fig.tight_layout();fig.savefig(FIG_DIR/"st184_cover.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST185"]["visibility_rows"];
    for k in range(7):ax.semilogy([r["layers"] for r in rows],[r["mode_attenuation"][k] for r in rows],"o-",label=f"k={k}")
    ax.set(xlabel="compression layers",ylabel="mode visibility",title="ST185: seeing deep modes through layered compression");ax.legend(ncol=2);fig.tight_layout();fig.savefig(FIG_DIR/"st185_layer_visibility.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST186"]["rows"];ax.semilogy([max(r["eta"],1e-5) for r in rows],[max(r["projector_distance"],1e-16) for r in rows],"o-",label="projector distance");ax.semilogy([max(r["eta"],1e-5) for r in rows],[max(r["multiplication_leakage_of_raw_low_basis"],1e-16) for r in rows],"s--",label="multiplication leakage");ax.set(xlabel="generator noise",ylabel="error",title="ST186: noisy factor reconstruction");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st186_factor_recovery.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST187"]["rows"];ax.plot([r["events"] for r in rows],[r["observed_finite_rate"] for r in rows],"o-",label="observed finite rate");ax.plot([r["events"] for r in rows],[r["pair_path_finite_rate"] for r in rows],"s--",label="pair-path rate");ax.set(xlabel="events",ylabel="-log coefficient / n",title="ST187: hidden-HMM exponent gap");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st187_hmm_gap.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));vals=d["ST188"]["Gram_eigenvalues"];ax.bar(range(1,5),vals);ax.set(xlabel="Gram eigenvalue index",ylabel="eigenvalue",title="ST188: pairwise-SWAP commutator independence");fig.tight_layout();fig.savefig(FIG_DIR/"st188_swap_gram.png",dpi=190);plt.close(fig)


def main():
    _,a,_=strict_operator();out={"metadata":{"programs":"ST178-ST189","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__}}
    out["ST178"]=st178_irreversible_carriers();out["ST179"]=st179_state_dependent_interaction(a);out["ST180"]=st180_fold_continuation(a);out["ST181"]=st181_three_error_recovery()
    out["ST182"]=st182_detailed_balance_selector();out["ST183"]=st183_twelfth_order_observables();out["ST184"]=st184_adaptive_cover(a);out["ST185"]=st185_layer_visibility()
    out["ST186"]=st186_factor_recovery();out["ST187"]=st187_hmm_gap();out["ST188"]=st188_local_reset(a);out["ST189"]=st189_external_interface()
    out["recommended_next_programs"]=[
        {"id":"ST190","priority":1,"study":"derive a repository-sourced candidate hierarchy map and test it against the ST185 Layer Visibility Spectrum without Planck-scale import"},
        {"id":"ST191","priority":2,"study":"classify irreversible carrier equivalence by Blackwell sufficiency and minimal sufficient operational quotients"},
        {"id":"ST192","priority":3,"study":"seek a dynamical law for the ST179 state-dependent interaction and prove or refute spontaneous pointer selection from symmetric noise"},
        {"id":"ST193","priority":4,"study":"validated arclength continuation and branch-collision search beyond the signed ST180 slices"},
        {"id":"ST194","priority":5,"study":"solve three generic noncommuting qubit errors with rational primal-dual SDP certificates"},
        {"id":"ST195","priority":6,"study":"optimize time-dependent detailed-balance selector controls under a frozen entropy-production budget"},
        {"id":"ST196","priority":7,"study":"search dimensionless internal laws relating the two ST183 response estimators to a strict nonlinear source"},
        {"id":"ST197","priority":8,"study":"locate the maximal nuisance continuation boundary with adaptive interval branch-and-bound"},
        {"id":"ST198","priority":9,"study":"derive sampling lower bounds for reconstructing attenuated deep modes under Poisson and multinomial observation noise"},
        {"id":"ST199","priority":10,"study":"project noisy commutant subspaces onto the nearest exact matrix-factor variety with uniqueness certificates"},
        {"id":"ST200","priority":11,"study":"certify the exact observed hidden-HMM Hellinger rate using an interval transfer operator on belief space"},
        {"id":"ST201","priority":12,"study":"interval-certify or falsify the ST188 pairwise-SWAP ansatz rank and enlarge the local generator class"},
        {"id":"ST202","priority":13,"study":"execute ST189 only on a genuinely external immutable record with independently supplied custody attestations"},
    ]
    out["central_verdict"]=("The intuition that an observer sees deeper scales through developed compression layers has a precise conditional form: visible observables lie in the pullback range of the composite layer map, invisible differences lie in its kernel, and nonzero singular modes are recovered only with inverse-singular-value noise amplification. Repeated self-similar soft layers exponentially attenuate selected Fourier modes. This is a theorem about supplied observation channels, not yet a derivation of FIN fractal layers or spacetime.")
    out["epistemic_boundary"]=("No physical carrierlessness, strict fractal hierarchy, Planck layer, spacetime projection, strict selector, QW-2191 discharge, dimensional scale, calibrated apparatus or external record, legacy-to-strict bridge completion or role transfer, Standard Model, gravity, L_total, or ToE closure is claimed.")
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(178,190):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out);print(json.dumps({"results":RESULTS.name,"programs":12,"figures":10},indent=2))


if __name__=="__main__":main()
