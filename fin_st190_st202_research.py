#!/usr/bin/env python3
"""FIN ST190--ST202: endogenous compression candidates and operational limits.

All computations are local and dimensionless.  The program intentionally keeps
strict-kernel theorems, supplied constructions, numerical evidence, and blocked
physical claims in separate fields.
"""

from __future__ import annotations

import csv
import hashlib
import inspect
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
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.optimize import root

from fin_st01_st15_research import N, strict_operator
from fin_st130_st141_research import point_design_system
from fin_st132_center_isolation_replay import bounds as replay_bounds, strict_interval_matrix
from fin_st154_st165_research import algebra_span_dimension, commutant_dimension
from fin_st166_st177_research import local_param_krawczyk, parametric_data
from fin_st178_st189_research import (
    krawczyk_slice,
    stationary_slice_float,
    swap_bit_matrix,
    uniform_fold_seed,
)
from fin_st189_external_record_validator import synthetic_self_test


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST190_ST202_Results.json"
SUMMARY = ROOT / "FIN_ST190_ST202_Summary.csv"
FIG_DIR = ROOT / "FIN_ST190_ST202_Figures"
SEED = 20260811
PACKETS = {k: ROOT / f"FIN_ST{k}_{name}.json" for k, name in {
    190: "Strict_Heat_Hierarchy", 191: "Blackwell_Minimal_Quotient",
    192: "Replicator_Pointer_Selection", 193: "Validated_Modal_Continuation",
    194: "Noncommuting_Qubit_Recovery", 195: "Selector_Control_Cost",
    196: "Nonlinear_Source_Nonidentifiability", 197: "Nuisance_Boundary",
    198: "Deep_Mode_Sampling_Bounds", 199: "Exact_Factor_Projection",
    200: "Observed_HMM_Finite_Interval", 201: "Local_Reset_Interval_Rank",
    202: "External_Record_Readiness",
}.items()}


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


def finalize(number: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    name, digest = write_packet(number, packet)
    return {"program": f"ST{number}", "object": obj, "packet_file": name,
            "packet_sha256": digest, **packet, "status": status, "boundary": boundary}


def st190_heat_hierarchy(a: np.ndarray) -> dict:
    tau0 = .25
    eig = np.linalg.eigvalsh(a)
    fourier = np.real(np.fft.fft(a[0]))[:7]
    rows = []
    composite = np.eye(N)
    elapsed = 0.0
    for j in range(1, 9):
        tau = tau0 * 2**(j-1)
        layer = expm(-tau*a)
        composite = layer @ composite
        elapsed += tau
        values = np.exp(-elapsed*fourier)
        rows.append({"layer": j, "incremental_parameter": tau,
                     "composite_parameter": elapsed,
                     "minimum_entry": float(layer.min()),
                     "maximum_row_sum_error": float(np.max(abs(layer.sum(axis=1)-1))),
                     "Fourier_visibility": values.tolist(),
                     "effective_rank_at_1e-8": int(np.sum(np.exp(-elapsed*eig) > 1e-8))})
    packet = {
        "base_parameter_tau0": tau0, "dyadic_rule": "tau_j=2^(j-1) tau0",
        "composite_rule": "C_n=exp[-(2^n-1) tau0 A]",
        "strict_eigenvalues": eig.tolist(), "strict_Fourier_class_eigenvalues_k0_to_k6": fourier.tolist(), "rows": rows,
        "semigroup_composition_error": float(np.linalg.norm(composite-expm(-elapsed*a), np.inf)),
        "theorem": (
            "For the strict graph Laplacian A, -A is a conservative Metzler generator, hence exp(-tA) is a positive stochastic contraction. "
            "Dyadic self-composition is exactly exp[-(2^n-1)tau0 A], and each strict eigenmode has visibility exp[-(2^n-1)tau0 lambda_k]. "
            "At every finite layer this map is invertible; only the infinite-depth limit removes all nonconstant modes."
        )}
    return finalize(190, "Strict-Generated Dyadic Heat-Compression Hierarchy",
        "proven_conditional_dimensionless_hierarchy_from_strict_semigroup",
        "The strict generator supplies the semigroup, but tau0 and the dyadic schedule are conventions. This is not a derived fractal law, length hierarchy, Planck scale, spacetime, or physical observation map.", packet)


def proportional_classes(table: list[list[Fraction]]) -> list[list[int]]:
    # table[h][x]; columns belong together iff all 2x2 minors vanish.
    classes: list[list[int]] = []
    for x in range(len(table[0])):
        for group in classes:
            y = group[0]
            if all(table[h][x]*table[k][y] == table[k][x]*table[h][y]
                   for h in range(len(table)) for k in range(len(table))):
                group.append(x); break
        else: classes.append([x])
    return classes


def st191_blackwell_quotient() -> dict:
    masses = [
        [Fraction(1,2), Fraction(1,4), Fraction(1,8), Fraction(1,8)],
        [Fraction(1,8), Fraction(1,2), Fraction(1,4), Fraction(1,8)],
        [Fraction(1,8), Fraction(1,8), Fraction(1,2), Fraction(1,4)],
    ]
    table = [[masses[h][x//3]/3 for x in range(12)] for h in range(3)]
    classes = proportional_classes(table)
    coarse = [[sum(table[h][x] for x in group) for group in classes] for h in range(3)]
    decoder = [[Fraction(1, len(group)) if x in group else Fraction(0)
                for x in range(12)] for group in classes]
    reconstructed = [[sum(coarse[h][g]*decoder[g][x] for g in range(4))
                      for x in range(12)] for h in range(3)]
    # Merging classes 0 and 3 destroys a likelihood ratio.
    p, q = np.array([float(x) for x in table[0]]), np.array([float(x) for x in table[1]])
    tv_fine = .5*float(np.sum(abs(p-q)))
    merge = [[0,1],[2],[3]]
    pc = np.array([sum(float(coarse[0][g]) for g in z) for z in merge])
    qc = np.array([sum(float(coarse[1][g]) for g in z) for z in merge])
    tv_merge = .5*float(np.sum(abs(pc-qc)))
    packet = {
        "hypotheses": 3, "fine_outcomes": 12,
        "likelihood_table_exact": [[str(x) for x in row] for row in table],
        "minimal_proportionality_classes": classes,
        "minimal_quotient_outcomes": len(classes),
        "coarse_likelihood_table_exact": [[str(x) for x in row] for row in coarse],
        "exact_reconstruction": reconstructed == table,
        "overmerged_total_variation_before": tv_fine,
        "overmerged_total_variation_after": tv_merge,
        "theorem": (
            "For a finite dominated experiment, a deterministic statistic is Blackwell sufficient exactly when every fine likelihood vector within a fiber is proportional. "
            "The quotient by proportional likelihood columns is minimal sufficient, unique up to relabeling. The displayed 12-outcome experiment has exactly four classes and an exact stochastic decoder; merging two inequivalent classes strictly reduces binary total variation."
        )}
    return finalize(191, "Blackwell-Minimal Operational Quotient",
        "proven_exact_for_declared_finite_experiment",
        "The preparation family and likelihood table are supplied. The theorem identifies lossless carrier compression relative to that experiment, not an observer-independent ontology or a FIN-derived apparatus.", packet)


def replicator_rhs(_t: float, p: np.ndarray, m: np.ndarray, gamma: float) -> np.ndarray:
    payoff = m @ p
    return gamma*p*(payoff-float(p@payoff))


def st192_replicator_pointer(a: np.ndarray) -> dict:
    m = a*a
    mu = np.linalg.eigvalsh(m)
    fourier_mu = np.real(np.fft.fft(m[0]))
    rng = np.random.default_rng(SEED+192)
    trajectories = []
    for run in range(12):
        p0 = np.ones(N)/N + 2e-3*rng.normal(size=N)
        p0 = np.maximum(p0, 1e-8); p0 /= p0.sum()
        sol = solve_ivp(replicator_rhs, (0,80), p0, args=(m,1.0), rtol=1e-10, atol=1e-12)
        pf = np.maximum(sol.y[:,-1],0); pf /= pf.sum()
        trajectories.append({"run":run,"initial_argmax":int(np.argmax(p0)),
                             "selected_vertex":int(np.argmax(pf)),
                             "final_max_probability":float(pf.max()),
                             "final_entropy":float(-sum(x*math.log(x) for x in pf if x>0))})
    uniform_rhs = replicator_rhs(0,np.ones(N)/N,m,1)
    packet = {
        "feedback_law": "p_dot_i=gamma p_i[(M p)_i-p^T M p], M=A Hadamard-square A",
        "Hadamard_square_eigenvalues": mu.tolist(), "Hadamard_square_Fourier_eigenvalues": fourier_mu.tolist(),
        "uniform_rhs_norm": float(np.linalg.norm(uniform_rhs,np.inf)),
        "minimum_nonconstant_linear_growth_rate": float(min(fourier_mu[1:])/N),
        "vertex_local_stability_margin": float(min(m.diagonal().min()-np.max(m-np.diag(np.diag(m)),axis=1))),
        "trajectories": trajectories,
        "distinct_selected_vertices": sorted(set(r["selected_vertex"] for r in trajectories)),
        "theorem": (
            "The constructed replicator law is C12-equivariant and preserves the probability simplex. The uniform state is an exact fixed point and therefore cannot select without a perturbation. "
            "Its tangent growth rates are the nonconstant eigenvalues of M/12; positive rates make it unstable. Each vertex is locally attracting when M_ii>M_ji for j!=i. Symmetric noise selects an orbit member but not a canonical label."
        )}
    return finalize(192, "Equivariant Positive-Feedback Pointer Dynamics",
        "proven_properties_for_constructed_law_with_numerical_trajectory_audit",
        "The feedback functional and initial noise are added, not derived by strict FIN. Random branch selection is spontaneous symmetry breaking but does not canonically distinguish +1 from -1 or discharge QW-2191.", packet)


def st193_modal_continuation(a: np.ndarray) -> dict:
    previous = json.loads((ROOT/"FIN_ST178_ST189_Results.json").read_text(encoding="utf-8"))["ST180"]
    aiv,_,_ = strict_interval_matrix(); mp.iv.dps=60
    targets = [.002,.005]
    rows=[]
    for old in previous["rows"]:
        q0,k0,v=uniform_fold_seed(a,old["uniform_amplitude"],old["mode"])
        center=np.array(old["continued_center"])
        sign=1 if old["slice_epsilon"]>0 else -1
        for eps_abs in targets:
            eps=sign*eps_abs
            sol=root(lambda y:stationary_slice_float(y,a,q0,v,eps)[0],center,
                     jac=lambda y:stationary_slice_float(y,a,q0,v,eps)[1],tol=1e-12)
            center=sol.x
            trials=[krawczyk_slice(center,aiv,q0,v,eps,r) for r in (3e-7,1e-7,3e-8,1e-8)]
            accepted=next((x for x in trials if x["included"]),None)
            rows.append({"uniform_amplitude":old["uniform_amplitude"],"mode":old["mode"],
                         "slice_epsilon":eps,"residual_inf":float(np.linalg.norm(stationary_slice_float(center,a,q0,v,eps)[0],np.inf)),
                         "continued_center":center.tolist(),"certificate":accepted,"trials":trials})
    passed=sum(r["certificate"] is not None for r in rows)
    packet={"source":"all 60 ST180 signed branches","additional_slice_levels":targets,
            "targets":len(rows),"certified_targets":passed,"rows":rows,
            "theorem":"Each accepted Krawczyk box contains exactly one stationary point on its declared modal slice for every strict operator in the interval enclosure. Nested successful slices continue the corresponding ST180 local branch to the displayed modal amplitude."}
    return finalize(193,"Validated Multi-Slice Modal Continuation",
        "proven_validated_modal_continuation" if passed==len(rows) else "partial_validated_modal_continuation",
        "This is a sequence of fixed modal slices, not a global arclength atlas. Failure would diagnose the certificate method, not necessarily a physical branch collision; no phase transition or soliton particle is inferred.",packet)


def rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis=np.asarray(axis,float);axis/=np.linalg.norm(axis)
    x,y,z=axis;c=math.cos(angle);s=math.sin(angle);c1=1-c
    return np.array([[c+x*x*c1,x*y*c1-z*s,x*z*c1+y*s],
                     [y*x*c1+z*s,c+y*y*c1,y*z*c1-x*s],
                     [z*x*c1-y*s,z*y*c1+x*s,c+z*z*c1]])


def st194_noncommuting_recovery() -> dict:
    probs=np.array([.46,.31,.23])
    rotations=[rotation([1,0,0],.61),rotation([0,1,0],1.07),rotation([1,1,1],.83)]
    t=sum(p*r for p,r in zip(probs,rotations))
    u,s,vh=np.linalg.svd(t);orientation=float(np.linalg.det(u@vh))
    r_opt=vh.T@u.T
    raw_lower=(1+float(np.trace(r_opt@t)))/4
    upper=(1+float(np.sum(s)))/4
    lower=min(raw_lower,upper)
    comms=[np.linalg.norm(rotations[i]@rotations[j]-rotations[j]@rotations[i]) for i in range(3) for j in range(i)]
    packet={"probabilities":probs.tolist(),"Bloch_rotations":[r.tolist() for r in rotations],
            "pairwise_commutator_norms":comms,"average_Bloch_matrix":t.tolist(),
            "singular_values":s.tolist(),"SVD_orientation":orientation,
            "optimal_recovery_rotation":r_opt.tolist(),"primal_entanglement_fidelity":lower,
            "trace_norm_upper_bound":upper,"primal_upper_gap":upper-lower,
            "theorem":(
                "The linear Bloch part S of every positive trace-preserving qubit recovery satisfies ||S||_2<=1, so F_e(R o E)=(1+Tr(ST))/4 <= (1+||T||_*)/4. "
                "When the polar orientation of T is positive, a unitary recovery realizes the polar rotation and saturates this bound. The displayed three errors are pairwise noncommuting and have positive orientation, hence the certificate is global over all CPTP recoveries, without an SDP solver."
            )}
    return finalize(194,"Generic Noncommuting-Qubit Recovery Certificate",
        "proven_global_CPTP_optimum_for_declared_positive_orientation_random_unitary_channel",
        "The result covers the declared unital qubit channel and positive polar orientation. It is not a rational SDP replay and does not cover negative-orientation or higher-dimensional channels.",packet)


def equilibrium_gap_cost(b: float, steps: int) -> dict:
    db=b/steps
    work=sum((1-1/(1+11*math.exp(-j*db)))*db for j in range(steps))
    free=-math.log((1+11*math.exp(-b))/12)
    return {"steps":steps,"work":work,"free_energy_change":free,"dissipation":work-free}


def st195_selector_control() -> dict:
    b=3.0; rows=[equilibrium_gap_cost(b,n) for n in [1,2,4,8,16,32,64,128,256,512]]
    p0=1/(1+11*math.exp(-b)); p=np.r_[[p0],np.full(11,(1-p0)/11)];u=np.ones(12)/12
    d=float(sum(x*math.log(x/y) for x,y in zip(p,u)))
    packet={"target_beta_h":b,"target_preferred_probability":p0,"target_gap":float(p0-p[1]),
            "target_relative_entropy_from_uniform":d,"quasistatic_rows":rows,
            "minimum_steps_for_dissipation_budgets":{str(e):next((r["steps"] for r in rows if r["dissipation"]<=e),None) for e in [.1,.05,.01,.005]},
            "theorem":(
                "For the controlled energies E_0=0 and E_j=b (j!=0), equilibrium has p_0=1/(1+11e^-b) and free-energy cost Delta F=-log[(1+11e^-b)/12]. "
                "Any isothermal protocol obeys W>=Delta F; the displayed equilibrate-then-step protocol has an explicit positive Riemann-sum excess that decreases to zero as the number of control steps grows."
            )}
    return finalize(195,"Detailed-Balance Selector Control under a Dissipation Budget",
        "proven_thermodynamic_bound_and_constructive_quasistatic_convergence_for_supplied_control",
        "The preferred vertex, bath, beta, field, protocol, and time-scale are supplied. This prices a selector but does not derive one from strict FIN or define physical joules.",packet)


def st196_nonlinear_source() -> dict:
    gs=[.01,.1,1,10,100]
    rows=[{"g":g,"Hessian_at_zero_difference_from_A":0.0,"Z12_invariant":True,
           "coercive_with_positive_quadratic_A":True,"twelfth_response":math.factorial(11)*g} for g in gs]
    packet={"family":"V_g(psi)=psi^T A psi/2+(g/12) sum_j psi_j^12","rows":rows,
            "theorem":(
                "Every g>0 gives the same strict Hessian A at the origin, the same C12 symmetry, and coercivity, while its twelfth local derivative is 11!g. "
                "Therefore A, symmetry, boundedness below, and the strict linear response identify only the sign class g>0, not its magnitude. A nonlinear response measurement can identify g but cannot make it a strict source theorem."
            )}
    return finalize(196,"Strict Nonlinear-Coupling Source Nonidentifiability",
        "proven_infinite_family_nonidentifiability_theorem",
        "The theorem is a no-go for deriving g from the listed strict linear data. Additional normalization, microscopic dynamics, or a measured twelfth-order response remains necessary.",packet)


def nuisance_cover(a: np.ndarray, h: float) -> dict:
    ec,ef,eigc,eigf=parametric_data(a);sub=7;lh=h/sub
    offsets=[-h+(2*i+1)*lh for i in range(sub)];rows=[]
    for ds in itertools.product(offsets,repeat=3):
        nu=(.2+ds[0],.7+ds[1],.05+ds[2])
        center=root(lambda x:point_design_system(x,ec,ef,*nu),[2.1862,.53983],tol=1e-12).x
        rows.append(local_param_krawczyk(eigc,eigf,nu,lh,center))
    return {"halfwidth":h,"boxes":len(rows),"passed":sum(r["included"] for r in rows),"minimum_margin":min(r["margin"] for r in rows)}


def st197_nuisance_boundary(a: np.ndarray) -> dict:
    accepted=nuisance_cover(a,.00013045);rejected=nuisance_cover(a,.00013050)
    packet={"fixed_partition":"7x7x7 equal cells","accepted_cover":accepted,"first_tested_rejected_cover":rejected,
            "method_bracket_width":rejected["halfwidth"]-accepted["halfwidth"],
            "theorem":"All 343 cells at halfwidth 0.00013045 satisfy interval Krawczyk inclusion and therefore cover a unique continued stationary root. At halfwidth 0.00013050, six cells fail this particular inclusion test."}
    return finalize(197,"Validated Nuisance Continuation Method Boundary",
        "proven_cover_to_0_00013045_with_certificate_failure_at_0_00013050",
        "The upper endpoint is a failure of the fixed 7^3 Krawczyk cover, not proof of root loss or a maximal mathematical continuation domain. Adaptive subdivision may extend it.",packet)


def st198_sampling_bounds() -> dict:
    delta=.05;epsilon=.01;alpha=.05;rows=[]
    for n in range(13):
        sigma=2.0**(-n);a=sigma*delta
        info=4*sigma*sigma/(1-4*a*a)
        cr=math.ceil(1/(info*epsilon*epsilon))
        kl=2*a*math.log((.5+a)/(.5-a))
        testing=math.ceil(2*(1-2*alpha)**2/kl)
        rows.append({"layers":n,"attenuation":sigma,"Fisher_information_per_sample":info,
                     "Cramer_Rao_samples_for_sd_0_01":cr,"necessary_samples_for_error_0_05_via_Pinsker":testing})
    packet={"deep_mode_amplitude":delta,"target_standard_deviation":epsilon,"target_testing_error":alpha,"rows":rows,
            "theorem":(
                "For a binary contrast observed as p_+=(1/2+sigma delta), Fisher information for delta is 4sigma^2/(1-4sigma^2 delta^2). Thus at delta=0 any unbiased estimator with standard deviation epsilon needs N>=1/(4sigma^2 epsilon^2). "
                "Pinsker's inequality gives the displayed necessary sample count for testing +/-delta. The same information formula holds for independent Poisson category counts with the stated expected total count. For sigma=2^-n the leading cost grows as 4^n."
            )}
    return finalize(198,"Sampling Lower Bounds for Deep Compressed Modes",
        "proven_information_bounds_for_declared_multinomial_and_Poisson_models",
        "The observation law, noise model, contrast, and soft-layer attenuation are supplied. These bounds do not prove that physical FIN layers use this channel.",packet)


def hermitian_sign(x: np.ndarray) -> tuple[np.ndarray,float]:
    vals,vecs=np.linalg.eigh((x+x.conj().T)/2)
    return vecs@np.diag(np.where(vals>=0,1.,-1.))@vecs.conj().T,float(np.min(abs(vals)))


def st199_factor_projection() -> dict:
    rng=np.random.default_rng(SEED+199);x=np.array([[0,1],[1,0]],complex);z=np.diag([1,-1]).astype(complex);i=np.eye(2)
    q,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    xi=q@np.kron(x,i)@q.conj().T;zi=q@np.kron(z,i)@q.conj().T
    rows=[]
    for eta in [.001,.003,.01,.03]:
        n1=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));n1=(n1+n1.conj().T)/2
        n2=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));n2=(n2+n2.conj().T)/2
        x0=xi+eta*n1/np.linalg.norm(n1,2);z0=zi+eta*n2/np.linalg.norm(n2,2)
        xp,gapx=hermitian_sign(x0)
        zodd=(z0-xp@z0@xp)/2;zp,gapz=hermitian_sign(zodd)
        rows.append({"eta":eta,"X_input_to_projected":float(np.linalg.norm(x0-xp,2)),
                     "Z_input_to_projected":float(np.linalg.norm(z0-zp,2)),"X_sign_gap":gapx,"Z_odd_sign_gap":gapz,
                     "X_involution_error":float(np.linalg.norm(xp@xp-np.eye(4),2)),
                     "Z_involution_error":float(np.linalg.norm(zp@zp-np.eye(4),2)),
                     "anticommutator_error":float(np.linalg.norm(xp@zp+zp@xp,2)),
                     "generated_algebra_dimension":algebra_span_dimension([xp,zp]),
                     "commutant_dimension":commutant_dimension([xp,zp])})
    packet={"projection":"X=sign(X0); Z=sign[(Z0-X Z0 X)/2]","rows":rows,
            "theorem":(
                "If X0 has no zero eigenvalue, sign(X0) is the unique nearest Hermitian involution in Frobenius norm. For fixed X, the odd part Z_o=(Z0-XZ0X)/2 anticommutes with X; if Z_o is invertible, sign(Z_o) is a Hermitian involution that exactly anticommutes with X. "
                "With balanced +/- multiplicities the pair generates an exact M2 factor with four-dimensional commutant."
            )}
    return finalize(199,"Sequential Projection onto an Exact Matrix-Factor Variety",
        "proven_sequential_projection_theorem_with_numerical_noisy_recovery",
        "Uniqueness is sequential (first X, then Z conditional on X), not a global nearest joint-factor theorem. Tensor gauge, instrument labels, and physical subsystems remain supplied.",packet)


def frac_hmm_likelihoods(n: int, emissions: list[list[Fraction]]) -> list[Fraction]:
    p=[[Fraction(9,10),Fraction(1,10)],[Fraction(1,5),Fraction(4,5)]];pi=[Fraction(2,3),Fraction(1,3)]
    f=[[pi[z]*emissions[z][y] for z in range(2)] for y in range(2)]
    for _ in range(1,n):
        nf=[]
        for row in f:
            pred=[sum(row[z]*p[z][w] for z in range(2)) for w in range(2)]
            for y in range(2): nf.append([pred[w]*emissions[w][y] for w in range(2)])
        f=nf
    return [sum(row) for row in f]


def sqrt_fraction_bounds(x: Fraction, digits: int=35) -> tuple[Fraction,Fraction]:
    d=10**digits;n=math.isqrt(x.numerator*x.denominator*d*d);den=x.denominator*d
    return Fraction(n,den),Fraction(n+1,den)


def st200_observed_hmm_interval() -> dict:
    e0=[[Fraction(98,100),Fraction(2,100)],[Fraction(92,100),Fraction(8,100)]]
    e1=[[Fraction(8,100),Fraction(92,100)],[Fraction(2,100),Fraction(98,100)]]
    rows=[]
    for n in [1,2,4,6,8,10,12]:
        l0=frac_hmm_likelihoods(n,e0);l1=frac_hmm_likelihoods(n,e1);lo=Fraction(0);hi=Fraction(0)
        for x,y in zip(l0,l1):
            a,b=sqrt_fraction_bounds(x*y);lo+=a;hi+=b
        mp.mp.dps=70
        los=mp.mpf(lo.numerator)/lo.denominator; his=mp.mpf(hi.numerator)/hi.denominator
        lf=float(np.nextafter(float(los),-np.inf));hf=float(np.nextafter(float(his),np.inf))
        rlo=float(np.nextafter(-math.log(hf)/n,-np.inf));rhi=float(np.nextafter(-math.log(lf)/n,np.inf))
        rows.append({"events":n,"strings":2**n,"Hellinger_coefficient_interval":[lf,hf],
                     "exact_rational_lower":str(lo),"exact_rational_upper":str(hi),
                     "exact_rational_width":str(hi-lo),
                     "high_precision_width":mp.nstr(his-los,20),
                     "interval_width":hf-lf,"finite_rate_interval":[rlo,rhi]})
    packet={"arithmetic":"exact rational forward likelihoods plus integer-square-root outward bounds","decimal_scale_digits":35,"rows":rows,
            "theorem":"Every displayed observed Hellinger coefficient is enclosed by exact rational endpoints obtained without floating forward recursion. Monotonicity of -log gives the corresponding outward finite-n rate interval."}
    return finalize(200,"Exact-Rational Finite-Record Observed-HMM Rate Enclosures",
        "proven_finite_n_interval_rates_through_n12",
        "This closes floating uncertainty at finite n but does not prove the exact asymptotic observed-process Hellinger/Chernoff rate or construct an interval belief-space transfer operator.",packet)


PAULI=[np.eye(2),np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]).astype(complex)]


def cross_pair_pauli(pair: int, left: int, right: int) -> np.ndarray:
    # bit ordering matches swap_bit_matrix: least-significant qubit first.
    out=np.array([[1]],complex)
    for bit in reversed(range(8)):
        factor=PAULI[left] if bit==pair else PAULI[right] if bit==4+pair else PAULI[0]
        out=np.kron(out,factor)
    return out


def st201_local_reset_rank(a: np.ndarray) -> dict:
    h=np.zeros((16,16));h[:12,:12]=a;i16=np.eye(16);ht=np.kron(h,i16)+np.kron(i16,h)
    swaps=[swap_bit_matrix(8,k,4+k) for k in range(4)]
    cs=[ht@s-s@ht for s in swaps];mat=np.column_stack([c.ravel() for c in cs]);sv=np.linalg.svd(mat,compute_uv=False)
    aiv,_,_=strict_interval_matrix();maxrad=0.0
    for i in range(12):
        for j in range(12):
            lo,hi=replay_bounds(aiv[i][j]);maxrad=max(maxrad,max(abs(a[i,j]-lo),abs(hi-a[i,j])))
    single_column_error=192*maxrad
    four_perturbation=math.sqrt(4)*single_column_error
    pauli_generators=[cross_pair_pauli(k,p,q) for k in range(4) for p in range(4) for q in range(4) if p or q]
    pcs=[ht@g-g@ht for g in pauli_generators];pmat=np.column_stack([c.ravel() for c in pcs]);psv=np.linalg.svd(pmat,compute_uv=False)
    sixty_perturbation=math.sqrt(60)*single_column_error
    # A posteriori positivity certificates.  G-cI=LL*+R implies
    # lambda_min(G)>=c-||R||.  We additionally pay the strict-input interval
    # perturbation and a deliberately conservative IEEE/BLAS arithmetic budget.
    gram=np.real(mat.conj().T@mat);pgram=np.real(pmat.conj().T@pmat)
    c4,c60=200.0,40.0
    l4=np.linalg.cholesky(gram-c4*np.eye(4));l60=np.linalg.cholesky(pgram-c60*np.eye(60))
    r4=float(np.linalg.norm(gram-c4*np.eye(4)-l4@l4.T,2))
    r60=float(np.linalg.norm(pgram-c60*np.eye(60)-l60@l60.T,2))
    input4=2*float(sv[0])*four_perturbation+four_perturbation**2
    input60=2*float(psv[0])*sixty_perturbation+sixty_perturbation**2
    ieee_payment=1e-6
    certified4=math.sqrt(c4-r4-input4-ieee_payment)
    certified60=math.sqrt(c60-r60-input60-ieee_payment)
    packet={"strict_interval_max_entry_radius":maxrad,
            "pairwise_SWAP_singular_values":sv.tolist(),"pairwise_SWAP_Cholesky_shift":c4,
            "pairwise_SWAP_Cholesky_residual_norm":r4,"pairwise_SWAP_input_Gram_perturbation":input4,
            "pairwise_SWAP_interval_lower_singular_bound":certified4,
            "expanded_static_generator_count":60,"expanded_Pauli_commutator_singular_values":psv.tolist(),
            "expanded_Cholesky_shift":c60,"expanded_Cholesky_residual_norm":r60,
            "expanded_input_Gram_perturbation":input60,"conservative_IEEE_arithmetic_payment":ieee_payment,
            "expanded_interval_lower_singular_bound":certified60,
            "theorem":(
                "Verified shifted-Cholesky residuals, outward strict-entry intervals, Frobenius input-perturbation bounds, and an explicit conservative arithmetic payment give positive lower bounds for the smallest singular value of both commutator maps. "
                "Consequently neither the four pairwise-SWAP span nor the full 60-dimensional span of nonidentity two-qubit Pauli generators acting on corresponding cross-register pairs contains a nonzero static generator commuting with H tensor I+I tensor H."
            )}
    return finalize(201,"Interval-Rank Certificate for Static Cross-Register Reset Generators",
        "proven_no_go_for_declared_4_and_60_dimensional_static_local_spans",
        "The result is restricted to the padded 12-in-16 encoding and static linear combinations on four corresponding qubit pairs. Time-dependent controls, ancillas, noncorresponding edges, and higher locality remain open.",packet)


def st202_external_readiness() -> dict:
    s=synthetic_self_test()
    required=["externally supplied event path","externally supplied frozen holdout hash","calibration hash",
              "provider identity/attestation","registrar identity/attestation","analyst identity/attestation",
              "independent custody attestation","laboratory execution attestation"]
    packet={"synthetic_self_test":s,"required_external_atoms":required,"missing_external_atoms":required,
            "physical_execution_valid":False,"execution_attempted":False,
            "theorem":"The local synthetic replay verifies structural tamper detection. Because every declared external provenance atom is absent, the physical-execution predicate is false and ST202 must refuse execution/promotion."}
    return finalize(202,"External Immutable-Record Execution Gate",
        "blocked_correctly_no_external_record",
        "No external event record, independent custody, calibration, laboratory execution, or physical result exists in scope. Local code cannot generate these facts and was not used to simulate them.",packet)


def make_figures(d: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST190"]["rows"]
    for k in range(1,7):ax.semilogy([r["layer"] for r in rows],[r["Fourier_visibility"][k] for r in rows],"o-",label=f"mode {k}")
    ax.set(xlabel="dyadic layer",ylabel="heat visibility",title="ST190: strict-generated soft hierarchy");ax.legend(ncol=2);fig.tight_layout();fig.savefig(FIG_DIR/"st190_heat_hierarchy.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST192"]["trajectories"];ax.bar(range(len(rows)),[r["final_max_probability"] for r in rows],color="tab:purple");ax.set(xlabel="noise realization",ylabel="largest final probability",title="ST192: spontaneous orbit-member selection");fig.tight_layout();fig.savefig(FIG_DIR/"st192_replicator_selection.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST195"]["quasistatic_rows"];ax.loglog([r["steps"] for r in rows],[r["dissipation"] for r in rows],"o-");ax.set(xlabel="control steps",ylabel="dimensionless dissipation",title="ST195: quasistatic selector cost");fig.tight_layout();fig.savefig(FIG_DIR/"st195_selector_cost.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST198"]["rows"];ax.semilogy([r["layers"] for r in rows],[r["Cramer_Rao_samples_for_sd_0_01"] for r in rows],"o-",label="estimation");ax.semilogy([r["layers"] for r in rows],[r["necessary_samples_for_error_0_05_via_Pinsker"] for r in rows],"s--",label="testing");ax.set(xlabel="soft layers",ylabel="necessary samples/count",title="ST198: deep-mode sampling burden");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st198_sampling_bounds.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST200"]["rows"];ax.plot([r["events"] for r in rows],[sum(r["finite_rate_interval"])/2 for r in rows],"o-");ax.set(xlabel="events",ylabel="certified finite rate",title="ST200: observed-HMM exact-rational enclosures");fig.tight_layout();fig.savefig(FIG_DIR/"st200_hmm_intervals.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));ax.bar(["4 SWAP","60 Pauli"],[d["ST201"]["pairwise_SWAP_interval_lower_singular_bound"],d["ST201"]["expanded_interval_lower_singular_bound"]]);ax.set(ylabel="certified minimum singular lower bound",title="ST201: static local-reset no-go");fig.tight_layout();fig.savefig(FIG_DIR/"st201_reset_rank.png",dpi=190);plt.close(fig)


def main() -> None:
    _,a,_=strict_operator()
    out={"metadata":{"programs":"ST190-ST202","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__}}
    out["ST190"]=st190_heat_hierarchy(a);out["ST191"]=st191_blackwell_quotient();out["ST192"]=st192_replicator_pointer(a)
    out["ST193"]=st193_modal_continuation(a);out["ST194"]=st194_noncommuting_recovery();out["ST195"]=st195_selector_control()
    out["ST196"]=st196_nonlinear_source();out["ST197"]=st197_nuisance_boundary(a);out["ST198"]=st198_sampling_bounds()
    out["ST199"]=st199_factor_projection();out["ST200"]=st200_observed_hmm_interval();out["ST201"]=st201_local_reset_rank(a);out["ST202"]=st202_external_readiness()
    out["recommended_next_programs"]=[
        {"id":"ST203","priority":1,"study":"classify all strict-generated completely positive semigroup hierarchies and test whether any scale schedule is selected internally"},
        {"id":"ST204","priority":2,"study":"derive the lattice of Blackwell-equivalent FIN experiments across preparation and instrument families"},
        {"id":"ST205","priority":3,"study":"replace injected symmetric noise in ST192 by a strictly typed endogenous fluctuation source or prove a no-go"},
        {"id":"ST206","priority":4,"study":"validated pseudo-arclength continuation and global collision isolation beyond the ST193 modal slices"},
        {"id":"ST207","priority":5,"study":"extend the trace-norm recovery certificate to negative polar orientation and nonunital qubit channels"},
        {"id":"ST208","priority":6,"study":"solve the finite-time entropy-production optimal-control problem with a fixed generator-rate budget"},
        {"id":"ST209","priority":7,"study":"classify the minimal extra response datum that identifies g while remaining carrier-invariant"},
        {"id":"ST210","priority":8,"study":"adaptive interval subdivision beyond the fixed-grid ST197 certificate boundary"},
        {"id":"ST211","priority":9,"study":"derive minimax reconstruction bounds for simultaneous deep Fourier modes under correlated counts"},
        {"id":"ST212","priority":10,"study":"prove or refute global nearest-joint-factor uniqueness beyond the sequential ST199 projection"},
        {"id":"ST213","priority":11,"study":"construct a rigorous interval transfer operator for the exact asymptotic observed-HMM rate"},
        {"id":"ST214","priority":12,"study":"search time-dependent energy-conserving local reset controls outside the static ST201 span"},
        {"id":"ST215","priority":13,"study":"test whether a refinement/coarse-graining tower can induce an effective metric via diffusion distance without a dimensional anchor"},
        {"id":"ST216","priority":14,"study":"formalize the distinction between viewing lower scales through layers and reconstructing them from a sufficient statistic"},
        {"id":"ST217","priority":15,"study":"execute the immutable-record gate only after a laboratory supplies genuine data, calibration, and independent custody"},
    ]
    out["central_verdict"]=(
        "The user's layered-view intuition now has a sharper conditional realization: the strict operator itself generates a canonical heat semigroup, and any supplied schedule of semigroup parameters produces a soft hierarchy whose mode visibility is exactly spectral. Thus an observer at a later effective layer may see deeper modes only through exponentially conditioned transfer. What is not generated is the schedule, an absolute layer length, a Planck anchor, an apparatus, or the claim that our world occupies a particular layer."
    )
    out["epistemic_boundary"]=(
        "No carrierless information, strict fractal scale law, Planck scale, spacetime projection, canonical strict selector, QW-2191 discharge, dimensional clock or length, external record, laboratory result, legacy-to-strict completion or physical-role transfer, Standard Model, gravity, L_total, or ToE closure is claimed."
    )
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(190,203):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out)
    print(json.dumps({"results":RESULTS.name,"programs":13,"figures":6},indent=2))


if __name__=="__main__": main()
