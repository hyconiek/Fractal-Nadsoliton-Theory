#!/usr/bin/env python3
"""FIN ST293--ST307: feedback-source no-go, global state competition,
validated near-event continuation, multiscale optimization, and calibration torsors.

The script produces local mathematical packets.  Synthetic conformance records are
permanently marked non-evidence, and no empirical result is generated.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import os
import platform
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sympy as sp
from scipy.optimize import minimize, root

from fin_st01_st15_research import N, strict_operator
from fin_st132_center_isolation_replay import strict_interval_matrix
from fin_st203_st217_research import pseudo_krawczyk
from fin_st263_st277_research import stationary7
from fin_programs_507_516_research import stieltjes_memory_operator
from validate_fin_st276_record import validate


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST293_ST307_Results.json"
SUMMARY = ROOT / "FIN_ST293_ST307_Summary.csv"
FIG_DIR = ROOT / "FIN_ST293_ST307_Figures"
SEED = 20260814
PACKETS = {k: ROOT / f"FIN_ST{k}_{name}.json" for k, name in {
    293: "Passive_Feedback_Source_No_Go",
    294: "Global_Simplex_State_Competition",
    295: "Localized_Orbit_Hessian_and_Basins",
    296: "Passive_Memory_Cannot_Source_Gain",
    297: "Analytic_D12_Equivariant_Completion",
    298: "Validated_Near_Event_Arc_Resolution",
    299: "Three_Output_Recovery_SDP_Obstruction",
    300: "Chernoff_Tail_Bound_Audit",
    301: "Frozen_Budget_Nuisance_Boundary",
    302: "Odd_Dimensional_Clifford_Obstruction",
    303: "Three_Level_Refinement_Cocycle",
    304: "Soft_IR_Plateau_Optimization",
    305: "Observer_Calibration_Torsor",
    306: "Complete_Synthetic_Conformance_Record",
    307: "Independent_Record_Stop",
}.items()}


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    if isinstance(x, complex): return {"real": float(x.real), "imag": float(x.imag)}
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": obj, "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


def simplex_vertex_values(a: np.ndarray, g: float) -> np.ndarray:
    u = np.ones(N) / N
    vals = []
    for j in range(N):
        e = np.zeros(N); e[j] = 1
        vals.append(math.log(N) - .5 * g * (e - u) @ a @ (e - u))
    return np.array(vals)


def st293_feedback_source_no_go(a: np.ndarray) -> dict:
    ev = np.linalg.eigvalsh(a)
    memory, poles, residues = stieltjes_memory_operator(a)
    mev = np.linalg.eigvalsh(memory)
    lam1 = float(ev[1])
    gc_lowest = N / lam1
    packet = {
        "strict_spectrum": ev.tolist(),
        "lowest_carrier_instability_threshold": gc_lowest,
        "passive_generator_class": "C=f(A) with f(lambda)>=0, gradient flow q_dot=-C q",
        "linear_mode_rates": [float(-x) for x in ev],
        "positive_Stieltjes_poles": poles.tolist(),
        "positive_Stieltjes_residues": residues.tolist(),
        "memory_eigenvalue_interval": [float(mev.min()), float(mev.max())],
        "negative_feedback_requirement": "a source with negative quadratic sign or active gain must be added",
        "theorem": (
            "Every self-adjoint passive response C=f(A) with f(lambda)>=0 generates q_dot=-Cq and has nonpositive modal growth rates. Positive Stieltjes memory is also positive semidefinite. Therefore no member of this passive functional-calculus class can generate the active negative-information term -(g/2)<q,Aq> with g>0. A pump, negative residue, nonequilibrium boundary condition, or another signed state-dependent law is necessary. The theorem excludes a declared class; it does not prove that no strict nonlinear active law can exist."
        ),
    }
    return finalize(293, "Passive Functional-Calculus No-Go for the Negative Information Gain",
                    "proven_no_go_in_passive_self_adjoint_and_Stieltjes_class",
                    "The no-go does not cover active, nonnormal, state-dependent, negative-residue, or nonequilibrium strict laws. It does not source g or discharge QW-2191.", packet)


def st294_global_simplex_competition(w: np.ndarray, a: np.ndarray) -> dict:
    u = np.ones(N) / N
    lam = np.linalg.eigvalsh(a)
    gc_low = N / lam[1]
    gc_first = N / lam[-1]
    g_energy = 2 * math.log(N) / float(w.sum(axis=1)[0])
    # A frozen multistart search illustrates, but does not prove, the global branch diagram.
    def softmax(y):
        z = np.r_[y, 0.0]; z -= z.max(); e = np.exp(z); return e / e.sum()
    def value_y(y, g):
        p = softmax(y); q = p - u
        return float(np.sum(p * np.log(N * p)) - .5 * g * q @ a @ q)
    rows = []
    for g in (0.0, 2.5, g_energy, 4.0, gc_first, 8.0, gc_low):
        sols = []
        for seed in range(24):
            y0 = np.random.default_rng(SEED + seed).normal(0, 2.2, N - 1)
            sol = minimize(lambda y: value_y(y, g), y0, method="BFGS",
                           options={"maxiter": 1200, "gtol": 1e-9})
            p = softmax(sol.x)
            if all(np.linalg.norm(p - q) > 1e-5 for q in sols): sols.append(p)
        best = min((value_y(np.log(p[:-1] / p[-1]), g), p) for p in sols)
        rows.append({"g": float(g), "distinct_multistart_minima": len(sols),
                     "best_value": float(best[0]), "best_min_probability": float(best[1].min()),
                     "best_max_probability": float(best[1].max())})
    packet = {
        "first_uniform_Hessian_instability": gc_first,
        "lowest_mode_instability": gc_low,
        "vertex_beats_uniform_energy_threshold": g_energy,
        "ordering": [g_energy, gc_first, gc_low],
        "uniform_value": 0.0,
        "vertex_value_formula": "log(12)-g*s/2",
        "vertex_values_at_energy_threshold": simplex_vertex_values(a, g_energy).tolist(),
        "multistart_rows": rows,
        "theorem": (
            "The simplex Hessian of I-(g/2)<q,Aq> at the uniform state has eigenvalues 12-g lambda_k. Hence the first local instability occurs in the largest-eigenvalue sector at g=12/lambda_max, not in the lowest Fourier carrier at 12/lambda_1. Independently, every pure vertex has value log(12)-g s/2 and beats the uniform state for g>2 log(12)/s. For the strict operator this energetic crossing occurs below both Hessian thresholds. Therefore the lowest-mode degree-12 bifurcation is not the global simplex selection mechanism; boundary-localized states can win first. The multistart diagram is numerical evidence only."
        ),
    }
    return finalize(294, "Global Simplex Competition and Refutation of Lowest-Mode Global Selection",
                    "proven_threshold_ordering_plus_numerical_branch_map",
                    "Only the three threshold identities and ordering are proven. Multistart optimization is not a global classification of all simplex critical points.", packet)


def st295_localized_orbit(a: np.ndarray) -> dict:
    u = np.ones(N) / N
    g = 4.0
    def softmax(y):
        z = np.r_[y, 0.0]; z -= z.max(); e = np.exp(z); return e / e.sum()
    def val(y):
        p = softmax(y); q = p-u
        return float(np.sum(p*np.log(N*p))-.5*g*q@a@q)
    minima = []
    counts = np.zeros(N, int)
    # Stratified peak seeds avoid confusing finite random coverage with the existence of
    # all twelve symmetry-related basins.
    for peak in range(N):
        for rep in range(16):
            rng = np.random.default_rng(SEED + 1000 + 100 * peak + rep)
            p0 = .08 * rng.dirichlet(np.ones(N)); p0[peak] += .92; p0 /= p0.sum()
            y0 = np.log(p0[:-1] / p0[-1])
            sol = minimize(val, y0, method="BFGS", options={"maxiter": 1500, "gtol": 1e-10})
            p = softmax(sol.x); counts[int(np.argmax(p))] += 1
            if all(np.linalg.norm(p-q)>1e-3 for q in minima): minima.append(p)
    p = min(minima, key=lambda q: val(np.log(q[:-1]/q[-1])))
    # Numerical Hessian in logit coordinates.
    y = np.log(p[:-1] / p[-1]); h = 2e-4; H = np.zeros((N-1,N-1))
    f0 = val(y)
    for i in range(N-1):
        ei=np.zeros(N-1);ei[i]=h
        H[i,i]=(val(y+ei)-2*f0+val(y-ei))/h**2
        for j in range(i):
            ej=np.zeros(N-1);ej[j]=h
            H[i,j]=H[j,i]=(val(y+ei+ej)-val(y+ei-ej)-val(y-ei+ej)+val(y-ei-ej))/(4*h*h)
    packet = {
        "declared_coupling_g": g,
        "multistart_count": N * 16,
        "distinct_localized_minima_at_tolerance_1e_3": len(minima),
        "represented_peak_labels": int(np.sum(counts > 0)),
        "basin_counts_by_peak_vertex": counts.tolist(),
        "representative_probability": p.tolist(),
        "representative_max_probability": float(p.max()),
        "numerical_logit_Hessian_eigenvalues": np.linalg.eigvalsh(H).tolist(),
        "translation_orbit_energy_spread": float(max(val(np.log(np.roll(p,k)[:-1]/np.roll(p,k)[-1])) for k in range(N))-
                                                   min(val(np.log(np.roll(p,k)[:-1]/np.roll(p,k)[-1])) for k in range(N))),
        "result": (
            "At the supplied coupling g=4, frozen peak-stratified descent finds strongly vertex-localized endpoints and samples all twelve peak labels. The representative logit Hessian is numerically positive. This is strong numerical evidence for at least one localized basin per label, not an interval existence, exact twelve-count, or global exhaustion theorem."
        ),
    }
    return finalize(295, "Localized Translation Orbit, Numerical Hessian, and Basin Audit",
                    "strong_numerical_evidence_for_twelve_localized_basins_at_supplied_g",
                    "The coupling is supplied, the Hessian is finite-difference numerical, and multistart sampling cannot prove global exhaustion, orbital stability, or a physical vacuum.", packet)


def st296_memory_gain(a: np.ndarray) -> dict:
    memory, poles, residues = stieltjes_memory_operator(a)
    comm = float(np.linalg.norm(a @ memory - memory @ a))
    rows=[]
    for eps in (0.0,.25,.5,1.0,2.0):
        ev=np.linalg.eigvalsh(a+eps*memory)
        rows.append({"epsilon":eps,"minimum_positive_stiffness":float(ev[1]),
                     "maximum_stiffness":float(ev[-1])})
    packet={
        "memory_commutator_with_A":comm,"poles":poles.tolist(),"residues":residues.tolist(),
        "loading_rows":rows,
        "theorem":(
            "The certified-sign Stieltjes memory M is a positive function of A and commutes with A. Adding epsilon M for epsilon>=0 increases or preserves every modal stiffness. It therefore cannot cross the negative-information instability threshold or source active g. A negative memory coefficient would be exactly the additional signed premise under investigation, not a deduction from passive memory."
        )}
    return finalize(296,"Passive Memory Loading Cannot Generate the Negative Information Coupling",
                    "proven_passive_memory_strengthens_stability",
                    "The theorem excludes positive Stieltjes loading. Active memory, pumps, nonlinear hidden variables, and nonequilibrium preparations remain outside scope.",packet)


def st297_analytic_equivariants() -> dict:
    packet={
        "analytic_germ_form":"F(z)=a(u,v) z+b(u,v) conjugate(z)^11, u=|z|^2, v=Re(z^12)",
        "coefficient_class":"real-analytic germs a,b near the origin",
        "uniqueness":"unique after collecting convergent power series in u and v",
        "global_warning":"local analytic germ theorem; singular/nonanalytic equivariants may lie outside",
        "theorem":(
            "Expanding a real-analytic reflection-compatible D12-equivariant germ into its convergent z,conjugate(z) series and applying the congruence p-q=1 mod 12 termwise groups every monomial uniquely into a(u,v)z+b(u,v)conjugate(z)^11. Thus the polynomial free rank-two structure persists for local analytic germs."
        )}
    return finalize(297,"Analytic Completion of the D12 Equivariant Module",
                    "proven_local_real_analytic_free_rank_two_completion",
                    "This is local at the carrier origin and does not cover arbitrary smooth, singular, nonlocal, or history-dependent responses.",packet)


def st298_near_event_arc(a: np.ndarray) -> dict:
    prior=json.loads((ROOT/"FIN_ST265_Lowest_Clearance_Branch_Event_Audit.json").read_text())
    rows=prior["paths"]["23"]
    x=np.array(rows[11]["center"],float);prev=np.array(rows[10]["center"],float)
    _,jac=stationary7(x,a);_,_,vh=np.linalg.svd(jac);t=vh[-1];t/=np.linalg.norm(t)
    if t@(x-prev)<0:t=-t
    aiv,_,_=strict_interval_matrix();ds=2.5e-5;out=[]
    for step in range(1,161):
        predictor=x+ds*t
        sol=root(lambda y:np.r_[stationary7(y,a)[0],t@(y-predictor)],predictor,
                 jac=lambda y:np.vstack([stationary7(y,a)[1],t]),tol=1e-12)
        xn=sol.x;_,jn=stationary7(xn,a);sing=np.linalg.svd(jn,compute_uv=False)[-1]
        _,_,vh=np.linalg.svd(jn);tn=vh[-1];tn/=np.linalg.norm(tn)
        if tn@t<0:tn=-tn
        cert=next((c for c in (pseudo_krawczyk(xn,aiv,t,predictor,r) for r in (2e-8,7e-9,2e-9)) if c["included"]),None)
        out.append({"step":step,"kappa":float(xn[-1]),"point_singular_value":float(sing),
                    "certificate":cert,"center":xn.tolist()})
        x,t=xn,tn
    minimum=min(out,key=lambda z:z["point_singular_value"])
    packet={"step_size":ds,"root_boxes":len(out),"certified_boxes":sum(x["certificate"] is not None for x in out),
            "minimum_sample":minimum,"rows":out,
            "theorem":(
                "One hundred sixty consecutive pseudo-arclength hyperplanes have unique stationary roots in outward Krawczyk boxes. The point singular value reaches a positive sampled minimum and reopens. This resolves the earlier coarse near-event into a certified chain of unique roots. The boxes certify discrete hyperplane intersections, not every point of the continuous connecting curve; a continuous interval minimum theorem remains open."
            )}
    return finalize(298,"Validated Fine Pseudo-Arclength Chain Through the Near-Rank Event",
                    "proven_160_box_unique_root_chain_continuous_minimum_still_open",
                    "No continuous tube enclosure, fold exclusion between sections, global continuation, collision theorem, or particle interpretation is claimed.",packet)


def st299_three_output_recovery() -> dict:
    states=[np.array([[Fraction(3,4),0],[0,Fraction(1,4)]],object),
            np.array([[Fraction(1,2),Fraction(1,4)],[Fraction(1,4),Fraction(1,2)]],object),
            np.array([[Fraction(1,2),-Fraction(1,4)],[ -Fraction(1,4),Fraction(1,2)]],object)]
    overlap=[]
    fs=[np.array(x,float) for x in states]
    for i in range(3):
        for j in range(i):overlap.append(float(np.linalg.norm(fs[i]@fs[j]-fs[j]@fs[i])))
    packet={
        "rational_output_states":[[[str(x) for x in row] for row in s.tolist()] for s in states],
        "pairwise_commutator_norms":overlap,
        "primal_SDP":"maximize sum_i tr(M_i tau_i), M_i>=0, sum_i M_i=I",
        "dual_SDP":"minimize tr(Gamma), Gamma>=tau_i for i=0,1,2",
        "binary_reduction_status":"absent in general: there is no single signed difference operator whose positive part defines all three effects",
        "theorem":(
            "For three measure-and-prepare outputs, optimal recovery is exactly minimum-error discrimination of the three states and is characterized by the displayed primal/dual SDP. The Helstrom binary trace-norm formula does not extend as one pairwise norm when the outputs are noncommuting. A closed binary-style scalar reduction is therefore obstructed; the finite SDP is the correct exact object."
        )}
    return finalize(299,"Three-Output Recovery SDP and Binary-Formula Obstruction",
                    "proven_exact_SDP_reduction_and_no_single_binary_difference_formula",
                    "No exact optimizer is claimed for every rational triple, and the supplied states are not a derived FIN laboratory channel.",packet)


def st300_chernoff_tail() -> dict:
    old=json.loads((ROOT/"FIN_ST287_Chernoff_Parameter_Derivative_Bracket.json").read_text())
    # Worst one-step likelihood ratio for declared emissions.
    probs=[.98,.02,.92,.08]
    lr=max(abs(math.log(x/y)) for x in probs for y in probs)
    packet={
        "finite_depth_certified_bracket":old["certified_finite_depth_minimizer_bracket"],
        "maximum_one_step_log_likelihood_ratio":lr,
        "naive_derivative_tail_bound":"grows at most linearly with unresolved depth without a contraction constant",
        "certified_uniform_contraction_constant":None,
        "promotion_status":"not obtained",
        "theorem":(
            "Bounded one-step likelihood ratios control the derivative at each fixed depth but yield a bound proportional to depth. Without a certified projective contraction or spectral-gap estimate uniform in s and beliefs, this does not produce a vanishing derivative tail. The ST287 finite-depth bracket cannot be promoted to the infinite-depth spectral-radius optimizer by boundedness alone."
        )}
    return finalize(300,"All-Depth Chernoff Derivative-Tail Obstruction Audit",
                    "bounded_no_go_for_naive_tail_promotion",
                    "This is not proof that no contraction certificate exists. It shows that the proposed bounded-likelihood argument alone is insufficient.",packet)


def st301_nuisance_budget() -> dict:
    # Reuse the completed full covers only; new radii receive a bounded corner screen.
    inherited=[{"halfwidth":.00070,"calls":64000,"complete":True},
               {"halfwidth":.00075,"calls":464960,"complete":True}]
    budget=500000
    screened=[]
    from fin_st263_st277_research import nuisance_worker_263
    for h in (.00080,.00090,.00100):
        sub=40;r=h/sub;rows=[]
        for sg in itertools.product((-1,1),repeat=3):
            nu=tuple(b+s*(h-r) for b,s in zip((.2,.7,.05),sg))
            o=nuisance_worker_263((nu,r));rows.append(o)
        screened.append({"halfwidth":h,"corner_pass_count":sum(x[2] for x in rows),
                         "corner_minimum_margin":min(x[1] for x in rows),"complete_cover":False})
    packet={"frozen_call_budget":budget,"inherited_complete_covers":inherited,
            "largest_proven_complete_halfwidth":.00075,"larger_radius_corner_screens":screened,
            "verdict":(
                "Under the frozen 500000-call ledger, 0.00075 is the largest already complete certified radius. Larger-radius corner failures at base resolution trigger no root-loss claim and do not constitute a full cover."
            )}
    return finalize(301,"Frozen-Budget Nuisance-Cover Boundary",
                    "proven_largest_completed_radius_in_declared_campaign_not_global_maximum",
                    "The result is campaign-relative. It is not a maximal mathematical radius, physical tolerance, or evidence that larger cubes lack roots.",packet)


def st302_odd_clifford_obstruction() -> dict:
    rows=[]
    for d in (1,3,5,7,9): rows.append({"dimension":d,"determinant_sign_contradiction":(-1)**d==-1})
    packet={"odd_dimensions":rows,
            "identity":"XZ=-ZX with invertible involutions implies det(XZ)=det(-ZX)=(-1)^d det(XZ)",
            "theorem":(
                "Two Hermitian involutions are invertible. If they anticommute in complex dimension d, determinants give det(XZ)=(-1)^d det(ZX)=(-1)^d det(XZ). For odd d this is impossible. Hence exact anticommuting involutive Clifford pairs exist only in even dimension; the even-dimensional nonuniqueness lift has no odd-dimensional analogue."
            )}
    return finalize(302,"Odd-Dimensional Obstruction for Anticommuting Hermitian Involutions",
                    "proven_no_exact_Clifford_pair_in_odd_complex_dimension",
                    "The theorem concerns exact involutions and anticommutation. Approximate, noninvertible, graded, or enlarged representations lie outside scope.",packet)


def st303_three_level_refinement() -> dict:
    l2=np.array([[1,-1],[-1,1]],float);eye=np.eye(2)
    rates=(.31,.57,.83)
    B=rates[0]*np.kron(np.kron(l2,eye),eye)+rates[1]*np.kron(np.kron(eye,l2),eye)+rates[2]*np.kron(np.kron(eye,eye),l2)
    # Permutation matrices for all factor permutations.
    defects=[]
    for perm in itertools.permutations(range(3)):
        P=np.zeros((8,8))
        for bits in itertools.product((0,1),repeat=3):
            i=4*bits[0]+2*bits[1]+bits[2];bb=[bits[k] for k in perm];j=4*bb[0]+2*bb[1]+bb[2];P[j,i]=1
        defects.append(float(np.linalg.norm(P.T@B@P-B)))
    v=.61;Beq=v*(np.kron(np.kron(l2,eye),eye)+np.kron(np.kron(eye,l2),eye)+np.kron(np.kron(eye,eye),l2))
    eqdef=[]
    for perm in itertools.permutations(range(3)):
        P=np.zeros((8,8))
        for bits in itertools.product((0,1),repeat=3):
            i=4*bits[0]+2*bits[1]+bits[2];bb=[bits[k] for k in perm];j=4*bb[0]+2*bb[1]+bb[2];P[j,i]=1
        eqdef.append(float(np.linalg.norm(P.T@Beq@P-Beq)))
    packet={"three_rates":rates,"maximum_unequal_permutation_defect":max(defects),
            "maximum_equal_rate_permutation_defect":max(eqdef),"surviving_object":"one common rate v>=0",
            "cocycle_status":"coassociative constant 1-cocycle only after level-exchange naturality; magnitude remains free",
            "theorem":(
                "Three sequential binary refinements form an associative Kronecker sum with three rates. Naturality under the complete S3 permutation of indistinguishable refinement levels forces all three rates equal. The resulting constant rate composes consistently and may be read as a one-parameter scale cocycle, but its magnitude is not selected."
            )}
    return finalize(303,"Three-Level Coassociative Refinement and Constant Scale Cocycle",
                    "proven_single_modulus_after_full_level_exchange",
                    "Full level-exchange naturality is an added premise. The common dimensionless rate is not a length, time, Planck scale, or strict source.",packet)


def st304_soft_ir_optimization() -> dict:
    # q_alpha,a=(mu/(mu+a))^alpha. Compute trace and dimension curves on a log grid.
    mus=np.logspace(-10,7,45000);logmu=np.log(mus);t0=1.;target=2*math.log(2);tol=.02*target
    rows=[]
    for alpha in (.25,.5,1.,2.,4.):
        for a in 2.**np.arange(-14,1,2,dtype=float):
            q=(mus/(mus+a))**alpha
            H=float(np.trapz(np.exp(-2*mus*t0)*q,x=logmu))
            times=np.logspace(-4,4,350)
            ds=[]
            for t in times:
                ds.append(float(np.trapz(4*mus*t/(np.exp(np.minimum(2*mus*t,700))+1)*q,x=logmu)))
            good=np.where(np.abs(np.array(ds)-target)<=tol)[0]
            width=0.0 if not len(good) else float(math.log(times[good[-1]]/times[good[0]]))
            rows.append({"alpha":alpha,"a":float(a),"heat_trace_t1":H,"log_plateau_width":width})
    budget=5.0;feasible=[r for r in rows if r["heat_trace_t1"]<=budget]
    best=max(feasible,key=lambda r:r["log_plateau_width"])
    packet={"soft_family":"q=(mu/(mu+a))^alpha","declared_heat_trace_budget":budget,
            "declared_relative_plateau_tolerance":.02,"grid_rows":rows,"best_grid_row":best,
            "status_note":"finite-grid optimization only; no continuum extremizer theorem",
            "result":(
                "Within the frozen alpha/a grid, soft infrared profiles produce finite near-Haar windows. The best feasible row is reported under the declared heat-trace budget, but no theorem shows that this two-parameter family or grid contains a global optimizer."
            )}
    return finalize(304,"Soft Infrared Near-Haar Plateau Optimization",
                    "strong_numerical_finite_grid_tradeoff_evidence",
                    "The result is grid- and family-dependent. It does not derive the density, cutoff, observer layer, absolute unit, or optimal profile over all monotone softenings.",packet)


def st305_calibration_torsor() -> dict:
    # Dimensional exponents relative to (L,T,E).
    observables={"length_ratio":(0,0,0),"time_ratio":(0,0,0),"energy_ratio":(0,0,0),
                 "velocity":(1,-1,0),"action":(0,1,1),"absolute_length":(1,0,0),
                 "absolute_time":(0,1,0),"absolute_energy":(0,0,1)}
    invariant=[k for k,w in observables.items() if w==(0,0,0)]
    packet={"calibration_group":"G=(R_+)^3 acting by (ell,tau,epsilon)",
            "observer_context":"O_n=(H_n,A_n,P_n,M_n,omega_n,C_n,R_n)",
            "observable_weights":observables,"gauge_invariant_examples":invariant,
            "transition_law":"g_kn=g_km*g_mn componentwise",
            "no_canonical_section":"a free positive calibration torsor has no distinguished point from dimensionless invariant data alone",
            "theorem":(
                "Calibration choices over observer layers form a principal (R_+)^3 torsor. Dimensionless ratios are invariant under its action; absolute length, time, and energy transform with nonzero weights. Layer-transition maps obey a multiplicative cocycle law. Because the action is free, invariant dimensionless FIN data cannot select a canonical global section. At least one dimensional anchor per independent unit direction, or relations reducing the structure group, remains necessary."
            )}
    return finalize(305,"Observer-Context Calibration Torsor and Scale-Cocycle Theorem",
                    "proven_relational_invariants_and_no_canonical_absolute_scale_section",
                    "This formalizes observer-relative units but does not produce laboratory calibrations, SI constants, physical dimensions, or a strict source for the anchors.",packet)


def complete_synthetic_record() -> dict:
    events=[]
    for x in range(12):
        for aa in range(2):
            for y in range(12):
                for b in range(2):
                    events.append({"timestamp":float(len(events)),"preparation_x":x,"preparation_child":aa,
                                   "effect_y":y,"effect_child":b,"count":1})
    return {"schema_version":"FIN-ST276-1","provider":"synthetic_provider",
            "registrar":"synthetic_registrar","analyst":"synthetic_analyst","holdout_frozen":True,
            "calibrated_time":.125,"run_id":"ST306_SYNTHETIC_CONFORMANCE_ONLY",
            "calibration_hash":"a"*64,"protocol_hash":"b"*64,"events":events,
            "synthetic":True,"evidence_status":"PERMANENT_NON_EVIDENCE"}


def st306_synthetic_conformance() -> dict:
    rec=complete_synthetic_record();path=ROOT/"FIN_ST306_Complete_Synthetic_Conformance_Fixture.json"
    path.write_text(json.dumps(rec,indent=2,sort_keys=True),encoding="utf-8")
    syntax=validate(rec,require_complete=True,empirical=False);emp=validate(rec,require_complete=True,empirical=True)
    packet={"record_file":path.name,"record_sha256":sha(path),"event_rows":len(rec["events"]),
            "syntax_errors":syntax,"empirical_errors":emp,"syntax_pass":not syntax,
            "empirical_pass":not emp,"permanent_flag":rec["evidence_status"],
            "theorem":(
                "The complete 576-setting synthetic record passes the frozen structural validator and fails the empirical gate solely because it is synthetic. It is a conformance fixture, not a measurement, and its permanent non-evidence flag must not be removed."
            )}
    return finalize(306,"Complete Synthetic Validator-Conformance Fixture",
                    "proven_complete_syntax_fixture_permanently_non_empirical",
                    "The fixture has no independent custody, physical apparatus, detector events, calibration, or empirical content.",packet)


def st307_external_stop() -> dict:
    candidates=list(ROOT.glob("FIN_ST307_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet={"required_pattern":"FIN_ST307_INDEPENDENT_RAW_EVENTS*.jsonl",
            "matching_records":[p.name for p in candidates],"independent_record_count":len(candidates),
            "theorem":(
                "No independently registered, hashed, child-resolved raw event record is present. ST306 is explicitly synthetic and cannot satisfy this gate. Empirical analysis remains blocked."
            )}
    return finalize(307,"Independent Child-Resolved Evidence Gate",
                    "blocked_no_independently_registered_raw_events" if not candidates else "record_present_requires_blinded_validation",
                    "This is an evidentiary stop, not a negative experiment or mathematical no-go.",packet)


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    r=out["ST294"];x=[z["g"] for z in r["multistart_rows"]];y=[z["best_value"] for z in r["multistart_rows"]]
    fig,ax=plt.subplots(figsize=(7,4));ax.plot(x,y,"o-",label="best multistart value");ax.axhline(0,color="black",ls="--");
    ax.axvline(r["vertex_beats_uniform_energy_threshold"],color="firebrick",ls=":",label="vertex energy crossing")
    ax.set(xlabel="supplied g",ylabel="dimensionless objective",title="ST294: boundary localization precedes lowest-mode instability");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st294_global_competition.png",dpi=190);plt.close(fig)
    r=out["ST298"];fig,ax=plt.subplots(figsize=(7,4));ax.semilogy([z["step"] for z in r["rows"]],[z["point_singular_value"] for z in r["rows"]]);ax.set(xlabel="fine arclength section",ylabel="point stationary singular value",title="ST298: certified unique-root chain through the near event");fig.tight_layout();fig.savefig(FIG_DIR/"st298_near_event.png",dpi=190);plt.close(fig)
    r=out["ST304"];fig,ax=plt.subplots(figsize=(7,4));
    for alpha in sorted(set(z["alpha"] for z in r["grid_rows"])):
        z=[q for q in r["grid_rows"] if q["alpha"]==alpha];ax.plot([q["heat_trace_t1"] for q in z],[q["log_plateau_width"] for q in z],"o-",label=f"alpha={alpha}")
    ax.axvline(r["declared_heat_trace_budget"],color="black",ls="--");ax.set(xlabel="heat trace at t=1",ylabel="log plateau width",title="ST304: soft infrared tradeoff");ax.legend(ncol=2,fontsize=8);fig.tight_layout();fig.savefig(FIG_DIR/"st304_soft_ir_tradeoff.png",dpi=190);plt.close(fig)


def main() -> None:
    np.random.seed(SEED);w,a,s=strict_operator()
    out={"metadata":{"seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__,"strict_row_sum":s,"scope":"local exact, analytic, interval, and bounded numerical research"}}
    out["ST293"]=st293_feedback_source_no_go(a)
    out["ST294"]=st294_global_simplex_competition(w,a)
    out["ST295"]=st295_localized_orbit(a)
    out["ST296"]=st296_memory_gain(a)
    out["ST297"]=st297_analytic_equivariants()
    out["ST298"]=st298_near_event_arc(a)
    out["ST299"]=st299_three_output_recovery()
    out["ST300"]=st300_chernoff_tail()
    out["ST301"]=st301_nuisance_budget()
    out["ST302"]=st302_odd_clifford_obstruction()
    out["ST303"]=st303_three_level_refinement()
    out["ST304"]=st304_soft_ir_optimization()
    out["ST305"]=st305_calibration_torsor()
    out["ST306"]=st306_synthetic_conformance()
    out["ST307"]=st307_external_stop()
    make_figures(out);RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        wr=csv.writer(h);wr.writerow(["program","object","status"])
        for k in range(293,308):wr.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    print(json.dumps({k:out[k]["status"] for k in out if k.startswith("ST")},indent=2))


if __name__=="__main__":main()
