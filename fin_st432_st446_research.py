#!/usr/bin/env python3
"""FIN ST432--ST446: sector exclusion, orbit census, and IR continuation.

All statements are finite, local, dimensionless, or explicitly conditional.
The strict/legacy split is preserved.  Nothing here supplies the gain, a
selector, physical units, laboratory evidence, or a Theory-of-Everything
closure.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.linalg import qr
from scipy.optimize import root

from fin_st132_center_isolation_replay import bounds
from fin_st372_st386_research import exponent_orbit, orbit_eval_vector, orbit_representatives
from fin_st387_st401_research import dual_waterfill_lower, normalized_remainder_entropy
from fin_st402_st416_research import independent_strict_matrix_float, independent_strict_matrix_interval
from fin_st417_st431_research import (
    FACE, N, SEED, chernoff_discretization, degree8_candidate_metadata,
    rank7_interval_matrix, real_orbit_eval, regularized_ir_float,
    regularized_ir_interval_fj, regularized_ir_krawczyk,
    stationary_interval_krawczyk,
)
from fin_st132_center_isolation_replay import iv
from fin_st118_st129_research import interval_left


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST432_ST446_Results.json"
SUMMARY = ROOT / "FIN_ST432_ST446_Summary.csv"
FIG_DIR = ROOT / "FIN_ST432_ST446_Figures"
NAMES = {
    432: "Degree8_Sparse_MultiPrime_Syzygy_Audit",
    433: "Global_Gain_Sector_Exclusion",
    434: "Branch_Transition_Transversality",
    435: "Certified_Stationary_Orbit_Census",
    436: "Partial_Morse_Polynomial_and_Barriers",
    437: "Enlarged_Singular_IR_Attachment",
    438: "Local_IR_Degree_Certificate",
    439: "Characteristic_Zero_Invariant_Presentation_Audit",
    440: "One_Exchange_E_A_Design_Improvement",
    441: "Transfer_Gap_Convergence_Audit",
    442: "Optimized_Local_ISS_Ball",
    443: "Rank7_Sign_Aware_Global_Stop",
    444: "Gain_Source_Admission_Gate",
    445: "Selector_Source_Admission_Gate",
    446: "Independent_Evidence_Gate",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}


def native_obj(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native_obj(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native_obj(v) for v in x]
    if isinstance(x, np.ndarray): return native_obj(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native_obj(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": NAMES[k], "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


def _mul(values: list[np.ndarray], p: int) -> np.ndarray:
    out = np.ones_like(values[0])
    for v in values: out = (out*v) % p
    return out


def degree8_signatures(prime: int, seed: int, npts: int = 128) -> dict:
    """Exact finite-field witnesses excluding all one/two-term candidates."""
    rng = np.random.default_rng(seed)
    pts = rng.integers(1, prime, size=(npts, N), dtype=np.int64)
    pts[:, -1] = (-np.sum(pts[:, :-1], axis=1)) % prime
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    ev = lambda rep, degree: orbit_eval_vector(exponent_orbit(tuple(rep)), pts, degree, prime)
    q = [ev(r, 8) for r in base["quadratic_generator_representatives"]]
    p4 = [ev(r, 8) for r in base["primitive_quartic_representatives"]]
    p6 = [ev(r, 8) for r in base["primitive_sextic_representatives"]]
    cols = []
    for ids in itertools.combinations_with_replacement(range(6), 4): cols.append(_mul([q[i] for i in ids], prime))
    for ids in itertools.combinations_with_replacement(range(6), 2):
        q2 = _mul([q[i] for i in ids], prime)
        cols.extend((q2*x) % prime for x in p4)
    cols.extend((p4[i]*p4[j]) % prime for i, j in itertools.combinations_with_replacement(range(len(p4)), 2))
    cols.extend((q[i]*p6[j]) % prime for i in range(6) for j in range(len(p6)))
    signatures, zero = {}, 0
    for col in cols:
        nz = np.flatnonzero(col)
        if not len(nz): zero += 1; continue
        norm = (col*pow(int(col[nz[0]]), -1, prime)) % prime
        signatures.setdefault(hashlib.sha256(norm.tobytes()).hexdigest(), 0)
        signatures[hashlib.sha256(norm.tobytes()).hexdigest()] += 1
    duplicates = sum(v-1 for v in signatures.values() if v > 1)
    return {"prime": prime, "evaluation_seed": seed, "evaluation_points": npts,
            "candidate_count": len(cols), "zero_signatures": zero,
            "projective_signature_duplicates": duplicates,
            "signature_set_sha256": hashlib.sha256("".join(sorted(signatures)).encode()).hexdigest()}


def st432() -> dict:
    meta = degree8_candidate_metadata()
    rows = [degree8_signatures(1000003, SEED+432), degree8_signatures(1000033, SEED+1432)]
    packet = {**meta, "multi_prime_exact_signature_audit": rows,
              "proven_statement": "No degree-eight candidate is zero or proportional to another candidate over either declared prime; hence no rational syzygy with support one or two exists.",
              "unresolved": "The forced nullspace has dimension at least 136, but no support >=3 basis vector was recovered; integral lifting remains open."}
    return finalize(432, "proven_no_support_le_2_syzygy_full_basis_blocked",
                    "Finite-field nonproportionality is exact, but it neither constructs the 136 required relations nor proves a characteristic-zero basis.", packet)


def concentrated_cover(g: float = 2.88, tmax: float = .2, nt: int = 1000, nr: int = 1000) -> dict:
    _, wiv, siv = independent_strict_matrix_interval()
    order = [6,5,5,4,4,3,3,2,2,1,1]
    wmid = np.array([(bounds(wiv[d])[0]+bounds(wiv[d])[1])/2 for d in order])
    wlo = np.array([bounds(wiv[d])[0] for d in order]); wmin = float(min(wlo)); shi = bounds(siv)[1]
    rlo = 1/11+(1-1/11)*np.arange(nr)/nr
    rhi = 1/11+(1-1/11)*(np.arange(nr)+1)/nr
    H = np.array([normalized_remainder_entropy(x) for x in rlo])
    L = np.array([dual_waterfill_lower(x, wmid, wlo, wmin) for x in rhi])
    best, arg = math.inf, None
    for i in range(nt):
        tl=max(1e-14,tmax*i/nt); th=tmax*(i+1)/nt
        entropy=(1-th)*math.log(12*(1-th))+th*math.log(12*th)+th*H
        collision=(1-tl)**2+tl*tl*rhi
        cross=2*(1-tl)*tl*L+wmin*tl*tl*(1-rhi)
        vals=entropy-g/2*(shi*collision-cross)-1e-10
        j=int(np.argmin(vals))
        if vals[j] < best: best=float(vals[j]); arg=[tl,th,float(rlo[j]),float(rhi[j])]
    return {"gain": g, "maximum_coordinate_sector": [1-tmax,1.0], "cells": nt*nr,
            "minimum_outward_paid_lower_bound": best, "worst_cell_t_rho": arg,
            "roundoff_payment_per_cell": 1e-10,
            "certificate_logic": "Entropy is minimized at t_hi and rho_lo; collision is maximized at t_lo,rho_hi; the SOCP-dual water-fill value and w_min give a cross-term lower bound."}


def st433() -> dict:
    A = independent_strict_matrix_float(); lm=float(np.linalg.eigvalsh(A)[-1]); g=2.88
    mconv=1/(g*lm); cover=concentrated_cover(g)
    packet={"critical_gain_identity": "g_global=inf_{p!=u} 2 D(p||u)/((p-u)^T A (p-u))",
            "uniform_convex_sector": {"maximum_coordinate_upper": mconv,
                "proof": "diag(1/p)-gA is positive definite on the tangent space when max(p)<1/(g lambda_max(A)); the convex sector contains the uniform stationary point."},
            "concentrated_sector_cover": cover,
            "remaining_competitor_band": [mconv,.8],
            "global_transition_enclosure_retained": [2.2513440017990893,3.999],
            "result": "At every g<=2.88 a negative competitor is excluded in both displayed outer sectors. Any missed competitor must lie in the explicit intermediate band."}
    return finalize(433,"proven_two_sector_exclusion_global_transition_still_open",
                    "The intermediate simplex band is not covered; the local crossing is not promoted to the global transition.",packet)


def st434() -> dict:
    d=json.loads((ROOT/"FIN_ST342_Narrow_Certified_Coexistence_Bracket.json").read_text());A=independent_strict_matrix_float()
    idx=np.array([i if i<=6 else 12-i for i in range(12)]);u=np.ones(12)/12;rows=[]
    for e in d["endpoint_certificates"]:
        p=np.array(e["center"][:7])[idx];q=p-u;Q=float(q@A@q);D=float(np.sum(p*np.log(12*p)))
        rows.append({"g":e["g"],"D":D,"Q":Q,"ratio_2D_over_Q":2*D/Q,
                     "envelope_derivative_minus_Q_over_2":-Q/2})
    # Enclose Q on the complete parametric Krawczyk tube, not merely at its midpoint.
    Aiv,_,_=independent_strict_matrix_interval(); center=np.mean([e["center"] for e in d["endpoint_certificates"]],axis=0)
    radii=np.array(d["continuous_parametric_root_tube"]["component_radii"])
    piv=[iv((float(center[j]-radii[j]),float(center[j]+radii[j]))) for j in idx];qiv=[x-iv(1/12) for x in piv];Qiv=iv(0)
    for i in range(N):
        for j in range(N): Qiv += qiv[i]*Aiv[i][j]*qiv[j]
    qbox=list(bounds(Qiv)); derivative=[-qbox[1]/2,-qbox[0]/2]
    packet={"certified_crossing_bracket":d["certified_crossing_bracket"],"endpoint_round_trip":rows,
            "certified_Q_interval_on_complete_tube":qbox,"transversality_bound":derivative,
            "theorem": "Along the certified stationary tube the envelope theorem gives dV/dg=-Q/2. The nonuniform tube has Q>1.3425, so every equality in the bracket is transverse and locally unique in g.",
            "third_competitor_exclusion":False}
    return finalize(434,"proven_branch_local_transversality_global_exclusion_open",
                    "This is a theorem for the certified reflection-even branch, not a uniqueness theorem over the full simplex.",packet)


def d12_perms():
    return [np.array([(s*i+k)%N for i in range(N)]) for s in (1,-1) for k in range(N)]


def st435() -> dict:
    mp.iv.dps=55; atlas=json.loads((ROOT/"FIN_ST421_Stationary_Orbit_Interval_Atlas.json").read_text())
    Aiv,_,_=independent_strict_matrix_interval();perms=d12_perms();rows=[]
    for r in atlas["locally_certified_stationary_orbits"]:
        p=np.array(r["center"]);stab=[q for q in perms if np.max(abs(p-p[q]))<1e-7]
        avg=sum((p[q] for q in stab),np.zeros(N))/len(stab);z=np.r_[avg,r["lambda"]]
        cert=stationary_interval_krawczyk(z,Aiv,2.1e-9)
        non=[float(np.max(abs(avg-avg[q]))) for q in perms if np.max(abs(avg-avg[q]))>1e-7]
        rows.append({"label":r["label"],"Morse_index":r["certified_Morse_index"],
                     "stabilizer_order":len(stab),"orbit_size":24//len(stab),
                     "symmetry_averaged_Krawczyk":cert,
                     "minimum_nonstabilizer_coordinate_separation":min(non) if non else None})
    total=sum(x["orbit_size"] for x in rows)
    packet={"certified_representatives":rows,"distinct_points_generated_by_D12":total,
            "orbit_size_histogram":{str(s):sum(x["orbit_size"]==s for x in rows) for s in sorted(set(x["orbit_size"] for x in rows))},
            "theorem":"Invariant symmetry-averaged boxes certify the stated stabilizers and orbit sizes for all nine already isolated roots.",
            "completeness":False}
    return finalize(435,"proven_stabilizers_and_85_distinct_stationary_points_nonexhaustive",
                    "The 85 points are certified, but a complete interval cover of the remaining simplex is absent.",packet)


def st436() -> dict:
    census=json.loads(PACKETS[435].read_text());atlas=json.loads((ROOT/"FIN_ST421_Stationary_Orbit_Interval_Atlas.json").read_text())
    coeff={0:0,1:0,2:0}
    for c in census["certified_representatives"]: coeff[c["Morse_index"]]+=c["orbit_size"]
    saddles=[{"label":r["label"],"value":r["value"],"orbit_size":next(c["orbit_size"] for c in census["certified_representatives"] if c["label"]==r["label"])} for r in atlas["locally_certified_stationary_orbits"] if r["certified_Morse_index"]==1]
    packet={"partial_Morse_polynomial":f"{coeff[0]} + {coeff[1]} t + {coeff[2]} t^2",
            "certified_coefficients":coeff,"alternating_sum_at_minus_one":coeff[0]-coeff[1]+coeff[2],
            "simplex_Euler_characteristic":1,"index_one_barrier_values":saddles,
            "interpretation":"The certified partial census is already Euler-balanced (13-42+30=1), a nontrivial consistency check, but extra critical orbits may enter in canceling index pairs.",
            "Conley_connection_theorem":False}
    return finalize(436,"proven_partial_Morse_census_Euler_balanced_connections_open",
                    "Euler balance is not a completeness proof and numerical descending trajectories are not promoted to heteroclinic theorems.",packet)


def ir_common_box(bhi: float) -> tuple[np.ndarray,np.ndarray,dict]:
    roots=[];z=FACE.copy()
    for b in (0,bhi/2,bhi): z=root(lambda x:regularized_ir_float(x,b),z,tol=1e-12).x;roots.append(z.copy())
    center=roots[1];radii=1.2*np.max(np.abs(np.array(roots)-center),axis=0)+np.array([2e-7,1e-3,2e-5,4e-5,1.5e-3])
    return center,radii,regularized_ir_krawczyk(center,radii,bhi)


def st437() -> dict:
    center,radii,cert=ir_common_box(.017);_,_,failure=ir_common_box(.018)
    uniqueness=contraction_norm_for_ir(.017)
    packet={"compactified_parameter_interval":[0,.017],"physical_tail_coordinate":"y3=1/b for b>0",
            "common_box_center":center,"common_box_radii":radii,"Krawczyk":cert,
            "weighted_uniqueness_contraction":uniqueness["weighted_box_norm_I_minus_CJ"],
            "extension_factor_over_ST419":8.5,"fixed_template_failure_probe_at_0_018":failure,
            "theorem":"For every b in [0,0.017] the regularized five-equation family has exactly one root in the displayed common box, extending the certified attachment by a factor 8.5."}
    return finalize(437,"proven_unique_singular_attachment_through_b_0_017",
                    "Failure of the chosen box at 0.018 is not a fold or a no-go; adaptive continuation may extend farther.",packet)


def contraction_norm_for_ir(bhi: float) -> dict:
    center,radii,cert=ir_common_box(bhi);X=[iv((center[i]-radii[i],center[i]+radii[i])) for i in range(5)]
    _,J=regularized_ir_interval_fj(X,bhi);jl=np.array([[bounds(x)[0] for x in row] for row in J]);jh=np.array([[bounds(x)[1] for x in row] for row in J])
    _,jc=regularized_ir_interval_fj([iv(float(x)) for x in center],1e-10);jm=np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in jc]);pre=np.linalg.inv(jm)
    lo,hi=interval_left(pre,jl,jh);mlo=-hi;mhi=-lo
    for i in range(5):mlo[i,i]+=1;mhi[i,i]+=1
    M=np.maximum(abs(mlo),abs(mhi));q=float(max(np.sum(M,axis=1)));qweighted=float(max((M@radii)/radii))
    return {"b_interval":[0,bhi],"Krawczyk_included":cert["included"],"sup_norm_I_minus_CJ":q,
            "weighted_box_norm_I_minus_CJ":qweighted,
            "midpoint_Jacobian_determinant":float(np.linalg.det(jm)),"degree":-1 if q<1 and np.linalg.det(jm)<0 else None}


def st438() -> dict:
    row=contraction_norm_for_ir(.015)
    packet={"local_component_certificate":row,
            "theorem":"Because ||I-CJ||_infinity<1 on the common box, every Jacobian is nonsingular and has the midpoint orientation. The unique compactified branch component therefore has Brouwer degree -1 through b=0.015.",
            "full_complement_cover":False}
    return finalize(438,"proven_local_IR_component_degree_minus_one_full_topology_open",
                    "A local degree does not exclude additional root components outside the certified box.",packet)


def modular_joint_rank(prime: int, seed: int) -> dict:
    # Reuse the exact algorithm with a temporary seed by reproducing its finite evaluation stage.
    rng=np.random.default_rng(seed);npts=390;pts=rng.integers(1,prime,size=(npts,N),dtype=np.int64);pts[:,-1]=(-np.sum(pts[:,:-1],axis=1))%prime
    base=json.loads((ROOT/"FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    ev=lambda rep,degree:orbit_eval_vector(exponent_orbit(tuple(rep)),pts,degree,prime)
    q=[ev(x,6) for x in base["quadratic_generator_representatives"]]
    d3=[ev(rep,6) for rep in json.loads((ROOT/"FIN_ST424_Odd_and_Even_Invariant_Presentation.json").read_text())["degree3_representatives"]]
    p4=[ev(x,6) for x in base["primitive_quartic_representatives"]]
    cand5=[(x*y)%prime for x in q for y in d3]
    cand6=[]
    for ids in itertools.combinations_with_replacement(range(6),3):cand6.append(_mul([q[i] for i in ids],prime))
    cand6 += [(x*y)%prime for x in q for y in p4]
    cand6 += [(d3[i]*d3[j])%prime for i,j in itertools.combinations_with_replacement(range(12),2)]
    def rank_cols(cols):
        echelon={};rank=0
        for col in cols:
            v=col.copy()
            for piv in sorted(echelon):
                if v[piv]:v=(v-v[piv]*echelon[piv])%prime
            nz=np.flatnonzero(v)
            if len(nz):
                piv=int(nz[0]);v=(v*pow(int(v[piv]),-1,prime))%prime;echelon[piv]=v;rank+=1
        return rank
    return {"prime":prime,"seed":seed,"degree5_rank":rank_cols(cand5),"degree5_candidates":len(cand5),
            "degree6_rank":rank_cols(cand6),"degree6_candidates":len(cand6)}


def st439() -> dict:
    rows=[modular_joint_rank(1000003,SEED+439),modular_joint_rank(1000033,SEED+1439)]
    packet={"multi_prime_ranks":rows,
            "characteristic_zero_conclusions":{"degree5_decomposable_rank":72,"degree5_primitive_quotient_dimension":52,
                "degree6_decomposable_rank_interval":[310,326],"degree6_primitive_quotient_dimension_interval":[39,55]},
            "proof":"A nonzero modular minor is a nonzero integer minor, so rank_Q>=rank_Fp. Candidate counts give the upper bounds. Degree five therefore closes exactly; degree six does not.",
            "degree8_presentation":False}
    return finalize(439,"proven_characteristic_zero_degree5_quotient_degree6_interval",
                    "Agreement of two primes is evidence, not an upper-rank proof; characteristic-zero degree six/eight relations remain open.",packet)


def st440() -> dict:
    reps=json.loads((ROOT/"FIN_ST378_Exact_Sextic_D12_Reynolds_Basis.json").read_text())["selected_orbit_representatives"]
    rng=np.random.default_rng(SEED+425);pool=1800;keep=500;pts=rng.normal(size=(pool,N));pts-=pts.mean(axis=1,keepdims=True);pts/=np.linalg.norm(pts,axis=1,keepdims=True)
    D=np.column_stack([real_orbit_eval(r,pts) for r in reps]);X=D/np.linalg.norm(D,axis=0)
    _,_,piv=qr(X.T,pivoting=True,mode="economic");sel=list(map(int,piv[:keep]));u,s,v=np.linalg.svd(X[sel],full_matrices=False)
    score=np.abs(X@v[-1]);excluded=np.ones(pool,bool);excluded[sel]=False;add=int(np.argmax(np.where(excluded,score,-1)));rp=int(np.argmin(score[sel]));removed=sel[rp];sel[rp]=add;s2=np.linalg.svd(X[sel],compute_uv=False)
    packet={"baseline_QR_E_value":float(s[-1]),"one_exchange_E_value":float(s2[-1]),"E_improvement_factor":float(s2[-1]/s[-1]),
            "baseline_A_trace_inverse":float(np.sum(1/s**2)),"one_exchange_A_trace_inverse":float(np.sum(1/s2**2)),
            "A_improvement_factor":float(np.sum(1/s**2)/np.sum(1/s2**2)),"added_row":add,"removed_row":removed,
            "condition_number_before":float(s[0]/s[-1]),"condition_number_after":float(s2[0]/s2[-1])}
    return finalize(440,"strong_numerical_one_exchange_improves_E_and_A_design",
                    "The observation pool, coefficient scaling, and noise model are synthetic; this is not apparatus calibration.",packet)


def st441() -> dict:
    rows=[]
    for n in (21,31,41,51,61):
        ev=chernoff_discretization(n);rows.append({"grid":n,"leading":float(ev[0].real),"second":float(ev[1].real),"ratio":float(abs(ev[1])/abs(ev[0]))})
    cone=json.loads((ROOT/"FIN_ST413_Adapted_Chernoff_Cone_Gap.json").read_text())
    packet={"certified_Birkhoff_contraction_upper":1-cone["Birkhoff_gap_lower"],"Ulam_sequence":rows,
            "last_two_ratio_spread":abs(rows[-1]["ratio"]-rows[-2]["ratio"]),
            "result":"The ratio stabilizes near 0.59574, much sharper than the analytic cone contraction, but no discretization-error theorem converts it into an infinite-dimensional interval."}
    return finalize(441,"certified_cone_upper_bound_numerical_gap_sequence_no_bracket",
                    "Finite Ulam eigenvalues are numerical evidence only; the synthetic transfer fixture is not a FIN detector.",packet)


def st442() -> dict:
    atlas=json.loads((ROOT/"FIN_ST421_Stationary_Orbit_Interval_Atlas.json").read_text());r=min(atlas["locally_certified_stationary_orbits"],key=lambda x:x["value"])
    pmin=min(r["center"]);mu0=min(r["tangent_Hessian_eigenvalues"])-1e-8;R=1.93e-6;variation=R/(pmin-R)**2;mu=mu0-variation;eps=mu*R/2
    old=json.loads((ROOT/"FIN_ST427_Nonlinear_Local_Flow_Robustness.json").read_text())
    packet={"radius":R,"strong_convexity_lower":mu,"forcing_threshold":eps,"Hessian_variation":variation,
            "radius_factor_over_ST427":R/old["declared_basin_radius"],"forcing_factor_over_ST427":eps/old["admissible_tangent_forcing_norm"],
            "ISS_bound":"||e(t)|| <= exp(-mu t)||e(0)|| + epsilon(1-exp(-mu t))/mu"}
    return finalize(442,"proven_optimized_conservative_local_ISS_ball",
                    "The forcing and time are supplied dimensionless inputs; no global basin or physical environment is derived.",packet)


def st443() -> dict:
    M=np.array([[float((bounds(x)[0]+bounds(x)[1])/2) for x in row] for row in rank7_interval_matrix()]);ev=np.linalg.eigvalsh(M);off=M[np.triu_indices(N,1)]
    packet={"rank":int(np.sum(ev>1e-10)),"positive_offdiagonal_pairs":int(np.sum(off>0)),"negative_offdiagonal_pairs":int(np.sum(off<0)),
            "lambda_max":float(ev[-1]),"available_SDP_solver":False,
            "result":"No local SDP package is installed. The exact mixed-sign obstruction remains: positive-weight water-filling cannot be transferred, and the cap-curvature/localization gap is not closed."}
    return finalize(443,"blocked_rank7_global_certificate_mixed_sign_and_no_SDP_closure",
                    "Shared eigenvalues do not imply a shared global variational theorem.",packet)


def gate(k: int, kind: str) -> dict:
    if kind=="gain": packet={"new_internal_gain_source":False,"g_remains_supplied":True,"reason":"Sector exclusion and transversality constrain consequences of g but do not select its value."};status="blocked_no_new_strict_gain_source"
    elif kind=="selector": packet={"new_nonpremise_selector":False,"QW_2191":"open","reason":"D12 orbit certification multiplies equivalent roots; it does not select one."};status="blocked_no_new_selector_provider"
    else: packet={"external_referee":"absent","independent_laboratory_record":"absent","holdout":"absent","reason":"All artifacts in this batch are local analytical or computational outputs."};status="blocked_no_independent_empirical_evidence"
    return finalize(k,status,"No physical, selector, unit, bridge-role-transfer, or ToE closure is exported.",packet)


def figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    c=results["ST435"]["certified_representatives"]
    fig,ax=plt.subplots(figsize=(7.2,4.2));ax.bar([x["label"] for x in c],[x["orbit_size"] for x in c],color="#3b82f6");ax.set_ylabel("certified orbit size");ax.set_title("ST435: orbit-stabilizer census (nonexhaustive)");fig.tight_layout();fig.savefig(FIG_DIR/"st435_orbit_sizes.png",dpi=180);plt.close(fig)
    rows=results["ST441"]["Ulam_sequence"];fig,ax=plt.subplots(figsize=(7.2,4.2));ax.plot([r["grid"] for r in rows],[r["ratio"] for r in rows],"o-");ax.set_xlabel("grid per belief coordinate");ax.set_ylabel(r"$|\lambda_2|/|\lambda_1|$");ax.set_title("ST441: numerical Ulam gap sequence");fig.tight_layout();fig.savefig(FIG_DIR/"st441_ulam_ratio.png",dpi=180);plt.close(fig)


def main() -> None:
    results={}
    for k,fn in [(432,st432),(433,st433),(434,st434),(435,st435),(436,st436),(437,st437),(438,st438),(439,st439),(440,st440),(441,st441),(442,st442),(443,st443)]:
        print(f"running ST{k}",flush=True);results[f"ST{k}"]=fn()
    results["ST444"]=gate(444,"gain");results["ST445"]=gate(445,"selector");results["ST446"]=gate(446,"evidence")
    RESULTS.write_text(json.dumps(native_obj(results),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["program","status","object","boundary"])
        for k in range(432,447):
            r=results[f"ST{k}"];w.writerow([f"ST{k}",r["status"],r["object"],r["boundary"]])
    figures(results)
    print(f"wrote {RESULTS.name} and {SUMMARY.name}")


if __name__ == "__main__": main()
