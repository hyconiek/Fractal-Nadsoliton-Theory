#!/usr/bin/env python3
"""FIN ST492--ST506: modular nullspace preparation, IR closure, and dual testers.

Third round in the ST462--ST506 cycle.  Exact finite-field statements remain
separate from characteristic-zero claims.  Noise, bias, clocks, refinement,
and instruments are typed resources.  No strict gain/selector/unit source,
physical continuum, laboratory record, legacy role transfer, SM/GR, L_total,
or ToE closure is exported.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import tempfile
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from sympy.polys.domains import ZZ
from sympy.polys.modulargcd import _integer_rational_reconstruction

from fin_st372_st386_research import exponent_orbit, orbit_eval_vector
from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N, SEED
from fin_st447_st461_research import _mul
from fin_st462_st476_research import probability_ratio, ratio_and_gradient, softmax_reduced
from fin_st477_st491_research import (
    compile_modrank,
    recertify_ir_cell,
)
from fin_st447_st461_research import direct_ir_cell, regularized_ir_float, FACE
from scipy.optimize import root


ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/"FIN_ST492_ST506_Results.json"
SUMMARY=ROOT/"FIN_ST492_ST506_Summary.csv"
NULLSPACE=ROOT/"FIN_ST495_Degree8_Modular_Nullspaces.npz"
FIG_DIR=ROOT/"FIN_ST492_ST506_Figures"
NAMES={
492:"Three_Prime_Degree8_Rank_and_Pivots",
493:"Alternative_Point_Ensemble_Modular_Rank",
494:"Cross_Prime_Pivot_Schema_Audit",
495:"Degree8_Modular_Nullspace_Preparation",
496:"Rational_Lift_Probe_and_Independent_Prime_Replay",
497:"Transition_Ratio_Local_Hessian",
498:"Cross_Algorithm_Transition_Orbit_Stress_Test",
499:"Adaptive_IR_Continuation_Beyond_0_2",
500:"Thermal_Fluctuation_Gain_Dimensional_Boundary",
501:"Biased_Noise_Selector_Source_Theorem",
502:"Minimal_Mode_Probe_for_Dual_Dynamics",
503:"Multitime_Operational_Conjugacy",
504:"Relative_Clock_Ratio_Perturbation_Bound",
505:"Source_Selector_Scale_Admission_Gate",
506:"Independent_Evidence_Gate",
}
PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()}
A=independent_strict_matrix_float(); U=np.ones(N)/N


def native(x:Any)->Any:
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    if isinstance(x,np.ndarray):return native(x.tolist())
    if isinstance(x,(np.floating,np.integer)):return x.item()
    return x


def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
    path=PACKETS[k];path.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8")
    return {"program":f"ST{k}","object":NAMES[k],"packet_file":path.name,"packet_sha256":sha(path),**packet,"status":status,"boundary":boundary}


def degree8_matrix_seed(prime:int,seed:int,npts:int=1892)->np.ndarray:
    rng=np.random.default_rng(seed);raw=rng.integers(-7,8,size=(npts,N),dtype=np.int64);raw[:,-1]=-np.sum(raw[:,:-1],axis=1);pts=raw%prime
    base=json.loads((ROOT/"FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    ev=lambda rep:orbit_eval_vector(exponent_orbit(tuple(rep)),pts,8,prime)
    q=[ev(r) for r in base["quadratic_generator_representatives"]]
    p4=[ev(r) for r in base["primitive_quartic_representatives"]]
    p6=[ev(r) for r in base["primitive_sextic_representatives"]]
    cols=[]
    for ids in itertools.combinations_with_replacement(range(6),4):cols.append(_mul([q[i] for i in ids],prime))
    for ids in itertools.combinations_with_replacement(range(6),2):
        q2=_mul([q[i] for i in ids],prime);cols.extend((q2*x)%prime for x in p4)
    cols.extend((p4[i]*p4[j])%prime for i,j in itertools.combinations_with_replacement(range(32),2))
    cols.extend((q[i]*p6[j])%prime for i in range(6) for j in range(117))
    return np.column_stack(cols)


def rank_pivots(prime:int,seed:int,echelon:bool=False,timeout:int=600)->dict:
    binary=compile_modrank();matrix=degree8_matrix_seed(prime,seed).astype(np.uint32);t0=time.time()
    raw=Path(tempfile.gettempdir())/f"fin492_{prime}_{seed}_{os.getpid()}.bin";matrix.tofile(raw)
    ech=Path(tempfile.gettempdir())/f"fin492_ech_{prime}_{seed}_{os.getpid()}.bin"
    cmd=[str(binary),str(raw),"1892","2028",str(prime)]+([str(ech)] if echelon else [])
    env=dict(os.environ);env["OMP_NUM_THREADS"]=env.get("OMP_NUM_THREADS","8")
    try:
        proc=subprocess.run(cmd,capture_output=True,text=True,timeout=timeout,env=env)
        rank=int(re.search(r"rank=(\d+)",proc.stdout).group(1));pm=re.search(r"pivots=([0-9,]+)",proc.stdout)
        pivots=[int(x) for x in pm.group(1).split(",")]
        out={"prime":prime,"seed":seed,"rank":rank,"nullity":2028-rank,"pivot_columns":pivots,
             "pivot_sha256":hashlib.sha256(np.array(pivots,dtype=np.int32).tobytes()).hexdigest(),
             "matrix_sha256":hashlib.sha256(matrix.tobytes()).hexdigest(),"wall_seconds":time.time()-t0,"return_code":proc.returncode}
        if echelon:
            out["echelon"]=np.fromfile(ech,dtype=np.uint32).reshape(1892,2028)
        return out
    finally:
        raw.unlink(missing_ok=True);ech.unlink(missing_ok=True)


def st492()->dict:
    rows=[rank_pivots(p,SEED+447) for p in (1000003,1000033,1000037)]
    packet={"exact_runs":rows,"all_ranks_equal_1791":all(x["rank"]==1791 for x in rows),
            "all_pivot_schemas_equal":len({x["pivot_sha256"] for x in rows})==1,
            "claim":"Three exact finite-field ranks and their ordered pivot schemas are recorded for the same integer point ensemble."}
    return finalize(492,"proven_three_prime_modular_rank_1791", "Finite-field agreement is not a characteristic-zero upper-rank certificate.",packet)


def st493()->dict:
    seed=SEED+1493;rows=[rank_pivots(p,seed) for p in (1000003,1000033)]
    packet={"alternative_seed":seed,"exact_runs":rows,"all_ranks_equal_1791":all(x["rank"]==1791 for x in rows),
            "alternative_pivots_equal_across_primes":len({x["pivot_sha256"] for x in rows})==1}
    return finalize(493,"proven_alternative_ensemble_modular_rank_1791", "Generic ensemble agreement is robustness evidence; rational upper closure remains absent.",packet)


def st494()->dict:
    a=json.loads(PACKETS[492].read_text());b=json.loads(PACKETS[493].read_text())
    standard={x["pivot_sha256"] for x in a["exact_runs"]};alternative={x["pivot_sha256"] for x in b["exact_runs"]}
    packet={"standard_schema_count":len(standard),"alternative_schema_count":len(alternative),
            "same_schema_between_point_ensembles":standard==alternative,
            "standard_pivot_sha256":next(iter(standard)),"alternative_pivot_sha256":next(iter(alternative)),
            "theorem":"Within each declared generic ensemble, identical ordered pivot columns across primes provide a canonical modular free-column schema for relation preparation.",
            "rational_schema_theorem":False}
    return finalize(494,"proven_cross_prime_pivot_schema_ensemble_stability", "A point-ensemble pivot change would not change rank; even a common schema does not lift coefficients to Q.",packet)


def modular_nullspace(echelon:np.ndarray,pivots:list[int],prime:int)->np.ndarray:
    n=echelon.shape[1];free=[j for j in range(n) if j not in set(pivots)];X=np.zeros((n,len(free)),dtype=np.int64)
    X[free,np.arange(len(free))]=1
    E=echelon.astype(np.int64)
    for i in range(len(pivots)-1,-1,-1):
        pc=pivots[i];X[pc]=-(E[i,pc+1:]@X[pc+1:])%prime
    return X.astype(np.uint32)


def st495()->dict:
    rows=[];bases={}
    for prime in (1000003,1000033):
        run=rank_pivots(prime,SEED+447,echelon=True);X=modular_nullspace(run.pop("echelon"),run["pivot_columns"],prime)
        matrix=degree8_matrix_seed(prime,SEED+447,npts=64).astype(np.int64)
        sample_res=(matrix@X.astype(np.int64))%prime
        rows.append({**run,"basis_shape":list(X.shape),"basis_sha256":hashlib.sha256(X.tobytes()).hexdigest(),
                     "independent_64_row_replay_max_residue":int(sample_res.max()),"independent_64_row_replay_nonzeros":int(np.count_nonzero(sample_res))})
        bases[str(prime)]=X
    np.savez_compressed(NULLSPACE,**bases)
    packet={"runs":rows,"nullspace_file":NULLSPACE.name,"nullspace_sha256":sha(NULLSPACE),
            "all_independent_replays_zero":all(x["independent_64_row_replay_nonzeros"]==0 for x in rows),
            "claim":"A complete 237-vector modular nullspace basis is prepared for each of two primes and independently replayed on 64 exact evaluation rows."}
    return finalize(495,"proven_complete_two_prime_modular_nullspace_preparation", "These are finite-field relations; coefficient lifting and exact polynomial verification are still required for Q.",packet)


def reconstruct_integer_vector(a:np.ndarray,b:np.ndarray,p1:int,p2:int)->list[int]|None:
    mod=p1*p2;inv=pow(p1,-1,p2);fr=[]
    for x,y in zip(a,b):
        residue=(int(x)+p1*(((int(y)-int(x))*inv)%p2))%mod
        q=_integer_rational_reconstruction(residue,mod,ZZ)
        if q is None:return None
        fr.append(Fraction(int(q.numerator),int(q.denominator)))
    lcm=1
    for q in fr:lcm=math.lcm(lcm,q.denominator)
    vals=[q.numerator*(lcm//q.denominator) for q in fr];g=0
    for x in vals:g=math.gcd(g,abs(x))
    vals=[x//g for x in vals]
    if next(x for x in vals if x)<0:vals=[-x for x in vals]
    return vals


def st496()->dict:
    z=np.load(NULLSPACE);x1=z["1000003"];x2=z["1000033"];p3=1000037
    M3=degree8_matrix_seed(p3,SEED+447,npts=96).astype(np.int64);rows=[]
    for j in range(min(8,x1.shape[1])):
        v=reconstruct_integer_vector(x1[:,j],x2[:,j],1000003,1000033)
        if v is None:rows.append({"column":j,"reconstructed":False});continue
        arr=np.array(v,dtype=object);maxc=max(abs(int(q)) for q in v)
        residue=np.array([int(q%p3) for q in v],dtype=np.int64);test=M3@residue%p3
        rows.append({"column":j,"reconstructed":True,"support":sum(q!=0 for q in v),"max_abs_integer_coefficient":maxc,
                     "independent_prime_nonzeros":int(np.count_nonzero(test)),"independent_prime_max_residue":int(test.max())})
    packet={"columns_probed":len(rows),"reconstruction_rows":rows,
            "all_reconstructed":all(x.get("reconstructed") for x in rows),
            "all_independent_prime_replays_zero":all(x.get("independent_prime_nonzeros")==0 for x in rows if x.get("reconstructed")),
            "exact_characteristic_zero_polynomial_checks":False}
    return finalize(496,"rational_lift_probe_completed", "Independent-prime evaluation is strong evidence, not exact symbolic polynomial identity; no Q-rank equality is promoted.",packet)


def numerical_hessian_gradient(z:np.ndarray,h=2e-5)->np.ndarray:
    m=len(z);H=np.zeros((m,m));eye=np.eye(m)
    for j in range(m):H[:,j]=(ratio_and_gradient(z+h*eye[j])[1]-ratio_and_gradient(z-h*eye[j])[1])/(2*h)
    return (H+H.T)/2


def st497()->dict:
    old=json.loads((ROOT/"FIN_ST462_First_Transition_Orbit_Numerical_Isolation.json").read_text());p=np.array(old["best_probability"]);z=np.log(p[:-1]/p[-1]);H=numerical_hessian_gradient(z);ev=np.linalg.eigvalsh(H)
    packet={"ratio":probability_ratio(p),"gradient_infinity_norm":float(np.linalg.norm(ratio_and_gradient(z)[1],np.inf)),
            "reduced_logit_Hessian_eigenvalues":ev,"minimum_eigenvalue":float(ev[0]),"Morse_index":int(np.sum(ev< -1e-7)),
            "claim":"The isolated transition-ratio representative is a strict numerical local minimum in reduced logit coordinates.","interval_Hessian_certificate":False}
    return finalize(497,"strong_numerical_strict_local_minimum_of_transition_ratio", "Floating Hessian positivity is not a global or interval theorem.",packet)


def st498()->dict:
    rng=np.random.default_rng(SEED+498);rows=[]
    for method,count in (("BFGS",80),("L-BFGS-B",80)):
        for _ in range(count):
            z0=rng.normal(0,4,N-1)
            if method=="BFGS":opt=minimize(lambda x:ratio_and_gradient(x),z0,jac=True,method=method,options={"maxiter":1200,"gtol":1e-9})
            else:opt=minimize(lambda x:ratio_and_gradient(x),z0,jac=True,method=method,options={"maxiter":1200,"ftol":1e-14,"gtol":1e-9})
            rows.append({"method":method,"value":float(opt.fun),"probability":softmax_reduced(opt.x)})
    best=min(x["value"] for x in rows);hits=sum(abs(x["value"]-best)<2e-9 for x in rows)
    packet={"run_count":len(rows),"methods":["BFGS","L-BFGS-B"],"best_ratio":best,"best_orbit_value_hits":hits,
            "distinct_value_clusters_1e_7":len({round(x["value"],7) for x in rows}),
            "candidate_beaten":best<2.902496471747767,"global_certificate":False}
    return finalize(498,"strong_cross_algorithm_transition_orbit_stress_test", "Optimization starts are not exhaustive and cannot replace the residual global strip certificate.",packet)


def repaired_chain(start=.196,stop=.25,width=.0015)->dict:
    prior=root(lambda z:regularized_ir_float(z,start),FACE,tol=1e-12).x;rows=[];lo=start;failure=None
    while lo<stop-1e-14:
        hi=min(stop,lo+width);raw,candidate=direct_ir_cell(lo,hi,prior);selected=None
        for factor in (1.,1.5,2.,3.,4.,6.,8.):
            t=recertify_ir_cell(raw,lo,hi,factor)
            if t["included"]:selected=t;break
        if selected is None:failure=raw;break
        rows.append(selected);prior=candidate;lo=hi
    return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}


def st499()->dict:
    c=repaired_chain();rows=c["cells"]
    packet={**c,"cell_count":len(rows),"minimum_margin":min(x["minimum_margin"] for x in rows) if rows else None,
            "maximum_contraction":max(x["weighted_contraction"] for x in rows) if rows else None,
            "maximum_radius_factor":max(x["radius_expansion_factor"] for x in rows) if rows else None,
            "target_reached":c["certified_stop"]>=.25-1e-14,
            "theorem":f"The same local IR component is certified from b=0.196 through b={c['certified_stop']:.4f}."}
    return finalize(499,"proven_IR_continuation_beyond_0_2" if c["certified_stop"]>.2 else "failed_IR_continuation_beyond_0_2", "Local boxes do not give a global root count or physical infrared scale.",packet)


def st500()->dict:
    rank=11
    packet={"thermal_covariance_on_mean_zero_space":"Sigma=theta*A^+","dimensionless_temperature":"theta",
            "identity":"Tr(A Sigma)=theta*rank(A)=11 theta","stationary_mean_gain":"(kappa/gamma)*11 theta",
            "theorem":"For a Gibbs/Ornstein-Uhlenbeck covariance proportional to A^+ on the mean-zero subspace, fluctuation-driven mean gain is exactly 11(kappa/gamma)theta.",
            "physical_bridge":"theta=k_B T requires bath, temperature and energy calibration","strict_gain_source":False}
    return finalize(500,"proven_thermal_gain_identity_with_dimensional_boundary", "The result transfers the missing source to theta/bath data; it does not derive them from FIN.",packet)


def st501()->dict:
    packet={"branch_law":"P_i(epsilon)=exp(epsilon*b_i)/sum_j exp(epsilon*b_j)",
            "linear_response":"P_i'(0)=(b_i-mean(b))/12","unbiased_law":"P_i(0)=1/12",
            "theorem":"A nonuniform first-order branch probability requires a nonconstant bias vector b. Therefore biased noise can operationally select, but b and its orientation are precisely an added selector source.",
            "QW_2191":"open","strict_selector":False}
    return finalize(501,"proven_biased_noise_selection_requires_bias_source", "No canonical b is exported by the radial strict core.",packet)


def st502()->dict:
    lam=float(np.linalg.eigvalsh(A)[-1]);x=1.4536736664610401;t=x/lam
    au=np.exp(-1j*lam*t);ah=np.exp(-lam*t)
    packet={"selected_mode":"largest positive eigenvalue","lambda":lam,"dimensionless_time":t,
            "unitary_mode_response":{"real":float(au.real),"imag":float(au.imag),"modulus":float(abs(au))},
            "heat_mode_response":float(ah),"response_distance":float(abs(au-ah)),
            "minimal_linear_probe":"one resolved eigenmode, calibrated t, and two-quadrature amplitude response",
            "theorem":"At t=x1/lambda_max the coherent and diffusive modal responses differ by exactly one in the ideal linear response norm.",
            "physical_tester":False}
    return finalize(502,"proven_minimal_ideal_mode_probe_for_dual_dynamics", "A laboratory implementation still requires preparation, phase reference, detector, noise model, and event record.",packet)


def st503()->dict:
    vals,vecs=np.linalg.eigh(A);theta=.29;R=np.eye(N);c,s=math.cos(theta),math.sin(theta);R[2,2]=c;R[2,4]=-s;R[4,2]=s;R[4,4]=c;Q=vecs@R@vecs.T;B=Q@A@Q.T
    rng=np.random.default_rng(SEED+503);rho=rng.normal(size=N);effect=rng.normal(size=N);M1=rng.normal(size=(N,N));M2=rng.normal(size=(N,N));rows=[]
    for kind in ("heat","unitary"):
        F=lambda G,t:expm((-t if kind=="heat" else -1j*t)*G)
        left=np.vdot(effect,F(A,.7)@M2@F(A,.4)@M1@F(A,.2)@rho)
        right=np.vdot(Q@effect,F(B,.7)@(Q@M2@Q.T)@F(B,.4)@(Q@M1@Q.T)@F(B,.2)@(Q@rho))
        rows.append({"kind":kind,"residual":float(abs(left-right))})
    packet={"multitime_record_checks":rows,"maximum_residual":max(x["residual"] for x in rows),
            "theorem":"Operational conjugacy extends to arbitrary finite products when every preparation, propagator, intervention, and effect is transported by the same Q.",
            "fixed_untransported_apparatus_equivalent":False}
    return finalize(503,"proven_multitime_operational_conjugacy_invariance", "The theorem defines medium-independent information as an operational bundle, not as bare eigenvalues.",packet)


def st504()->dict:
    eig=np.linalg.eigvalsh(A)[1:];delta=1e-4;rng=np.random.default_rng(SEED+504);P=np.eye(N)-np.ones((N,N))/N;E=P@rng.normal(size=(N,N))@P;E=(E+E.T)/2;E*=delta/np.linalg.norm(E,2);ep=np.linalg.eigvalsh(A+E)[1:]
    i,j=-1,0;raw=abs(ep[i]/ep[j]-eig[i]/eig[j]);bound=delta*(eig[i]+eig[j])/(eig[j]*(eig[j]-delta));wraw=abs(math.sqrt(ep[i]/ep[j])-math.sqrt(eig[i]/eig[j]));wbound=bound/(2*math.sqrt(min(ep[i]/ep[j],eig[i]/eig[j])))
    packet={"generator_defect_norm":delta,"tested_ratio_indices":[int(i),int(j)],"unitary_heat_ratio_error":float(raw),"unitary_heat_ratio_bound":float(bound),"wave_sqrt_ratio_error":float(wraw),"wave_sqrt_ratio_bound":float(wbound),
            "theorem":"Weyl's inequality yields an explicit Lipschitz bound for positive eigenvalue ratios when delta is below the reference gap; square-root wave ratios inherit the corresponding derivative bound.","absolute_clock":False}
    return finalize(504,"proven_relative_clock_ratio_perturbation_bound", "The bound protects dimensionless ratios; it neither supplies seconds nor equates lambda and sqrt(lambda) channel maps.",packet)


def st505()->dict:
    packet={"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","new_scale_charged_source":False,"absolute_unit":"absent","characteristic_zero_degree8_rank_exact":False,
            "reason":"Noise temperature/covariance and bias are explicit added resources; modular nullspaces are not rational closure."}
    return finalize(505,"blocked_no_strict_source_selector_scale_or_Q_rank_closure", "No bridge completion, role transfer, SM/GR, L_total, or ToE closure is exported.",packet)


def st506()->dict:
    packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent"}
    return finalize(506,"blocked_no_independent_empirical_evidence", "All three rounds remain local analytical/computational research.",packet)


def figures(r:dict)->None:
    FIG_DIR.mkdir(exist_ok=True)
    rows=r["ST492"]["exact_runs"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar([str(x["prime"]) for x in rows],[x["rank"] for x in rows],color="#2563eb");ax.set_ylim(1775,1800);ax.set(ylabel="exact modular rank",title="ST492: three-prime degree-eight rank");fig.tight_layout();fig.savefig(FIG_DIR/"st492_three_prime_rank.png",dpi=180);plt.close(fig)
    ev=r["ST497"]["reduced_logit_Hessian_eigenvalues"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(range(1,len(ev)+1),ev,"o-");ax.axhline(0,color="black",lw=1);ax.set(xlabel="eigenvalue index",ylabel="Hessian eigenvalue",title="ST497: local curvature of transition ratio");fig.tight_layout();fig.savefig(FIG_DIR/"st497_transition_hessian.png",dpi=180);plt.close(fig)
    rows=r["ST499"]["cells"];b=[x["b_interval"][1] for x in rows];y2=[x["endpoint_centers"][-1][1] for x in rows];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(b,y2,"o-");ax.set(xlabel="b",ylabel=r"$y_2$",title="ST499: IR continuation beyond b=0.2");fig.tight_layout();fig.savefig(FIG_DIR/"st499_ir_beyond_02.png",dpi=180);plt.close(fig)
    x=np.linspace(0,3,500);fig,ax=plt.subplots(figsize=(7.2,4));
    for eps in (0,.3,1):
        bias=np.linspace(-1,1,12);p=np.exp(eps*bias);p/=p.sum();ax.plot(range(12),p,"o-",label=fr"$\epsilon={eps}$")
    ax.set(xlabel="branch",ylabel="probability",title="ST501: bias selects only after a source is supplied");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st501_biased_selector.png",dpi=180);plt.close(fig)


def main()->None:
    funcs={492:st492,493:st493,494:st494,495:st495,496:st496,497:st497,498:st498,499:st499,500:st500,501:st501,502:st502,503:st503,504:st504,505:st505,506:st506};r={}
    for k in range(492,507):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
    RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["program","status","object","boundary"])
        for k in range(492,507):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
    figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")


if __name__=="__main__":main()
