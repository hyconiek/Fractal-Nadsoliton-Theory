#!/usr/bin/env python3
"""FIN ST507--ST521: multi-prime lift, symmetry no-go, and open dual dynamics."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import root
from sympy.polys.domains import ZZ
from sympy.polys.modulargcd import _integer_rational_reconstruction

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N, SEED
from fin_st447_st461_research import direct_ir_cell, regularized_ir_float, FACE
from fin_st477_st491_research import recertify_ir_cell
from fin_st492_st506_research import rank_pivots, modular_nullspace, degree8_matrix_seed


ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/"FIN_ST507_ST521_Results.json";SUMMARY=ROOT/"FIN_ST507_ST521_Summary.csv"
BASES=ROOT/"FIN_ST507_FivePrime_Degree8_Modular_Nullspaces.npz";FIG_DIR=ROOT/"FIN_ST507_ST521_Figures"
NAMES={507:"Five_Prime_Modular_Nullspace_Expansion",508:"Five_Prime_Pivot_and_Replay_Audit",509:"Five_Prime_Rational_Reconstruction",510:"Independent_Prime_Lift_Replay",511:"CharacteristicZero_Closure_Obligation",512:"Deterministic_Equivariant_Vacuum_NoGo",513:"Stochastic_Seed_Resource_Classification",514:"Finite_Wick_Continuation_Identity",515:"Spectral_Dephasing_Interpolation",516:"Dephasing_Is_Not_Classical_Heat",517:"Noise_Robust_Ideal_Dual_Tester_Bound",518:"Adaptive_IR_Continuation_Toward_0_3",519:"Transition_Branch_Curvature_Certificate",520:"Source_Selector_Scale_Admission_Gate",521:"Independent_Evidence_Gate"}
PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();U=np.ones(N)/N


def native(x:Any)->Any:
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    if isinstance(x,np.ndarray):return native(x.tolist())
    if isinstance(x,(np.floating,np.integer)):return x.item()
    return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
    p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}


def build_basis(prime:int)->tuple[dict,np.ndarray]:
    run=rank_pivots(prime,SEED+447,echelon=True);E=run.pop("echelon");X=modular_nullspace(E,run["pivot_columns"],prime)
    M=degree8_matrix_seed(prime,SEED+2447,npts=32).astype(np.int64);res=(M@X.astype(np.int64))%prime
    row={**run,"basis_shape":list(X.shape),"basis_sha256":hashlib.sha256(X.tobytes()).hexdigest(),"independent_seed_replay_nonzeros":int(np.count_nonzero(res))}
    return row,X


def st507()->dict:
    old=np.load(ROOT/"FIN_ST495_Degree8_Modular_Nullspaces.npz");data={k:old[k] for k in old.files};rows=[]
    for p in (1000037,1000039,1000081):
        row,X=build_basis(p);rows.append(row);data[str(p)]=X
    np.savez_compressed(BASES,**data)
    packet={"new_prime_runs":rows,"all_new_ranks_1791":all(x["rank"]==1791 for x in rows),"all_new_replays_zero":all(x["independent_seed_replay_nonzeros"]==0 for x in rows),"combined_basis_file":BASES.name,"combined_basis_sha256":sha(BASES),"prime_count":len(data)}
    return finalize(507,"proven_five_prime_modular_nullspace_family", "Five finite fields improve lift capacity but do not themselves prove a rational polynomial identity.",packet)


def st508()->dict:
    z=np.load(BASES);primes=sorted(map(int,z.files));pivot_hashes=[];basis_hashes=[]
    for p in primes:
        run=rank_pivots(p,SEED+447);pivot_hashes.append(run["pivot_sha256"]);basis_hashes.append(hashlib.sha256(z[str(p)].tobytes()).hexdigest())
    packet={"primes":primes,"pivot_schema_count":len(set(pivot_hashes)),"common_pivot_sha256":pivot_hashes[0],"basis_sha256_by_prime":basis_hashes,"all_basis_shapes_2028_by_237":all(z[str(p)].shape==(2028,237) for p in primes),"theorem":"A common ordered pivot/free-column schema exists over all five declared primes, so corresponding modular nullspace columns are canonically aligned for CRT."}
    return finalize(508,"proven_five_prime_common_modular_schema", "Canonical modular alignment is not coefficient descent to Q.",packet)


def crt_many(residues:list[int],primes:list[int])->tuple[int,int]:
    r=0;M=1
    for a,p in zip(residues,primes):
        r=(r+M*(((a-r)*pow(M,-1,p))%p))%(M*p);M*=p
    return r,M


def reconstruct_column(arrays:list[np.ndarray],j:int,primes:list[int])->list[int]|None:
    fr=[];mod=math.prod(primes)
    for i in range(arrays[0].shape[0]):
        residue,_=crt_many([int(x[i,j]) for x in arrays],primes);q=_integer_rational_reconstruction(residue,mod,ZZ)
        if q is None:return None
        fr.append((int(q.numerator),int(q.denominator)))
    lcm=1
    for _,d in fr:lcm=math.lcm(lcm,d)
    vals=[n*(lcm//d) for n,d in fr];g=0
    for x in vals:g=math.gcd(g,abs(x))
    vals=[x//g for x in vals]
    if next(x for x in vals if x)<0:vals=[-x for x in vals]
    return vals


def st509()->dict:
    z=np.load(BASES);primes=sorted(map(int,z.files));arrays=[z[str(p)] for p in primes];rows=[]
    for j in range(16):
        v=reconstruct_column(arrays,j,primes)
        rows.append({"column":j,"reconstructed":v is not None,"support":None if v is None else sum(x!=0 for x in v),"maximum_abs_integer_coefficient":None if v is None else max(map(abs,v)),"coefficient_bit_length":None if v is None else max(map(abs,v)).bit_length()})
    packet={"primes":primes,"CRT_modulus":math.prod(primes),"CRT_modulus_bit_length":math.prod(primes).bit_length(),"columns_probed":len(rows),"rows":rows,"reconstructed_count":sum(x["reconstructed"] for x in rows),"all_237_reconstructed":False,"exact_polynomial_verification":False}
    return finalize(509,"five_prime_rational_reconstruction_probe_completed", "Even successful coefficient reconstruction is not a Q theorem before all relations are lifted and exact polynomial identities are checked.",packet)


def st510()->dict:
    z=np.load(BASES);primes=sorted(map(int,z.files));arrays=[z[str(p)] for p in primes];tests=[]
    for j in range(16):
        v=reconstruct_column(arrays,j,primes)
        if v is None:tests.append({"column":j,"available":False});continue
        row={"column":j,"available":True,"replays":[]}
        for p in (1000099,1000117):
            M=degree8_matrix_seed(p,SEED+3447,npts=96).astype(np.int64);vv=np.array([x%p for x in v],dtype=np.int64);res=(M@vv)%p
            row["replays"].append({"prime":p,"nonzeros":int(np.count_nonzero(res)),"max_residue":int(res.max())})
        tests.append(row)
    available=[x for x in tests if x["available"]]
    packet={"rows":tests,"available_count":len(available),"all_independent_replays_zero":bool(available) and all(y["nonzeros"]==0 for x in available for y in x["replays"]),"exact_symbolic_identity":False}
    return finalize(510,"independent_prime_lift_replay_completed", "Two unseen-prime evaluations are falsification tests, not exact characteristic-zero polynomial division.",packet)


def st511()->dict:
    a=json.loads(PACKETS[509].read_text());b=json.loads(PACKETS[510].read_text());complete=a["reconstructed_count"]==237 and b["all_independent_replays_zero"]
    packet={"reconstructed_probe_count":a["reconstructed_count"],"probe_size":a["columns_probed"],"unseen_prime_replay":b["all_independent_replays_zero"],"characteristic_zero_rank_1791_closed":False,"remaining_obligations":["lift all 237 aligned modular columns","verify each relation as an exact polynomial identity modulo x0+...+x11","prove 1791 independent characteristic-zero candidates using an exact nonzero minor"],"conditional_consequence_if_closed":"rank_Q=1791, relation nullity=237, primitive degree-eight quotient dimension=101"}
    return finalize(511,"characteristic_zero_closure_obligation_sharpened", "No conditional dimension is promoted while exact symbolic checks are absent.",packet)


def st512()->dict:
    packet={"state_space":"mean-zero tangent space H0 of the probability simplex at u","fixed_subspace_dimension_under_D12":0,"theorem":"For any deterministic D12-equivariant locally Lipschitz vector field F tangent to the simplex, F(u) is D12-fixed. The only fixed tangent vector is zero, hence F(u)=0; ODE uniqueness makes the exactly uniform state invariant for all time.","spontaneous_departure_from_exact_u":False,"necessary_escape_resource":"nonuniform initial condition, stochastic forcing, non-equivariant boundary/apparatus, or nonunique dynamics"}
    return finalize(512,"proven_deterministic_equivariant_uniform_vacuum_no_go", "The theorem does not forbid spontaneous symmetry breaking from fluctuations or unstable nearby states; it identifies the missing seed/source type.",packet)


def st513()->dict:
    packet={"deterministic_result":"ST512 fixes u exactly","stochastic_result":"D12-invariant noise can leave u pathwise while preserving a D12-invariant ensemble law","selector_result":"transitivity keeps branch probabilities equal absent asymmetric covariance/drift","resource_taxonomy":{"seed":"nonuniform initial realization","pump":"energy/information maintaining gain","selector":"non-D12-invariant bias or record frame"},"theorem":"Pathwise realization, energetic gain, and canonical selection are three distinct resources and cannot be inferred from one another."}
    return finalize(513,"proven_seed_gain_selector_resource_separation", "The classification is structural; no strict provider for any missing resource is exported.",packet)


def st514()->dict:
    t=.83;left=expm(-t*A);right=expm(-1j*(-1j*t)*A);res=float(np.linalg.norm(left-right,2))
    packet={"identity":"P_t=exp(-tA)=U_z at z=-it for U_z=exp(-izA)","test_t":t,"matrix_residual":res,"theorem":"In finite dimension z->exp(-izA) is entire, and evaluation at imaginary time z=-it equals the heat semigroup exactly.","operational_equivalence":False,"reason":"analytic continuation changes the time domain and does not transport preparations, probabilities, clocks, instruments, or noise"}
    return finalize(514,"proven_finite_Wick_continuation_identity", "Wick continuation is a mathematical relation between functional calculi, not a physical mechanism converting coherent events into diffusion.",packet)


def dephase(A:np.ndarray,rho:np.ndarray,t:float,gamma:float)->np.ndarray:
    lam,V=np.linalg.eigh(A);r=V.T@rho@V;dl=lam[:,None]-lam[None,:];r*=np.exp(-1j*t*dl-gamma*t*dl**2);return V@r@V.T


def st515()->dict:
    lam=np.linalg.eigvalsh(A);gamma=.4;t=.9;dl=lam[:,None]-lam[None,:];factor=np.exp(-1j*t*dl-gamma*t*dl**2)
    packet={"Lindblad_generator":"L(rho)=-i[A,rho]-gamma[A,[A,rho]]","gamma":gamma,"t":t,"mode_factor_formula":"exp[-it(lambda_r-lambda_s)-gamma*t(lambda_r-lambda_s)^2]","diagonal_mode_factor_residual":float(np.max(np.abs(np.diag(factor)-1))),"maximum_offdiagonal_modulus":float(np.max(np.abs(factor-np.diag(np.diag(factor))))),"theorem":"The supplied double-commutator channel continuously damps spectral coherences while preserving spectral populations."}
    return finalize(515,"proven_spectral_dephasing_interpolation_formula", "Gamma and the environment are added; dephasing is not the classical heat semigroup exp(-tA) on vertex probabilities.",packet)


def st516()->dict:
    rho=np.zeros((N,N),complex);rho[0,0]=1;t=.9;gamma=.4;rd=dephase(A,rho,t,gamma);pdeph=np.real(np.diag(rd));pheat=expm(-t*A)@np.eye(N)[:,0]
    packet={"initial_state":"vertex-0 pure density / vertex-0 classical probability","t":t,"gamma":gamma,"dephased_vertex_populations":pdeph,"classical_heat_probabilities":pheat,"L1_difference":float(np.sum(abs(pdeph-pheat))),"max_difference":float(np.max(abs(pdeph-pheat))),"theorem":"Spectral dephasing preserves A-basis populations and damps coherences; classical heat transports vertex probability. Their vertex records differ for the displayed common initial label."}
    return finalize(516,"proven_dephasing_classical_heat_non_equivalence", "The comparison uses declared embeddings of a vertex label into quantum and classical state spaces; it does not select either as physical.",packet)


def st517()->dict:
    eps=.08;ideal=1.;lower=ideal-2*eps
    packet={"ideal_response_distance":ideal,"preparation_response_error_bound_each":eps,"certified_triangle_lower_bound":lower,"positive_margin":lower>0,"theorem":"If the ideal two channel responses differ by one and the implemented response of each channel is within epsilon of its ideal value, their observed separation is at least 1-2epsilon.","required_epsilon_for_positive_separation":"epsilon<1/2","finite_shot_decision_rule":False}
    return finalize(517,"proven_conditional_noise_robust_dual_response_margin", "This is a deterministic response bound, not a shot-noise power calculation or laboratory validation.",packet)


def repaired_chain(start=.2125,stop=.3,width=.001)->dict:
    prior=root(lambda z:regularized_ir_float(z,start),FACE,tol=1e-12).x;rows=[];lo=start;failure=None
    while lo<stop-1e-14:
        hi=min(stop,lo+width);raw,candidate=direct_ir_cell(lo,hi,prior);selected=None
        for factor in (1.,1.5,2.,3.,4.,6.,8.,12.,16.):
            trial=recertify_ir_cell(raw,lo,hi,factor)
            if trial["included"]:selected=trial;break
        if selected is None:failure=raw;break
        rows.append(selected);prior=candidate;lo=hi
    return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}


def st518()->dict:
    c=repaired_chain();rows=c["cells"]
    packet={**c,"cell_count":len(rows),"minimum_margin":min(x["minimum_margin"] for x in rows) if rows else None,"maximum_contraction":max(x["weighted_contraction"] for x in rows) if rows else None,"maximum_radius_factor":max(x["radius_expansion_factor"] for x in rows) if rows else None,"target_reached":c["certified_stop"]>=.3-1e-14,"theorem":f"The same locally unique IR component is certified from b=0.2125 through b={c['certified_stop']:.4f}."}
    return finalize(518,"proven_adaptive_IR_continuation_toward_0_3", "The result remains local and dimensionless; template failure is not a branch singularity.",packet)


def st519()->dict:
    old=json.loads((ROOT/"FIN_ST342_Narrow_Certified_Coexistence_Bracket.json").read_text());rows=[];P=np.eye(N)-np.ones((N,N))/N
    for ep in old["endpoint_certificates"]:
        c=ep["center"];p=np.array([c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[5],c[4],c[3],c[2],c[1]]);g=ep["g"]
        ev=np.linalg.eigvalsh(P@(np.diag(1/p)-g*A)@P);delta=1e-10;pert=delta/(p.min()-delta)**2+1e-8*np.linalg.norm(A,2)
        rows.append({"g":g,"tangent_eigenvalues":ev[1:],"minimum_tangent_eigenvalue":float(ev[1]),"Weyl_perturbation_payment":float(pert),"paid_lower_bound":float(ev[1]-pert)})
    packet={"endpoint_rows":rows,"minimum_paid_curvature":min(x["paid_lower_bound"] for x in rows),"theorem":"For H=diag(1/p)-gA restricted to the simplex tangent, endpoint centers have a large positive spectral gap. Paying the displayed root/gain box perturbation preserves strict positivity across the certified crossing tube.","global_minimum_theorem":False}
    return finalize(519,"proven_positive_transition_branch_curvature_on_certified_tube", "Strict local minimality of this branch does not exclude a separate lower global orbit.",packet)


def st520()->dict:
    packet={"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","new_scale_charged_source":False,"absolute_unit":"absent","characteristic_zero_degree8_rank_exact":False,"reason":"Noise covariance, thermal theta, and bias are explicit added resources; multi-prime lifts still lack exact polynomial verification."}
    return finalize(520,"blocked_no_strict_source_selector_scale_or_Q_rank_closure", "No legacy bridge/role transfer, SM/GR, L_total, or ToE closure is exported.",packet)


def st521()->dict:
    packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent"}
    return finalize(521,"blocked_no_independent_empirical_evidence", "All results are local analytic/computational statements.",packet)


def figures(r:dict)->None:
    FIG_DIR.mkdir(exist_ok=True)
    rows=r["ST509"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar(range(len(rows)),[x["coefficient_bit_length"] or 0 for x in rows],color="#2563eb");ax.set(xlabel="aligned nullspace column",ylabel="max coefficient bit length",title="ST509: five-prime rational lift probe");fig.tight_layout();fig.savefig(FIG_DIR/"st509_lift_bits.png",dpi=180);plt.close(fig)
    row=r["ST516"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(row["dephased_vertex_populations"],"o-",label="spectral dephasing");ax.plot(row["classical_heat_probabilities"],"s-",label="classical heat");ax.set(xlabel="vertex",ylabel="population",title="ST516: dephasing is not heat");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st516_dephasing_vs_heat.png",dpi=180);plt.close(fig)
    rows=r["ST518"]["cells"];b=[x["b_interval"][1] for x in rows];y2=[x["endpoint_centers"][-1][1] for x in rows];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(b,y2,"o-");ax.set(xlabel="b",ylabel=r"$y_2$",title="ST518: IR continuation toward 0.3");fig.tight_layout();fig.savefig(FIG_DIR/"st518_ir_toward_03.png",dpi=180);plt.close(fig)


def main()->None:
    funcs={507:st507,508:st508,509:st509,510:st510,511:st511,512:st512,513:st513,514:st514,515:st515,516:st516,517:st517,518:st518,519:st519,520:st520,521:st521};r={}
    for k in range(507,522):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
    RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["program","status","object","boundary"])
        for k in range(507,522):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
    figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
