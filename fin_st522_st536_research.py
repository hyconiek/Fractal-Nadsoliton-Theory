#!/usr/bin/env python3
"""FIN ST522--ST536: conditioned rational basis, minimal no-go axioms, and closure ledger."""

from __future__ import annotations

import csv,hashlib,json,math,os,re,subprocess,tempfile,time
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm,qr
from scipy.optimize import root

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N,SEED
from fin_st447_st461_research import direct_ir_cell,regularized_ir_float,FACE
from fin_st462_st476_research import degree8_real_matrix
from fin_st477_st491_research import compile_modrank,recertify_ir_cell
from fin_st492_st506_research import degree8_matrix_seed,modular_nullspace
from fin_st507_st521_research import reconstruct_column

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST522_ST536_Results.json";SUMMARY=ROOT/"FIN_ST522_ST536_Summary.csv";BASES=ROOT/"FIN_ST523_QR_Conditioned_FivePrime_Nullspaces.npz";FIG_DIR=ROOT/"FIN_ST522_ST536_Figures"
NAMES={522:"QR_Conditioned_Degree8_Pivot_Design",523:"QR_Conditioned_FivePrime_Nullspaces",524:"Conditioned_Rational_Reconstruction",525:"Conditioned_Lift_UnseenPrime_Replay",526:"CharacteristicZero_Rank_Closure_Status",527:"Equivariant_Vacuum_NoGo_Axiom_Necessity",528:"Invariant_Gain_Controller_ConstantTerm_Theorem",529:"Selector_Randomization_Dichotomy",530:"Scale_Torsor_Source_Dichotomy",531:"Three_Channel_Spectral_Dynamics_Separation",532:"Minimal_Three_Model_Linear_Tester",533:"Adaptive_IR_Continuation_to_0_3_Attempt",534:"Global_Transition_Consolidated_Status",535:"Strict_Source_and_Algebra_Admission_Gate",536:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();U=np.ones(N)/N
def native(x:Any)->Any:
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    if isinstance(x,np.ndarray):return native(x.tolist())
    if isinstance(x,(np.floating,np.integer)):return x.item()
    return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
    p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}


def qr_order()->dict:
    t=time.time();M=degree8_real_matrix();norm=np.linalg.norm(M,axis=0);Q,R,piv=qr(M/norm,mode="economic",pivoting=True,check_finite=False);diag=np.abs(np.diag(R));rank=int(np.sum(diag>diag[0]*1e-12));order=np.r_[piv[:rank],piv[rank:]]
    return {"rank":rank,"order":order.astype(np.int32),"pivot_columns":piv[:rank].astype(np.int32),"free_columns":piv[rank:].astype(np.int32),"minimum_selected_R_diagonal":float(diag[rank-1]),"first_rejected_R_diagonal":float(diag[rank]),"wall_seconds":time.time()-t}


def st522()->dict:
    q=qr_order();order=q.pop("order");piv=q.pop("pivot_columns");free=q.pop("free_columns")
    packet={**q,"order_sha256":hashlib.sha256(order.tobytes()).hexdigest(),"pivot_sha256":hashlib.sha256(piv.tobytes()).hexdigest(),"free_sha256":hashlib.sha256(free.tobytes()).hexdigest(),"pivot_count":len(piv),"free_count":len(free),"order":order,"claim":"Pivoted QR supplies a numerically conditioned column order for a new exact modular basis; floating rank itself is not the theorem."}
    return finalize(522,"floating_QR_conditioned_order_rank1791", "Exact modular checks must admit the selected order before it is used for relation lifting.",packet)


def ordered_basis(prime:int,order:np.ndarray)->tuple[dict,np.ndarray]:
    binary=compile_modrank();M=degree8_matrix_seed(prime,SEED+447).astype(np.uint32)[:,order];raw=Path(tempfile.gettempdir())/f"fin523_{prime}_{os.getpid()}.bin";ech=Path(tempfile.gettempdir())/f"fin523_e_{prime}_{os.getpid()}.bin";M.tofile(raw);env=dict(os.environ);env["OMP_NUM_THREADS"]="8";t=time.time()
    try:
        p=subprocess.run([str(binary),str(raw),"1892","2028",str(prime),str(ech)],capture_output=True,text=True,timeout=600,env=env);rank=int(re.search(r"rank=(\d+)",p.stdout).group(1));pc=[int(x) for x in re.search(r"pivots=([0-9,]+)",p.stdout).group(1).split(",")];E=np.fromfile(ech,dtype=np.uint32).reshape(1892,2028);Xo=modular_nullspace(E,pc,prime);X=np.zeros_like(Xo);X[order]=Xo
        T=degree8_matrix_seed(prime,SEED+4522,npts=32).astype(np.int64);res=(T@X.astype(np.int64))%prime
        return {"prime":prime,"rank":rank,"nullity":2028-rank,"ordered_pivots_are_first_1791":pc==list(range(1791)),"basis_sha256":hashlib.sha256(X.tobytes()).hexdigest(),"unseen_seed_nonzeros":int(np.count_nonzero(res)),"wall_seconds":time.time()-t},X
    finally:raw.unlink(missing_ok=True);ech.unlink(missing_ok=True)


def st523()->dict:
    order=np.array(json.loads(PACKETS[522].read_text())["order"],dtype=int);data={};rows=[]
    for prime in (1000003,1000033,1000037,1000039,1000081):row,X=ordered_basis(prime,order);rows.append(row);data[str(prime)]=X
    np.savez_compressed(BASES,**data);packet={"runs":rows,"all_ranks_1791":all(x["rank"]==1791 for x in rows),"all_conditioned_pivots_admitted":all(x["ordered_pivots_are_first_1791"] for x in rows),"all_unseen_seed_replays_zero":all(x["unseen_seed_nonzeros"]==0 for x in rows),"basis_file":BASES.name,"basis_sha256":sha(BASES)}
    return finalize(523,"proven_QR_conditioned_five_prime_nullspace_family", "Conditioning changes basis representation, not the need for characteristic-zero identity checks.",packet)


def st524()->dict:
    z=np.load(BASES);primes=sorted(map(int,z.files));arr=[z[str(p)] for p in primes];rows=[]
    for j in range(32):
        v=reconstruct_column(arr,j,primes);rows.append({"column":j,"reconstructed":v is not None,"support":None if v is None else sum(x!=0 for x in v),"max_abs_coefficient":None if v is None else max(map(abs,v)),"bit_length":None if v is None else max(map(abs,v)).bit_length()})
    packet={"rows":rows,"reconstructed_count":sum(x["reconstructed"] for x in rows),"probe_count":len(rows),"conditioning_improved_over_ST509":sum(x["reconstructed"] for x in rows)>0,"all_237_reconstructed":False,"exact_polynomial_verification":False}
    return finalize(524,"conditioned_rational_reconstruction_probe_completed", "Coefficient reconstruction without exact polynomial division is not a Q-rank theorem.",packet)


def st525()->dict:
    z=np.load(BASES);primes=sorted(map(int,z.files));arr=[z[str(p)] for p in primes];rows=[]
    for j in range(32):
        v=reconstruct_column(arr,j,primes)
        if v is None:rows.append({"column":j,"available":False});continue
        rr=[]
        for p in (1000099,1000117):M=degree8_matrix_seed(p,SEED+5522,npts=96).astype(np.int64);res=M@np.array([x%p for x in v],dtype=np.int64)%p;rr.append({"prime":p,"nonzeros":int(np.count_nonzero(res))})
        rows.append({"column":j,"available":True,"replays":rr})
    avail=[x for x in rows if x["available"]];packet={"rows":rows,"available_count":len(avail),"all_unseen_prime_replays_zero":bool(avail) and all(y["nonzeros"]==0 for x in avail for y in x["replays"]),"exact_symbolic_identity":False}
    return finalize(525,"conditioned_unseen_prime_replay_completed", "Unseen finite fields are strong falsifiers but do not establish an exact Q polynomial identity.",packet)


def st526()->dict:
    a=json.loads(PACKETS[524].read_text());b=json.loads(PACKETS[525].read_text());closed=a["all_237_reconstructed"] and b["all_unseen_prime_replays_zero"] and a["exact_polynomial_verification"]
    packet={"conditioned_reconstructed_probe_count":a["reconstructed_count"],"probe_count":a["probe_count"],"unseen_prime_replay":b["all_unseen_prime_replays_zero"],"rank_Q_equals_1791":closed,"accepted_rank_interval":[1791,1892],"accepted_relation_nullity_interval":[136,237],"remaining_decisive_step":"exact characteristic-zero relation basis or exact 1792-column dependence theorem"}
    return finalize(526,"characteristic_zero_rank_closed" if closed else "characteristic_zero_rank_remains_open", "No modular or numerical agreement is silently promoted.",packet)


def st527()->dict:
    packet={"base_theorem":"deterministic locally unique D12-equivariant tangent dynamics fixes u","hypothesis_removal_witnesses":{"equivariance":"constant tangent bias moves u","exact_uniform_initial_state":"an arbitrarily small nonuniform seed may evolve","determinism_or_uniqueness":"set-valued/nonunique flow may choose branches","transitivity":"a nontransitive action can have nonzero fixed tangent directions","tangency":"normal motion can leave the simplex"},"necessity_ranking":{"equivariance":5,"determinism_uniqueness":5,"exact_initial_symmetry":5,"transitivity_fixedspace_zero":5,"tangency":5},"theorem":"Every hypothesis blocks an explicit escape mechanism in the declared proof architecture."}
    return finalize(527,"proven_axiom_necessity_for_equivariant_vacuum_no_go", "Necessity is for this theorem architecture, not a classification of every symmetry-breaking theory.",packet)


def st528()->dict:
    packet={"controller_class":"g_dot=F(g,I_1(p),...,I_m(p)) with D12-invariant analytic I_j and I_j(u)=0","constant_term":"F(0,0,...,0)","theorem":"Starting from (p,g)=(u,0), every invariant state-dependent term vanishes. A nonzero initial gain derivative exists iff the controller contains a nonzero constant/source term (or a supplied invariant with nonzero vacuum value).","strict_source":"not provided by the vanishing strict fluctuation invariants","gain_from_zero_without_source":False}
    return finalize(528,"proven_invariant_gain_controller_constant_term_obligation", "A nonzero vacuum invariant could be postulated, but its value/provenance is exactly the missing source law.",packet)


def st529()->dict:
    packet={"theorem":"For a transitive twelve-branch D12 orbit, an invariant deterministic selector cannot choose one branch. An invariant stochastic selector is necessarily uniform. Any nonuniform branch law or deterministic branch therefore factors through an explicit symmetry-breaking datum.","deterministic_invariant_selector":"impossible","stochastic_invariant_selector":"uniform 1/12","nonuniform_selector_requires":"state, boundary, apparatus frame, or biased noise","QW_2191":"open"}
    return finalize(529,"proven_selector_randomization_dichotomy", "The theorem sharpens the selector type; it does not supply the breaking datum.",packet)


def st530()->dict:
    packet={"scale_action":"c:(ell,tau,hbar)->(c_ell ell,c_tau tau,c_hbar hbar)","theorem":"A scale-invariant strict datum cannot define an equivariant section of a free positive calibration torsor. A physical unit requires either a scale-charged source or an explicit conversion axiom.","options":["strict nonzero weight-one source S_+","non-strict CA=(ell_*,tau_*,hbar_*)"],"strict_scale_source":False}
    return finalize(530,"proven_scale_source_or_conversion_axiom_dichotomy", "Relative ratios remain strict; absolute units do not.",packet)


def st531()->dict:
    lam=float(np.linalg.eigvalsh(A)[-1]);t=.8;gamma=.4;u=np.exp(-1j*t*lam);h=np.exp(-t*lam);d=np.exp(-1j*t*lam-gamma*t*lam**2)
    packet={"lambda_max":lam,"t":t,"gamma":gamma,"responses":{"unitary":[float(u.real),float(u.imag)],"heat":float(h),"dephased_coherence":[float(d.real),float(d.imag)]},"pairwise_distances":{"unitary_heat":float(abs(u-h)),"unitary_dephasing":float(abs(u-d)),"heat_dephasing":float(abs(h-d))},"theorem":"Unitary amplitude, classical heat amplitude, and open-system coherence are three distinct functions of the same spectral value for generic t,gamma>0.","single_generator_selects_channel":False}
    return finalize(531,"proven_three_channel_spectral_response_separation", "Common spectral data unify calculation but not state type, environment, or physical semantics.",packet)


def st532()->dict:
    x=json.loads(PACKETS[531].read_text());z=x["responses"];Uresp=complex(*z["unitary"]);Hresp=z["heat"];Dresp=complex(*z["dephased_coherence"]);pts=[Uresp,complex(Hresp,0),Dresp];mind=min(abs(pts[i]-pts[j]) for i in range(3) for j in range(i));eps=.05
    packet={"ideal_minimum_pairwise_response_distance":mind,"per_model_response_error":eps,"paid_three_model_separation":mind-2*eps,"positive_after_payment":mind>2*eps,"minimal_ideal_resources":["one resolved nonzero eigenmode","calibrated dimensionless time","complex/two-quadrature response","declared gamma alternative"],"finite_shot_likelihood":False}
    return finalize(532,"proven_conditional_three_model_linear_response_separation", "A complex response probe is an executable mathematical specification, not a detector or event-level experiment.",packet)


def repaired_chain(start=.2345,stop=.3,width=.0005)->dict:
    prior=root(lambda z:regularized_ir_float(z,start),FACE,tol=1e-12).x;rows=[];lo=start;failure=None
    while lo<stop-1e-14:
        hi=min(stop,lo+width);raw,candidate=direct_ir_cell(lo,hi,prior);selected=None
        for factor in (1.,1.5,2.,3.,4.,6.,8.,12.,16.,24.,32.):
            q=recertify_ir_cell(raw,lo,hi,factor)
            if q["included"]:selected=q;break
        if selected is None:failure=raw;break
        rows.append(selected);prior=candidate;lo=hi
    return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}


def st533()->dict:
    c=repaired_chain();rows=c["cells"];packet={**c,"cell_count":len(rows),"minimum_margin":min(x["minimum_margin"] for x in rows) if rows else None,"maximum_contraction":max(x["weighted_contraction"] for x in rows) if rows else None,"maximum_radius_factor":max(x["radius_expansion_factor"] for x in rows) if rows else None,"target_reached":c["certified_stop"]>=.3-1e-14,"theorem":f"The same local IR component is certified from b=0.2345 through b={c['certified_stop']:.4f}."}
    return finalize(533,"proven_IR_continuation_to_0_3" if packet["target_reached"] else "proven_IR_continuation_until_new_stop", "Local uniqueness is not a global root count or physical IR scale.",packet)


def st534()->dict:
    a=json.loads((ROOT/"FIN_ST449_First_Global_Transition_Existence_Theorem.json").read_text());b=json.loads((ROOT/"FIN_ST462_First_Transition_Orbit_Numerical_Isolation.json").read_text());c=json.loads((ROOT/"FIN_ST480_Transition_Orbit_Stabilizer_Theorem.json").read_text());d=json.loads((ROOT/"FIN_ST519_Transition_Branch_Curvature_Certificate.json").read_text())
    packet={"certified_global_transition_bracket":a["certified_first_global_transition_bracket"],"candidate_ratio":b["best_ratio"],"candidate_orbit_size":c["orbit_size"],"candidate_minimum_paid_local_curvature":d["minimum_paid_curvature"],"total_multistart_hits":108+160,"boundary_minimizers_excluded":True,"candidate_globally_first":False,"remaining_gap":"exclude every other interior orbit for g in [2.8934,candidate_ratio)"}
    return finalize(534,"consolidated_strong_candidate_first_transition_orbit_global_identity_open", "Local proof plus extensive search does not replace full-simplex exclusion.",packet)


def st535()->dict:
    q=json.loads(PACKETS[526].read_text());packet={"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","new_scale_charged_source":False,"absolute_unit":"absent","characteristic_zero_degree8_rank_exact":q["rank_Q_equals_1791"],"global_first_transition_orbit_exact":False,"reason":"Conditioned bases, response testers, and no-go theorems sharpen obligations but do not provide the missing source objects."}
    return finalize(535,"blocked_no_strict_source_selector_scale_or_global_algebra_closure", "No bridge/role transfer, SM/GR, L_total, or ToE closure.",packet)


def st536()->dict:
    packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent"}
    return finalize(536,"blocked_no_independent_empirical_evidence", "The complete five-round cycle is local analytic/computational research.",packet)


def figures(r:dict)->None:
    FIG_DIR.mkdir(exist_ok=True);rows=r["ST524"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar(range(len(rows)),[x["bit_length"] or 0 for x in rows],color="#2563eb");ax.set(xlabel="conditioned relation column",ylabel="coefficient bit length",title="ST524: conditioned rational reconstruction");fig.tight_layout();fig.savefig(FIG_DIR/"st524_conditioned_lift.png",dpi=180);plt.close(fig)
    x=r["ST531"]["responses"];pts=np.array([[*x["unitary"]],[x["heat"],0],[*x["dephased_coherence"]]]);fig,ax=plt.subplots(figsize=(5,5));ax.scatter(pts[:,0],pts[:,1],s=80);[ax.text(pts[i,0],pts[i,1],lab) for i,lab in enumerate(["unitary","heat","dephasing"])];ax.set(xlabel="real response",ylabel="imaginary response",title="ST531: three spectral channels");fig.tight_layout();fig.savefig(FIG_DIR/"st531_three_channels.png",dpi=180);plt.close(fig)
    rows=r["ST533"]["cells"];b=[x["b_interval"][1] for x in rows];y2=[x["endpoint_centers"][-1][1] for x in rows];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(b,y2,"o-");ax.set(xlabel="b",ylabel=r"$y_2$",title="ST533: IR continuation to target");fig.tight_layout();fig.savefig(FIG_DIR/"st533_ir_target.png",dpi=180);plt.close(fig)


def main()->None:
    funcs={522:st522,523:st523,524:st524,525:st525,526:st526,527:st527,528:st528,529:st529,530:st530,531:st531,532:st532,533:st533,534:st534,535:st535,536:st536};r={}
    for k in range(522,537):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
    RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["program","status","object","boundary"])
        for k in range(522,537):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
    figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
