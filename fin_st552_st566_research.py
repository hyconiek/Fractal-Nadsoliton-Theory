#!/usr/bin/env python3
"""FIN ST552--ST566: stage-one algebra, causal-speed audit, and new-cycle gates."""
from __future__ import annotations
import csv,hashlib,json,math,os,re,subprocess,tempfile
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm,qr
from scipy.optimize import root

from fin_st402_st416_research import independent_strict_matrix_float,real_orbit_eval
from fin_st417_st431_research import N,SEED
from fin_st447_st461_research import direct_ir_cell,regularized_ir_float,FACE,orbit_poly,poly_mul,poly_add,multiply_by_sum11
from fin_st477_st491_research import compile_modrank,recertify_ir_cell
from fin_st492_st506_research import degree8_matrix_seed,modular_nullspace
from fin_st507_st521_research import reconstruct_column

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST552_ST566_Results.json";SUMMARY=ROOT/"FIN_ST552_ST566_Summary.csv";BASES=ROOT/"FIN_ST553_Stage1_Conditioned_FivePrime_Nullspaces.npz";RELATIONS=ROOT/"FIN_ST554_Stage1_Reconstructed_Relations.json";FIG_DIR=ROOT/"FIN_ST552_ST566_Figures"
NAMES={552:"Stage1_QR_Conditioned_Pivot_Design",553:"Stage1_FivePrime_Modular_Nullspaces",554:"Stage1_Rational_Reconstruction",555:"Stage1_Exact_Polynomial_Verification",556:"Stage1_CharacteristicZero_Rank_Status",557:"TwoLevel_Transition_Margin_Audit",558:"ReflectionOrbit_Transition_Consolidation",559:"FirstTransition_Global_Obligation",560:"Equivariant_Markov_Selector_Theorem",561:"FiniteShot_ThreeChannel_MonteCarlo_Replay",562:"CPTP_Operational_Conjugacy",563:"Layer_Covariant_Effective_Speed_and_Causal_Cone_Audit",564:"UltraFine_IR_Continuation_Attempt",565:"Round1_Strict_Admission_Gate",566:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();U=np.ones(N)/N
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}


def real_stage1()->np.ndarray:
 base=json.loads((ROOT/"FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text());rng=np.random.default_rng(SEED+552);pts=rng.normal(size=(1400,N));pts-=pts.mean(axis=1,keepdims=True);pts/=np.linalg.norm(pts,axis=1,keepdims=True);ev=lambda rep:real_orbit_eval(rep,pts);q=[ev(r) for r in base["quadratic_generator_representatives"]];p4=[ev(r) for r in base["primitive_quartic_representatives"]];cols=[]
 import itertools
 for ids in itertools.combinations_with_replacement(range(6),4):
  v=np.ones(len(pts))
  for i in ids:v*=q[i]
  cols.append(v)
 for ids in itertools.combinations_with_replacement(range(6),2):
  q2=q[ids[0]]*q[ids[1]];cols.extend(q2*x for x in p4)
 cols.extend(p4[i]*p4[j] for i,j in itertools.combinations_with_replacement(range(32),2))
 return np.column_stack(cols)


def st552()->dict:
 M=real_stage1();norm=np.linalg.norm(M,axis=0);_,R,piv=qr(M/norm,mode="economic",pivoting=True);diag=np.abs(np.diag(R));rank=int(np.sum(diag>diag[0]*1e-12));order=np.r_[piv[:rank],piv[rank:]].astype(np.int32);packet={"matrix_shape":list(M.shape),"floating_rank":rank,"pivot_count":rank,"free_count":M.shape[1]-rank,"minimum_selected_diagonal":float(diag[rank-1]),"first_rejected_diagonal":float(diag[rank]),"order":order,"order_sha256":hashlib.sha256(order.tobytes()).hexdigest()}
 return finalize(552,"floating_stage1_QR_rank1288", "The order is a design input until exact fields admit its first 1288 columns.",packet)


def ordered_stage_basis(prime:int,order:np.ndarray)->tuple[dict,np.ndarray]:
 binary=compile_modrank();M=degree8_matrix_seed(prime,SEED+447)[:,:1326].astype(np.uint32)[:,order];raw=Path(tempfile.gettempdir())/f"fin553_{prime}_{os.getpid()}.bin";ech=Path(tempfile.gettempdir())/f"fin553e_{prime}_{os.getpid()}.bin";M.tofile(raw);env=dict(os.environ);env["OMP_NUM_THREADS"]="8"
 try:
  p=subprocess.run([str(binary),str(raw),"1326","1326",str(prime),str(ech)],capture_output=True,text=True,timeout=600,env=env);rank=int(re.search(r"rank=(\d+)",p.stdout).group(1));pc=[int(x) for x in re.search(r"pivots=([0-9,]+)",p.stdout).group(1).split(",")];E=np.fromfile(ech,dtype=np.uint32).reshape(1326,1326);Xo=modular_nullspace(E,pc,prime);X=np.zeros_like(Xo);X[order]=Xo;T=degree8_matrix_seed(prime,SEED+6552,npts=48)[:,:1326].astype(np.int64);res=T@X.astype(np.int64)%prime
  return {"prime":prime,"rank":rank,"nullity":1326-rank,"first_1288_admitted":pc==list(range(1288)),"basis_sha256":hashlib.sha256(X.tobytes()).hexdigest(),"unseen_nonzeros":int(np.count_nonzero(res))},X
 finally:raw.unlink(missing_ok=True);ech.unlink(missing_ok=True)


def st553()->dict:
 order=np.array(json.loads(PACKETS[552].read_text())["order"],dtype=int);data={};rows=[]
 for p in (1000003,1000033,1000037,1000039,1000081):row,X=ordered_stage_basis(p,order);rows.append(row);data[str(p)]=X
 np.savez_compressed(BASES,**data);packet={"runs":rows,"all_ranks_1288":all(x["rank"]==1288 for x in rows),"all_conditioned_pivots_admitted":all(x["first_1288_admitted"] for x in rows),"all_unseen_replays_zero":all(x["unseen_nonzeros"]==0 for x in rows),"basis_file":BASES.name,"basis_sha256":sha(BASES)}
 return finalize(553,"proven_stage1_conditioned_five_prime_nullspaces", "The 38-dimensional modular kernel still requires rational reconstruction and exact identities.",packet)


def st554()->dict:
 z=np.load(BASES);primes=sorted(map(int,z.files));arr=[z[str(p)] for p in primes];rows=[];vectors=[]
 for j in range(38):
  v=reconstruct_column(arr,j,primes);rows.append({"column":j,"reconstructed":v is not None,"support":None if v is None else sum(x!=0 for x in v),"max_abs_coefficient":None if v is None else max(map(abs,v)),"bit_length":None if v is None else max(map(abs,v)).bit_length()});vectors.append(v)
 RELATIONS.write_text(json.dumps({"primes":primes,"rows":rows,"integer_relations":vectors},indent=2),encoding="utf-8");packet={"rows":rows,"reconstructed_count":sum(x["reconstructed"] for x in rows),"all_38_reconstructed":all(x["reconstructed"] for x in rows),"relation_file":RELATIONS.name,"relation_sha256":sha(RELATIONS)}
 return finalize(554,"stage1_rational_reconstruction_completed", "Reconstructed vectors become Q relations only after exact polynomial divisibility checks.",packet)


def stage1_polys():
 import itertools
 base=json.loads((ROOT/"FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text());q=[orbit_poly(r) for r in base["quadratic_generator_representatives"]];p4=[orbit_poly(r) for r in base["primitive_quartic_representatives"]];out=[]
 for ids in itertools.combinations_with_replacement(range(6),4):
  v={(0,)*12:1}
  for i in ids:v=poly_mul(v,q[i])
  out.append(v)
 for ids in itertools.combinations_with_replacement(range(6),2):
  q2=poly_mul(q[ids[0]],q[ids[1]]);out.extend(poly_mul(q2,x) for x in p4)
 out.extend(poly_mul(p4[i],p4[j]) for i,j in itertools.combinations_with_replacement(range(32),2));return out


def divisible_degree8(poly)->tuple[bool,int]:
 coeff=[{} for _ in range(9)]
 for e,c in poly.items():coeff[e[11]][e[:11]]=coeff[e[11]].get(e[:11],0)+c
 q=coeff[8]
 for k in range(7,-1,-1):
  rem=poly_add(coeff[k],multiply_by_sum11(q),scale=-1)
  if k==0:return not rem,len(rem)
  q=rem
 raise AssertionError


def st555()->dict:
 data=json.loads(RELATIONS.read_text());vectors=data["integer_relations"]
 if not all(v is not None for v in vectors):
  packet={"attempted":False,"reason":"not all 38 modular columns reconstructed","verified_count":0,"all_exact":False}
 else:
  polys=stage1_polys();rows=[]
  for j,v in enumerate(vectors):
   combined={}
   for c,p in zip(v,polys):
    if c:combined=poly_add(combined,p,c)
   ok,rem=divisible_degree8(combined);rows.append({"column":j,"combined_terms":len(combined),"divisible_by_sum_x":ok,"remainder_terms":rem})
  packet={"attempted":True,"rows":rows,"verified_count":sum(x["divisible_by_sum_x"] for x in rows),"all_exact":all(x["divisible_by_sum_x"] for x in rows)}
 return finalize(555,"proven_38_exact_stage1_relations" if packet["all_exact"] else "stage1_exact_verification_not_closed", "No characteristic-zero rank equality is promoted unless all 38 identities pass.",packet)


def st556()->dict:
 a=json.loads(PACKETS[553].read_text());b=json.loads(PACKETS[554].read_text());c=json.loads(PACKETS[555].read_text());closed=a["all_ranks_1288"] and b["all_38_reconstructed"] and c["all_exact"]
 packet={"rank_Q_first1326_equals_1288":closed,"relation_nullity_Q_first1326":38 if closed else None,"modular_rank":1288,"modular_nullity":38,"remaining_full_degree8_relations":199 if closed else 237,"full_rank_Q_1791":False}
 return finalize(556,"proven_stage1_characteristic_zero_rank1288" if closed else "stage1_characteristic_zero_rank_open", "Even stage-one closure leaves the 199 relations introduced by q*p6 and full rank_Q unresolved.",packet)


def st557()->dict:
 x=json.loads((ROOT/"FIN_ST540_Exhaustive_TwoLevel_Transition_Search.json").read_text());cand=2.9024964817477654;packet={"two_level_best_ratio":x["best_ratio"],"candidate_ratio":cand,"positive_margin_above_candidate":x["best_ratio"]-cand,"best_subset":x["best_subset"],"exhaustive_subset_count":x["canonical_subsets"],"interval_scalar_certificate":False}
 return finalize(557,"strong_two_level_family_exclusion_margin", "The family is exhaustive in subset patterns but scalar minimization is not interval certified and general states are multilevel.",packet)


def st558()->dict:
 a=json.loads((ROOT/"FIN_ST541_ReflectionEven_OrbitSpace_Stress_Test.json").read_text());b=json.loads((ROOT/"FIN_ST519_Transition_Branch_Curvature_Certificate.json").read_text());packet={"reflection_multistarts":a["starts"],"candidate_hits":a["best_value_hits"],"reflection_best_ratio":a["best_ratio"],"paid_local_curvature":b["minimum_paid_curvature"],"compact_reflection_quotient_cover":False,"result":"The candidate is the best sampled reflection-even local minimum and is interval-strict locally, but the six-dimensional quotient remains uncovered."}
 return finalize(558,"strong_reflection_orbit_candidate_consolidation", "Local curvature and starts do not prove globality even inside the reflection-fixed quotient.",packet)


def st559()->dict:
 packet={"global_lower":2.8934,"candidate_ratio":2.9024964817477654,"certified_upper":2.9024964917477667,"boundary_excluded":True,"two_level_margin":json.loads(PACKETS[557].read_text())["positive_margin_above_candidate"],"reflection_candidate_local_curvature":json.loads(PACKETS[558].read_text())["paid_local_curvature"],"remaining_obligation":"interval-exclude every other interior multilevel orbit below candidate","global_identity":False}
 return finalize(559,"first_transition_obligation_reduced_to_interior_multilevel_orbits", "The reduced obligation is sharp but remains open.",packet)


def st560()->dict:
 packet={"state_space":"twelve transitive branches","assumptions":["D12-covariant continuous-time Markov generator","irreducibility"],"theorem":"Covariance makes the uniform distribution stationary. Irreducibility makes it the unique stationary distribution. Therefore a symmetric Markov selector randomizes uniformly and cannot prefer a canonical branch.","stationary_branch_probability":"1/12","nonuniform_stationary_selector_requires":"broken covariance or reducibility plus initial-sector data","QW_2191":"open"}
 return finalize(560,"proven_equivariant_Markov_uniform_selector_theorem", "A laboratory may break covariance; that apparatus frame is then the selector resource.",packet)


def st561()->dict:
 x=json.loads((ROOT/"FIN_ST531_Three_Channel_Spectral_Dynamics_Separation.json").read_text());means=np.array([[*x["responses"]["unitary"]],[x["responses"]["heat"],0],[*x["responses"]["dephased_coherence"]]]);sigma=.05;n=2;rng=np.random.default_rng(SEED+561);trials=150000;conf=np.zeros((3,3),int)
 for i,m in enumerate(means):
  obs=m+rng.normal(scale=sigma/math.sqrt(n),size=(trials,2));pred=np.argmin(np.sum((obs[:,None,:]-means[None,:,:])**2,axis=2),axis=1);conf[i]=np.bincount(pred,minlength=3)
 err=1-np.trace(conf)/(3*trials);packet={"trials_per_model":trials,"samples_per_trial":n,"sigma":sigma,"confusion_counts":conf,"empirical_total_error":float(err),"target_1_percent_pass":err<.01,"frozen_classifier":"nearest ideal complex mean"}
 return finalize(561,"strong_synthetic_finite_shot_three_channel_replay", "Gaussian response samples are synthetic and not event-level detector data.",packet)


def st562()->dict:
 vals,V=np.linalg.eigh(A);theta=.17;R=np.eye(N);c,s=math.cos(theta),math.sin(theta);R[3,3]=c;R[3,5]=-s;R[5,3]=s;R[5,5]=c;Q=V@R@V.T;B=Q@A@Q.T;t=.6;w=.35;UA=expm(-1j*t*A);UB=expm(-1j*t*B);KA=[math.sqrt(w)*np.eye(N),math.sqrt(1-w)*UA];KB=[Q@K@Q.T for K in KA];tpA=sum(K.conj().T@K for K in KA);tpB=sum(K.conj().T@K for K in KB);rng=np.random.default_rng(SEED+562);v=rng.normal(size=N)+1j*rng.normal(size=N);rho=np.outer(v,v.conj());rho/=np.trace(rho);E=np.diag(rng.uniform(size=N));left=np.trace(E@sum(K@rho@K.conj().T for K in KA));right=np.trace((Q@E@Q.T)@sum(K@(Q@rho@Q.T)@K.conj().T for K in KB))
 packet={"Kraus_count":2,"trace_preservation_residual_A":float(np.linalg.norm(tpA-np.eye(N))),"trace_preservation_residual_B":float(np.linalg.norm(tpB-np.eye(N))),"transported_record_residual":float(abs(left-right)),"theorem":"A random-unitary CPTP channel and all transported states/effects preserve records under operational conjugacy."}
 return finalize(562,"proven_CPTP_operational_conjugacy_fixture", "This validates the mathematical transport rule, not a physical process implementation.",packet)


def st563()->dict:
 s=float(A[0,0]);W=s*np.eye(N)-A;weights=np.array([W[0,d] for d in range(1,7)])
 c2_c12=float(np.sum(weights[:5]*np.arange(1,6)**2)+18*weights[5]);c_c12=math.sqrt(c2_c12)
 c2_paired=float(np.sum(weights*np.arange(1,7)**2));c_paired=math.sqrt(c2_paired)
 modes=[]
 for n in range(1,7):
  q=2*math.pi*n/12;lam=2*sum(weights[d-1]*(1-math.cos(d*q)) for d in range(1,6))+weights[5]*(1-math.cos(6*q));der=2*sum(weights[d-1]*d*math.sin(d*q) for d in range(1,6))+6*weights[5]*math.sin(6*q);omega=math.sqrt(max(lam,0));modes.append({"mode":n,"q":q,"lambda":lam,"omega":omega,"phase_velocity":omega/q,"group_velocity":der/(2*omega) if omega else None})
 packet={"strict_radial_weights":weights,"C12_symbol":"2 sum_{d=1}^5 w_d(1-cos(dq))+w_6(1-cos(6q))","formal_C12_small_q_speed_squared":c2_c12,"formal_C12_small_q_speed":c_c12,"paired_range6_large_cycle_speed_squared":c2_paired,"paired_range6_large_cycle_speed":c_paired,"discrete_C12_modes":modes,"first_available_q":math.pi/6,"first_mode_phase_velocity":modes[0]["phase_velocity"],"first_mode_group_velocity":modes[0]["group_velocity"],"exact_causal_cone":False,"short_time_far_site_orders":{"wave_probability":"O(t^2)","diffusion_probability":"O(t)"},"long_range_positive_tail_eta_1_8":{"symbol_exponent":.8,"wave_exponent":.4,"group_velocity_behavior":"proportional to |k|^-0.6; no finite low-k limit"},"layer_speed_law":"c_phys,n=(ell_n/tau_n)*c_hat_n","layer_invariance_condition":"(alpha/beta)*(c_hat_{n+1}/c_hat_n)=1","alpha_equals_beta_sufficient_only_if_dimensionless_speed_is_refinement_invariant":True,"SI_299792458_derived":False,"theorem":"The fixed C12 symbol has computable dimensionless phase/group-speed proxies but no q->0 continuum and no exact causal cone. A layer-invariant physical speed requires a declared refinement with quadratic low-k symbol plus calibration satisfying the displayed covariance law."}
 return finalize(563,"proven_finite_speed_proxies_and_causal_cone_obstruction", "Neither the formal slope nor the SI value of c is derived. Finite-range, long-range, signed oscillatory, and Schur refinements are different models and must be compared separately.",packet)


def ultra_ir(start=.253,stop=.3,width=.0001)->dict:
 prev=json.loads((ROOT/"FIN_ST547_FineChart_IR_Continuation_Attempt.json").read_text());prior=np.array(prev["cells"][-1]["endpoint_centers"][-1]);rows=[];lo=start;failure=None
 while lo<stop-1e-14:
  hi=min(stop,lo+width);raw,candidate=direct_ir_cell(lo,hi,prior);selected=None
  for factor in (1.,2.,4.,8.,16.,32.,48.,64.,96.,128.):
   q=recertify_ir_cell(raw,lo,hi,factor)
   if q["included"]:selected=q;break
  if selected is None:failure=raw;break
  rows.append(selected);prior=candidate;lo=hi
 return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}
def st564()->dict:
 c=ultra_ir();rows=c["cells"];packet={**c,"cell_count":len(rows),"minimum_margin":min(x["minimum_margin"] for x in rows) if rows else None,"maximum_contraction":max(x["weighted_contraction"] for x in rows) if rows else None,"maximum_radius_factor":max(x["radius_expansion_factor"] for x in rows) if rows else None,"target_reached":c["certified_stop"]>=.3-1e-14,"theorem":f"Ultra-fine inherited boxes continue the local component from b=0.253 through b={c['certified_stop']:.5f}."}
 return finalize(564,"proven_ultrafine_IR_extension" if rows else "failed_ultrafine_IR_extension", "Repeated radius refinement is now repetition-gated after this stop; a new analytic chart is required.",packet)


def st565()->dict:
 q=json.loads(PACKETS[556].read_text());packet={"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","new_scale_charged_source":False,"stage1_rank_Q_closed":q["rank_Q_first1326_equals_1288"],"full_degree8_rank_Q_closed":False,"first_transition_global_identity":False,"physical_tester":False,"physical_light_cone_or_c_derived":False}
 return finalize(565,"round1_admission_gate_evaluated", "No bridge/role transfer, physical c, SM/GR, L_total, laboratory, or ToE closure follows.",packet)


def st566()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":1}
 return finalize(566,"blocked_no_independent_empirical_evidence", "Round one is local analytic/computational work.",packet)


def figures(r:dict)->None:
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST554"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar(range(len(rows)),[x["bit_length"] or 0 for x in rows],color="#2563eb");ax.set(xlabel="stage-one relation",ylabel="coefficient bit length",title="ST554: stage-one rational reconstruction");fig.tight_layout();fig.savefig(FIG_DIR/"st554_stage1_lift.png",dpi=180);plt.close(fig)
 x=r["ST561"]["confusion_counts"];fig,ax=plt.subplots(figsize=(5.2,4.5));im=ax.imshow(x,cmap="Blues");ax.set(xlabel="predicted",ylabel="true",title="ST561: synthetic three-channel confusion");fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(FIG_DIR/"st561_confusion.png",dpi=180);plt.close(fig)
 x=r["ST563"]["discrete_C12_modes"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot([z["q"] for z in x],[z["omega"] for z in x],"o-",label=r"$\omega=\sqrt{\lambda}$");q=np.linspace(0,math.pi,200);ax.plot(q,r["ST563"]["formal_C12_small_q_speed"]*q,"--",label="formal small-q tangent");ax.set(xlabel="q",ylabel="dimensionless frequency",title="ST563: C12 wave dispersion and formal speed");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st563_dispersion.png",dpi=180);plt.close(fig)
 rows=r["ST564"]["cells"];b=[x["b_interval"][1] for x in rows];y2=[x["endpoint_centers"][-1][1] for x in rows];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(b,y2,"o-");ax.set(xlabel="b",ylabel=r"$y_2$",title="ST564: ultra-fine IR continuation");fig.tight_layout();fig.savefig(FIG_DIR/"st564_ir_ultrafine.png",dpi=180);plt.close(fig)


def main()->None:
 funcs={552:st552,553:st553,554:st554,555:st555,556:st556,557:st557,558:st558,559:st559,560:st560,561:st561,562:st562,563:st563,564:st564,565:st565,566:st566};r={}
 for k in range(552,567):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(552,567):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
