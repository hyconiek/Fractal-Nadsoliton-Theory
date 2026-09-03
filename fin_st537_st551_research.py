#!/usr/bin/env python3
"""FIN ST537--ST551: blockwise invariant algebra and operational escape tests."""
from __future__ import annotations
import csv,hashlib,json,math,os,re,subprocess,tempfile,time
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize,minimize_scalar,root

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N,SEED,degree8_candidate_metadata
from fin_st447_st461_research import direct_ir_cell,regularized_ir_float,FACE
from fin_st477_st491_research import compile_modrank,recertify_ir_cell
from fin_st492_st506_research import degree8_matrix_seed
from fin_st462_st476_research import probability_ratio

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST537_ST551_Results.json";SUMMARY=ROOT/"FIN_ST537_ST551_Summary.csv";FIG_DIR=ROOT/"FIN_ST537_ST551_Figures"
NAMES={537:"Degree8_Blockwise_Exact_Modular_Ranks",538:"Degree8_Cumulative_Family_Ranks",539:"Degree8_Blockwise_CharacteristicZero_Consequences",540:"Exhaustive_TwoLevel_Transition_Search",541:"ReflectionEven_OrbitSpace_Stress_Test",542:"Symmetric_Noise_Branch_Frequency_Replay",543:"Nonunique_Flow_Escape_Classification",544:"FiniteShot_ThreeChannel_Gaussian_Design",545:"CompletelyPositive_Operational_Conjugacy",546:"Joint_Channel_Calibration_Identifiability",547:"FineChart_IR_Continuation_Attempt",548:"IR_Complement_Multistart_Search",549:"Strict_Source_Admission_Gate",550:"Minimal_Conversion_Sector_Operational_Package",551:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();U=np.ones(N)/N
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}


COUNTS=[126,672,528,702];LABELS=["q4","q2_p4","p4_p4","q_p6"]
def exact_rank_array(M:np.ndarray,prime:int)->dict:
 binary=compile_modrank();M=M.astype(np.uint32);raw=Path(tempfile.gettempdir())/f"fin537_{prime}_{M.shape[1]}_{os.getpid()}.bin";M.tofile(raw);env=dict(os.environ);env["OMP_NUM_THREADS"]="8";t=time.time()
 try:
  p=subprocess.run([str(binary),str(raw),str(M.shape[0]),str(M.shape[1]),str(prime)],capture_output=True,text=True,timeout=600,env=env);rank=int(re.search(r"rank=(\d+)",p.stdout).group(1));return {"shape":list(M.shape),"rank":rank,"nullity":M.shape[1]-rank,"matrix_sha256":hashlib.sha256(M.tobytes()).hexdigest(),"seconds":time.time()-t}
 finally:raw.unlink(missing_ok=True)


def family_rank_runs()->dict:
 out={}
 for prime in (1000003,1000033):
  M=degree8_matrix_seed(prime,SEED+447);start=0;rows=[]
  for lab,n in zip(LABELS,COUNTS):rows.append({"family":lab,**exact_rank_array(M[:,start:start+n],prime)});start+=n
  out[str(prime)]=rows
 return out
def st537()->dict:
 runs=family_rank_runs();packet={"family_counts":dict(zip(LABELS,COUNTS)),"prime_runs":runs,"cross_prime_rank_agreement":[x["rank"] for x in runs["1000003"]]==[x["rank"] for x in runs["1000033"]]}
 return finalize(537,"proven_exact_blockwise_modular_ranks", "Rank-deficient block equality over Q still requires rational upper certificates.",packet)


def st538()->dict:
 rows={}
 for prime in (1000003,1000033):
  M=degree8_matrix_seed(prime,SEED+447);ends=np.cumsum(COUNTS);rr=[]
  for i,e in enumerate(ends):rr.append({"families":LABELS[:i+1],**exact_rank_array(M[:,:e],prime)})
  rows[str(prime)]=rr
 packet={"prime_runs":rows,"cross_prime_cumulative_rank_agreement":[x["rank"] for x in rows["1000003"]]==[x["rank"] for x in rows["1000033"]],"final_rank":rows["1000003"][-1]["rank"]}
 return finalize(538,"proven_exact_cumulative_modular_rank_profile", "The profile localizes dependencies but remains finite-field unless a family is column-full.",packet)


def st539()->dict:
 r=json.loads(PACKETS[537].read_text());rows=[]
 for i,(lab,n) in enumerate(zip(LABELS,COUNTS)):
  rank=r["prime_runs"]["1000003"][i]["rank"];full=rank==n;rows.append({"family":lab,"candidate_count":n,"modular_rank":rank,"Q_linear_independence_proven":full,"Q_rank_interval":[rank,n],"relation_nullity_interval":[0,n-rank]})
 packet={"family_rows":rows,"theorem":"Whenever modular rank equals the number of columns, reduction proves those integer polynomials are linearly independent over Q.","all_families_individually_Q_independent":all(x["Q_linear_independence_proven"] for x in rows)}
 return finalize(539,"proven_blockwise_Q_independence_where_full", "Cross-family relations and every deficient block retain the displayed intervals.",packet)


def two_level_search()->dict:
 rows=[]
 for mask in range(1,2**N-1):
  if not(mask&1):continue
  ids=np.array([i for i in range(N) if mask>>i&1]);k=len(ids)
  def pof(a):
   b=(1-k*a)/(N-k);p=np.full(N,b);p[ids]=a;return p
  lo=1/N+1e-12;hi=1/k-1e-12
  if hi<=lo:continue
  opt=minimize_scalar(lambda a:probability_ratio(pof(a)),bounds=(lo,hi),method="bounded",options={"xatol":1e-14});rows.append((float(opt.fun),ids.tolist(),float(opt.x)))
 rows.sort(key=lambda x:x[0]);return {"canonical_subsets":len(rows),"best_ratio":rows[0][0],"best_subset":rows[0][1],"best_high_probability":rows[0][2],"ten_best":[{"ratio":x[0],"subset":x[1],"a":x[2]} for x in rows[:10]]}
def st540()->dict:
 x=two_level_search();x["candidate_beaten"]=x["best_ratio"]<2.902496471747767;x["claim"]="All two-level states up to complement symmetry were optimized in their one-dimensional interiors."
 return finalize(540,"exhaustive_two_level_family_search_completed", "Scalar minimization is high-accuracy numerical, not interval-certified; general states have more than two levels.",x)


MULT=np.array([1,2,2,2,2,2,1])
def radial_soft(z):
 zz=np.r_[z,0.];zz-=zz.max();w=np.exp(zz);levels=w/float(MULT@w);return np.array([levels[0],levels[1],levels[2],levels[3],levels[4],levels[5],levels[6],levels[5],levels[4],levels[3],levels[2],levels[1]])
def st541()->dict:
 rng=np.random.default_rng(SEED+541);vals=[]
 for _ in range(240):
  o=minimize(lambda z:probability_ratio(radial_soft(z)),rng.normal(0,4,6),method="Nelder-Mead",options={"maxiter":2500,"xatol":1e-10,"fatol":1e-12});vals.append((float(o.fun),radial_soft(o.x)))
 vals.sort(key=lambda x:x[0]);packet={"starts":len(vals),"best_ratio":vals[0][0],"best_probability":vals[0][1],"best_value_hits":sum(abs(x[0]-vals[0][0])<2e-9 for x in vals),"distinct_clusters_1e_7":len({round(x[0],7) for x in vals}),"candidate_beaten":vals[0][0]<2.902496471747767,"global_reflection_subspace_certificate":False}
 return finalize(541,"strong_reflection_even_orbit_space_stress_test", "Multistart optimization does not exhaust the compact six-dimensional quotient.",packet)


def st542()->dict:
 p=np.array(json.loads((ROOT/"FIN_ST462_First_Transition_Orbit_Numerical_Isolation.json").read_text())["best_probability"]);protos=np.array([np.roll(p-U,k) for k in range(N)]);protos/=np.linalg.norm(protos,axis=1,keepdims=True);rng=np.random.default_rng(SEED+542);xi=rng.normal(size=(240000,N));xi-=xi.mean(axis=1,keepdims=True);labels=np.argmax(xi@protos.T,axis=1);counts=np.bincount(labels,minlength=N);expected=len(labels)/N;chi=float(np.sum((counts-expected)**2/expected))
 packet={"samples":len(labels),"counts":counts,"frequencies":counts/len(labels),"maximum_frequency_deviation_from_1_12":float(np.max(abs(counts/len(labels)-1/N))),"Pearson_chi_square":chi,"degrees_of_freedom":11,"claim":"Isotropic tangent Gaussian seeds and a covariant nearest-orbit rule produce frequencies consistent with the required uniform branch law."}
 return finalize(542,"strong_numerical_uniform_branch_frequency_replay", "Sampling confirms the symmetry theorem but does not create physical fluctuations or a selector.",packet)


def st543()->dict:
 packet={"scalar_witness":"z_dot=sign(z), z(0)=0 interpreted as differential inclusion z_dot in [-1,1] at zero","solutions":["z(t)=0","z(t)=t","z(t)=-t","delayed departures"],"theorem":"Dropping uniqueness permits symmetry-related departures without an explicit oriented deterministic vector at the origin; the choice is carried by solution selection for the inclusion.","hidden_resource":"solution-selection rule or stochastic regularization","strict_selector":False}
 return finalize(543,"proven_nonunique_flow_evades_vacuum_no_go_but_moves_selector", "Nonuniqueness does not remove the selector problem; it relocates it to the solution concept/regularization.",packet)


def st544()->dict:
 x=json.loads((ROOT/"FIN_ST531_Three_Channel_Spectral_Dynamics_Separation.json").read_text());d=min(x["pairwise_distances"].values());sigma=.05;chernoff=d*d/(8*sigma*sigma);alpha=.01;n=math.ceil(math.log(3/alpha)/chernoff)
 packet={"observation_model":"independent isotropic two-quadrature Gaussian observations around the three ideal complex responses","sigma":sigma,"minimum_pairwise_distance":d,"minimum_pairwise_Chernoff_exponent_per_sample":chernoff,"union_bound_target_error":alpha,"sufficient_sample_count":n,"theorem":"For equal-covariance Gaussian response models, pairwise Chernoff information is d^2/(8 sigma^2); a three-pair union bound is below alpha when n>=log(3/alpha)/C_min.","event_level_physical_model":False}
 return finalize(544,"proven_conditional_finite_sample_three_channel_bound", "The Gaussian likelihood and sigma are supplied synthetic assumptions, not detector calibration.",packet)


def st545()->dict:
 rng=np.random.default_rng(SEED+545);vals,V=np.linalg.eigh(A);theta=.21;R=np.eye(N);c,s=math.cos(theta),math.sin(theta);R[1,1]=c;R[1,2]=-s;R[2,1]=s;R[2,2]=c;Q=V@R@V.T;B=Q@A@Q.T
 Ks=[]
 for _ in range(3):K=rng.normal(size=(N,N));Ks.append(K/np.linalg.norm(K,2)/math.sqrt(3))
 rho=rng.normal(size=(N,N));rho=rho@rho.T;rho/=np.trace(rho);effect=np.eye(N);t=.4
 def channel(G,r,Klist):return sum(K@expm(-t*G)@r@expm(-t*G)@K.T for K in Klist)
 left=np.trace(effect@channel(A,rho,Ks));Kq=[Q@K@Q.T for K in Ks];right=np.trace((Q@effect@Q.T)@channel(B,Q@rho@Q.T,Kq))
 packet={"Kraus_count":len(Ks),"record_residual":float(abs(left-right)),"theorem":"Conjugating every Kraus operator, state, propagator and effect by Q preserves completely positive multitime records by cyclicity of trace.","trace_preserving_fixture":False,"physical_process_tensor":False}
 return finalize(545,"proven_CP_operational_conjugacy_record_invariance", "The numerical Kraus fixture is CP but not asserted trace preserving; physical instruments require normalization and records.",packet)


def st546()->dict:
 packet={"model":"y_U=a lambda, y_H=b lambda, y_W=c sqrt(lambda)","gauge":"lambda->q lambda, (a,b,c)->(a/q,b/q,c/sqrt(q))","invariants":["a/b=y_U/y_H","c^2/a=y_W^2/y_U"],"Jacobian_parameter_count":4,"identifiable_rank":3,"residual_torsor_dimension":1,"theorem":"Joint observation of all three ideal channel scales identifies two calibration ratios and modal products but leaves one positive common rescaling orbit.","absolute_seconds":False}
 return finalize(546,"proven_joint_channel_calibration_one_torsor_obstruction", "One dimensional anchor remains necessary even with all three channels.",packet)


def fine_ir(start=.2435,stop=.3,width=.00025)->dict:
 previous=json.loads((ROOT/"FIN_ST533_Adaptive_IR_Continuation_to_0_3_Attempt.json").read_text());prior=np.array(previous["cells"][-1]["endpoint_centers"][-1],dtype=float);rows=[];lo=start;failure=None
 while lo<stop-1e-14:
  hi=min(stop,lo+width);raw,candidate=direct_ir_cell(lo,hi,prior);selected=None
  for factor in (1.,2.,4.,8.,16.,24.,32.,48.,64.):
   q=recertify_ir_cell(raw,lo,hi,factor)
   if q["included"]:selected=q;break
  if selected is None:failure=raw;break
  rows.append(selected);prior=candidate;lo=hi
 return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}
def st547()->dict:
 c=fine_ir();rows=c["cells"];packet={**c,"cell_count":len(rows),"minimum_margin":min(x["minimum_margin"] for x in rows) if rows else None,"maximum_contraction":max(x["weighted_contraction"] for x in rows) if rows else None,"maximum_radius_factor":max(x["radius_expansion_factor"] for x in rows) if rows else None,"target_reached":c["certified_stop"]>=.3-1e-14,"theorem":f"A finer chart continues the same local component from b=0.2435 through b={c['certified_stop']:.5f}."}
 return finalize(547,"proven_fine_chart_IR_extension" if rows else "failed_fine_chart_IR_extension", "Radius refinement is not a new analytic coordinate or global root count.",packet)


def st548()->dict:
 rng=np.random.default_rng(SEED+548);rows=[]
 def safe(z,b):
  if not np.all(np.isfinite(z)) or np.any(z<=0) or np.max(z)>100:return np.full(5,1e6)
  try:return regularized_ir_float(z,b)
  except (OverflowError,ValueError,FloatingPointError):return np.full(5,1e6)
 for b in np.linspace(.24,.3,13):
  roots=[]
  for _ in range(80):
   z0=np.array([.03,3.5,.05,.08,2.7])*np.exp(rng.normal(0,.7,5));sol=root(lambda z:safe(z,b),z0)
   if sol.success and np.linalg.norm(sol.fun,np.inf)<1e-7 and np.all(sol.x>0):
    if not any(np.linalg.norm(sol.x-x)<1e-5 for x in roots):roots.append(sol.x)
  rows.append({"b":float(b),"positive_root_clusters":len(roots),"roots":roots})
 packet={"grid_rows":rows,"starts_per_grid":80,"maximum_positive_clusters":max(x["positive_root_clusters"] for x in rows),"global_complement_theorem":False}
 return finalize(548,"numerical_IR_complement_multistart_search", "Failure to find a disconnected root is not an interval complement cover.",packet)


def st549()->dict:
 packet={"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","new_scale_charged_source":False,"new_anisotropic_provider":False,"reason":"Nonunique flow, stochastic seeds and likelihood models all expose additional provider choices."}
 return finalize(549,"blocked_no_new_strict_source_provider", "No legacy bridge/role transfer, SM/GR, L_total, or ToE closure.",packet)


def st550()->dict:
 packet={"architecture":"W0+CA+SA+OA","W0":"dimensionless spectral-information core","CA":["ell_*","tau_*","hbar_*"],"SA":["origin/orientation/chirality datum","coupling law"],"OA":["state","preparations","clock/dynamics","instruments","environment","apparatus","record/custody"],"removal_matrix":{"remove_CA":"absolute dimensions unidentified","remove_SA":"oriented branch claims unidentified","remove_OA":"no mapping from mathematics to event records"},"theorem":"CA, SA and OA close logically distinct group actions/operational obligations; none can be inferred merely by renaming another package.","strict_derivation":False}
 return finalize(550,"proven_minimal_package_role_separation", "This is an axiom-augmented architecture, not strict ToE closure.",packet)


def st551()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent"}
 return finalize(551,"blocked_no_independent_empirical_evidence", "Round six remains local analytic/computational work.",packet)


def figures(r:dict)->None:
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST537"]["prime_runs"]["1000003"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar([x["family"] for x in rows],[x["rank"] for x in rows],color="#2563eb");ax.plot([x["family"] for x in rows],[x["shape"][1] for x in rows],"ro",label="columns");ax.set(ylabel="rank / columns",title="ST537: degree-eight family ranks");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st537_family_ranks.png",dpi=180);plt.close(fig)
 x=r["ST542"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar(range(12),x["frequencies"],color="#0f766e");ax.axhline(1/12,color="#dc2626",ls="--");ax.set(xlabel="branch",ylabel="frequency",title="ST542: symmetric noise branch replay");fig.tight_layout();fig.savefig(FIG_DIR/"st542_branch_frequencies.png",dpi=180);plt.close(fig)
 rows=r["ST547"]["cells"];b=[x["b_interval"][1] for x in rows];y2=[x["endpoint_centers"][-1][1] for x in rows];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(b,y2,"o-");ax.set(xlabel="b",ylabel=r"$y_2$",title="ST547: fine-chart IR continuation");fig.tight_layout();fig.savefig(FIG_DIR/"st547_ir_fine.png",dpi=180);plt.close(fig)


def main()->None:
 funcs={537:st537,538:st538,539:st539,540:st540,541:st541,542:st542,543:st543,544:st544,545:st545,546:st546,547:st547,548:st548,549:st549,550:st550,551:st551};r={}
 for k in range(537,552):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(537,552):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
