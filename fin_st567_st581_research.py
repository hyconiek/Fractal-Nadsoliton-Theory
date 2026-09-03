#!/usr/bin/env python3
"""FIN ST567--ST581: refinement-speed classes, causal tails, and sparse algebra audit."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar,root

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N,SEED
from fin_st447_st461_research import regularized_ir_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST567_ST581_Results.json";SUMMARY=ROOT/"FIN_ST567_ST581_Summary.csv";FIG_DIR=ROOT/"FIN_ST567_ST581_Figures"
NAMES={567:"FiniteRange_Refinement_Dispersion_Limit",568:"Positive_LongRange_Fractional_Dispersion",569:"Signed_Oscillatory_Refinement_Positivity_NoGo",570:"Layer_Covariant_Physical_Speed_Condition",571:"FiniteRange_LiebRobinson_Velocity_Proxy",572:"Exact_Refinement_WaveSpeed_Scaling",573:"ShortTime_CausalTail_Operational_Audit",574:"TwoLevel_Interval_Certificate_Preparation",575:"ReflectionQuotient_Cover_Stop",576:"Stage1_Modular_Sparse_Circuit_Search",577:"Stage1_CharacteristicZero_Method_Gate",578:"ThreeChannel_Likelihood_Robustness",579:"LogChart_IR_Conditioning_Audit",580:"Speed_Source_Selector_Admission_Gate",581:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();U=np.ones(N)/N
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}

s=float(A[0,0]);W=s*np.eye(N)-A;STRICT_W=np.array([W[0,d] for d in range(1,7)])
def paired_symbol(q,w):return 2*sum(float(w[d-1])*(1-math.cos(d*q)) for d in range(1,len(w)+1))

def st567()->dict:
 Ns=[24,48,96,192,384,768,1536];rows=[];exact_c=math.sqrt(float(np.sum(STRICT_W*np.arange(1,7)**2)))
 for n in Ns:
  q=2*math.pi/n;lam=paired_symbol(q,STRICT_W);rows.append({"N":n,"q1":q,"lambda1":lam,"phase_speed":math.sqrt(lam)/q,"relative_error_to_limit":abs(math.sqrt(lam)/q-exact_c)/exact_c})
 packet={"refinement":"paired finite range d<=6 with fixed weights","exact_small_q_speed":exact_c,"rows":rows,"convergence_ratio_last_two":rows[-1]["relative_error_to_limit"]/rows[-2]["relative_error_to_limit"],"theorem":"Taylor expansion gives lambda(q)=q^2 sum_d w_d d^2+O(q^4), so omega(q)/|q| converges to sqrt(sum_d w_d d^2).","strict_refinement_provenance":False,"exact_lattice_causal_cone":False}
 return finalize(567,"proven_finite_range_linear_dispersion_limit", "The paired finite-range continuation is an added refinement model; continuous-time lattice propagation still has analytic tails rather than an exact cone.",packet)


def longrange_rows(eta=1.8):
 Ns=[64,128,256,512,1024,2048,4096];rows=[]
 for n in Ns:
  d=np.arange(1,n//2+1);w=1/(1+d**eta);q=2*np.pi/n;lam=2*np.sum(w[:-1]*(1-np.cos(d[:-1]*q)))+w[-1]*(1-np.cos(d[-1]*q));rows.append({"N":n,"q1":q,"lambda1":float(lam),"wave_phase_speed":float(math.sqrt(lam)/q)})
 return rows
def st568()->dict:
 rows=longrange_rows();x=np.log([r["q1"] for r in rows[-5:]]);y=np.log([r["lambda1"] for r in rows[-5:]]);nu=float(np.polyfit(x,y,1)[0]);vgexp=nu/2-1
 packet={"positive_envelope":"w_d=1/(1+d^1.8)","rows":rows,"fitted_low_q_symbol_exponent_last5":nu,"Tauberian_expected_exponent":.8,"implied_wave_exponent":nu/2,"implied_group_velocity_exponent":vgexp,"finite_nonzero_low_q_speed":False,"theorem_scope":"numerical asymptotic audit consistent with fractional order eta-1"}
 return finalize(568,"strong_fractional_long_range_dispersion_evidence", "The positive envelope is not the signed strict kernel; a theorem requires Tauberian upper/lower bounds. The observed exponent is incompatible with relativistic linear wave dispersion.",packet)


def signed_spectrum(n):
 d=np.arange(1,n//2+1);w=np.cos(.18575*d+.1625)/(1+d**1.8);lam=[]
 for m in range(n):
  q=2*np.pi*m/n;v=2*np.sum(w[:-1]*(1-np.cos(d[:-1]*q)))+w[-1]*(1-np.cos(d[-1]*q));lam.append(float(v))
 return w,np.array(lam)
def st569()->dict:
 rows=[]
 for n in (12,24,48,96,192,384,768):
  w,l=signed_spectrum(n);rows.append({"N":n,"negative_weights":int(np.sum(w<0)),"minimum_laplacian_eigenvalue":float(l.min()),"negative_modes":int(np.sum(l<-1e-10)),"lambda_mode1":float(l[1])})
 first=next(r["N"] for r in rows if r["negative_modes"])
 packet={"literal_extension":"cos(0.18575d+0.1625)/(1+d^1.8) to all cyclic distances","rows":rows,"first_tested_indefinite_carrier":first,"theorem":"The literal sampled extension ceases to be a positive graph Laplacian in the tested family; already N=48 has negative modes.","universal_all_N_no_go":False}
 return finalize(569,"proven_finite_counterexample_to_naive_signed_strict_refinement", "A finite counterexample rejects this naive refinement as a universal positive generator; it does not exclude a different completion map or signed/Krein model.",packet)


def st570()->dict:
 c12=json.loads((ROOT/"FIN_ST563_Layer_Covariant_Effective_Speed_and_Causal_Cone_Audit.json").read_text());finite=json.loads(PACKETS[567].read_text());frac=json.loads(PACKETS[568].read_text())
 packet={"physical_speed_law":"c_phys,n=(ell_n/tau_n)*c_hat_n","invariance_law":"(alpha_n/beta_n)*(c_hat_{n+1}/c_hat_n)=1","C12_formal_speed":c12["formal_C12_small_q_speed"],"finite_range_limit_speed":finite["exact_small_q_speed"],"fractional_speed_limit":"diverges or has no finite nonzero limit","theorem":"Equal length/time scaling alpha=beta preserves physical speed only if the dimensionless dispersion slope is itself refinement invariant. Different admitted refinement classes fail that premise.","universal_layer_c_from_A12":False}
 return finalize(570,"proven_layer_speed_covariance_condition_and_refinement_dependence", "The condition is exact; FIN has not selected a refinement or dimensional anchor and therefore has not derived physical c.",packet)


def st571()->dict:
 def proxy(mu):return 2*sum(STRICT_W[d-1]*math.sinh(mu*d) for d in range(1,7))/mu
 opt=minimize_scalar(proxy,bounds=(1e-5,1.5),method="bounded");q=np.linspace(1e-6,math.pi,20000);der=np.array([2*sum(STRICT_W[d-1]*d*math.sin(d*x) for d in range(1,7)) for x in q]);lam=np.array([paired_symbol(x,STRICT_W) for x in q]);wave=np.abs(der)/(2*np.sqrt(np.maximum(lam,1e-300)))
 packet={"exponential_moment_proxy":"v_mu=2 sum_d w_d sinh(mu d)/mu","optimal_mu_in_scan":float(opt.x),"optimized_velocity_proxy":float(opt.fun),"maximum_sampled_unitary_symbol_velocity":float(np.max(abs(der))),"maximum_sampled_wave_group_velocity":float(np.max(wave)),"wave_low_q_limit":json.loads(PACKETS[567].read_text())["exact_small_q_speed"],"exact_Lieb_Robinson_theorem":False}
 return finalize(571,"computed_finite_range_LR_velocity_proxy", "The proxy is a declared exponential-moment bound candidate, not a fully proved many-body Lieb--Robinson constant or Lorentz cone.",packet)


def st572()->dict:
 packet={"exact_intertwining":"A_tilde R=R A","coarse_mode_eigenvalue":"lambda preserved","wave_frequency":"sqrt(lambda) preserved","label_wavenumber_under_N_to_2N":"q->q/2","dimensionless_phase_speed_scaling":"c_hat->2 c_hat under the raw mode label","physical_invariance_condition":"(ell_child/tau_child)=0.5*(ell_parent/tau_parent) for this labeling","theorem":"Exact spectral refinement naturality does not by itself preserve a dimensionless phase velocity because the wavenumber coordinate changes. Physical speed is preserved only after a compatible length/time transition cocycle is supplied.","canonical_transition_cocycle":False}
 return finalize(572,"proven_refinement_spectrum_speed_coordinate_distinction", "The factor two belongs to the declared Fourier labeling; alternative physical embeddings must state their coordinate map explicitly.",packet)


def st573()->dict:
 tvals=[1e-4,3e-4,1e-3,3e-3,1e-2];rows=[]
 from scipy.linalg import expm
 for t in tvals:
  Uop=expm(-1j*t*A);Hop=expm(-t*A);rows.append({"t":t,"far_wave_probability":float(abs(Uop[6,0])**2),"far_heat_probability":float(Hop[6,0])})
 sw=np.polyfit(np.log(tvals),np.log([x["far_wave_probability"] for x in rows]),1)[0];sh=np.polyfit(np.log(tvals),np.log([x["far_heat_probability"] for x in rows]),1)[0]
 packet={"source_vertex":0,"opposite_vertex":6,"rows":rows,"fitted_wave_power":float(sw),"fitted_heat_power":float(sh),"expected_orders":{"wave_probability":2,"heat_probability":1},"exact_zero_outside_cone":False}
 return finalize(573,"strong_short_time_far_tail_order_replay", "Floating power fits illustrate the analytic leading orders; the nonzero opposite coupling already proves absence of an exact finite cone.",packet)


def st574()->dict:
 old=json.loads((ROOT/"FIN_ST540_Exhaustive_TwoLevel_Transition_Search.json").read_text());margin=old["best_ratio"]-2.9024964817477654
 packet={"subset_orbits":old["canonical_subsets"],"numerical_best_ratio":old["best_ratio"],"candidate_margin":margin,"proposed_interval_partition":"adaptive scalar intervals around each subset optimum plus analytic endpoint entropy bounds","minimum_required_total_error_budget":margin/2,"interval_certificate_completed":False}
 return finalize(574,"two_level_interval_certificate_prepared_not_closed", "The prior exhaustive patterns remain numerical until all scalar intervals and endpoint regions are outward certified.",packet)


def st575()->dict:
 old=json.loads((ROOT/"FIN_ST541_ReflectionEven_OrbitSpace_Stress_Test.json").read_text());packet={"dimension":6,"multistarts":old["starts"],"stationary_clusters":old["distinct_clusters_1e_7"],"best_ratio":old["best_ratio"],"candidate_hits":old["best_value_hits"],"interval_cover_cells":0,"result":"No tractable six-dimensional outward cover is exported in this round.","global_reflection_quotient_theorem":False}
 return finalize(575,"bounded_stop_no_reflection_quotient_cover", "This is a method/resource stop, not evidence for another orbit or against candidate globality.",packet)


def st576()->dict:
 z=np.load(ROOT/"FIN_ST553_Stage1_Conditioned_FivePrime_Nullspaces.npz");X=z["1000003"].astype(np.int64);p=1000003;rng=np.random.default_rng(SEED+576);best=X.shape[0]+1;best_coeff=None;trials=160000
 for _ in range(trials//1000):
  C=rng.integers(-1,2,size=(X.shape[1],1000),dtype=np.int64);mask=np.any(C!=0,axis=0);Y=X@C[:,mask]%p;support=np.count_nonzero(Y,axis=0);j=int(np.argmin(support))
  if int(support[j])<best:best=int(support[j]);best_coeff=C[:,mask][:,j].copy()
 packet={"prime":p,"trials":trials,"coefficient_alphabet":[-1,0,1],"minimum_support_found":best,"coefficient_support":int(np.count_nonzero(best_coeff)),"best_combination_sha256":hashlib.sha256(best_coeff.astype(np.int8).tobytes()).hexdigest(),"exact_rational_relation":False}
 return finalize(576,"exact_modular_sparse_combination_search", "Small modular support under one prime may be a bad-prime artifact and is not a rational relation without multi-prime alignment and exact polynomial verification.",packet)


def st577()->dict:
 recon=json.loads((ROOT/"FIN_ST554_Stage1_Rational_Reconstruction.json").read_text());sparse=json.loads(PACKETS[576].read_text());packet={"five_prime_reconstructed_relations":recon["reconstructed_count"],"modular_minimum_support_search":sparse["minimum_support_found"],"new_structure_aware_Q_method_available":False,"decision":"Freeze direct CRT and random-combination escalation; next admissible algebra move requires representation-theoretic coordinates, p-adic lifting with height control, or a characteristic-zero CAS checkpoint."}
 return finalize(577,"stage1_characteristic_zero_lane_method_gated", "Failure of current methods is not nonexistence of 38 rational relations.",packet)


def st578()->dict:
 base=json.loads((ROOT/"FIN_ST561_FiniteShot_ThreeChannel_MonteCarlo_Replay.json").read_text());rows=[]
 for factor in (.5,1,1.5,2,3):
  sigma=.05*factor;d=.26321132940158;C=d*d/(8*sigma*sigma);n=math.ceil(math.log(300)/C);rows.append({"sigma":sigma,"Chernoff_per_sample":C,"sufficient_n_union_bound_1pct":n})
 packet={"rows":rows,"theorem":"Under the supplied equal-covariance Gaussian response model, required sample size scales quadratically with sigma and inversely with squared minimum channel separation.","physical_noise_calibration":False}
 return finalize(578,"proven_conditional_three_channel_noise_scaling", "The response likelihood remains synthetic; no photon/count detector model is supplied.",packet)


def finite_jac(fun,x,eps=1e-6):
 e=np.eye(len(x));return np.column_stack([(fun(x+eps*e[j])-fun(x-eps*e[j]))/(2*eps) for j in range(len(x))])
def st579()->dict:
 prev=json.loads((ROOT/"FIN_ST564_UltraFine_IR_Continuation_Attempt.json").read_text());seed=np.array(prev["cells"][-1]["endpoint_centers"][-1]);rows=[]
 for b in np.linspace(.2592,.3,18):
  sx=root(lambda z:regularized_ir_float(z,b),seed);sz=root(lambda z:regularized_ir_float(np.exp(z),b),np.log(np.maximum(seed,1e-12)));xlog=np.exp(sz.x);seed=xlog if sz.success else sx.x
  Jx=finite_jac(lambda x:regularized_ir_float(x,b),seed);Jz=finite_jac(lambda z:regularized_ir_float(np.exp(z),b),np.log(seed));rows.append({"b":float(b),"x_success":bool(sx.success),"log_success":bool(sz.success),"x_residual":float(np.linalg.norm(regularized_ir_float(sx.x,b),np.inf)) if sx.success else None,"log_residual":float(np.linalg.norm(regularized_ir_float(xlog,b),np.inf)) if sz.success else None,"x_condition":float(np.linalg.cond(Jx)),"log_condition":float(np.linalg.cond(Jz))})
 packet={"rows":rows,"all_log_chart_roots":all(x["log_success"] and x["log_residual"]<1e-8 for x in rows),"median_condition_improvement":float(np.median([x["x_condition"]/x["log_condition"] for x in rows])),"interval_log_chart_certificate":False}
 return finalize(579,"strong_log_chart_IR_conditioning_audit", "Numerical log-chart continuation is not an interval Krawczyk chain; it only selects the next analytic chart to implement.",packet)


def st580()->dict:
 packet={"physical_c_derived":False,"exact_causal_cone":False,"Lorentz_symmetry":False,"canonical_refinement":False,"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","new_scale_charged_source":False,"reason":"Finite-range, positive-tail and signed strict continuations have inequivalent dispersion/locality behavior."}
 return finalize(580,"blocked_no_physical_speed_cone_or_source_closure", "No bridge/role transfer, SM/GR, L_total, laboratory, or ToE closure follows.",packet)


def st581()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":2}
 return finalize(581,"blocked_no_independent_empirical_evidence", "Round two is local analytic/computational work.",packet)


def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST567"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.semilogx([x["N"] for x in rows],[x["phase_speed"] for x in rows],"o-");ax.axhline(r["ST567"]["exact_small_q_speed"],ls="--",color="#dc2626");ax.set(xlabel="N",ylabel="dimensionless phase speed",title="ST567: finite-range speed convergence");fig.tight_layout();fig.savefig(FIG_DIR/"st567_finite_speed.png",dpi=180);plt.close(fig)
 rows=r["ST568"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["q1"] for x in rows],[x["lambda1"] for x in rows],"o-");ax.set(xlabel="q1",ylabel="lambda(q1)",title="ST568: positive long-range fractional scaling");fig.tight_layout();fig.savefig(FIG_DIR/"st568_fractional.png",dpi=180);plt.close(fig)
 rows=r["ST569"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot([x["N"] for x in rows],[x["minimum_laplacian_eigenvalue"] for x in rows],"o-");ax.axhline(0,color="black");ax.set(xlabel="N",ylabel="minimum eigenvalue",title="ST569: naive signed refinement loses positivity");fig.tight_layout();fig.savefig(FIG_DIR/"st569_signed_no_go.png",dpi=180);plt.close(fig)
 rows=r["ST579"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.semilogy([x["b"] for x in rows],[x["x_condition"] for x in rows],"o-",label="raw x");ax.semilogy([x["b"] for x in rows],[x["log_condition"] for x in rows],"s-",label="log chart");ax.legend();ax.set(xlabel="b",ylabel="Jacobian condition",title="ST579: IR chart conditioning");fig.tight_layout();fig.savefig(FIG_DIR/"st579_log_chart.png",dpi=180);plt.close(fig)


def main():
 funcs={567:st567,568:st568,569:st569,570:st570,571:st571,572:st572,573:st573,574:st574,575:st575,576:st576,577:st577,578:st578,579:st579,580:st580,581:st581};r={}
 for k in range(567,582):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(567,582):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
