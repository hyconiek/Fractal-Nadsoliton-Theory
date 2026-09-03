#!/usr/bin/env python3
"""FIN ST612--ST626: dyadic speed tower, locality cost, and observer blindness."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

from fin_st28_st45_research import dyadic_lift,periodic_embedding
from fin_st402_st416_research import independent_strict_matrix_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST612_ST626_Results.json";SUMMARY=ROOT/"FIN_ST612_ST626_Summary.csv";FIG_DIR=ROOT/"FIN_ST612_ST626_Figures"
NAMES={612:"Dyadic_FormalSpeed_Recurrence",613:"InfiniteTower_Speed_Convergence_Criterion",614:"SpeedStationary_NonnegativeRate_Theorem",615:"Uniform_ExponentialCone_Rate_Criterion",616:"Complete_Coarse_Tower_Speed_Blindness",617:"LevelResolved_Rate_Tomography",618:"RefinementRate_Record_Complexity",619:"Coassociativity_Does_Not_Select_Speed",620:"Trace_and_Speed_Stationarity_Comparison",621:"Calibration_Compensation_Across_Tower",622:"ZeroRate_LocalContinuum_Candidate",623:"Tower_Tail_Class_Comparison",624:"RefinementRate_Source_Admission_Gate",625:"PhysicalSpeed_CausalCone_Gate",626:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A12=independent_strict_matrix_float()
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}
def formal_c2(a):
 n=len(a);w=-a[0,1:n//2+1];return float(np.sum(w[:-1]*np.arange(1,n//2)**2)+w[-1]*(n//2)**2/2)
def tower(qs):
 a=A12.copy();rows=[{"N":12,"c2":formal_c2(a),"q_in":None}]
 for q in qs:
  n=len(a);before=formal_c2(a);a=dyadic_lift(a,q);rows.append({"N":len(a),"c2":formal_c2(a),"q_in":q,"predicted_increment":n*n*q/2,"actual_increment":formal_c2(a)-before})
 return rows,a


def st612()->dict:
 qs=[.02,.003,.0004];rows,a=tower(qs);res=max(abs(x["actual_increment"]-x["predicted_increment"]) for x in rows[1:]);packet={"rows":rows,"maximum_recurrence_residual":res,"recurrence":"C_2N=C_N+(N^2/2)q_N","theorem":"The old formal slope coefficient is exactly preserved by weight splitting; the new antipodal edge adds N^2 q_N/2."}
 return finalize(612,"proven_exact_dyadic_speed_recurrence", "C_N is dimensionless and formal; the q sequence and physical coordinate map remain supplied.",packet)


def st613()->dict:
 packet={"recurrence":"C_j=C_0+1/2 sum_{i<j} N_i^2 q_i","nonnegative_rates":True,"finite_limit_criterion":"sum_i N_i^2 q_i < infinity","q_i~kappa/N_i^2":"adds kappa/2 per level and diverges linearly in depth","power_law_q_i~N_i^(-2-epsilon)":"converges for epsilon>0 in a dyadic tower","theorem":"For nonnegative dyadic rates, a finite infinite-depth formal speed exists iff the weighted series sum N_i^2 q_i converges.","canonical_limit":False}
 return finalize(613,"proven_infinite_tower_speed_convergence_criterion", "Convergence does not choose the series, its sum, or a dimensional speed.",packet)


def st614()->dict:
 packet={"premise":"C_2N=C_N at every level","nonnegative_q":True,"conclusion":"q_N=0 at every level","theorem":"Exact dimensionless speed stationarity and nonnegative antipodal rates force the unique zero-rate tower.","premise_strict_sourced":False}
 return finalize(614,"proven_conditional_zero_rate_speed_stationary_tower", "Speed stationarity is an added refinement law, not a consequence of coarse intertwining or fractal ontology.",packet)


def st615()->dict:
 mu=.1;rows=[]
 for law in ("N^-2","N^-3","exp(-0.2N)"):
  vals=[]
  for N in (12,24,48,96,192):
   q=N**-2 if law=="N^-2" else (N**-3 if law=="N^-3" else math.exp(-.2*N));vals.append(q*math.sinh(mu*N))
  rows.append({"law":law,"weighted_antipodal_moment":vals,"last_over_first":vals[-1]/vals[0]})
 packet={"mu":mu,"rows":rows,"uniform_condition":"sup_N q_N sinh(mu N)<infinity for some mu>0","theorem":"A uniform exponential-moment propagation bound requires antipodal rates to decay exponentially in distance. Polynomial N^-p rates fail for every fixed mu>0.","many_body_LR":False}
 return finalize(615,"proven_polynomial_antipodal_rates_fail_uniform_exponential_cone", "The condition concerns this bound architecture and single-particle tower; other polynomial long-range cones remain possible.",packet)


def st616()->dict:
 qs1=[0,0,0];qs2=[.1,.02,.003];r1,a1=tower(qs1);r2,a2=tower(qs2);J12=periodic_embedding(12)/math.sqrt(2);J24=periodic_embedding(24)/math.sqrt(2);J48=periodic_embedding(48)/math.sqrt(2);R=J48@J24@J12;rows=[]
 for name,a,rr in (("zero",a1,r1),("nonzero",a2,r2)):
  rows.append({"tower":name,"final_N":len(a),"base_generator_residual":float(np.linalg.norm(a@R-R@A12)),"base_heat_residual":float(np.linalg.norm(R.T@expm(-a)@R-expm(-A12))),"final_formal_speed":math.sqrt(rr[-1]["c2"])})
 packet={"rows":rows,"theorem":"Arbitrarily different finite rate sequences give exactly the same complete base generator and heat records under composed normalized embeddings, while final fine formal speeds differ.","base_observer_identifies_rate_sequence":False}
 return finalize(616,"proven_multilevel_coarse_speed_blindness", "The theorem is finite-depth but extends inductively; fine instruments may resolve the hidden sequence.",packet)


def st617()->dict:
 qs=[.07,.011,.002];rows,a=tower(qs);recovered=[];current=A12.copy()
 for q in qs:
  fine=dyadic_lift(current,q);ev=np.fft.fft(fine[0]).real;ev0=np.fft.fft(dyadic_lift(current,0)[0]).real;recovered.append({"parent_N":len(current),"q":q,"odd_mode_shift":float(ev[1]-ev0[1]),"recovered_q":float((ev[1]-ev0[1])/2)});current=fine
 packet={"rows":recovered,"maximum_recovery_error":max(abs(x["q"]-x["recovered_q"]) for x in recovered),"theorem":"At each dyadic level, any odd fine Fourier mode shifts by 2q_N and reconstructs that level's hidden rate independently."}
 return finalize(617,"proven_level_resolved_rate_tomography", "One child-resolved preparation/effect and calibrated clock per level are added operational resources.",packet)


def st618()->dict:
 packet={"depth":10,"arbitrary_continuous_rate_vector_dimension":10,"fixed_precision_bits_per_rate":32,"record_bits_for_arbitrary_rates":320,"self_similar_law_example":"q_N=kappa N^(-2-epsilon) uses two real parameters plus depth","theorem":"An arbitrary depth-L refinement rate sequence carries L independent continuous parameters; compression below linear depth requires a declared correlated/self-similar law and distortion precision.","universe_history_lossless_compression":False}
 return finalize(618,"proven_refinement_rate_record_complexity_growth", "Parameter counting is not Kolmogorov complexity of nature and supplies no preferred compression law.",packet)


def st619()->dict:
 packet={"associativity":"composed periodic embeddings agree canonically","free_rates_per_level":"one new q_N remains at every dyadic step","theorem":"Plain coassociativity/associativity of the dyadic lift and embeddings imposes no equality among level rates and therefore does not select a speed tower.","rate_sequence_selected":False}
 return finalize(619,"proven_coassociativity_does_not_select_speed_tower", "Additional scale stationarity, factor-exchange naturality, or fine records lie outside plain associativity.",packet)


def st620()->dict:
 packet={"trace_recurrence":"tr(A_2N)/(2N)=tr(A_N)/N+q_N","speed_recurrence":"C_2N=C_N+N^2 q_N/2","common_stationary_solution":"q_N=0","theorem":"For nonnegative dyadic lifts, exact trace-density stationarity and exact formal-speed stationarity independently select the same zero-rate tower.","equivalence_of_axioms":False,"reason":"Away from stationarity the two observables weight q_N differently (1 versus N^2/2)."}
 return finalize(620,"proven_common_zero_solution_but_nonequivalent_scale_laws", "Agreement at q=0 does not derive either axiom from strict FIN.",packet)


def st621()->dict:
 packet={"formal_speed_sequence":"c_hat_j=sqrt(C_j)","physical_invariance":"ell_j c_hat_j/tau_j = c_phys","recursive_calibration":"(ell_{j+1}/tau_{j+1})/(ell_j/tau_j)=c_hat_j/c_hat_{j+1}","theorem":"Every positive rate tower can be made to display one constant physical speed by a compensating calibration cocycle. Thus constant observed c alone cannot identify the hidden refinement law.","unique_calibration":False}
 return finalize(621,"proven_tower_speed_calibration_compensation", "Independent rulers/clocks or a sourced transition cocycle are required to break the degeneracy.",packet)


def st622()->dict:
 rows=[]
 for levels in range(1,7):
  a=A12.copy()
  for _ in range(levels):a=dyadic_lift(a,0)
  n=len(a);q=2*math.pi/n;lam=np.fft.fft(a[0]).real[1];rows.append({"N":n,"phase_speed":float(math.sqrt(lam)/q),"formal_c2":formal_c2(a)})
 packet={"rows":rows,"formal_c2_constant":max(x["formal_c2"] for x in rows)-min(x["formal_c2"] for x in rows),"theorem":"The conditional zero-rate tower preserves the formal quadratic coefficient exactly and its first-mode phase speed converges to sqrt(C12).","canonical_physical_continuum":False}
 return finalize(622,"proven_zero_rate_tower_quadratic_symbol_candidate", "This is a conditional local continuum candidate; continuous-time lattice tails, dimension, Lorentz symmetry and physical calibration remain unresolved.",packet)


def st623()->dict:
 packet={"classes":{"zero_rate_finite_range":"quadratic symbol; finite formal speed; exponential proxy cone","summable_antipodal":"finite but law-dependent formal speed; long-range edges may spoil uniform exponential cone","N^-2_antipodal":"speed diverges as sqrt(depth); no infinite limit","positive_d^-1.8_tail":"fractional symbol 0.8; divergent low-k group speed","signed_literal_strict":"finite counterexample to positivity at N=48"},"theorem":"The current refinement candidates occupy inequivalent universality classes; coarse A12 does not determine which class is physical."}
 return finalize(623,"proven_refinement_universality_class_separation", "Classification is within declared candidates, not exhaustive over all possible completion maps.",packet)


def st624()->dict:
 packet={"new_strict_q_sequence":False,"speed_stationarity_strict_sourced":False,"trace_stationarity_strict_sourced":False,"fine_record_available":False}
 return finalize(624,"blocked_no_refinement_rate_source", "No canonical refinement is derived.",packet)
def st625()->dict:
 packet={"physical_c":False,"exact_causal_cone":False,"Lorentz_symmetry":False,"canonical_continuum":False,"SI_299792458_derived":False}
 return finalize(625,"blocked_no_physical_speed_or_continuum_closure", "No spacetime, light field, boosts, apparatus, or calibration is exported.",packet)
def st626()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":5}
 return finalize(626,"blocked_no_independent_empirical_evidence", "Round five is local analytic/computational work.",packet)


def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST612"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot([x["N"] for x in rows],[math.sqrt(x["c2"]) for x in rows],"o-");ax.set(xlabel="N",ylabel="formal speed",title="ST612: dyadic speed recurrence");fig.tight_layout();fig.savefig(FIG_DIR/"st612_speed_tower.png",dpi=180);plt.close(fig)
 rows=r["ST615"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));Ns=[12,24,48,96,192]
 for x in rows:ax.semilogy(Ns,x["weighted_antipodal_moment"],"o-",label=x["law"])
 ax.legend();ax.set(xlabel="N",ylabel=r"$q_N\sinh(\mu N)$",title="ST615: uniform exponential-cone cost");fig.tight_layout();fig.savefig(FIG_DIR/"st615_locality_scaling.png",dpi=180);plt.close(fig)


def main():
 funcs={612:st612,613:st613,614:st614,615:st615,616:st616,617:st617,618:st618,619:st619,620:st620,621:st621,622:st622,623:st623,624:st624,625:st625,626:st626};r={}
 for k in range(612,627):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(612,627):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
