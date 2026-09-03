#!/usr/bin/env python3
"""FIN ST627--ST641: continuum scaling, one-clock no-go, and conditional wave cone."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST627_ST641_Results.json";SUMMARY=ROOT/"FIN_ST627_ST641_Summary.csv";FIG_DIR=ROOT/"FIN_ST627_ST641_Figures"
NAMES={627:"Wave_Heat_Unitary_Refinement_TimeScaling",628:"SingleClock_Continuum_NoGo",629:"Minimal_ChannelSpecific_Clock_Maps",630:"ZeroRate_Rescaled_Dispersion_Convergence",631:"ZeroRate_WaveCone_Tail_Convergence",632:"Conditional_Continuum_Wave_Operator_Theorem",633:"Conditional_Heat_Unitary_Continuum_Limits",634:"DualDynamics_Continuum_Scaling_Incompatibility",635:"PhysicalSpeed_Anchor_and_SI_Boundary",636:"Conditional_LightCone_Emergence_Theorem",637:"Lorentz_and_LightField_Obligation_Audit",638:"RefinementCalibration_NoGo_Consolidation",639:"Transition_and_Algebra_Frontier_Gate",640:"Strict_Physics_Source_Admission_Gate",641:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();s=float(A[0,0]);W=s*np.eye(12)-A;WEIGHTS=np.array([W[0,d] for d in range(1,7)]);C=float(np.sum(WEIGHTS*np.arange(1,7)**2));CHAT=math.sqrt(C)
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}
def lam(q):return float(2*np.sum(WEIGHTS*(1-np.cos(np.arange(1,7)*q))))


def st627()->dict:
 packet={"lattice_spacing":"a_N=L/N","low_mode_eigenvalue":"lambda_N(a_N k)~C a_N^2 k^2","unitary_nontrivial_limit":"dimensionless time t_U~T/a_N^2","heat_nontrivial_limit":"dimensionless time t_H~T/a_N^2","wave_nontrivial_limit":"dimensionless time t_W~T/a_N","theorem":"Under a fixed-physical-length finite-range refinement, unitary and heat channels require diffusive a^-2 time acceleration, while the second-order wave channel requires ballistic a^-1 acceleration.","common_raw_time_scaling":False}
 return finalize(627,"proven_channel_refinement_time_scaling_exponents", "The result is asymptotic and conditional on the zero-rate finite-range tower and physical coordinate embedding.",packet)


def st628()->dict:
 packet={"required_exponents":{"unitary":2,"heat":2,"wave":1},"theorem":"No single power-law clock t_dimless~a^-z yields simultaneous nontrivial unitary/heat and wave continuum dynamics: z must equal 2 and 1 at once.","one_universal_clock":False,"escape_options":["channel-specific conversion maps","rescale the generator differently by sector","drop one continuum channel"]}
 return finalize(628,"proven_single_clock_three_channel_continuum_no_go", "The no-go concerns one raw refinement clock; it does not contradict common finite spectral projectors or dimensionless duality.",packet)


def st629()->dict:
 packet={"minimal_maps":{"unitary_heat":"tau_D,N proportional to a_N^2","wave":"tau_W,N proportional to a_N"},"independent_exponent_count":2,"theorem":"Two channel classes of clock conversion are sufficient and necessary for nontrivial limits: one diffusive/Schrodinger map shared by U and P, and one ballistic map for the second-order wave equation.","absolute_normalization_count_after_exponents":2,"strict_clock_source":False}
 return finalize(629,"proven_minimal_two_clock_class_continuum_architecture", "Proportionality constants and physical seconds remain calibration inputs.",packet)


def st630()->dict:
 L=12.;ks=[2*math.pi*m/L for m in (1,2,3)];rows=[]
 for N in (48,96,192,384,768,1536,3072):
  a=L/N;rr=[]
  for k in ks:
   q=a*k;scaled=lam(q)/a**2;rr.append({"k":k,"scaled_lambda":scaled,"continuum_Ck2":C*k*k,"relative_error":abs(scaled-C*k*k)/(C*k*k)})
  rows.append({"N":N,"a":a,"modes":rr})
 packet={"C":C,"c_hat":CHAT,"rows":rows,"maximum_last_level_relative_error":max(x["relative_error"] for x in rows[-1]["modes"]),"theorem":"For every fixed physical k, a^-2 lambda(a k)->C k^2 with O(a^2 k^4) error because the finite-range fourth moment is finite."}
 return finalize(630,"proven_rescaled_zero_rate_dispersion_convergence", "The zero-rate tower and fixed circumference are added premises; convergence of symbols is not laboratory evidence.",packet)


def wave_tail(N,L=12.,t=.5,margin=1.05):
 a=L/N;q=2*np.pi*np.arange(N)/N;lv=np.zeros(N)
 for d in range(1,7):lv+=2*WEIGHTS[d-1]*(1-np.cos(d*q))
 amp=np.fft.ifft(np.cos(t*np.sqrt(np.maximum(lv,0))/a)).real;dist=np.minimum(np.arange(N),N-np.arange(N))*a;mask=dist>CHAT*t*margin;return float(np.sum(amp[mask]**2)/np.sum(amp**2)),float(np.max(np.abs(amp[mask])))
def st631()->dict:
 rows=[]
 for N in (48,96,192,384,768,1536,3072,6144,12288):
  tail,mx=wave_tail(N);rows.append({"N":N,"outside_1_05ct_relative_L2":tail,"maximum_outside_amplitude":mx})
 packet={"physical_time":.5,"cone_margin":1.05,"rows":rows,"last_tail":rows[-1]["outside_1_05ct_relative_L2"],"monotone_last5":all(rows[i+1]["outside_1_05ct_relative_L2"]<rows[i]["outside_1_05ct_relative_L2"] for i in range(len(rows)-5,len(rows)-1)),"theorem_status":"strong numerical continuum-cone diagnostic"}
 return finalize(631,"strong_zero_rate_wave_tail_convergence", "Finite grids retain analytic tails; the numerical sequence is not an outward uniform convergence proof.",packet)


def st632()->dict:
 fourth=float(np.sum(WEIGHTS*np.arange(1,7)**4));packet={"discrete_operator":"a^-2 A_N on the paired finite-range zero-rate tower","second_moment":C,"fourth_moment":fourth,"continuum_operator":"-C partial_x^2","Taylor_remainder_coefficient_bound":"(a^2/12) sum_d w_d d^4 times ||f''''||","theorem":"On periodic C^4 functions sampled on a fixed circle, the rescaled finite-range generator converges uniformly to -C d^2/dx^2 with an explicit O(a^2) consistency bound.","physical_field_theory":False}
 return finalize(632,"proven_conditional_continuum_wave_operator_consistency", "Consistency plus stability is a route to convergence; no strict selection of the tower, field ontology, or apparatus follows.",packet)


def st633()->dict:
 packet={"rescaled_generator":"A_N/a_N^2","heat_limit":"exp(T C partial_x^2)","unitary_limit":"exp(i T C partial_x^2) up to sign convention","wave_limit":"cos(T sqrt(-C partial_x^2))","theorem":"The same rescaled spatial stencil supports heat and Schrodinger-type first-order limits when their physical generator is A_N/a^2, and a wave limit when A_N/a^2 is interpreted as stiffness in a second-order equation.","shared_physical_time_dimension":False}
 return finalize(633,"proven_conditional_three_functional_continuum_calculi", "The channel equations use different dimensional roles and operational state spaces despite the shared spatial symbol.",packet)


def st634()->dict:
 packet={"finite_level_identity":"U_t=exp(-itA), P_t=exp(-tA)","continuum_scaling":"U/P use t~a^-2; wave cos(t sqrt(A)) uses t~a^-1","theorem":"Finite spectral duality commutes with diffusive scaling for U/P, but the wave square-root calculus belongs to a different dynamic scaling exponent. No single refinement-time functor preserves all three raw parameterizations.","dual_dynamics_refuted":False}
 return finalize(634,"proven_duality_survives_but_wave_scaling_separates", "This sharpens rather than weakens the common spectral theorem; channel-specific time maps remain mandatory.",packet)


def st635()->dict:
 packet={"dimensionless_continuum_speed":CHAT,"physical_speed":"c_phys=(ell_*/tau_W,*) c_hat","one_anchor_after_ratio":"one transverse length/time calibration fixes the numeric speed","SI_value":"299792458 m/s is exact by SI definition","theorem":"The continuum slope determines only a dimensionless factor. A physical number requires one independent length/time ratio after the wave clock class is chosen.","FIN_prediction_of_SI_c":False}
 return finalize(635,"proven_dimensionless_speed_plus_anchor_boundary", "Defining the meter by light propagation would make agreement circular; independent clock/ruler records are required for a predictive test.",packet)


def st636()->dict:
 packet={"premises":["zero-rate finite-range tower","fixed physical circle embedding","ballistic wave time scaling","continuum convergence/stability"],"limit_equation":"partial_T^2 phi-C partial_x^2 phi=0","causal_speed":"sqrt(C) in internal units","theorem":"Under the declared premises the limiting 1+1 wave equation has an exact domain of dependence with speed sqrt(C), even though every finite lattice has exponentially small analytic tails.","strict_FIN_consequence":False}
 return finalize(636,"proven_conditional_continuum_lightcone_mechanism", "The premises are not strict-sourced, and a scalar 1+1 wave cone is not observed 3+1 Lorentzian spacetime or electromagnetism.",packet)


def st637()->dict:
 packet={"available":"conditional scalar 1+1 wave equation and speed","missing":["3+1 dimensional refinement","isotropy under spatial rotations","Lorentz boost representation on states/observables","Maxwell or gauge connection","massless light sector source","matter coupling","operational rods/clocks"],"theorem":"Linear low-k dispersion alone does not establish Lorentz symmetry or identify the mode as physical light."}
 return finalize(637,"proven_Lorentz_lightfield_obligation_separation", "No physical photon/light field is derived.",packet)


def st638()->dict:
 packet={"coarse_A12_allows_speed_fiber":True,"zero_rate_selected_only_by_added_stationarity":True,"one_clock_three_channel_no_go":True,"conditional_wave_cone_available":True,"physical_c_derived":False,"conclusion":"Current FIN permits a mathematically coherent emergent-cone completion but does not make it inevitable. Refinement, channel clocks and calibration remain independent bridge data."}
 return finalize(638,"consolidated_refinement_calibration_no_go_with_conditional_bridge", "This is a no-go for inevitability under current axioms, not for all future FIN completions.",packet)


def st639()->dict:
 packet={"first_transition_global_identity":False,"rank_Q_first1326_closed":False,"full_degree8_rank_Q_closed":False,"lanes":"method gated"}
 return finalize(639,"transition_and_algebra_frontiers_remain_open", "Propagation progress does not close independent mathematical frontiers.",packet)
def st640()->dict:
 packet={"canonical_zero_rate_tower":False,"wave_clock_source":False,"length_time_anchor":False,"light_sector_source":False,"seed_gain_selector":False,"QW_2191":"open"}
 return finalize(640,"blocked_no_strict_physics_source_package", "No bridge/role transfer, SM/GR, L_total, apparatus, or ToE closure.",packet)
def st641()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":6,"six_round_cycle":"complete"}
 return finalize(641,"blocked_no_independent_empirical_evidence", "Round six and the complete cycle are local analytic/computational work.",packet)


def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST630"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["N"] for x in rows],[max(m["relative_error"] for m in x["modes"]) for x in rows],"o-");ax.set(xlabel="N",ylabel="max relative dispersion error",title="ST630: rescaled symbol convergence");fig.tight_layout();fig.savefig(FIG_DIR/"st630_dispersion_convergence.png",dpi=180);plt.close(fig)
 rows=r["ST631"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["N"] for x in rows],[x["outside_1_05ct_relative_L2"] for x in rows],"o-");ax.set(xlabel="N",ylabel="relative tail outside 1.05 c t",title="ST631: conditional wave-cone tail");fig.tight_layout();fig.savefig(FIG_DIR/"st631_wave_tail.png",dpi=180);plt.close(fig)


def main():
 funcs={627:st627,628:st628,629:st629,630:st630,631:st631,632:st632,633:st633,634:st634,635:st635,636:st636,637:st637,638:st638,639:st639,640:st640,641:st641};r={}
 for k in range(627,642):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(627,642):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
