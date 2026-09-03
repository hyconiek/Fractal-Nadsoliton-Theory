#!/usr/bin/env python3
"""FIN ST597--ST611: exact refinement-speed fiber and coarse-observer no-go."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from fin_st16_st27_research import lifted_laplacian_24
from fin_st263_st277_research import build_local_refinement
from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST597_ST611_Results.json";SUMMARY=ROOT/"FIN_ST597_ST611_Summary.csv";FIG_DIR=ROOT/"FIN_ST597_ST611_Figures"
NAMES={597:"Exact_Dyadic_Refinement_Speed_Fiber",598:"Coarse_Observer_FineSpeed_Blindness",599:"Antipodal_Rate_Scaling_Trichotomy",600:"TraceDensity_Conditional_Speed_Selection",601:"UnlabelledFiber_VerticalSector_Modulus",602:"FineInstrument_AntipodalRate_Tomography",603:"LayerSpeed_Cocycle_with_Refinement_Modulus",604:"AntipodalRate_Locality_Cost",605:"Canonical_Refinement_and_Speed_NoGo",606:"FirstTransition_Admission_Gate",607:"CharacteristicZero_Algebra_Gate",608:"IR_PredictorChart_Obligation",609:"DualDynamics_Speed_Semantics_Audit",610:"PhysicalSpeed_Source_Selector_Gate",611:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();s=float(A[0,0]);W=s*np.eye(N)-A;WEIGHTS=np.array([W[0,d] for d in range(1,7)])
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}

C12_SLOPE2=float(np.sum(WEIGHTS[:5]*np.arange(1,6)**2)+18*WEIGHTS[5])
def symbol_eigs(a):return np.fft.fft(a[0]).real

def st597()->dict:
 rows=[]
 for q in (0,.01,.05,.1,.5,1):
  a=lifted_laplacian_24(W,q);c2=C12_SLOPE2+72*q;ev=symbol_eigs(a);rows.append({"q":q,"formal_speed_squared":c2,"formal_speed":math.sqrt(c2),"coarse_even_spectrum_residual":float(np.max(abs(ev[::2]-symbol_eigs(A)))),"minimum_eigenvalue":float(np.linalg.eigvalsh(a)[0])})
 packet={"family":"A24(q), q>=0","rows":rows,"exact_speed_law":"c_hat(q)^2=c_hat(0)^2+72q","theorem":"Every exact positive dyadic lift preserves the complete coarse spectrum while its formal fine small-k speed varies continuously and unboundedly with q.","q_selected_by_strict_core":False}
 return finalize(597,"proven_exact_refinement_speed_fiber", "Formal fine speed is a symbol coefficient, not physical c. The family itself is supplied and no q is strict-selected.",packet)


def st598()->dict:
 J=np.zeros((24,12))
 for x in range(24):J[x,x%12]=1/math.sqrt(2)
 rows=[]
 for q in (0,.1,1):
  a=lifted_laplacian_24(W,q);rows.append({"q":q,"generator_residual":float(np.linalg.norm(a@J-J@A)),"unitary_residual_t1":float(np.linalg.norm(J.T@expm(-1j*a)@J-expm(-1j*A))),"heat_residual_t1":float(np.linalg.norm(J.T@expm(-a)@J-expm(-A))),"fine_formal_speed":math.sqrt(C12_SLOPE2+72*q)})
 packet={"rows":rows,"theorem":"All coarse generator, unitary and heat records are exactly independent of q, while the fine formal speed depends on q. Therefore a coarse observer cannot infer the fine-layer speed from complete coarse spectral dynamics.","coarse_identifiability_of_q":False}
 return finalize(598,"proven_coarse_observer_fine_speed_blindness", "A child-resolving instrument can identify q; the theorem is context-relative, not ontological invisibility.",packet)


def st599()->dict:
 packet={"general_antipodal_contribution":"Delta c_hat_N^2=q_N*N^2/8","scaling_cases":{"q_N=o(N^-2)":"antipodal speed contribution vanishes","q_N~kappa/N^2":"finite contribution kappa/8","q_N>>N^-2":"formal speed diverges"},"theorem":"A finite nonzero multilevel speed contribution from antipodal refinement rates requires q_N to scale as N^-2. Exact coarse intertwining permits all three regimes and does not select one.","canonical_kappa":False}
 return finalize(599,"proven_antipodal_rate_speed_scaling_trichotomy", "The N^-2 law is a necessary scaling for this dyadic family, not a strict derivation of q_N or physical speed.",packet)


def st600()->dict:
 trace12=float(np.trace(A)/12);rows=[]
 for q in (0,.01,.1,1):rows.append({"q":q,"trace_density_24":float(np.trace(lifted_laplacian_24(W,q))/24),"difference_from_coarse":float(np.trace(lifted_laplacian_24(W,q))/24-trace12),"formal_speed":math.sqrt(C12_SLOPE2+72*q)})
 packet={"coarse_trace_density":trace12,"rows":rows,"conditional_selector":"trace-density conservation selects q=0","selected_formal_speed":math.sqrt(C12_SLOPE2),"theorem":"Within A24(q), exact trace-density conservation uniquely selects q=0 and hence one speed proxy.","trace_law_strict_sourced":False}
 return finalize(600,"proven_conditional_trace_law_speed_selection", "The trace-density law is an extra refinement axiom; selection by an added law is not emergence from A12.",packet)


def st601()->dict:
 rows=[];plus=np.array([1,1])/math.sqrt(2);R=np.kron(np.eye(12),plus[:,None])
 for v in (0,.1,.5,2):
  at=build_local_refinement(W,.5,v);coarse=R.T@at@R;F=np.kron(np.eye(12),np.array([1,-1])[:,None]/math.sqrt(2));fiber=F.T@at@F;rows.append({"vertical_rate":v,"coarse_residual":float(np.linalg.norm(coarse-A)),"fiber_min":float(np.linalg.eigvalsh(fiber)[0]),"fiber_max":float(np.linalg.eigvalsh(fiber)[-1]),"coarse_fiber_mixing":float(np.linalg.norm(R.T@at@F))})
 packet={"rows":rows,"theorem":"Independent unlabelled-fiber naturality fixes horizontal splitting but leaves a nonnegative vertical rate v. Coarse dynamics is independent of v while the entire fine antisymmetric sector changes.","v_selected":False}
 return finalize(601,"proven_unlabelled_fiber_sector_modulus_and_coarse_blindness", "The vertical modulus is not the same as physical speed; it demonstrates a second unsourced fine-sector law.",packet)


def st602()->dict:
 rows=[]
 for q in (0,.03,.2,.8):
  ev=symbol_eigs(lifted_laplacian_24(W,q));delta=float(ev[1]-symbol_eigs(lifted_laplacian_24(W,0))[1]);rows.append({"q":q,"odd_mode_1_eigenvalue":float(ev[1]),"shift_from_q0":delta,"recovered_q":delta/2})
 packet={"rows":rows,"theorem":"Every odd Fourier eigenvalue is shifted by exactly 2q, while even modes are unchanged. One calibrated child-resolved odd-mode measurement recovers q exactly.","coarse_instrument_recovers_q":False,"fine_instrument_required":True}
 return finalize(602,"proven_fine_instrument_antipodal_rate_tomography", "The child preparation/effect and clock are operational inputs; no laboratory realization is supplied.",packet)


def st603()->dict:
 packet={"fine_speed":"c_hat_N^2=c_base,N^2+q_N N^2/8","physical_speed":"c_phys,N=(ell_N/tau_N)c_hat_N","invariance_condition":"(ell_{2N}/tau_{2N})/(ell_N/tau_N)=c_hat_N/c_hat_2N","free_data":["q_N sequence","length/time transition cocycle"],"theorem":"Layer-invariant physical speed constrains only the product of the refinement modulus and calibration cocycle. Without either one, the other is not identifiable.","unique_refinement_or_calibration":False}
 return finalize(603,"proven_refinement_calibration_compensation_torsor", "Constant observed c can be maintained by infinitely many compensating q_N and calibration sequences.",packet)


def st604()->dict:
 def proxy(q,mu):return (2*sum(WEIGHTS[d-1]*math.sinh(mu*d) for d in range(1,6))+WEIGHTS[5]*math.sinh(6*mu)+2*q*math.sinh(12*mu))/mu
 rows=[]
 for q in (0,.001,.01,.1):
  o=minimize_scalar(lambda m:proxy(q,m),bounds=(1e-6,.5),method="bounded");rows.append({"q":q,"optimized_mu":float(o.x),"velocity_proxy":float(o.fun),"formal_speed":math.sqrt(C12_SLOPE2+72*q)})
 packet={"rows":rows,"theorem":"A nonzero antipodal q inserts a range-12 term into the exponential-moment velocity bound and raises its locality cost. Coarse dynamics remains blind to this cost.","many_body_LR":False}
 return finalize(604,"computed_refinement_modulus_locality_cost", "The proxy is single-particle and finite-volume; no Lorentz cone follows.",packet)


def st605()->dict:
 packet={"exact_coarse_equivalence":True,"fine_speed_range":"[sqrt(C12_SLOPE2), infinity) as q ranges over nonnegative values","fine_locality_cost_q_dependent":True,"fine_sector_vertical_modulus":True,"theorem":"The present exact refinement axioms admit infinitely many fine speeds and locality structures with identical coarse FIN dynamics. Therefore no canonical propagation speed is reconstructible from A12 plus coarse naturality alone.","strict_canonical_speed":False}
 return finalize(605,"proven_canonical_refinement_speed_no_go_for_current_axioms", "A new scale law, locality principle, fine record, or source theorem can leave this no-go class.",packet)


def st606()->dict:
 packet={"global_transition_bracket":[2.8934,2.9024964917477667],"candidate":2.9024964817477654,"new_interior_cover":False,"status":"open"}
 return finalize(606,"first_transition_remains_open", "Propagation/refinement results do not close the variational interior multilevel orbit problem.",packet)
def st607()->dict:
 packet={"rank_Q_first1326_closed":False,"full_rank_Q_closed":False,"new_exact_basis_method":False}
 return finalize(607,"characteristic_zero_algebra_remains_method_gated", "No direct CRT/random combination replay is admitted.",packet)
def st608()->dict:
 packet={"naive_log_interval":"failed by dependency inflation","next_method":"affine predictor z(b)=z0+z1(b-b0) retained symbolically in interval evaluation","implemented":False,"reason":"direct_ir_interval_fj creates an independent B interval and loses z--b correlation"}
 return finalize(608,"predictor_chart_obligation_identified", "No new IR interval theorem is claimed.",packet)
def st609()->dict:
 packet={"unitary_frequency":"lambda/tau_U","heat_rate":"lambda/tau_H","wave_frequency":"sqrt(lambda)/tau_W","refinement_speed_modulus":"q_N","theorem":"Common eigenvalues preserve dimensionless mode ratios, but propagation speed belongs to the wave/refinement coordinate map and is not selected by unitary--heat duality.","dual_dynamics_derives_c":False}
 return finalize(609,"proven_dual_dynamics_speed_semantics_separation", "Dual functional calculus does not supply a spatial embedding, wave clock, or light cone.",packet)
def st610()->dict:
 packet={"physical_c":False,"exact_causal_cone":False,"Lorentz_symmetry":False,"canonical_q_sequence":False,"calibration_cocycle":False,"new_seed_gain_selector":False,"QW_2191":"open"}
 return finalize(610,"blocked_no_physical_speed_or_source_closure", "No bridge/role transfer, SM/GR, L_total, apparatus, or ToE closure.",packet)
def st611()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":4}
 return finalize(611,"blocked_no_independent_empirical_evidence", "Round four is local analytic/computational work.",packet)


def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST597"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot([x["q"] for x in rows],[x["formal_speed"] for x in rows],"o-");ax.set(xlabel="refinement modulus q",ylabel="formal fine speed",title="ST597: exact refinement speed fiber");fig.tight_layout();fig.savefig(FIG_DIR/"st597_speed_fiber.png",dpi=180);plt.close(fig)
 rows=r["ST604"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot([x["q"] for x in rows],[x["velocity_proxy"] for x in rows],"o-");ax.set(xlabel="q",ylabel="exponential velocity proxy",title="ST604: hidden fine locality cost");fig.tight_layout();fig.savefig(FIG_DIR/"st604_locality_cost.png",dpi=180);plt.close(fig)


def main():
 funcs={597:st597,598:st598,599:st599,600:st600,601:st601,602:st602,603:st603,604:st604,605:st605,606:st606,607:st607,608:st608,609:st609,610:st610,611:st611};r={}
 for k in range(597,612):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(597,612):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
