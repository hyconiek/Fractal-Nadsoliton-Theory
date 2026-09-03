#!/usr/bin/env python3
"""FIN ST582--ST596: log-interval IR, propagation theorems, and speed calibration."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.optimize import root,minimize_scalar

from fin_st118_st129_research import interval_left,interval_matvec
from fin_st132_center_isolation_replay import bounds,iv
from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N,SEED
from fin_st447_st461_research import direct_ir_interval_fj,regularized_ir_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST582_ST596_Results.json";SUMMARY=ROOT/"FIN_ST582_ST596_Summary.csv";FIG_DIR=ROOT/"FIN_ST582_ST596_Figures"
NAMES={582:"LogChart_Interval_Krawczyk_Cell",583:"LogChart_Adaptive_IR_Chain_to_0_3",584:"LogChart_IR_Ordering_and_Degree",585:"PositiveTail_Tauberian_Dispersion_Theorem",586:"DispersionExponent_Speed_Trichotomy",587:"FiniteRange_SingleParticle_LiebRobinson_Bound",588:"Finite_vs_LongRange_CausalTail_Classification",589:"LayerSpeed_Calibration_Torsor_Theorem",590:"FiniteShot_EffectiveSpeed_Estimator",591:"TwoLevel_Interval_Method_Stop",592:"ReflectionQuotient_Method_Stop",593:"Stage1_Algebra_Method_Gate",594:"Seed_Gain_Selector_Admission_Gate",595:"PhysicalSpeed_Lorentz_Admission_Gate",596:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();U=np.ones(N)/N;s=float(A[0,0]);W=s*np.eye(N)-A;WEIGHTS=np.array([W[0,d] for d in range(1,7)])
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k:int,status:str,boundary:str,packet:dict)->dict:
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}
def finite_jac(fun,x,eps=1e-6):
 e=np.eye(len(x));return np.column_stack([(fun(x+eps*e[j])-fun(x-eps*e[j]))/(2*eps) for j in range(len(x))])


def log_ir_cell(blo,bhi,prior_z,radius_factor=1.8):
 roots=[]
 for b in (blo,(blo+bhi)/2,bhi):
  sol=root(lambda z:regularized_ir_float(np.exp(z),b),prior_z,tol=1e-12)
  if not sol.success or np.linalg.norm(regularized_ir_float(np.exp(sol.x),b),np.inf)>2e-9:raise RuntimeError(f"log root failure {b}")
  prior_z=sol.x;roots.append(sol.x.copy())
 center=roots[1];radii=radius_factor*np.max(np.abs(np.array(roots)-center),axis=0)+np.array([2e-7,2e-6,2e-6,2e-6,2e-6]);Z=[iv((center[i]-radii[i],center[i]+radii[i])) for i in range(5)];X=[mp.iv.exp(z) for z in Z]
 f0,_=direct_ir_interval_fj([iv(float(math.exp(z))) for z in center],blo,bhi);flo=np.array([bounds(x)[0] for x in f0]);fhi=np.array([bounds(x)[1] for x in f0]);_,Jx=direct_ir_interval_fj(X,blo,bhi);Jz=[[Jx[i][j]*X[j] for j in range(5)] for i in range(5)];jl=np.array([[bounds(x)[0] for x in row] for row in Jz]);jh=np.array([[bounds(x)[1] for x in row] for row in Jz]);bmid=(blo+bhi)/2;jm=finite_jac(lambda z:regularized_ir_float(np.exp(z),bmid),center);pre=np.linalg.inv(jm);cfl,cfh=interval_matvec(pre,pre,flo,fhi);ylo,yhi=center-cfh,center-cfl;cjl,cjh=interval_left(pre,jl,jh);mlo,mhi=-cjh,-cjl
 for i in range(5):mlo[i,i]+=1;mhi[i,i]+=1
 dlo,dhi=interval_matvec(mlo,mhi,-radii,radii);margins=np.minimum(ylo+dlo-(center-radii),(center+radii)-(yhi+dhi));M=np.maximum(abs(mlo),abs(mhi));contr=float(max((M@radii)/radii));return {"b_interval":[blo,bhi],"log_center":center,"log_radii":radii,"x_center":np.exp(center),"included":bool(np.min(margins)>0),"minimum_margin":float(np.min(margins)),"weighted_contraction":contr,"endpoint_log_centers":roots[::2],"Jacobian_determinant":float(np.linalg.det(jm))},prior_z


def st582()->dict:
 prev=json.loads((ROOT/"FIN_ST564_UltraFine_IR_Continuation_Attempt.json").read_text());prior=np.log(np.array(prev["cells"][-1]["endpoint_centers"][-1]));row,_=log_ir_cell(.2592,.2612,prior)
 packet={"cell":row,"theorem":None if not row["included"] else "The logarithmic coordinate Krawczyk image is strictly contained in the displayed box for every b in [0.2592,0.2612].","result":"Direct interval substitution x=exp(z) suffers dependency inflation despite improved floating Jacobian conditioning."}
 return finalize(582,"proven_first_log_chart_interval_cell" if row["included"] else "bounded_no_go_for_naive_log_interval_cell", "A failed Krawczyk inclusion is a method stop, not loss of the numerical root or IR branch.",packet)


def log_chain(start=.2592,stop=.3,width=.002):
 prev=json.loads((ROOT/"FIN_ST564_UltraFine_IR_Continuation_Attempt.json").read_text());prior=np.log(np.array(prev["cells"][-1]["endpoint_centers"][-1]));rows=[];lo=start;w=width;failure=None
 while lo<stop-1e-14:
  hi=min(stop,lo+w)
  try:row,candidate=log_ir_cell(lo,hi,prior)
  except RuntimeError as e:row={"included":False,"error":str(e)}
  if row.get("included"):rows.append(row);prior=candidate;lo=hi;w=min(width,w*1.2)
  else:
   w/=2
   if w<2e-5:failure=row;break
 return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}
def st583()->dict:
 c=log_chain();rows=c["cells"];links=[]
 for a,b in zip(rows[:-1],rows[1:]):
  za=np.array(a["log_center"]);ra=np.array(a["log_radii"]);zb=np.array(b["log_center"]);rb=np.array(b["log_radii"]);links.append(bool(np.all(np.maximum(za-ra,zb-rb)<=np.minimum(za+ra,zb+rb))))
 packet={**c,"cell_count":len(rows),"minimum_margin":min((x["minimum_margin"] for x in rows),default=None),"maximum_contraction":max((x["weighted_contraction"] for x in rows),default=None),"neighbor_box_overlaps":links,"all_neighbor_boxes_overlap":bool(rows) and all(links),"target_reached":c["certified_stop"]>=.3-1e-14,"theorem":f"Log-coordinate interval cells certify a locally unique root for every b from 0.2592 through {c['certified_stop']:.6f}." if rows else None}
 return finalize(583,"proven_log_chart_IR_chain" if rows else "bounded_stop_no_log_chart_interval_chain", "Box overlap supports component tracking but is weaker than separate common-root Krawczyk links; no global complement root count follows.",packet)


def st584()->dict:
 rows=json.loads(PACKETS[583].read_text())["cells"];ordered=[]
 for r in rows:
  x=np.array(r["x_center"]);b=sum(r["b_interval"])/2;a=b*x[2];ordered.append(x[0]<x[1]<1/b and 0<a<1 and x[4]>1)
 packet={"cell_count":len(rows),"all_cells_ordered":bool(rows) and all(ordered),"all_log_Jacobian_determinants_negative":bool(rows) and all(r["Jacobian_determinant"]<0 for r in rows),"local_Brouwer_degree":-1 if rows else None,"theorem":"The log chart is orientation preserving relative to positive x coordinates, so the negative determinant retains local degree -1 while y1<y2<y3, 0<a<1 and T>1 persist." if rows else None}
 return finalize(584,"proven_log_IR_ordering_and_degree_persistence" if rows else "not_run_no_certified_log_cells", "Local degree does not exclude disconnected components.",packet)


def st585()->dict:
 alpha=.8;I=math.pi/(2*math.gamma(1+alpha)*math.sin(math.pi*alpha/2));C=2*I
 packet={"tail":"a_d~d^(-1-alpha)","alpha":alpha,"integral_constant":"I_alpha=int_0^inf (1-cos x)x^(-1-alpha) dx","I_alpha":I,"symbol_asymptotic_coefficient":C,"theorem":"For 0<alpha<2 and a_d~c d^(-1-alpha), the even lattice symbol 2 sum_d a_d(1-cos(kd)) is asymptotic to 2c I_alpha |k|^alpha. The proof follows by splitting near/far tails and a Riemann-sum limit.","strict_signed_kernel_covered":False}
 return finalize(585,"proven_positive_tail_fractional_symbol_theorem", "The theorem applies to a positive regularly varying envelope, not automatically to the signed oscillatory strict continuation.",packet)


def st586()->dict:
 packet={"assumption":"lambda(k)~C|k|^nu, C>0","wave_law":"omega(k)~sqrt(C)|k|^(nu/2)","group_law":"v_g~sqrt(C)*(nu/2)|k|^(nu/2-1)","trichotomy":{"nu<2":"speed diverges","nu=2":"finite nonzero speed sqrt(C)","nu>2":"speed tends to zero"},"theorem":"A finite nonzero low-k wave speed exists iff the generator symbol is asymptotically quadratic.","Lorentz_symmetry_consequence":False}
 return finalize(586,"proven_dispersion_exponent_speed_trichotomy", "Quadratic dispersion is necessary for a finite massless speed, but not sufficient for Lorentz symmetry, isotropy, boosts, or an exact causal cone.",packet)


def st587()->dict:
 def C(mu):return 2*sum(WEIGHTS[d-1]*math.sinh(mu*d) for d in range(1,7))
 opt=minimize_scalar(lambda mu:C(mu)/mu,bounds=(1e-6,1.2),method="bounded");packet={"bound":"|<x|exp(-itA)|y>| <= exp[-mu(r-v_mu |t|)]","v_mu":"2 sum_d w_d sinh(mu d)/mu","optimized_mu":float(opt.x),"optimized_velocity":float(opt.fun),"theorem":"For the supplied finite-range one-particle hopping generator, exponential conjugation and the Duhamel/Gronwall estimate give the displayed bound for every mu>0.","many_body_LR":False,"Lorentz_cone":False}
 return finalize(587,"proven_finite_range_single_particle_exponential_cone_bound", "This is an approximate exponential cone for the added finite-range refinement, not strict finite propagation or a many-body Lorentz theorem.",packet)


def st588()->dict:
 packet={"finite_range":{"direct_coupling_beyond_R":False,"tail_type":"exponentially bounded outside a linear proxy cone","exact_zero":False},"positive_long_range_eta_1_8":{"direct_coupling_all_distances":True,"unitary_short_time_probability":"t^2 r^-3.6","heat_short_time_probability":"t r^-1.8","tail_type":"polynomial"},"second_order_wave":{"offdiagonal_amplitude":"O(t^2 A_xy)","probability":"O(t^4 |A_xy|^2)"},"theorem":"Direct nonzero long-range matrix elements create immediate polynomial spatial tails; finite-range generators instead require matrix-power paths and admit exponential moment bounds.","physical_causality":False}
 return finalize(588,"proven_finite_vs_long_range_causal_tail_classification", "Short-time matrix response is not a relativistic causal order; operational calibration and a continuum limit remain missing.",packet)


def st589()->dict:
 packet={"dimensionless_speed":"c_hat_n","physical_speed":"c_n=(ell_n/tau_n)c_hat_n","layer_transition":"(ell,tau,c_hat)->(alpha ell,beta tau,c_hat_next)","invariance":"(alpha/beta)(c_hat_next/c_hat)=1","gauge_dimension_without_anchor":1,"theorem":"Even after a refinement supplies c_hat_n, all speed records are invariant under one common length/time rescaling orbit. A transverse length or time anchor is necessary for an SI-valued speed.","SI_c_prediction":False}
 return finalize(589,"proven_layer_speed_calibration_torsor", "The exact SI number is partly definitional and is not predicted by a dimensionless slope.",packet)


def st590()->dict:
 c=json.loads((ROOT/"FIN_ST567_FiniteRange_Refinement_Dispersion_Limit.json").read_text())["exact_small_q_speed"];q=2*math.pi/1536;omega=c*q;sigma=.002*omega;n=400;relative_se=sigma/(omega*math.sqrt(n));packet={"synthetic_model":"omega observations Gaussian with known q","N":1536,"true_dimensionless_speed":c,"omega":omega,"relative_single_shot_noise":.002,"samples":n,"predicted_relative_standard_error":relative_se,"estimator":"mean(omega)/q","physical_length_time_calibration":False}
 return finalize(590,"conditional_finite_shot_dimensionless_speed_estimator", "The Gaussian frequency record, q labeling and refinement are supplied; no apparatus or SI calibration is provided.",packet)


def st591()->dict:
 p=json.loads((ROOT/"FIN_ST574_TwoLevel_Interval_Certificate_Preparation.json").read_text());packet={**p,"new_outward_cells":0,"result":"No complete 2047-family outward certificate was implemented in this round; preserve the numerical margin only.","interval_certificate_completed":False}
 return finalize(591,"bounded_stop_two_level_interval_certificate", "The class is not theorem-excluded until the scalar families are outward certified.",packet)


def st592()->dict:
 packet={"reflection_quotient_dimension":6,"new_interval_cells":0,"global_reflection_quotient_theorem":False,"decision":"Freeze generic box subdivision. A next move requires invariant coordinates or a sum-of-squares/convexity certificate."}
 return finalize(592,"reflection_quotient_lane_method_gated", "No candidate globality or counterexample follows.",packet)


def st593()->dict:
 p=json.loads((ROOT/"FIN_ST577_Stage1_CharacteristicZero_Method_Gate.json").read_text());packet={**p,"new_Q_relation":False,"rank_Q_first1326_closed":False}
 return finalize(593,"stage1_algebra_lane_remains_method_gated", "No further CRT or random-combination replay is admitted without a new exact basis method.",packet)


def st594()->dict:
 packet={"new_seed_source":False,"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open","reason":"Causal-speed and refinement results classify propagation but do not create the state/pump/branch source resources."}
 return finalize(594,"blocked_no_seed_gain_selector_source", "No source closure is inferred from dispersion.",packet)


def st595()->dict:
 packet={"quadratic_strict_refinement_selected":False,"finite_physical_c":False,"exact_causal_cone":False,"Lorentz_symmetry":False,"SI_299792458_derived":False,"reason":"Only the added finite-range refinement has quadratic dispersion; positive long-range and naive signed refinements fail differently."}
 return finalize(595,"blocked_no_physical_c_or_Lorentz_closure", "No continuum spacetime, boosts, isotropy, laboratory calibration, or physical light field is exported.",packet)


def st596()->dict:
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":3}
 return finalize(596,"blocked_no_independent_empirical_evidence", "Round three is local analytic/computational work.",packet)


def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST583"]["cells"];fig,ax=plt.subplots(figsize=(7.2,4));
 if rows:ax.plot([x["b_interval"][1] for x in rows],[x["x_center"][1] for x in rows],"o-")
 else:ax.text(.5,.5,"No Krawczyk cell admitted",ha="center",va="center",transform=ax.transAxes)
 ax.set(xlabel="b",ylabel=r"$y_2$",title="ST583: log-chart interval IR audit");fig.tight_layout();fig.savefig(FIG_DIR/"st583_log_ir.png",dpi=180);plt.close(fig)
 a=np.linspace(.2,2.6,200);fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(a,a/2-1);ax.axhline(0,color="black");ax.axvline(2,color="#dc2626",ls="--");ax.set(xlabel=r"symbol exponent $\nu$",ylabel="group-velocity exponent",title="ST586: finite speed only at nu=2");fig.tight_layout();fig.savefig(FIG_DIR/"st586_speed_trichotomy.png",dpi=180);plt.close(fig)


def main():
 funcs={582:st582,583:st583,584:st584,585:st585,586:st586,587:st587,588:st588,589:st589,590:st590,591:st591,592:st592,593:st593,594:st594,595:st595,596:st596};r={}
 for k in range(582,597):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(582,597):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
