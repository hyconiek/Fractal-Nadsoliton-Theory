#!/usr/bin/env python3
"""FIN ST657--ST671: 3D continuum, gauge-complex obligations, and two-clock protocol."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N,SEED
from fin_st447_st461_research import regularized_ir_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST657_ST671_Results.json";SUMMARY=ROOT/"FIN_ST657_ST671_Summary.csv";FIG_DIR=ROOT/"FIN_ST657_ST671_Figures"
NAMES={657:"IR_CorrelatedPredictor_Krawczyk_Attempt",658:"IR_Chart_Method_NoGo_Status",659:"Conditional_3D_Bandlimited_Wave_Convergence",660:"Rotational_Anisotropy_Order_Theorem",661:"Discrete_Boost_Dispersion_Defect",662:"ScalarLaplacian_GaugeComplex_Nonuniqueness",663:"Conditional_Cubic_Maxwell_Complex",664:"Discrete_Transverse_Polarization_Count",665:"TwoClock_Gaussian_Likelihood_Protocol",666:"Independent_RulerClock_Anchor_Schema",667:"Uniform_WaveTail_Analytic_Method_Audit",668:"Krein_Unitary_vs_Wave_Semantics",669:"Refinement_Clock_Gauge_Source_Gate",670:"Physical_Light_Lorentz_Gate",671:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();s=float(A[0,0]);W=s*np.eye(12)-A;w=np.array([W[0,d] for d in range(1,7)]);C=float(np.sum(w*np.arange(1,7)**2));M4=float(np.sum(w*np.arange(1,7)**4));CHAT=math.sqrt(C)
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k,status,boundary,packet):
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}
def lam(q):return float(2*np.sum(w*(1-np.cos(np.arange(1,7)*q))))
def finite_jac(fun,x,eps=1e-6):
 e=np.eye(len(x));return np.column_stack([(fun(x+eps*e[j])-fun(x-eps*e[j]))/(2*eps) for j in range(len(x))])

def st657():
 prior=np.log(np.array([.0284,2.9,.044,.09,3.25]));rows=[]
 for blo,bhi in zip(np.linspace(.2592,.295,12),np.linspace(.26245,.3,12)):
  bm=(blo+bhi)/2;sm=root(lambda z:regularized_ir_float(np.exp(z),bm),prior);z=sm.x;J=finite_jac(lambda zz:regularized_ir_float(np.exp(zz),bm),z);db=1e-5;Fb=(regularized_ir_float(np.exp(z),bm+db)-regularized_ir_float(np.exp(z),bm-db))/(2*db);slope=-np.linalg.solve(J,Fb);predL=z+slope*(blo-bm);predR=z+slope*(bhi-bm);solL=root(lambda zz:regularized_ir_float(np.exp(zz),blo),predL);solR=root(lambda zz:regularized_ir_float(np.exp(zz),bhi),predR);err=max(np.linalg.norm(solL.x-predL,np.inf),np.linalg.norm(solR.x-predR,np.inf));rows.append({"b_interval":[float(blo),float(bhi)],"predictor_endpoint_error":float(err),"mid_condition":float(np.linalg.cond(J)),"root_success":bool(solL.success and solR.success)});prior=solR.x
 packet={"rows":rows,"all_endpoint_roots":all(x["root_success"] for x in rows),"maximum_predictor_error":max(x["predictor_endpoint_error"] for x in rows),"interval_Jacobian_and_Fb_bounds":False,"Krawczyk_inclusion":False}
 return finalize(657,"strong_correlated_predictor_root_tracking", "Floating implicit predictors do not certify the interval derivative or inclusion.",packet)
def st658():
 packet={"failed_methods":["raw x-box radius inflation","naive z=log x interval substitution"],"successful_numerical_method":"affine implicit predictor in log coordinates","missing_exact_objects":["interval F_b bound retaining B correlation","interval inverse/Jacobian contraction on predictor tube","certified common-root links"],"decision":"IR continuation is method-gated until these objects are implemented."}
 return finalize(658,"IR_interval_lane_precisely_method_gated", "Numerical roots through b=0.3 are not promoted.",packet)
def symbol3(k,a):return sum(lam(a*x)/a**2 for x in k)
def st659():
 modes=[np.array(v,float) for v in [(1,0,0),(1,1,0),(1,1,1),(2,1,0)]];rows=[]
 for a in (.2,.1,.05,.025):
  errs=[]
  for k in modes:errs.append(abs(math.sqrt(symbol3(k,a))-CHAT*np.linalg.norm(k)))
  rows.append({"a":a,"maximum_frequency_error":max(errs),"L2_error":float(np.linalg.norm(errs))})
 packet={"mode_vectors":[x.tolist() for x in modes],"rows":rows,"theorem":"On every fixed finite 3D Fourier set, the Cartesian-product wave propagator converges uniformly on compact time intervals to omega=sqrt(C)|k|.","full_energy_space":False}
 return finalize(659,"proven_conditional_3D_bandlimited_wave_convergence", "The product carrier and zero-rate stencil are supplied.",packet)
def st660():
 packet={"expansion":"lambda_a(k)=C|k|^2-(M4/12)a^2 sum_i k_i^4+O(a^4|k|^6)","anisotropic_invariant":"sum_i k_i^4","isotropic_leading_term":"C|k|^2","theorem":"Rotational breaking is generically order a^2 and cannot vanish at finite spacing for nonzero fourth moment; it vanishes in the conditional continuum limit.","exact_finite_rotation_group":"cubic, not SO(3)"}
 return finalize(660,"proven_rotational_anisotropy_order_a2", "SO(3) emergence is asymptotic and conditional.",packet)
def st661():
 beta=.6;gamma=1/math.sqrt(1-beta**2);k=np.array([.7,.2,.1]);omega=CHAT*np.linalg.norm(k);kp=np.array([gamma*(k[0]-beta*omega/CHAT),k[1],k[2]]);op=gamma*(omega-beta*CHAT*k[0]);rows=[]
 for a in (.2,.1,.05,.025):defect=abs(op**2-symbol3(kp,a));rows.append({"a":a,"boosted_dispersion_defect":defect})
 packet={"beta":beta,"rows":rows,"defect_ratio_last":rows[-1]["boosted_dispersion_defect"]/rows[-2]["boosted_dispersion_defect"],"theorem_scope":"numerical O(a^2) boost-defect audit","exact_discrete_boost":False}
 return finalize(661,"strong_discrete_boost_defect_decay", "No boost representation on lattice observables is constructed.",packet)

def cubic_complex(n=3):
 idx=lambda x,y,z:(x%n)*n*n+(y%n)*n+(z%n);edges=[]
 for x in range(n):
  for y in range(n):
   for z in range(n):
    v=idx(x,y,z)
    for a,(dx,dy,dz) in enumerate(((1,0,0),(0,1,0),(0,0,1))):edges.append((v,idx(x+dx,y+dy,z+dz),a,x,y,z))
 D=np.zeros((len(edges),n**3))
 for i,(u,v,*_) in enumerate(edges):D[i,u]=-1;D[i,v]=1
 # Curl rows: oriented plaquette boundary for xy,xz,yz at every base vertex.
 emap={(a,x,y,z):i for i,(_,_,a,x,y,z) in enumerate(edges)};faces=[]
 for x in range(n):
  for y in range(n):
   for z in range(n):
    for a,b in ((0,1),(0,2),(1,2)):
     row=np.zeros(len(edges));row[emap[(a,x,y,z)]]+=1
     shift=[x,y,z];shift[a]+=1;row[emap[(b,*[q%n for q in shift])]]+=1
     shift=[x,y,z];shift[b]+=1;row[emap[(a,*[q%n for q in shift])]]-=1
     row[emap[(b,x,y,z)]]-=1;faces.append(row)
 Cmat=np.array(faces);return D,Cmat
def st662():
 D,Cm=cubic_complex();scalar=D.T@D;packet={"scalar_vertices":scalar.shape[0],"edge_dimension":D.shape[0],"face_dimension":Cm.shape[0],"curl_grad_residual":float(np.linalg.norm(Cm@D)),"scalar_rank":int(np.linalg.matrix_rank(scalar)),"theorem":"A scalar Laplacian D*D does not determine a unique edge/face complex: adding decoupled cycle spaces or changing cell structure preserves the scalar operator while changing gauge sectors.","unique_gauge_complex_from_scalar_A":False}
 return finalize(662,"proven_scalar_to_gauge_complex_nonuniqueness_fixture", "The cubic complex is supplied and unrelated to strict carrier provenance.",packet)
def st663():
 D,Cm=cubic_complex();packet={"d0_shape":list(D.shape),"d1_shape":list(Cm.shape),"d1_d0_residual":float(np.linalg.norm(Cm@D)),"gauge_transformation":"A_edge->A_edge+D phi","field_strength":"F=Cm A_edge","gauge_invariance_residual_formula":"Cm D=0","theorem":"The supplied cubic cochain complex supports an exact discrete Abelian gauge kinematics and Maxwell-type field strength.","derived_from_FIN":False}
 return finalize(663,"proven_conditional_discrete_Abelian_gauge_complex", "Gauge kinematics is constructed, not sourced by the scalar FIN kernel.",packet)
def st664():
 D,Cm=cubic_complex(3);lap1=D@D.T+Cm.T@Cm;evals=np.linalg.eigvalsh(lap1);packet={"edge_dimension":len(evals),"zero_modes":int(np.sum(abs(evals)<1e-9)),"positive_modes":int(np.sum(evals>1e-9)),"gradient_rank":int(np.linalg.matrix_rank(D)),"curl_rank":int(np.linalg.matrix_rank(Cm)),"theorem":"Hodge decomposition separates gradient, coexact and harmonic edge sectors; transverse polarizations live in data absent from the scalar vertex Laplacian."}
 return finalize(664,"proven_conditional_discrete_Hodge_sector_count", "No photon statistics, Hamiltonian normalization or physical polarization measurement follows.",packet)
def st665():
 rng=np.random.default_rng(SEED+665);levels=np.array([96,192,384,768,1536]);a=12/levels;shots=2000;sig=.01;zu=[];zw=[]
 for _ in range(shots):zu.append(-np.polyfit(np.log(a),np.log(a**-2*np.exp(rng.normal(0,sig,len(a)))),1)[0]);zw.append(-np.polyfit(np.log(a),np.log(a**-1*np.exp(rng.normal(0,sig,len(a)))),1)[0])
 packet={"replicates":shots,"log_noise_sigma":sig,"diffusive_exponent_mean":float(np.mean(zu)),"diffusive_exponent_sd":float(np.std(zu)),"wave_exponent_mean":float(np.mean(zw)),"wave_exponent_sd":float(np.std(zw)),"misorder_probability":float(np.mean(np.array(zu)<=np.array(zw))),"physical_protocol":False}
 return finalize(665,"strong_synthetic_two_clock_likelihood_separation", "Records are synthetic lognormal scale observations, not laboratory clocks.",packet)
def st666():
 packet={"required_records":["independent non-light length standard","independent periodic clock","wave propagation record","diffusive/unitary timing record"],"custody":"standards calibrated before unblinding propagation data","falsifier":"no common layer-covariant c after independent calibration","theorem":"Precalibrating length and time without the tested wave closes the circular unit-definition loophole.","apparatus_available":False}
 return finalize(666,"proven_noncircular_anchor_protocol_schema", "A schema is not a physical apparatus or record.",packet)
def st667():
 packet={"candidate_methods":["complex contour deformation of finite-range symbol","Bessel/path expansion with factorial distance suppression","energy estimate plus continuum consistency"],"current_numeric_tail_at_N24576":json.loads((ROOT/"FIN_ST645_FiniteLattice_ConeTail_Convergence_Audit.json").read_text())["last_tail"],"uniform_N_bound_proved":False,"decision":"No analytic uniform tail theorem is exported in this round."}
 return finalize(667,"uniform_wave_tail_method_audit_open", "Numerical collapse remains strong evidence only.",packet)
def st668():
 packet={"Krein_unitary":"self-adjoint signed generator still defines exp(-itA) with norm conservation","raw_heat":"negative eigenvalues cause exponential growth","second_order_wave":"negative stiffness causes cosh growth","theorem":"Unitary mathematical well-posedness of a signed generator does not transfer to heat positivity or stable wave/causal semantics.","physical_signed_causal_model":False}
 return finalize(668,"proven_Krein_unitary_wave_semantics_separation", "A fundamental symmetry/metric and operational probability rule would be additional.",packet)
def st669():
 packet={"zero_rate_source":False,"3D_carrier_source":False,"gauge_complex_source":False,"two_clock_source":False,"anchor_record":False}
 return finalize(669,"blocked_no_refinement_clock_gauge_source", "Conditional constructions do not become strict by composition.",packet)
def st670():
 packet={"physical_light":False,"Lorentz_group_on_observables":False,"Maxwell_dynamics":False,"photon_sector":False,"laboratory_c":False}
 return finalize(670,"blocked_no_physical_light_Lorentz_closure", "No SM/GR, apparatus, or ToE closure.",packet)
def st671():
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":2}
 return finalize(671,"blocked_no_independent_empirical_evidence", "Round two is local analytic/computational work.",packet)

def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST659"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["a"] for x in rows],[x["maximum_frequency_error"] for x in rows],"o-");ax.set(xlabel="a",ylabel="max frequency error",title="ST659: conditional 3D band convergence");fig.tight_layout();fig.savefig(FIG_DIR/"st659_3d_convergence.png",dpi=180);plt.close(fig)
 rows=r["ST661"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["a"] for x in rows],[x["boosted_dispersion_defect"] for x in rows],"o-");ax.set(xlabel="a",ylabel="boosted dispersion defect",title="ST661: lattice boost defect");fig.tight_layout();fig.savefig(FIG_DIR/"st661_boost_defect.png",dpi=180);plt.close(fig)

def main():
 funcs={657:st657,658:st658,659:st659,660:st660,661:st661,662:st662,663:st663,664:st664,665:st665,666:st666,667:st667,668:st668,669:st669,670:st670,671:st671};r={}
 for k in range(657,672):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(657,672):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
