#!/usr/bin/env python3
"""FIN ST642--ST656: locality selection, strong continuum, and light-sector obligations."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import root

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N,SEED
from fin_st447_st461_research import regularized_ir_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST642_ST656_Results.json";SUMMARY=ROOT/"FIN_ST642_ST656_Summary.csv";FIG_DIR=ROOT/"FIN_ST642_ST656_Figures"
NAMES={642:"FixedRange_Locality_ZeroRate_Selection",643:"Coarse_Green_Schur_Selector_NoGo",644:"BandLimited_Strong_Wave_Convergence",645:"FiniteLattice_ConeTail_Convergence_Audit",646:"Conditional_3D_Isotropic_Symbol",647:"Discrete_LorentzBreaking_Correction",648:"ScalarKernel_GaugeLight_NoGo",649:"TwoClock_Multilevel_Identifiability",650:"Noncircular_LengthTime_Anchor_Theorem",651:"LongRange_DirectTail_LowerBound",652:"SignedRefinement_WaveInstability",653:"IR_AffinePredictor_Correlation_Audit",654:"Transition_Algebra_StateMap_Gate",655:"Strict_Physics_Source_Gate",656:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();s=float(A[0,0]);W=s*np.eye(12)-A;WEIGHTS=np.array([W[0,d] for d in range(1,7)]);C=float(np.sum(WEIGHTS*np.arange(1,7)**2));M4=float(np.sum(WEIGHTS*np.arange(1,7)**4));CHAT=math.sqrt(C)
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_,complex)):return x.item() if not isinstance(x,complex) else {"real":x.real,"imag":x.imag}
 return x
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k,status,boundary,packet):
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}
def lam(q):return float(2*np.sum(WEIGHTS*(1-np.cos(np.arange(1,7)*q))))

def st642():
 packet={"class":"nonnegative dyadic antipodal rates q_N","fixed_range_premise":"there exists R independent of N with no edges longer than R","conclusion":"q_N=0 whenever N>R","self_similar_scale_shift_premise":"the same local rule applies at every dyadic level","strengthened_conclusion":"q_N=0 at every level","theorem":"A new antipodal q_N edge has range N. Uniform fixed-range locality forbids it at all sufficiently fine levels; exact scale-shift self-similarity propagates the zero choice to the complete tower.","strict_locality_premise_sourced":False}
 return finalize(642,"proven_conditional_fixed_range_zero_rate_selection", "The locality/self-similarity premises are new; early finite-level q values remain free without scale-shift naturality.",packet)
def st643():
 packet={"inputs":"all coarse Green resolvents, unitary/heat functions, and coarse preparations/effects","hidden_family":"A_2N(q) with exact intertwining","theorem":"Every functional of complete coarse spectral/Green data is constant on the q fiber and therefore cannot select a refinement rate or fine speed.","coarse_Green_selector_exists":False,"required_escape":"fine instrument, noncoarse locality axiom, or source law"}
 return finalize(643,"proven_coarse_Green_Schur_refinement_selector_no_go", "Green--Schur response preservation strengthens blindness; it does not make fine locality observable.",packet)
def st644():
 L=12.;modes=range(1,9);rows=[]
 for NN in (96,192,384,768,1536,3072):
  a=L/NN;errs=[]
  for m in modes:
   k=2*math.pi*m/L;om=math.sqrt(lam(a*k))/a;errs.append(abs(om-CHAT*abs(k)))
  rows.append({"N":NN,"maximum_frequency_error_first8":max(errs),"L2_operator_error_band":math.sqrt(sum(e*e for e in errs))})
 packet={"fixed_mode_band":list(modes),"rows":rows,"theorem":"On every fixed finite Fourier band, rescaled wave frequencies converge uniformly to sqrt(C)|k|; hence the cosine and sine wave propagators converge strongly and uniformly on compact time intervals in that band.","full_Sobolev_space_theorem":False}
 return finalize(644,"proven_bandlimited_strong_wave_propagator_convergence", "Extension to full energy/Sobolev spaces needs a uniform high-mode stability/density argument.",packet)
def wave_tail(NN,L=12.,t=.5,margin=1.05):
 a=L/NN;q=2*np.pi*np.arange(NN)/NN;lv=sum(2*WEIGHTS[d-1]*(1-np.cos(d*q)) for d in range(1,7));amp=np.fft.ifft(np.cos(t*np.sqrt(lv)/a)).real;dist=np.minimum(np.arange(NN),NN-np.arange(NN))*a;mask=dist>CHAT*t*margin;return float(np.sum(amp[mask]**2)/np.sum(amp**2))
def st645():
 rows=[{"N":n,"tail":wave_tail(n)} for n in (384,768,1536,3072,6144,12288,24576)];rat=[rows[i+1]["tail"]/rows[i]["tail"] for i in range(len(rows)-1)]
 packet={"rows":rows,"successive_ratios":rat,"last_tail":rows[-1]["tail"],"analytic_uniform_bound":False,"result":"The tail continues to collapse, but nonmonotone preasymptotic ratios prevent promotion to a simple geometric bound."}
 return finalize(645,"strong_extended_wave_cone_tail_evidence", "A stationary-phase/contour theorem uniform in N is still required.",packet)
def symbol3(k,a):return sum(lam(a*x)/a**2 for x in k)
def st646():
 rows=[]
 for a in (.25,.125,.0625,.03125):
  axis=np.array([1.,0,0]);diag=np.ones(3)/math.sqrt(3);la=symbol3(axis,a);ld=symbol3(diag,a);rows.append({"a":a,"axis":la,"diagonal":ld,"relative_anisotropy":abs(la-ld)/(C)})
 packet={"construction":"3D Kronecker sum of three identical zero-rate 1D stencils","symbol":"sum_i a^-2 lambda(a k_i)","rows":rows,"theorem":"The leading symbol is C(kx^2+ky^2+kz^2); rotational anisotropy first appears at O(a^2 sum k_i^4).","strict_3D_carrier_source":False}
 return finalize(646,"proven_conditional_3D_isotropic_leading_symbol", "The Cartesian product carrier is added and does not derive observed spatial dimension or rotations from FIN.",packet)
def st647():
 packet={"dispersion_expansion":"omega^2=C|k|^2-(M4/12)a^2 sum_i k_i^4+O(a^4|k|^6)","C":C,"M4":M4,"Lorentz_breaking_order":"O(a^2 k^4)","theorem":"The conditional 3D product refinement recovers rotational/Lorentz scalar dispersion only at leading order; lattice corrections select axes and break boosts at finite spacing.","exact_finite_lattice_Lorentz":False}
 return finalize(647,"proven_leading_Lorentz_dispersion_with_lattice_breaking", "Convergence of dispersion does not construct a boost action on the full state/observable algebra.",packet)
def st648():
 packet={"available":"one real scalar Laplacian/wave field","missing":["vector or one-form module","exterior derivative complex","U(1) connection","Gauss constraint","transverse polarization sectors","gauge-invariant field-strength observables"],"theorem":"A scalar radial kernel and its functional calculus do not uniquely generate Maxwell data; many inequivalent gauge complexes share the same scalar Laplacian spectrum.","physical_light_sector":False}
 return finalize(648,"proven_scalar_kernel_gauge_light_nonuniqueness", "A supplied discrete differential complex can build a gauge model, but it is additional structure.",packet)
def st649():
 rng=np.random.default_rng(SEED+649);a=12/np.array([96,192,384,768,1536],float);TU=a**-2*(1+rng.normal(0,2e-3,len(a)));TW=a**-1*(1+rng.normal(0,2e-3,len(a)));zu=-np.polyfit(np.log(a),np.log(TU),1)[0];zw=-np.polyfit(np.log(a),np.log(TW),1)[0]
 packet={"levels":len(a),"relative_noise":.002,"estimated_diffusive_exponent":float(zu),"estimated_wave_exponent":float(zw),"true_exponents":[2,1],"separation":float(abs(zu-zw)),"theorem_scope":"synthetic log-regression identifiability"}
 return finalize(649,"strong_synthetic_two_clock_exponent_identifiability", "No physical clock apparatus or timing record is supplied.",packet)
def st650():
 packet={"circular_protocol":"define ell_*=c_target tau_* and then report c_target","rejected":True,"noncircular_requirement":"ell_* and tau_* must be calibrated by records not using the tested propagation speed","theorem":"A speed prediction is falsifiable only when at least one length/time ratio is fixed independently of that speed. Otherwise any dimensionless c_hat can be mapped to the target SI number.","independent_anchor_available":False}
 return finalize(650,"proven_noncircular_speed_anchor_requirement", "The theorem identifies an experimental obligation; it does not provide the anchor.",packet)
def st651():
 packet={"tail":"w_r~r^-1.8","heat":"P_t(x,y)=t w_r+O(t^2)","unitary_Born":"|U_t(x,y)|^2=t^2 w_r^2+O(t^3)","consequences":{"heat_spatial_lower_order":"r^-1.8","unitary_probability_lower_order":"r^-3.6"},"theorem":"For every directly coupled distance, the first nonzero short-time coefficient fixes a polynomial tail; no uniform spatial bound with faster asymptotic decay can hold at sufficiently small fixed t without cancellation/removal of the edge.","strict_positive_all_distance_extension":False}
 return finalize(651,"proven_direct_edge_polynomial_tail_lower_order", "The theorem applies to a positive long-range class; oscillatory signed cancellation needs separate analysis.",packet)
def signed_min(n=48):
 d=np.arange(1,n//2+1);w=np.cos(.18575*d+.1625)/(1+d**1.8);vals=[]
 for m in range(n):q=2*np.pi*m/n;vals.append(float(2*np.sum(w[:-1]*(1-np.cos(d[:-1]*q)))+w[-1]*(1-np.cos(d[-1]*q))))
 return min(vals)
def st652():
 lmin=signed_min();times=(1,2,5,10);rows=[{"t":t,"unstable_wave_factor":math.cosh(t*math.sqrt(-lmin)),"heat_growth_factor":math.exp(t*(-lmin))} for t in times]
 packet={"N":48,"minimum_eigenvalue":lmin,"rows":rows,"theorem":"A negative stiffness eigenvalue produces cosh(t sqrt(-lambda)) growth in the second-order wave channel and exp(t|lambda|) growth in raw heat evolution.","stable_causal_wave":False}
 return finalize(652,"proven_naive_signed_refinement_dynamic_instability", "Krein/unitary reformulations may remain mathematically valid but change positivity, probability and causal interpretation.",packet)
def st653():
 prev=json.loads((ROOT/"FIN_ST579_LogChart_IR_Conditioning_Audit.json").read_text());rows=prev["rows"];z0=None;out=[]
 for rr in rows:
  b=rr["b"];sol=root(lambda z:regularized_ir_float(np.exp(z),b),np.log(np.array([.0285,3.,.043,.09,3.3])) if z0 is None else z0);z=sol.x
  if z0 is not None:
   pred=zprev+slope*(b-bprev);out.append({"b":b,"constant_residual":float(np.linalg.norm(regularized_ir_float(np.exp(zprev),b),np.inf)),"affine_predictor_residual":float(np.linalg.norm(regularized_ir_float(np.exp(pred),b),np.inf)),"actual_log_error":float(np.linalg.norm(z-pred,np.inf))})
  if z0 is not None:slope=(z-zprev)/(b-bprev)
  else:slope=np.zeros(5)
  zprev=z;bprev=b;z0=z
 packet={"rows":out,"median_residual_improvement":float(np.median([x["constant_residual"]/max(x["affine_predictor_residual"],1e-300) for x in out])),"interval_correlated_predictor":False}
 return finalize(653,"strong_affine_predictor_IR_correlation_evidence", "Floating predictor improvement does not retain interval z--b correlation; no Krawczyk theorem follows.",packet)
def st654():
 packet={"transition_global_identity":False,"degree8_rank_Q_closed":False,"new_method":False,"status":"gated"}
 return finalize(654,"transition_algebra_state_map_remains_gated", "No propagation result closes these independent lanes.",packet)
def st655():
 packet={"zero_rate_source":False,"3D_carrier_source":False,"two_clock_source":False,"independent_anchor":False,"gauge_light_source":False,"QW_2191":"open"}
 return finalize(655,"blocked_no_strict_continuum_light_source_package", "No physical light, SM/GR, L_total, apparatus, or ToE closure.",packet)
def st656():
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":1}
 return finalize(656,"blocked_no_independent_empirical_evidence", "Round one is local analytic/computational work.",packet)

def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST644"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["N"] for x in rows],[x["maximum_frequency_error_first8"] for x in rows],"o-");ax.set(xlabel="N",ylabel="max frequency error",title="ST644: band-limited wave convergence");fig.tight_layout();fig.savefig(FIG_DIR/"st644_wave_convergence.png",dpi=180);plt.close(fig)
 rows=r["ST646"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["a"] for x in rows],[x["relative_anisotropy"] for x in rows],"o-");ax.set(xlabel="lattice spacing a",ylabel="axis/diagonal anisotropy",title="ST646: conditional 3D isotropy");fig.tight_layout();fig.savefig(FIG_DIR/"st646_isotropy.png",dpi=180);plt.close(fig)
 rows=r["ST652"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.semilogy([x["t"] for x in rows],[x["unstable_wave_factor"] for x in rows],"o-",label="wave cosh");ax.semilogy([x["t"] for x in rows],[x["heat_growth_factor"] for x in rows],"s-",label="heat growth");ax.legend();ax.set(xlabel="t",ylabel="growth",title="ST652: signed-refinement instability");fig.tight_layout();fig.savefig(FIG_DIR/"st652_instability.png",dpi=180);plt.close(fig)

def main():
 funcs={642:st642,643:st643,644:st644,645:st645,646:st646,647:st647,648:st648,649:st649,650:st650,651:st651,652:st652,653:st653,654:st654,655:st655,656:st656};r={}
 for k in range(642,657):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(642,657):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
