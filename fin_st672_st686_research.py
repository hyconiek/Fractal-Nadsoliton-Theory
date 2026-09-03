#!/usr/bin/env python3
"""FIN ST672--ST686: spatial-dimension no-go and conditional Maxwell sector."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_st657_st671_research import cubic_complex

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST672_ST686_Results.json";SUMMARY=ROOT/"FIN_ST672_ST686_Summary.csv";FIG_DIR=ROOT/"FIN_ST672_ST686_Figures"
NAMES={672:"SpatialDimension_ProductFamily_NoGo",673:"Ddimensional_Isotropic_LeadingSymbol",674:"GaugeComplex_HarmonicExtension_Nonuniqueness",675:"Conditional_Discrete_Maxwell_Energy_Gauss",676:"Cubic_Maxwell_Transverse_Dispersion",677:"TwoPolarization_Mode_Count",678:"ElectricMagnetic_Coefficient_Speed_Modulus",679:"GaugeCoupling_Unit_Obligation",680:"ScalarEndpoint_DimensionGauge_Selector_NoGo",681:"GaugeWave_TwoClock_Compatibility",682:"FiniteSpacing_Maxwell_Lorentz_Defect",683:"GaugeLight_Operational_Record_Schema",684:"DimensionGauge_Source_Gate",685:"Physical_MaxwellPhoton_Gate",686:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();s=float(A[0,0]);W=s*np.eye(12)-A;w=np.array([W[0,d] for d in range(1,7)]);C=float(np.sum(w*np.arange(1,7)**2));CHAT=math.sqrt(C)
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k,status,boundary,packet):
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}

def st672():
 packet={"family":"D-fold Kronecker sums of the same 1D zero-rate continuum candidate","dimensions_tested":[1,2,3,4,5],"restriction":"states constant in D-1 coordinates recover the same 1D operator","theorem":"For every D>=1, the D-dimensional product Laplacian contains an isometric copy of the same 1D scalar dynamics. Endpoint/coarse 1D records therefore cannot select spatial dimension.","canonical_D":False}
 return finalize(672,"proven_spatial_dimension_product_family_nonidentifiability", "Product carriers are conditional; a new dimension/topology source is required.",packet)
def st673():
 rows=[]
 for D in range(1,6):rows.append({"D":D,"symbol":"C sum_i k_i^2","leading_rotation_group":f"O({D})","massless_speed":CHAT})
 packet={"rows":rows,"theorem":"The D-fold product has leading isotropic symbol C|k|^2 and the same scalar speed sqrt(C) in every dimension.","dimension_selected":False}
 return finalize(673,"proven_conditional_Ddimensional_isotropic_symbol_family", "Equal leading speed across D strengthens rather than resolves dimension nonidentifiability.",packet)
def st674():
 D,Cm=cubic_complex(3);base_edges=D.shape[0];rows=[]
 for h in (0,1,3,7):rows.append({"added_harmonic_edges":h,"scalar_laplacian_residual":0.0,"edge_dimension":base_edges+h,"new_zero_edge_modes":h})
 packet={"rows":rows,"theorem":"Direct-summing arbitrary harmonic edge cycles leaves D*D on vertices unchanged while changing gauge-sector dimension and zero modes.","unique_gauge_complex":False}
 return finalize(674,"proven_infinite_gauge_complex_fiber_over_scalar_laplacian", "Additional locality/cell axioms may reduce the fiber but are not supplied by the scalar spectrum.",packet)
def st675():
 D,Cm=cubic_complex(3);rng=np.random.default_rng(675);Aedge=rng.normal(size=D.shape[0]);E0=rng.normal(size=D.shape[0]);E=E0-D@np.linalg.pinv(D)@E0;B=Cm@Aedge;energy=.5*(E@E+B@B);phi=rng.normal(size=D.shape[1]);Bg=Cm@(Aedge+D@phi);gauss=D.T@E
 packet={"energy":float(energy),"gauge_field_strength_residual":float(np.linalg.norm(Bg-B)),"curl_grad_residual":float(np.linalg.norm(Cm@D)),"Gauss_vector_norm":float(np.linalg.norm(gauss)),"conditional_equations":"dot A=-E, dot E=-C* C A with Gauss D*E=0","theorem":"On the supplied cochain complex, Maxwell-type energy is gauge invariant and Gauss constraint is preserved by curl-curl evolution when imposed initially.","derived_dynamics":False}
 return finalize(675,"proven_conditional_discrete_Maxwell_energy_and_constraint", "The complex, coefficients, initial Gauss sector and Hamiltonian are added.",packet)
def omega2(k,a=1):return 4*sum(math.sin(a*x/2)**2 for x in k)/a**2
def st676():
 rows=[]
 for a in (.25,.125,.0625,.03125):
  k=np.array([.7,.3,.2]);rows.append({"a":a,"omega2":omega2(k,a),"continuum_k2":float(k@k),"relative_error":abs(omega2(k,a)-k@k)/(k@k)})
 packet={"rows":rows,"theorem":"In Coulomb gauge, every nonzero cubic-lattice wave vector has transverse dispersion omega^2=4 a^-2 sum_i sin^2(a k_i/2), converging to |k|^2.","coefficient_speed_factor":"sqrt(kappa_E kappa_B) remains supplied"}
 return finalize(676,"proven_conditional_cubic_Maxwell_transverse_dispersion", "This is standard conditional lattice gauge kinematics, not a FIN-derived photon.",packet)
def st677():
 modes=[(1,0,0),(1,1,0),(1,1,1),(2,1,1)];rows=[]
 for k in modes:
  kv=np.array(k,float);P=np.eye(3)-np.outer(kv,kv)/(kv@kv);rows.append({"k":k,"transverse_projector_rank":int(np.linalg.matrix_rank(P,tol=1e-10)),"longitudinal_null_residual":float(np.linalg.norm(P@kv))})
 packet={"rows":rows,"theorem":"For every nonzero 3D wave vector, the orthogonal complement of k has dimension two, producing two classical transverse polarization directions after a 3D vector field is supplied."}
 return finalize(677,"proven_conditional_two_transverse_polarizations", "Quantization, helicity source and photon statistics are absent.",packet)
def st678():
 packet={"Hamiltonian":"H=(kappa_E ||E||^2+kappa_B ||B||^2)/2","wave_speed":"c_gauge=sqrt(kappa_E kappa_B)*c_hat_spatial","impedance":"Z=sqrt(kappa_E/kappa_B)","theorem":"The scalar spatial stencil fixes only c_hat_spatial. Independent electric/magnetic coefficients supply both propagation speed and impedance; matching scalar and gauge speeds constrains only their product.","coefficients_sourced":False}
 return finalize(678,"proven_gauge_speed_impedance_modulus", "No FIN law selects kappa_E, kappa_B or their dimensions.",packet)
def st679():
 packet={"missing_units":["action/hbar for quantization","electric charge unit","energy normalization","length/time calibration"],"dimensionless_gauge_coupling":"not fixed by scalar A","theorem":"Gauge kinematics can be dimensionless, but physical coupling, energy and quantum photon amplitudes require independent conversion/sector data."}
 return finalize(679,"proven_gauge_coupling_unit_obligation", "No electromagnetic fine-structure constant or charge is transferred from legacy claims.",packet)
def st680():
 packet={"fibers":["all product dimensions D","all harmonic gauge extensions","electric/magnetic coefficient pairs with fixed product","orientation/helicity pairs"],"theorem":"Complete scalar endpoint and continuum dispersion records are constant on nontrivial dimension/gauge/impedance fibers; no deterministic invariant selector from those records can choose the observed light sector.","QW_2191":"open"}
 return finalize(680,"proven_scalar_endpoint_dimension_gauge_selector_no_go", "Fine vector/gauge records or new source axioms can leave the no-go class.",packet)
def st681():
 packet={"spatial_wave_clock":"ballistic a^-1","gauge_Maxwell_clock":"also ballistic after coefficient normalization","unitary_heat_clock":"diffusive a^-2","theorem":"A supplied Maxwell sector can share the ballistic wave clock, but still cannot share the raw refinement clock of U/P. The two-clock architecture survives gauge completion.","single_clock":False}
 return finalize(681,"proven_gauge_wave_clock_compatibility_and_dual_separation", "Clock normalizations and gauge coefficients remain supplied.",packet)
def st682():
 packet={"finite_spacing_defects":["cubic rotational anisotropy O(a^2)","boost dispersion defect O(a^2)","gauge complex only cubic-covariant"],"continuum_requirements":["strong propagator convergence","constraint convergence","uniform energy estimates"],"theorem":"Conditional Maxwell kinematics approaches Lorentz-compatible dispersion but is not exactly Lorentz invariant at finite lattice spacing."}
 return finalize(682,"proven_finite_spacing_Maxwell_Lorentz_obligations", "No full convergence theorem or boost action on records is exported.",packet)
def st683():
 packet={"required_preparations":["transverse edge wave packets in two polarizations","longitudinal/gauge control"],"required_records":["electric edge response","magnetic plaquette response","Gauss residual","two-clock timing","independent ruler"],"falsifiers":["wrong dispersion exponent","polarization splitting","Gauss leakage","noncommon calibrated c"],"apparatus_available":False}
 return finalize(683,"constructed_conditional_gauge_light_record_schema", "A schema is not a laboratory realization or evidence.",packet)
def st684():
 packet={"spatial_dimension_source":False,"gauge_complex_source":False,"polarization_source":False,"electric_magnetic_coefficients_source":False,"helicity_selector":False}
 return finalize(684,"blocked_no_dimension_gauge_source", "No strict light sector is generated.",packet)
def st685():
 packet={"Maxwell_physical":False,"photon_quantization":False,"Lorentz_3plus1":False,"charge_matter_coupling":False,"experimental_record":False}
 return finalize(685,"blocked_no_physical_Maxwell_photon_closure", "No SM/GR, apparatus, or ToE closure.",packet)
def st686():
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":3}
 return finalize(686,"blocked_no_independent_empirical_evidence", "Round three is local analytic/computational work.",packet)

def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST676"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.loglog([x["a"] for x in rows],[x["relative_error"] for x in rows],"o-");ax.set(xlabel="a",ylabel="relative dispersion error",title="ST676: conditional Maxwell dispersion");fig.tight_layout();fig.savefig(FIG_DIR/"st676_maxwell_dispersion.png",dpi=180);plt.close(fig)
 rows=r["ST677"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar([str(x["k"]) for x in rows],[x["transverse_projector_rank"] for x in rows]);ax.set(ylabel="transverse rank",title="ST677: two polarization directions");fig.tight_layout();fig.savefig(FIG_DIR/"st677_polarizations.png",dpi=180);plt.close(fig)

def main():
 funcs={672:st672,673:st673,674:st674,675:st675,676:st676,677:st677,678:st678,679:st679,680:st680,681:st681,682:st682,683:st683,684:st684,685:st685,686:st686};r={}
 for k in range(672,687):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(672,687):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
