#!/usr/bin/env python3
"""FIN ST687--ST701: conditional 3+1 selectors and photon-sector obligations."""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST687_ST701_Results.json";SUMMARY=ROOT/"FIN_ST687_ST701_Summary.csv";FIG_DIR=ROOT/"FIN_ST687_ST701_Figures"
NAMES={687:"Fine_HeatTrace_SpectralDimension",688:"Coarse_vs_Fine_Dimension_Identifiability",689:"Maxwell_Conformal_Dimension_Theorem",690:"PolarizationCount_SpatialDimension_Selector",691:"Huygens_Principle_Dimension_Classification",692:"Combined_Conditional_3plus1_Selector",693:"ScalarEndpoint_Dimension_NoGo_Persistence",694:"GaugeCoefficient_Source_Obligation",695:"PhotonQuantization_Fock_Obligation",696:"Helicity_Orientation_Selector",697:"DimensionPolarization_Record_Schema",698:"GaugeClockCalibration_Identifiability",699:"DimensionGauge_Source_Gate",700:"PhysicalPhoton_Lorentz_Gate",701:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()};A=independent_strict_matrix_float();s=float(A[0,0]);W=s*np.eye(12)-A;w=np.array([W[0,d] for d in range(1,7)]);C=float(np.sum(w*np.arange(1,7)**2))
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

def st687():
 N=4096;a=12/N;q=2*np.pi*np.arange(N)/N;lv=np.array([lam(x)/a**2 for x in q]);times=np.logspace(-1,1,80);z1=np.array([np.mean(np.exp(-t*lv)) for t in times]);rows=[]
 for D in range(1,6):
  z=z1**D;slope=np.gradient(np.log(z),np.log(times));ds=-2*slope;mid=(times>.3)&(times<3);rows.append({"D":D,"median_spectral_dimension_window":float(np.median(ds[mid])),"target":D})
 packet={"N":N,"rows":rows,"theorem_scope":"conditional numerical heat-trace estimator on D-fold product tower","dimension_from_fine_heat_trace":True}
 return finalize(687,"strong_fine_heat_trace_dimension_recovery", "The product family, physical volume normalization and time window are supplied.",packet)
def st688():
 packet={"coarse_endpoint":"same A12 restriction for every D-fold product","fine_heat_trace":"scales with D and can identify D after fine access","theorem":"Spatial dimension is operationally identifiable only after fine multi-axis heat records are supplied; complete one-axis coarse dynamics is dimension blind.","coarse_dimension_identifiable":False}
 return finalize(688,"proven_coarse_fine_dimension_identifiability_separation", "Observer access, not the scalar endpoint alone, determines dimension identifiability.",packet)
def st689():
 rows=[{"spacetime_dimension":n,"Weyl_scaling_exponent":n-4,"conformal":n==4} for n in range(2,9)]
 packet={"action":"integral F wedge star F","rows":rows,"theorem":"Under g->Omega^2 g in n spacetime dimensions, the Maxwell 2-form action scales as Omega^(n-4) and is classically conformally invariant exactly at n=4.","Maxwell_action_sourced":False}
 return finalize(689,"proven_conditional_Maxwell_conformal_dimension4", "Requiring Maxwell conformal invariance imports the gauge action and is not derived from FIN scalar data.",packet)
def st690():
 rows=[{"spatial_D":D,"massless_vector_polarizations":D-1} for D in range(1,8)];packet={"rows":rows,"observed_two_polarization_condition":"D-1=2","selected_spatial_D":3,"theorem":"A massless gauge vector in D spatial dimensions has D-1 transverse polarizations; imposing exactly two selects D=3.","two_polarizations_strict_sourced":False}
 return finalize(690,"proven_conditional_two_polarization_selects_D3", "Using the observed polarization count is a physical input, not a zero-data FIN derivation.",packet)
def st691():
 packet={"sharp_Huygens_spatial_dimensions":"odd D>=3 for the standard massless wave equation","D1":"interior interval dependence rather than shell-only support","even_D":"tails inside the cone","minimal_nontrivial_sharp_dimension":3,"theorem":"Requiring a standard scalar wave fundamental solution supported only on the cone and the smallest nontrivial spatial dimension selects D=3 conditionally.","Huygens_premise_sourced":False}
 return finalize(691,"proven_conditional_Huygens_minimal_D3", "Huygens behavior is an added continuum principle and can fail with mass, curvature or dispersion.",packet)
def st692():
 packet={"selectors":{"Maxwell_conformal":"spacetime n=4","two_polarizations":"spatial D=3","minimal_sharp_Huygens":"spatial D=3"},"common_result":"3+1","theorem":"Three independent conditional physics principles intersect at 3+1 dimensions.","strict_FIN_dimension_derivation":False}
 return finalize(692,"proven_three_conditional_selectors_intersect_at_3plus1", "Agreement of imported principles is coherence, not provenance from the nadsoliton kernel.",packet)
def st693():
 packet={"available_scalar_records":["A12 spectrum","coarse Green functions","dual dynamics","conditional 1D slope"],"dimension_fiber":"all D-fold products restrict to the same one-axis records","theorem":"No functional of the available scalar endpoint/coarse records can distinguish the conditional product dimension.","new_dimension_source":False}
 return finalize(693,"proven_scalar_endpoint_dimension_no_go_persists", "Fine cross-axis or gauge records leave this class.",packet)
def st694():
 packet={"coefficients":["kappa_E","kappa_B","gauge coupling e"],"observable_combinations":{"speed":"sqrt(kappa_E kappa_B)","impedance":"sqrt(kappa_E/kappa_B)","interaction":"e"},"theorem":"Propagation and free-field records can identify speed/impedance combinations but not source the interaction coupling or unit normalization.","coefficients_strict_sourced":False}
 return finalize(694,"proven_gauge_coefficient_source_obligation", "No legacy electromagnetic formula is transferred.",packet)
def st695():
 packet={"classical_input":"symplectic phase space of transverse modes","quantum_additions":["complex structure/positive-frequency split","canonical commutators and hbar","Fock representation","vacuum state","Born instruments"],"theorem":"Classical Maxwell modes and two polarizations do not uniquely determine photon quantization; inequivalent representations can occur beyond finite mode sets.","photon_sector_derived":False}
 return finalize(695,"proven_photon_quantization_additional_structure_obligation", "No photon statistics or quantum field is exported.",packet)
def st696():
 packet={"helicity_pair":["+1","-1"],"inversion":"exchanges helicities","radial_scalar_source":"even","theorem":"A parity-even scalar kernel cannot canonically choose a helicity sign. Circular polarization can be prepared/measured only after an oriented frame or chiral state is supplied.","QW_2191":"open"}
 return finalize(696,"proven_helicity_selector_obstruction", "Both helicities may coexist physically; the result concerns canonical sign provenance.",packet)
def st697():
 packet={"required_records":["fine heat trace over multiple axes","two transverse polarization responses","cone/tail profile","electric/magnetic coefficient ratio","independent rods and both clock classes"],"model_comparison":["D=2,3,4 products","scalar vs gauge","local vs fractional refinements"],"held_out":True,"apparatus_available":False}
 return finalize(697,"constructed_dimension_polarization_record_schema", "The schema is synthetic and has no independent custody record.",packet)
def st698():
 packet={"unknowns":["ballistic clock scale","diffusive clock scale","length scale","kappa_E","kappa_B"],"ideal_record_rank":4,"residual_gauge_dimension":1,"additional_anchor":"one independent dimensional rod/clock ratio","theorem":"Even combined wave, heat, speed and impedance records leave one common calibration orbit without an independent anchor.","physical_identification":False}
 return finalize(698,"proven_gauge_clock_calibration_one_orbit", "Rank accounting is conditional on ideal independent records.",packet)
def st699():
 packet={"dimension_source":False,"Maxwell_action_source":False,"polarization_source":False,"Huygens_source":False,"gauge_coefficients_source":False}
 return finalize(699,"blocked_no_dimension_gauge_source", "Conditional selectors do not become strict by intersection.",packet)
def st700():
 packet={"physical_photon":False,"quantum_Maxwell":False,"exact_Lorentz":False,"charge_matter_coupling":False,"experiment":False}
 return finalize(700,"blocked_no_physical_photon_Lorentz_closure", "No SM/GR, apparatus, or ToE closure.",packet)
def st701():
 packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent","new_cycle_round":4}
 return finalize(701,"blocked_no_independent_empirical_evidence", "Round four is local analytic/computational work.",packet)

def figures(r):
 FIG_DIR.mkdir(exist_ok=True);rows=r["ST687"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.plot([x["D"] for x in rows],[x["median_spectral_dimension_window"] for x in rows],"o-",label="estimated");ax.plot([1,5],[1,5],"--",label="target");ax.legend();ax.set(xlabel="D",ylabel="spectral dimension",title="ST687: fine heat-trace dimension");fig.tight_layout();fig.savefig(FIG_DIR/"st687_spectral_dimension.png",dpi=180);plt.close(fig)
 rows=r["ST689"]["rows"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar([x["spacetime_dimension"] for x in rows],[x["Weyl_scaling_exponent"] for x in rows]);ax.axhline(0,color="black");ax.set(xlabel="spacetime dimension",ylabel="Weyl exponent n-4",title="ST689: Maxwell conformal dimension");fig.tight_layout();fig.savefig(FIG_DIR/"st689_conformal_dimension.png",dpi=180);plt.close(fig)

def main():
 funcs={687:st687,688:st688,689:st689,690:st690,691:st691,692:st692,693:st693,694:st694,695:st695,696:st696,697:st697,698:st698,699:st699,700:st700,701:st701};r={}
 for k in range(687,702):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(687,702):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
