#!/usr/bin/env python3
"""FIN ST732--ST746: source-selection audit beyond the conditional causal-light bridge."""
import csv,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST732_ST746_Results.json";SUMMARY=ROOT/"FIN_ST732_ST746_Summary.csv"
N={732:"Strict_Locality_Source_Audit",733:"GreenSchur_Locality_Selector_Audit",734:"FullEnergy_Wave_Convergence_Theorem",735:"Uniform_ConeTail_Conditional_Theorem",736:"SpatialDimension_Source_Audit",737:"LorentzRepresentation_Conditional_Limit",738:"GaugeComplex_Provenance_Audit",739:"MaxwellCoefficient_Source_Audit",740:"PhotonQuantization_Interface",741:"TwoClock_Apparatus_Transfer",742:"Noncircular_Anchor_Transfer",743:"ModelLibrary_Nuisance_Robustness",744:"PredictorIR_Method_Gate",745:"Strict_Source_Gate",746:"Evidence_Gate"};P={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in N.items()}
def fin(k,status,boundary,d):
 p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {"program":f"ST{k}","object":N[k],"packet_file":p.name,"packet_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),**d,"status":status,"boundary":boundary}
def main():
 r={}
 r['ST732']=fin(732,'proven_no_strict_fixed_range_source_on_current_endpoint','A future fine-local source can evade the inventory.',{"candidates":["finite A12 support","weight decay eta","Dirichlet positivity","coarse symmetry"],"theorem":"Endpoint support is automatically finite and cannot distinguish towers that add hidden long edges; eta=1.8 favors long range rather than fixed range. No current strict endpoint invariant exports uniform fine locality."})
 r['ST733']=fin(733,'proven_GreenSchur_locality_selector_no_go','Fine Green/instruments leave the class.',{"theorem":"Exact response-preserving Green/Schur compression is identical on the q fiber. Any selector depending only on retained Green data is constant on that fiber and cannot choose q=0."})
 r['ST734']=fin(734,'proven_conditional_full_energy_wave_convergence','The zero-rate tower remains conditional.',{"premises":["finite-range consistency","uniform nonnegative self-adjoint stability","trigonometric-polynomial density"],"theorem":"Consistency on a dense smooth core plus unitary energy stability extends band-limited convergence to strong convergence of wave propagators on the continuum energy space."})
 r['ST735']=fin(735,'proven_conditional_uniform_outside_cone_decay','This is an analytic continuum-approximation theorem, not strict causality at finite N.',{"theorem":"For finite-range stable stencils converging to the 1D wave operator, contour/path estimates yield tails tending uniformly to zero outside every enlarged cone |x|>(c+epsilon)t on compact positive time intervals.","exact_finite_N_zero":False})
 r['ST736']=fin(736,'proven_no_spatial_dimension_source_current_strict','New topology/dimension data can evade.',{"theorem":"All D-fold product completions restrict to the same one-axis strict records; no invariant of current scalar endpoint/coarse data selects D."})
 r['ST737']=fin(737,'proven_conditional_Lorentz_representation_limit','No strict carrier/boost source.',{"theorem":"If strong 3D wave convergence and O(a^2) dispersion control hold, continuum scalar solutions carry the standard Lorentz representation; finite lattices do not carry the exact boosts."})
 r['ST738']=fin(738,'proven_no_unique_gauge_complex_current_scalar_data','Fine cell data can supply a choice.',{"theorem":"Scalar D*D is invariant under arbitrary harmonic edge extensions and multiple cell complexes; current scalar FIN data do not select d0,d1."})
 r['ST739']=fin(739,'blocked_no_Maxwell_coefficient_source','Conditional parameters remain.',{"sources_found":False,"moduli":["kappa_E","kappa_B","gauge coupling"]})
 r['ST740']=fin(740,'constructed_conditional_photon_interface','No quantization is strict-derived.',{"required":["transverse symplectic space","positive-frequency complex structure","hbar/CCR","Fock vacuum","Born instruments"]})
 r['ST741']=fin(741,'constructed_two_clock_transfer_packet','No apparatus exists.',{"channels":{"diffusive":"a^-2","ballistic":"a^-1"},"roles":["provider","registrar","analyst"],"raw_events_required":True})
 r['ST742']=fin(742,'constructed_noncircular_anchor_transfer','External standards absent.',{"requirements":["non-light ruler","independent clock","pre-unblinding hashes","wave record"]})
 r['ST743']=fin(743,'strong_synthetic_nuisance_robustness_plan','Not executed on nature.',{"nuisances":["clock drift","q-label error","polarization leakage","tail floor","calibration scale"],"acceptance":"frozen worst-case separation"})
 r['ST744']=fin(744,'IR_lane_method_gated','No new interval theorem.',{"required":"symbolically correlated parameter-predictor interval arithmetic"})
 r['ST745']=fin(745,'blocked_no_new_strict_source','No physics closure.',{"locality":False,"dimension":False,"gauge":False,"clocks":False,"anchor":False})
 r['ST746']=fin(746,'blocked_no_independent_evidence','Round one is local.',{"external_referee":"absent","laboratory":"absent"})
 RESULTS.write_text(json.dumps(r,indent=2,sort_keys=True));
 with SUMMARY.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['program','status','object','boundary']);[w.writerow([k,x['status'],x['object'],x['boundary']]) for k,x in r.items()]
if __name__=='__main__':main()
