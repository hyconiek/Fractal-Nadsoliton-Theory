#!/usr/bin/env python3
"""FIN ST1887--ST1901: source-action gate update and final synthesis."""
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np

from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1887,1901
NAMES=["GaugeAttemptLedger","PartialGateScore","GaugeMatterPositiveResult","UniquenessNoGo",
 "CurvatureTopologyNoGo","RefinementContinuumNoGo","FermionGravityNoGo","UpdatedDecisiveGate",
 "RouteRanking","HighestValueNextTheorem","FalsificationOrder","StrictConditionalBoundary",
 "FinalInterpretation","StopRules","CycleSixReportTrigger"]
FIGDIR=ROOT/'FIN_ST1812_ST1901_Figures';FIG=FIGDIR/'gauge_source_action_audit.png';GATE=ROOT/'FIN_ST1887_ST1901_GaugeSourceActionGate.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def make_figure():
 FIGDIR.mkdir(exist_ok=True);fig,ax=plt.subplots(1,3,figsize=(11,3.7))
 ax[0].bar(['G1 regular','G2 scale','G3 gauge','G4 EOM','G5 OA','G6 predict','G7 evidence'],[.5,0,.4,.3,0,0,0]);ax[0].set(title='Partial source-action gate readiness',ylabel='qualitative readiness',ylim=(0,1));ax[0].tick_params(axis='x',rotation=45);ax[0].grid(axis='y',alpha=.25)
 r=np.arange(1,7);edges=np.where(r<6,12*r,66);ax[1].plot(r,edges,'o-',label='edges');ax[1].plot(r,edges-11,'s--',label='cycle rank');ax[1].set(title='Band-threshold topology',xlabel='max cyclic distance');ax[1].legend(fontsize=7);ax[1].grid(alpha=.25)
 dims=[12,78];zeros=[1,56];ax[2].bar(['vertex A','incidence Dirac'],dims,label='Hilbert dim');ax[2].bar(['vertex A','incidence Dirac'],zeros,label='zero modes');ax[2].set(title='First-order enlargement',ylabel='dimension/count');ax[2].legend(fontsize=7);ax[2].grid(axis='y',alpha=.25);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def main():
 make_figure();x={}
 x["ST1887"]=packet(1887,"constructed_direct_gauge_source_action_attempt_ledger","Immediate construction and falsification.",{
  "positive":["A=d*d","globally regular U1 matter action","literal gauge-covariant EOM","exact coarse refinement"],"negative":["group/charge/source nonunique","no canonical faces","no continuum","spectral triple incomplete"]})
 x["ST1888"]=packet(1888,"scored_decisive_gate_as_partial_progress_but_zero_complete_strict_rows","Conditional subrows do not count as strict pass.",{
  "complete_strict_rows":0,"conditional_partial_rows":["G1","G3","G4"]})
 x["ST1889"]=packet(1889,"proven_minimal_conditional_gauge_matter_action_repairs_prior_seagull_and_variation_defects","Once local U1 complex scalar and links are chosen.",{
  "globally_regular":True,"gauge_invariant":True,"single_action_EOM":True})
 x["ST1890"]=packet(1890,"proven_A_does_not_select_gauge_group_charge_complexification_or_connection_state","No unique physics.",{})
 x["ST1891"]=packet(1891,"proven_complete_support_does_not_select_unique_2complex_face_measure_or_Maxwell_action","K12 clique vs band choices.",{})
 x["ST1892"]=packet(1892,"proven_exact_refinement_transport_does_not_select_q_locality_dimension_or_continuum","Moduli grow.",{})
 x["ST1893"]=packet(1893,"proven_incidence_Dirac_does_not_automatically_supply_chiral_fermion_SM_or_gravity_data","56 zero modes and representation gaps.",{})
 gate={"schema":"FIN-GAUGE-SOURCE-ACTION-GATE-v2","rows":[
  {"id":"H1","name":"strict source of state category/gauge group/charge","strict":False},
  {"id":"H2","name":"strict source of connection/holonomy state","strict":False},
  {"id":"H3","name":"canonical 2-complex, face measure and gauge kinetic action","strict":False},
  {"id":"H4","name":"refinement scaling to local continuum dimension/causality","strict":False},
  {"id":"H5","name":"real even spectral triple/chiral fermion representation","strict":False},
  {"id":"H6","name":"dimensional action and nonabsorbable prediction","strict":False},
  {"id":"H7","name":"OA platform and independent record","strict":False}]}
 GATE.write_text(json.dumps(gate,indent=2,sort_keys=True)+'\n')
 x["ST1894"]=packet(1894,"constructed_updated_gauge_source_action_gate","All seven strict rows open.",{
  "file":GATE.name,"rows":7,"strict_passes":0})
 ranking=[
  {"route":"derive canonical 2-complex/face measure from FIN refinement","rigorous_success":.28,"physics_decision":1.0},
  {"route":"publish conditional U1 Dirichlet lattice matter model","rigorous_success":.8,"physics_decision":.45},
  {"route":"complete real even spectral triple from incidence complex","rigorous_success":.3,"physics_decision":.75},
  {"route":"continue APD/moment Ltotal repair","rigorous_success":.01,"physics_decision":.5},
  {"route":"OA laboratory test of already frozen channels","rigorous_success":.45,"physics_decision":.55}]
 x["ST1895"]=packet(1895,"constructed_post_attempt_route_ranking","Strategic estimates only.",{"ranking":ranking})
 x["ST1896"]=packet(1896,"selected_highest_value_next_theorem_as_canonical_cell_complex_face_measure_or_no_go_from_refinement_data","This decides whether gauge field dynamics can be strict-sourced.",{
  "target":"functor (A,refinement)->oriented weighted 2-complex with natural Wilson/spectral gauge action"})
 x["ST1897"]=packet(1897,"constructed_immediate_falsification_order_for_cell_complex_candidate","Stop early.",{
  "tests":["automorphism naturality","threshold independence","refinement functoriality","locality/dimension","positive face weights","continuum Maxwell scaling"]})
 x["ST1898"]=packet(1898,"strict_core_remains_finite_Dirichlet_geometry_conditional_gauge_model_must_be_labelled","No new physics by assumption.",{})
 x["ST1899"]=packet(1899,"final_interpretation_direct_gauge_completion_is_best_constructive_bridge_but_not_unique_or_physical_without_cell_complex_and_sources","Survives all six rounds.",{
  "fundamental_physics":False,"conditional_gauge_matter":True})
 x["ST1900"]=packet(1900,"installed_release1104_stop_rules","Prevent replay/promotion.",{
  "stop_if":["choice of U1 called derivation","K12 clique called spacetime","thresholded C12 called canonical","incidence grading called SM chirality","refinement transport called continuum","conditional action called strict physics"]})
 x["ST1901"]=packet(1901,"six_round_cycle_complete_report_required","Release 11.04.",{
  "programs":90,"rounds":6,"release":"11.04","main_result":"minimal gauge matter action exists conditionally; unique curvature/continuum/chiral physics no-go; canonical 2-complex selected as decisive next theorem","strict_ToE_closure":False,"gate_sha256":hashlib.sha256(GATE.read_bytes()).hexdigest(),"figure":str(FIG.relative_to(ROOT))})
 write_round(LO,HI,x)
if __name__=="__main__":main()
