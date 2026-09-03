#!/usr/bin/env python3
"""FIN ST1617--ST1631: Release 11.01 synthesis."""
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np

from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1617,1631
NAMES=["RobustnessHierarchy","NoncirculantVerdict","ContinuousMinimaxVerdict","ContinuousPriorVerdict",
 "IndependentCalibrationVerdict","FiberQVerdict","Protocol1101Freeze","ClaimLattice","EvidenceBoundary",
 "AnnihilationScope","RemainingBlockers","RecommendedNextProgrammes","StopRules","FinalScientificVerdict","CycleSix_ReportTrigger"]
FIGDIR=ROOT/'FIN_ST1542_ST1631_Figures';FIG=FIGDIR/'continuum_and_fiber_summary.png';PROTO=ROOT/'FIN_ST1617_ST1631_Protocol_11_01.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def h(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def figure():
 FIGDIR.mkdir(exist_ok=True);fig,ax=plt.subplots(1,3,figsize=(11,3.7))
 e=np.logspace(-14,-7,200);m=3.089572412462438e-6-3*(.5945308+10)*e;ax[0].semilogx(e,m);ax[0].axhline(0,color='k',lw=.8);ax[0].set(title='Noncirculant location margin',xlabel='$\\|B-A\\|$',ylabel='certified margin');ax[0].grid(alpha=.25)
 ax[1].bar(['lower','upper'],[.05,.05370588868619479]);ax[1].set(title='Continuous minimax value bracket',ylabel='bidirectional KL value',ylim=(0,.06));ax[1].grid(axis='y',alpha=.25)
 q=np.linspace(.1,2,400);f1=np.cos(.6*q)**2-(1+np.exp(-1.2*q))/2;f2=np.cos(q)**2-(1+np.exp(-2*q))/2;ax[2].plot(q,f1,label='$t=0.6$');ax[2].plot(q,f2,label='$t=1.0$');ax[2].axhline(0,color='k',lw=.8);ax[2].set(title='Fine heat-unitary residuals',xlabel='$q$',ylabel='$R_q-H_q$');ax[2].legend(fontsize=7);ax[2].grid(alpha=.25);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def main():
 figure();x={}
 x["ST1617"]=packet(1617,"constructed_extended_six_layer_robustness_hierarchy","From point model toward operations.",{
  "layers":["decimal spectrum","circulant box","noncirculant graph norm ball","continuous nuisance","independent calibration","exact refinement/fine fiber"]})
 x["ST1618"]=packet(1618,"proven_noncirculant_graph_perturbations_preserve_global_maximum_location_below_norm_threshold","Uniqueness not preserved.",{
  "epsilon_threshold":9.720652006795929e-08,"target":[.59,.60]})
 x["ST1619"]=packet(1619,"proven_continuous_minimax_value_narrow_bracket_exact_weights_open","Frozen metric/model.",{
  "value":[.05,.05370588868619479],"exact_weights":False})
 x["ST1620"]=packet(1620,"continuous_prior_eprocess_theorem_exact_positive_quadrature_executable","Both one-sided.",{
  "continuous_theorem":True,"quadrature_components":63,"exact_integral_code":False})
 x["ST1621"]=packet(1621,"proven_independent_5pct_2pct_affine_calibration_composite_separation","Sharper than triangle bound.",{
  "RSS_lower":.022937735365240054,"one_time_residual_lower":.07572604466965123})
 x["ST1622"]=packet(1622,"proven_fine_q_identifiability_split_by_channel_and_range","Conditional result.",{
  "classical_one_time_global":True,"quantum_one_time_first_branch":True,"two_time_channel_RSS_lower":.009270580748993435,"q_zero_obstruction":True})
 protocol={"schema":"FIN-OA-11.01-v1","parents":{"11.00":h('FIN_ST1527_ST1541_Protocol_11_00.json'),"continuous_quadrature":h('FIN_ST1572_ST1586_ContinuousPriorQuadrature.json'),"independent_calibration":h('FIN_ST1587_ST1601_IndependentCalibrationPolytope.json')},"noncirculant":{"graph_laplacian":True,"operator_norm_threshold":9.720652006795929e-08,"global_location":[.59,.60],"unique_not_claimed":True},"continuous_minimax":{"value_bracket":[.05,.05370588868619479],"candidate_weights":[.31753172,0,0,.68246828]},"fine_fiber":{"q_range":[.1,2],"times":[.6,1.0],"RSS_lower":.009270580748993435,"q_zero_invalid":True},"roles":{"external_review":False,"nature_events":False}}
 PROTO.write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n')
 x["ST1623"]=packet(1623,"constructed_protocol_11_01_hash_freeze","Composes new artifacts.",{
  "file":PROTO.name,"sha256":h(PROTO.name),"figure":str(FIG.relative_to(ROOT))})
 x["ST1624"]=packet(1624,"constructed_proven_strong_numerical_conditional_blocked_claim_lattice","Status discipline.",{
  "proven":["Duhamel thresholds","continuous value bracket","continuous-prior theorem","independent calibration RSS","q interval cover"],"strong_numerical":["candidate weights","quadrature power"],"blocked":["platform","events","strict OA"]})
 x["ST1625"]=packet(1625,"blocked_no_platform_calibration_raw_events_or_independent_custody","Physical boundary.",{
  "platform":False,"calibration_record":False,"events":False,"custody":False})
 x["ST1626"]=packet(1626,"protocol_remains_internal_channel_and_fiber_test_not_total_annihilation_test","Ontology guardrail.",{
  "literal_nonexistence":False})
 x["ST1627"]=packet(1627,"listed_remaining_frontier_after_continuum_upgrades","No replay.",{
  "blockers":["noncirculant composite dephasing sensitivity","symbolic kernel box","exact minimax weights","symmetric composite decision",
              "vertex detector matrix","platform forward law"]})
 x["ST1628"]=packet(1628,"recommended_next_programmes","Next nonreplay batch.",{
  "programmes":["Duhamel bounds for Lindblad/dephasing generators","interval symbolic kernel parameter propagation",
   "KKT/dual isolation of exact minimax weights","uniformly valid composite acceptance or confidence set","detector matrix polytope",
   "fine-fiber q clock-scale joint identifiability","multi-level q consistency","optimal continuous prior for power",
   "platform transfer function","external calibration form","independent audit","raw events only after freeze"]})
 x["ST1629"]=packet(1629,"installed_release1101_stop_rules","Prevent overpromotion.",{
  "stop_if":["location robustness called uniqueness","value bracket called exact optimizer","quadrature called exact integral",
             "affine detector called arbitrary matrix","q interval result used at q=0 or outside range","synthetic power called evidence"]})
 x["ST1630"]=packet(1630,"final_verdict_conditional_OA_now_has_noncirculant_location_continuum_value_and_independent_calibration_certificates","Strict physics unchanged.",{
  "conditional_math_strengthened":True,"strict_state_category":False,"laboratory":False})
 x["ST1631"]=packet(1631,"six_round_cycle_complete_report_required","Release 11.01.",{
  "programs":90,"rounds":6,"release":"11.01","breakthrough":"noncirculant location, continuous minimax bracket, independent calibration and q certificates","strict_ToE_closure":False})
 write_round(LO,HI,x)
if __name__=="__main__":main()
