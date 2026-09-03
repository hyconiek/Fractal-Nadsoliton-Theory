#!/usr/bin/env python3
"""FIN ST1527--ST1541: Release 11.00 synthesis and robustness ladder."""
import hashlib
import json
import math

import matplotlib.pyplot as plt
import numpy as np

from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1527,1541
NAMES=["RobustnessLadder","StructuredUncertaintyVerdict","AllocationVerdict","MixtureEProcessVerdict",
 "CalibrationVerdict","RefinementVerdict","Protocol1100Freeze","EvidenceBoundary","AnnihilationScope",
 "StrictCoreStatus","RemainingBlockers","RecommendedNextProgrammes","StopRules","FinalScientificVerdict","CycleSix_ReportTrigger"]
FIGDIR=ROOT/'FIN_ST1452_ST1541_Figures';FIG=FIGDIR/'robust_refinement_summary.png';PROTO=ROOT/'FIN_ST1527_ST1541_Protocol_11_00.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def make_figure():
 FIGDIR.mkdir(exist_ok=True);fig,ax=plt.subplots(1,3,figsize=(11,3.7))
 ax[0].bar(['0.3','0.6','1.2','2.0'],[.315539497,0,0,.684460503]);ax[0].set(title='Grid-minimax allocation',xlabel='time',ylabel='weight',ylim=(0,.75));ax[0].grid(axis='y',alpha=.25)
 factors=np.linspace(0,1,101);ax[1].plot(factors,.41121822546536*factors);ax[1].axvline(.912,color='crimson',ls='--',label='common calibration lower');ax[1].set(title='Gap under common affine factor',xlabel='factor',ylabel='return gap');ax[1].legend(fontsize=7);ax[1].grid(alpha=.25)
 q=np.linspace(0,2,301);t=.6;h=.5*(1+np.exp(-2*q*t));u=np.cos(q*t)**2;ax[2].plot(q,h,label='fiber heat');ax[2].plot(q,u,label='fiber unitary');ax[2].set(title='Fine fiber dual dynamics',xlabel='$q$',ylabel='fiber return');ax[2].legend(fontsize=7);ax[2].grid(alpha=.25);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def h(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def main():
 make_figure();x={}
 x["ST1527"]=packet(1527,"constructed_four_level_certificate_robustness_ladder","Every level declares extra premises.",{
  "levels":["point decimal spectrum","structured circulant rounding box","common calibration polytope","exact coarse refinement"]})
 x["ST1528"]=packet(1528,"proven_key_certificates_survive_structured_entry_radius_5e_minus15","Not arbitrary matrices.",{
  "root_common_bracket":[.5945307,.5945309],"RSS_lower":0.04669201121378897})
 x["ST1529"]=packet(1529,"allocation_candidate_improves_grid_metric_but_does_not_replace_certified_four_time_default","Numerical/interval split.",{
  "weights":[.3155394971590464,0,0,.6844605028409537],"heldout_metric":.0537030927793231,"continuous_KL_optimality":False})
 x["ST1530"]=packet(1530,"finite_mixture_eprocess_provides_valid_anytime_C_rejection_only","Grid-alternative scope.",{
  "components":18,"alpha":.01,"symmetric_composite_acceptance":False})
 x["ST1531"]=packet(1531,"common_calibration_polytope_preserves_certified_separation","Independent maps need tighter bounds.",{
  "common_RSS_lower":.03883580017500168,"common_one_time_residual_lower":.09853400450479226,"independent_5pct_2pct_passes":False,"independent_3pct_2pct_passes":True})
 x["ST1532"]=packet(1532,"coarse_protocol_transports_exactly_through_two_fiber_refinement","Fine test optional.",{
  "coarse_q_blind":True,"fine_heat":"(1+exp(-2qt))/2","fine_unitary":"cos^2(qt)","q_sourced":False})
 protocol={"schema":"FIN-OA-11.00-v1","parents":{"10.99":h('FIN_ST1422_ST1436_Protocol_10_99.json'),"mixture":h('FIN_ST1482_ST1496_MixtureGrid.json'),"calibration":h('FIN_ST1497_ST1511_CalibrationPolytope.json')},"structured_uncertainty":{"entry_radius":5e-15,"circulant":True,"symmetric":True,"row_sum_zero":True},"certificates":{"root_bracket":[.5945307,.5945309],"composite_RSS_lower":.04669201121378897,"common_calibration_RSS_lower":.03883580017500168},"allocation":{"certified_default_times":[.3,.6,1.2,2.0],"numerical_candidate_weights":[.3155394971590464,0,0,.6844605028409537]},"refinement":{"A24":"A12 tensor I2 + I12 tensor Bq","coarse_transport":True,"q":"unsourced"},"roles":{"independent_custody_required":True,"nature_events_present":False}}
 PROTO.write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n');ph=h(PROTO.name)
 x["ST1533"]=packet(1533,"constructed_protocol_11_00_hash_freeze","Composes parent artifacts without altering their claims.",{
  "file":PROTO.name,"sha256":ph,"figure":str(FIG.relative_to(ROOT))})
 x["ST1534"]=packet(1534,"blocked_no_symbolic_kernel_box_platform_calibration_or_nature_events","Transfer boundary.",{
  "symbolic_kernel":False,"platform":False,"external_calibration":False,"raw_events":False})
 x["ST1535"]=packet(1535,"protocol_still_tests_internal_channel_models_not_literal_total_annihilation","Ontology guardrail.",{
  "whole_nadsoliton_nonexistence_tested":False})
 x["ST1536"]=packet(1536,"strict_FIN_physics_status_unchanged","Conditional OA mathematics strengthened only.",{
  "state_category_strict":False,"selector":False,"units":False,"laboratory":False})
 x["ST1537"]=packet(1537,"listed_remaining_mathematical_blockers","Next proof frontier.",{
  "blockers":["noncirculant interval matrix","symbolic kernel parameter box","continuous minimax KL","continuous-mixture eprocess",
              "full detector polytope","platform forward map"]})
 x["ST1538"]=packet(1538,"recommended_next_programmes","Nonreplay continuation.",{
  "programmes":["interval-eigenvector bounds for noncirculant A boxes","symbolic omega-phi-beta-eta uncertainty propagation",
   "continuous KL minimax certificate","continuous prior mixture e-process","independent calibration polytope optimization",
   "vertex-dependent detector matrix certificate","optimal shots with integer allocation","fine-fiber q-identifiability",
   "multi-level refinement consistency","one platform-specific transfer function","independent protocol review","raw data only after custody freeze"]})
 x["ST1539"]=packet(1539,"installed_release1100_stop_rules","No silent scope extension.",{
  "stop_if":["noncirculant error treated as circulant","grid optimum called continuum theorem","mixture Q rejection called symmetric decision",
             "common calibration used for independent devices","coarse transport called q derivation","synthetic fixture called evidence"]})
 x["ST1540"]=packet(1540,"final_verdict_conditional_OA_discriminator_is_now_robust_across_four_explicit_layers","Physical realization remains missing.",{
  "theorem_grade":["structured uncertainty","common calibration","coarse refinement"],"numerical":["KL allocation"],"operational":["finite mixture eprocess"]})
 x["ST1541"]=packet(1541,"six_round_cycle_complete_report_required","Release 11.00.",{
  "programs":90,"rounds":6,"release":"11.00","breakthrough":"structured uncertainty and calibration/refinement robustness ladder","strict_ToE_closure":False})
 write_round(LO,HI,x)
if __name__=="__main__":main()
