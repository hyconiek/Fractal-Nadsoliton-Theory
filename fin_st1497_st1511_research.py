#!/usr/bin/env python3
"""FIN ST1497--ST1511: calibration polytopes and no-click accounting."""
import hashlib
import json
import math

from fin_oa_discrimination_common import classical_return, quantum_return
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1497,1511
NAMES=["CommonCalibrationPolytope","CommonMapFormula","CommonGapFactor","CommonCompositeRSS",
 "IndependentCalibrationBound","IndependentDefaultBudgetFailure","IndependentSafeBudget","KnownLossEfficiencyBox",
 "AttemptClickAccounting","ExpectedAttemptInflation","VertexDependentConfusionBoundary","CalibrationArtifact",
 "CalibrationStatus","RoundFour_Verdict","RoundFour_Recommendation"]
CAL=ROOT/'FIN_ST1497_ST1511_CalibrationPolytope.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};epsmax=.05;amax=bmax=.02;etamin=.8;base=quantum_return(.6)-classical_return(.6);rss=.04669201121378897
 cal={"schema":"FIN-OA-CAL-11.00-v1","common":{"uniform_preparation_epsilon":[0,epsmax],"false_positive_a":[0,amax],"false_negative_b":[0,bmax],"uniform_click_efficiency_eta":[etamin,1]},"independent_safe_example":{"preparation_epsilon_each":.03,"detector_probability_error_each":.02},"requires":["attempts","click/no-click","vertex if click","calibration run id"]};CAL.write_text(json.dumps(cal,indent=2,sort_keys=True)+'\n')
 x["ST1497"]=packet(1497,"constructed_common_calibration_polytope","Common nuisance map applied to both hypotheses.",cal['common'])
 x["ST1498"]=packet(1498,"proven_common_preparation_detector_map_is_affine","Binary return channel.",{
  "formula":"T(p)=a+(1-a-b)[(1-epsilon)p+epsilon/12]"})
 factor=(1-amax-bmax)*(1-epsmax)
 x["ST1499"]=packet(1499,"proven_common_polytope_preserves_positive_base_gap","Worst factor at upper nuisance corner.",{
  "factor_lower":factor,"base_gap":base,"observed_gap_lower":factor*base})
 x["ST1500"]=packet(1500,"proven_common_affine_map_preserves_composite_RSS_certificate_up_to_scale","Same common map at all four times.",{
  "RSS_lower":rss*factor**2,"one_time_residual_lower":math.sqrt(rss*factor**2/4)})
 delta=math.sqrt(rss/4);formula="11 epsilon_p/12 + delta_detector"
 x["ST1501"]=packet(1501,"proven_independent_calibration_error_condition","Each hypothesis may have a different map.",{
  "per_model_error_bound":formula,"required_less_than":delta/2})
 default=.05*11/12+.02
 x["ST1502"]=packet(1502,"proven_default_independent_5pct_plus_2pct_budget_not_certified_by_current_residual_bound","Does not prove actual nonidentifiability.",{
  "per_model_bound":default,"threshold":delta/2,"passes":default<delta/2})
 safe=.03*11/12+.02
 x["ST1503"]=packet(1503,"proven_independent_3pct_plus_2pct_example_passes_bound","Conservative sufficient condition.",{
  "per_model_bound":safe,"threshold":delta/2,"passes":safe<delta/2,"remaining_one_time_gap":delta-2*safe})
 x["ST1504"]=packet(1504,"proven_known_uniform_loss_does_not_change_click_conditioned_model_gap","Attempts remain essential.",{
  "eta_band":[etamin,1],"conditional_gap_factor":1})
 x["ST1505"]=packet(1505,"constructed_attempt_click_vertex_record_chain","Fail-closed ordering.",{
  "constraints":"0<=returns<=clicks<=attempts","no_clicks":"attempts-clicks"})
 x["ST1506"]=packet(1506,"constructed_worst_mean_attempt_inflation_for_target_clicks","Expectation only, not high-probability guarantee.",{
  "eta_min":etamin,"mean_factor":1/etamin,"1146_clicks_mean_attempts":1146/etamin})
 x["ST1507"]=packet(1507,"vertex_dependent_efficiency_or_confusion_requires_matrix_polytope","Scalar certificate does not cover it.",{
  "needed":"stochastic/substochastic detector matrix box","implemented":False})
 ch=hashlib.sha256(CAL.read_bytes()).hexdigest()
 x["ST1508"]=packet(1508,"constructed_calibration_polytope_hash_freeze","Reproducibility.",{
  "file":CAL.name,"sha256":ch})
 x["ST1509"]=packet(1509,"common_calibration_robustness_proven_independent_calibration_only_conditionally_bounded","Scope split.",{
  "common_polytope":True,"independent_safe_example":True,"full_detector_polytope":False})
 x["ST1510"]=packet(1510,"calibration_is_an_OA_premise_not_a_strict_FIN_prediction","No physical calibration record.",{
  "laboratory_calibration":False})
 x["ST1511"]=packet(1511,"recommended_ST1512_ST1526","Transport coarse protocol through exact 12-to-24 refinement.",{
  "next":["Kronecker factorization","coarse return invariance","fine fiber dynamics","q blindness","detector expansion"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
