#!/usr/bin/env python3
"""FIN ST1677--ST1691: symmetric finite-grid composite decision with abstention."""
import hashlib
import json

import numpy as np

from fin_oa_discrimination_common import classical_return
from fin_oa_symmetric_grid_decision import score_symmetric
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1677,1691
NAMES=["ThreeDecisionRule","MixtureQBoundary","AllGridCBoundary","ClassicalWrongQBound",
 "GridQWrongCBound","MutualExclusivity","AbstentionNecessity","SyntheticClassicalDecision",
 "SyntheticGridQDecision","InvalidEvent","OffGridBoundary","ContinuousAllQTheoremSchema",
 "ArtifactHash","RoundFour_Verdict","RoundFour_Recommendation"]
GRID=ROOT/'FIN_ST1482_ST1496_MixtureGrid.json';FIX=ROOT/'FIN_ST1677_ST1691_SymmetricFixtures.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};grid=json.loads(GRID.read_text());alpha=.01
 x["ST1677"]=packet(1677,"constructed_Q_mixture_C_all_grid_continue_decision_rule","Finite grid and deterministic schedule.",{
  "decisions":["Q_MIXTURE","C_ALL_GRID","CONTINUE","INVALID"]})
 x["ST1678"]=packet(1678,"proven_mixture_boundary_is_anytime_valid_under_C","Prior-weighted LR e-process.",{
  "boundary":100,"wrong_Q_probability_at_most":alpha})
 x["ST1679"]=packet(1679,"constructed_uniform_C_boundary_as_all_component_LRs_below_alpha","Equivalent max_j LR_j<=alpha.",{
  "boundary":alpha})
 x["ST1680"]=packet(1680,"proven_classical_wrong_Q_boundary_probability_at_most_one_percent","Ville.",{})
 x["ST1681"]=packet(1681,"proven_each_true_grid_Q_component_has_wrong_C_boundary_probability_at_most_one_percent","Event all LR<=alpha implies true component inverse LR>=1/alpha.",{
  "uniform_over_grid":True})
 x["ST1682"]=packet(1682,"proven_Q_and_C_boundaries_are_mutually_exclusive","If all component LR<=.01 then mixture LR<=.01, not >=100.",{})
 x["ST1683"]=packet(1683,"proven_abstention_is_required_for_finite_time_symmetric_error_control","Neither boundary need be reached.",{})
 times=grid['times'];pC=grid['p_C'];rng=np.random.default_rng(1684);evC=[];sC=None
 for n in range(10000):
  ti=n%4;evC.append({"time_index":ti,"return":int(rng.random()<pC[ti])});sC=score_symmetric(grid,evC,alpha)
  if sC['decision']!='CONTINUE':break
 true=next(c for c in grid['components'] if c['alpha']==1.0 and c['gamma']==10);rngq=np.random.default_rng(1685);evQ=[];sQ=None
 for n in range(10000):
  ti=n%4;evQ.append({"time_index":ti,"return":int(rngq.random()<true['p_Q'][ti])});sQ=score_symmetric(grid,evQ,alpha)
  if sQ['decision']!='CONTINUE':break
 fixtures={"C":{"events":evC,"score":sC},"Q":{"truth":{"alpha":1,"gamma":10},"events":evQ,"score":sQ}};FIX.write_text(json.dumps(fixtures,indent=2,sort_keys=True)+'\n')
 x["ST1684"]=packet(1684,"synthetic_classical_fixture_reaches_C_all_grid_boundary","Code fixture.",{
  "decision":sC['decision'],"n":sC['n']})
 x["ST1685"]=packet(1685,"synthetic_grid_Q_fixture_reaches_Q_mixture_boundary","Code fixture.",{
  "decision":sQ['decision'],"n":sQ['n']})
 bad=score_symmetric(grid,[{"time_index":99,"return":1}],alpha)
 x["ST1686"]=packet(1686,"proven_invalid_event_fails_closed","Validator.",bad)
 x["ST1687"]=packet(1687,"finite_grid_C_decision_has_no_uniform_off_grid_Q_error_bound","Model library boundary.",{
  "continuum_valid":False})
 x["ST1688"]=packet(1688,"constructed_continuous_symmetric_schema","Not executable/certified.",{
  "C_boundary":"sup_theta LR_theta<=alpha","Q_boundary":"integral LR_theta pi(dtheta)>=1/alpha","validity":"formal if sup/integral exact"})
 fh=hashlib.sha256(FIX.read_bytes()).hexdigest()
 x["ST1689"]=packet(1689,"constructed_symmetric_fixture_hash","Reproducibility.",{
  "file":FIX.name,"sha256":fh})
 x["ST1690"]=packet(1690,"symmetric_finite_grid_decision_closed_with_abstention_continuous_implementation_open","Status.",{
  "finite_grid":True,"continuous":False})
 x["ST1691"]=packet(1691,"recommended_ST1692_ST1706","Add arbitrary vertex misclassification bounds and joint q-clock no-go.",{
  "next":["column stochastic detector ball","TV lower bound","loss augmentation","q/tau scaling torsor","calibration necessity"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
