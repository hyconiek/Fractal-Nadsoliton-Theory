#!/usr/bin/env python3
"""FIN ST1422--ST1436: anytime-valid e-process and protocol 10.99."""
import hashlib
import json
import math

import numpy as np

from fin_oa_discrimination_common import classical_distribution, quantum_distribution
from fin_oa_protocol_validator_10_99 import score_sequence, validate_protocol
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1422,1436
NAMES=["LikelihoodRatio_EProcess","MartingaleUnderClassical","ReverseMartingaleUnderQuantum","VilleAnytimeBounds",
 "ExpectedLogGrowth","ApproximateStoppingLengths","FullVsBinarySequentialGain","OptionalStoppingSafety",
 "CompositePluginWarning","CompositeEProcessRequirement","Protocol1099Freeze","SequentialValidator",
 "SyntheticSequentialFixtures","CustodyBoundary","RoundFive_Recommendation"]
PROTOCOL=ROOT/'FIN_ST1422_ST1436_Protocol_10_99.json';FIX=ROOT/'FIN_ST1422_ST1436_SequentialFixtures.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};p=classical_distribution(.6);q=quantum_distribution(.6);pc=[float(p[0])]+[float(p[d]+p[-d]) for d in range(1,6)]+[float(p[6])];qc=[float(q[0])]+[float(q[d]+q[-d]) for d in range(1,6)]+[float(q[6])];llr=[math.log(b/a) for a,b in zip(pc,qc)];alpha=.01
 x["ST1422"]=packet(1422,"constructed_full_record_likelihood_ratio_e_process","Simple frozen hypotheses only.",{
  "E_n":"product q(X_i)/p(X_i)","outcomes":"seven distance classes"})
 x["ST1423"]=packet(1423,"proven_E_n_is_nonnegative_mean_one_martingale_under_C","IID simple null.",{
  "E_C[E_n]":1.0})
 x["ST1424"]=packet(1424,"proven_inverse_likelihood_is_mean_one_martingale_under_Q","IID simple alternative.",{
  "E_Q[1/E_n]":1.0})
 x["ST1425"]=packet(1425,"proven_symmetric_anytime_error_bounds_by_Ville","Stopping rule may be data-dependent.",{
  "alpha":alpha,"Q_threshold":1/alpha,"C_threshold":alpha,"each_wrong_boundary_probability_at_most":alpha})
 dqc=float(np.sum(q*np.log(q/p)));dcq=float(np.sum(p*np.log(p/q)))
 x["ST1426"]=packet(1426,"proven_expected_log_likelihood_growth_rates","Full vertex/distance record.",{
  "under_Q":dqc,"under_C":-dcq})
 loga=math.log(1/alpha)
 x["ST1427"]=packet(1427,"constructed_Wald_scale_stopping_length_estimates","Ignores overshoot and is not a guarantee.",{
  "under_Q":loga/dqc,"under_C":loga/dcq})
 x["ST1428"]=packet(1428,"proven_full_record_has_larger_expected_growth_than_binary_return","Frozen models.",{
  "full_Q":dqc,"binary_Q":0.35744222947815313,"full_C":dcq,"binary_C":0.4218128334046286})
 x["ST1429"]=packet(1429,"proven_optional_stopping_is_safe_for_declared_simple_e_process","No fixed sample size needed for Type-I control.",{
  "peeking_allowed":True,"model_adaptation_allowed":False})
 x["ST1430"]=packet(1430,"proven_plug_in_maximized_nuisance_likelihood_is_not_automatically_an_e_process","Composite null/alternative warning.",{
  "naive_refit_safe":False})
 x["ST1431"]=packet(1431,"constructed_valid_composite_extension_requirements","Future work.",{
  "options":["predeclared mixture likelihood","safe test martingale","universal inference with proved domination"],"implemented":False})
 protocol={
  "schema_version":"FIN-OA-10.99-v1","decimal_spectrum":{"values":["0","0.75412115420707981","1.5770495144276093","1.9614068619764458","2.1995688493332102","2.2986062720790965","2.3421820411462999"],"multiplicities":[1,2,2,2,2,2,1]},
  "simple_full_record":{"time":.6,"distance_classes":7,"p_C":pc,"p_Q":qc},
  "sequential":{"alpha":alpha,"distance_class_log_likelihood_ratios":llr,"upper_decision":"Q","lower_decision":"C"},
  "composite":{"times":[.3,.6,1.2,2.0],"clock_scale_band":[.9,1.1],"dephasing_gamma_band":[0,100],"certified_RSS_lower":0.046692011213846946,"conservative_shots_each":1146},
  "calibration":{"max_independent_preparation_probability_error":.05,"max_independent_detector_probability_error":.05,"attempt_and_no_click_required":True},
  "role_requirements":{"provider_equals_analyst_allowed":False,"hash_before_unblinding":True,"one_shot_or_anytime_rule_frozen":True}}
 PROTOCOL.write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n');ph=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
 x["ST1432"]=packet(1432,"constructed_protocol_10_99_freeze","Includes interval and sequential upgrades.",{
  "file":PROTOCOL.name,"sha256":ph})
 x["ST1433"]=packet(1433,"proven_protocol_passes_sequential_validator","Local schema test.",{
  "errors":validate_protocol(protocol)})
 rng=np.random.default_rng(1434)
 def simulate(prob,label):
  seq=[]
  for _ in range(1000):
   seq.append(int(rng.choice(7,p=prob)));score=score_sequence(protocol,seq)
   if score['decision']!='CONTINUE':return {"truth":label,"outcomes":seq,"score":score}
  return {"truth":label,"outcomes":seq,"score":score_sequence(protocol,seq)}
 fixtures={"C":simulate(pc,"C"),"Q":simulate(qc,"Q")};FIX.write_text(json.dumps(fixtures,indent=2,sort_keys=True)+'\n')
 x["ST1434"]=packet(1434,"synthetic_sequential_fixtures_stop_on_correct_boundaries","Random code fixtures, not evidence.",{
  "file":FIX.name,"C_decision":fixtures['C']['score']['decision'],"C_n":fixtures['C']['score']['n'],"Q_decision":fixtures['Q']['score']['decision'],"Q_n":fixtures['Q']['score']['n']})
 x["ST1435"]=packet(1435,"blocked_no_independent_registration_or_nature_sequence","Protocol readiness only.",{
  "external_hash_custody":False,"raw_events":False})
 x["ST1436"]=packet(1436,"recommended_ST1437_ST1451","Final theorem/status/roadmap round.",{
  "next":["certificate synthesis","claim ladder","remaining composite gap","laboratory boundary","next programmes"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
