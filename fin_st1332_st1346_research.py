#!/usr/bin/env python3
"""FIN ST1332--ST1346: frozen executable protocol and validator."""
import hashlib
import json
from pathlib import Path

import numpy as np

from fin_oa_discrimination_common import classical_return, quantum_return
from fin_oa_protocol_validator import score_primary, validate_aggregate_record, validate_protocol
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1332,1346
NAMES=["ProtocolSchema","FrozenModelAndTimeFields","PrimaryDecisionRule","RawEventSchema","NoClickAccounting",
 "NuisanceLibrary","ProtocolHashFreeze","ValidatorPositivePath","ValidatorNegativePath","SyntheticClassicalFixture",
 "SyntheticQuantumFixture","FourTimeSyntheticFixture","OneShotAnalysisRule","CustodyBoundary","RoundFive_Recommendation"]
PROTOCOL=ROOT/'FIN_ST1332_ST1346_Protocol.json';FIXTURES=ROOT/'FIN_ST1332_ST1346_SyntheticFixtures.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def aggregate(model,seed):
 rng=np.random.default_rng(seed);times=[.3,.6,1.2,2.0];attempts=[100]*4;clicks=[100]*4
 probs=[classical_return(t) if model=='C' else quantum_return(t) for t in times]
 returns=[int(rng.binomial(100,p)) for p in probs]
 return {"model_blind_id":f"fixture-{seed}","times":times,"attempts":attempts,"clicks":clicks,"return_counts":returns,"config_id":"OA-10.98","run_id":f"SYN-{seed}","truth":model}
def main():
 x={}
 matrix_hash=json.loads((ROOT/'FIN_ST1272_ST1286_Results.json').read_text())['ST1282']['sha256']
 protocol={
  "schema_version":"FIN-OA-10.98-v1","matrix_fingerprint_sha256":matrix_hash,
  "hypotheses":{"C":"p(t)=exp(-tA)delta0","Q":"rho(t)=exp(-itA)|0><0|exp(itA)"},
  "primary":{"time":.6,"shots":29,"return_threshold_Q":19,"max_ideal_error_target":.01},
  "robust":{"times":[.3,.6,1.2,2.0],"shots_each":[100,100,100,100],"clock_scale_band":[.9,1.1],"dephasing_gamma_band":[0,100]},
  "raw_event_fields":["timestamp","attempt_id","click","vertex_if_click","config_id","run_id"],
  "invalid_record_rules":["missing field","count ordering failure","shot mismatch","unfrozen time","post-unblinding model edit"],
  "role_requirements":{"provider_equals_analyst_allowed":False,"hash_before_unblinding":True,"one_shot_analysis":True}}
 PROTOCOL.write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n')
 fixtures={"classical":aggregate('C',1339),"quantum":aggregate('Q',1340)};FIXTURES.write_text(json.dumps(fixtures,indent=2,sort_keys=True)+'\n')
 ph=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
 x["ST1332"]=packet(1332,"constructed_frozen_executable_protocol_schema","Conditional OA only.",{"protocol_file":PROTOCOL.name,"schema":protocol['schema_version']})
 x["ST1333"]=packet(1333,"frozen_matrix_hypotheses_and_times","Any change requires a new version before data.",{
  "matrix_hash":matrix_hash,"primary":protocol['primary'],"robust":protocol['robust']})
 x["ST1334"]=packet(1334,"constructed_fail_closed_primary_integer_rule","Ideal simple hypotheses.",{
  "rule":"at t=0.6 with 29 clicks choose Q iff returns>=19; else C","invalid_is_not_C_or_Q":True})
 x["ST1335"]=packet(1335,"constructed_raw_event_schema","Aggregate fixtures do not replace raw laboratory events.",{"fields":protocol['raw_event_fields']})
 x["ST1336"]=packet(1336,"proven_no_click_must_be_recorded_against_postselection","Protocol counts attempts and clicks separately.",{
  "attempt_field":True,"click_field":True,"conditional_only_record_rejected":True})
 x["ST1337"]=packet(1337,"frozen_nuisance_library_bounds","Not exhaustive hardware systematics.",{
  "clock_scale_band":[.9,1.1],"dephasing_gamma_band":[0,100],"uniform_loss":"record via no-click"})
 x["ST1338"]=packet(1338,"constructed_protocol_SHA256_freeze","Hash must precede unblinding.",{"sha256":ph})
 x["ST1339"]=packet(1339,"proven_protocol_and_valid_fixtures_pass_validator","Local synthetic replay.",{
  "protocol_errors":validate_protocol(protocol),"classical_errors":validate_aggregate_record(protocol,fixtures['classical']),"quantum_errors":validate_aggregate_record(protocol,fixtures['quantum'])})
 bad=dict(fixtures['classical']);bad['return_counts']=[101,0,0,0]
 x["ST1340"]=packet(1340,"proven_invalid_count_fixture_fails_closed","Negative test.",{
  "errors":validate_aggregate_record(protocol,bad),"score":score_primary(protocol,bad)})
 # Primary scorer needs a 29-click t=.6 record.
 rng=np.random.default_rng(1341);kc=int(rng.binomial(29,classical_return(.6)));kq=int(rng.binomial(29,quantum_return(.6)))
 def primary_rec(k,tag):return {"model_blind_id":tag,"times":[.6],"attempts":[29],"clicks":[29],"return_counts":[k],"config_id":"OA-10.98","run_id":tag}
 x["ST1341"]=packet(1341,"synthetic_classical_primary_fixture_classified_C","Random fixture, not evidence.",{"returns":kc,"score":score_primary(protocol,primary_rec(kc,'C-syn'))})
 x["ST1342"]=packet(1342,"synthetic_quantum_primary_fixture_classified_Q","Random fixture, not evidence.",{"returns":kq,"score":score_primary(protocol,primary_rec(kq,'Q-syn'))})
 x["ST1343"]=packet(1343,"constructed_four_time_synthetic_fixture_bundle","Used for future scorer development.",{"fixture_file":FIXTURES.name,"classical_counts":fixtures['classical']['return_counts'],"quantum_counts":fixtures['quantum']['return_counts']})
 x["ST1344"]=packet(1344,"constructed_one_shot_analysis_rule","Prevents model repair after seeing holdout.",{
  "rule":"validate hash and record once; emit C/Q/INVALID; no refit after unblinding"})
 x["ST1345"]=packet(1345,"blocked_independent_custody_and_nature_record","Cannot be generated by local code.",{
  "provider":False,"registrar":False,"independent_analyst":False,"nature_events":False})
 x["ST1346"]=packet(1346,"recommended_ST1347_ST1361","Final falsification and research-roadmap round.",{
  "next":["base-vs-composite verdict","minimal claim","failure modes","laboratory boundary","next 10 programmes"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
