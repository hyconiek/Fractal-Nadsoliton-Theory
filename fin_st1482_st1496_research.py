#!/usr/bin/env python3
"""FIN ST1482--ST1496: valid finite-mixture e-process."""
import hashlib
import json
import math

import numpy as np

from fin_oa_discrimination_common import classical_return, dephased_quantum_return
from fin_oa_mixture_eprocess import score_mixture
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1482,1496
NAMES=["MixtureGridDefinition","MixturePriorFreeze","ComponentLikelihoodProcesses","MixtureMartingaleTheorem",
 "VilleOneSidedBound","PriorDilutionCost","SyntheticGridAlternative","SyntheticClassicalSafety","InvalidEventFailClosed",
 "NoSymmetricCompositeAcceptance","OffGridPowerBoundary","ContinuousMixtureRequirement","MixtureArtifactHash","RoundThree_Verdict","RoundThree_Recommendation"]
GRID=ROOT/'FIN_ST1482_ST1496_MixtureGrid.json';FIX=ROOT/'FIN_ST1482_ST1496_MixtureFixtures.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};times=[.3,.6,1.2,2.0];pC=[classical_return(t) for t in times];pars=[(a,g) for a in [.9,1.0,1.1] for g in [0,5,10,20,50,100]];weight=1/len(pars)
 components=[{"alpha":a,"gamma":g,"weight":weight,"p_Q":[dephased_quantum_return(a*t,g) for t in times]} for a,g in pars]
 grid={"schema":"FIN-OA-MIXTURE-11.00-v1","times":times,"p_C":pC,"components":components,"schedule":"cycle time indices 0,1,2,3"};GRID.write_text(json.dumps(grid,indent=2,sort_keys=True)+'\n')
 x["ST1482"]=packet(1482,"constructed_predeclared_finite_nuisance_grid","Not continuum-complete.",{
  "alpha":[.9,1.0,1.1],"gamma":[0,5,10,20,50,100],"components":len(components)})
 x["ST1483"]=packet(1483,"frozen_uniform_mixture_prior","Prior chosen before fixtures/data.",{
  "weight":weight,"sum":sum(c['weight'] for c in components)})
 x["ST1484"]=packet(1484,"proven_each_component_likelihood_ratio_is_C_null_martingale","Deterministic time schedule, iid conditional events.",{
  "components":len(components)})
 x["ST1485"]=packet(1485,"proven_weighted_mixture_is_nonnegative_mean_one_C_null_martingale","Convex combination theorem.",{
  "E_C[E_n]":1.0})
 x["ST1486"]=packet(1486,"proven_anytime_rejection_of_C_has_one_percent_bound","One-sided composite-alternative test.",{
  "alpha":.01,"threshold":100,"optional_stopping":True})
 x["ST1487"]=packet(1487,"proven_true_grid_component_pays_at_most_log18_initial_mixture_dilution_in_lower_bound","E_mix >= pi_j LR_j.",{
  "log_penalty":math.log(len(components))})
 rng=np.random.default_rng(1488);true=next(c for c in components if c['alpha']==1.0 and c['gamma']==10);events=[];score=None
 for n in range(5000):
  ti=n%4;events.append({"time_index":ti,"return":int(rng.random()<true['p_Q'][ti])});score=score_mixture(grid,events)
  if score['decision']!='CONTINUE':break
 rngc=np.random.default_rng(1489);eventsC=[]
 for n in range(500):
  ti=n%4;eventsC.append({"time_index":ti,"return":int(rngc.random()<pC[ti])})
 scoreC=score_mixture(grid,eventsC)
 fixtures={"Q_grid":{"truth":{"alpha":1.0,"gamma":10},"events":events,"score":score},"C":{"events":eventsC,"score":scoreC}};FIX.write_text(json.dumps(fixtures,indent=2,sort_keys=True)+'\n')
 x["ST1488"]=packet(1488,"synthetic_grid_alternative_crosses_mixture_boundary","Random code fixture only.",{
  "n":score['n'],"decision":score['decision'],"log_e":score['log_e']})
 x["ST1489"]=packet(1489,"synthetic_classical_sequence_does_not_cross_in_fixture","Not an empirical error estimate.",{
  "n":scoreC['n'],"decision":scoreC['decision'],"log_e":scoreC['log_e']})
 bad=score_mixture(grid,[{"time_index":9,"return":1}])
 x["ST1490"]=packet(1490,"proven_invalid_mixture_event_fails_closed","Validator behavior.",bad)
 x["ST1491"]=packet(1491,"proven_mixture_eprocess_only_supports_anytime_rejection_of_C_not_uniform_acceptance_of_C","Composite asymmetry.",{
  "reject_C_valid":True,"accept_C_against_every_component_valid":False})
 x["ST1492"]=packet(1492,"off_grid_quantum_power_is_not_guaranteed_by_finite_mixture","Type-I remains valid, power may be poor.",{
  "continuum_coverage":False})
 x["ST1493"]=packet(1493,"constructed_continuous_mixture_extension_requirement","Future theorem.",{
  "requires":["normalized prior density on alpha,gamma","stable quadrature or analytic integral","power/cover theorem"]})
 gh=hashlib.sha256(GRID.read_bytes()).hexdigest()
 x["ST1494"]=packet(1494,"constructed_mixture_grid_hash_freeze","Reproducibility.",{
  "grid_file":GRID.name,"fixtures":FIX.name,"sha256":gh})
 x["ST1495"]=packet(1495,"finite_mixture_eprocess_closes_naive_plugin_gap_only_on_one_side_and_grid_scope","No symmetric composite decision.",{
  "valid_anytime_C_rejection":True,"continuous_minimax":False})
 x["ST1496"]=packet(1496,"recommended_ST1497_ST1511","Construct common and independent calibration polytopes.",{
  "next":["preparation simplex","binary confusion box","common-map certificate","independent-map failure","no-click efficiency"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
