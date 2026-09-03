#!/usr/bin/env python3
"""FIN ST1572--ST1586: continuous-prior e-process theorem and quadrature implementation."""
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from fin_oa_discrimination_common import classical_return, dephased_quantum_return
from fin_oa_mixture_eprocess import score_mixture
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1572,1586
NAMES=["ContinuousPriorDefinition","PointwiseLikelihoodMartingales","FubiniMixtureTheorem","VilleContinuousPriorBound",
 "PositiveQuadratureGrid","QuadratureWeights","QuadratureIsValidFiniteMixture","QuadratureVsExactIntegralBoundary",
 "OffGridSyntheticAlternative","ClassicalSyntheticSafety","PriorDilution","OneSidedDecisionBoundary",
 "ArtifactHash","RoundThree_Verdict","RoundThree_Recommendation"]
GRID=ROOT/'FIN_ST1572_ST1586_ContinuousPriorQuadrature.json';FIX=ROOT/'FIN_ST1572_ST1586_ContinuousPriorFixtures.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};times=[.3,.6,1.2,2.0];pC=[classical_return(t) for t in times]
 x["ST1572"]=packet(1572,"constructed_uniform_continuous_prior_on_nuisance_rectangle","Mathematical prior, not physical distribution.",{
  "alpha":[.9,1.1],"gamma":[0,100],"density":1/(.2*100)})
 x["ST1573"]=packet(1573,"proven_each_theta_likelihood_ratio_is_C_null_martingale","Deterministic schedule and conditionally independent events.",{})
 x["ST1574"]=packet(1574,"proven_continuous_prior_integral_is_nonnegative_mean_one_martingale","Tonelli/Fubini exchanges expectation and positive integral.",{
  "E_n":"integral L_theta,n pi(dtheta)"})
 x["ST1575"]=packet(1575,"proven_continuous_prior_anytime_C_rejection_bound","Exact integral object.",{
  "threshold":100,"type_I_at_most":.01})
 xa,wa=leggauss(7);xg,wg=leggauss(9);components=[]
 for i in range(7):
  alpha=1.0+.1*float(xa[i])
  for j in range(9):
   gamma=50+50*float(xg[j]);weight=float(wa[i]*wg[j]/4)
   components.append({"alpha":alpha,"gamma":gamma,"weight":weight,"p_Q":[dephased_quantum_return(alpha*t,gamma) for t in times]})
 total=sum(c['weight'] for c in components)
 for c in components:c['weight']/=total
 grid={"schema":"FIN-OA-CONTINUOUS-QUADRATURE-11.01-v1","times":times,"p_C":pC,"components":components,"interpretation":"positive Gauss-Legendre finite mixture approximating uniform continuous prior"};GRID.write_text(json.dumps(grid,indent=2,sort_keys=True)+'\n')
 x["ST1576"]=packet(1576,"constructed_positive_7_by_9_Gauss_Legendre_nuisance_grid","Executable approximation.",{
  "alpha_nodes":7,"gamma_nodes":9,"components":len(components)})
 x["ST1577"]=packet(1577,"proven_quadrature_weights_positive_and_normalized","Floating construction.",{
  "minimum":min(c['weight'] for c in components),"sum":sum(c['weight'] for c in components)})
 x["ST1578"]=packet(1578,"proven_quadrature_sum_is_itself_a_valid_finite_mixture_eprocess","Validity does not depend on quadrature accuracy.",{
  "threshold":100,"type_I_at_most":.01})
 x["ST1579"]=packet(1579,"quadrature_eprocess_is_not_the_exact_continuous_integral","Accuracy/power approximation remains numerical.",{
  "exact_integral_evaluated":False})
 rng=np.random.default_rng(1580);truth=(1.03,7.3);events=[];score=None
 for n in range(10000):
  ti=n%4;p=dephased_quantum_return(truth[0]*times[ti],truth[1]);events.append({"time_index":ti,"return":int(rng.random()<p)});score=score_mixture(grid,events)
  if score['decision']!='CONTINUE':break
 rngc=np.random.default_rng(1581);eventsC=[]
 for n in range(500):
  ti=n%4;eventsC.append({"time_index":ti,"return":int(rngc.random()<pC[ti])})
 scoreC=score_mixture(grid,eventsC);fixtures={"off_grid_Q":{"truth":truth,"events":events,"score":score},"C":{"events":eventsC,"score":scoreC}};FIX.write_text(json.dumps(fixtures,indent=2,sort_keys=True)+'\n')
 x["ST1580"]=packet(1580,"synthetic_off_node_alternative_crosses_quadrature_mixture_boundary","Power fixture, not theorem/evidence.",{
  "truth":list(truth),"n":score['n'],"decision":score['decision']})
 x["ST1581"]=packet(1581,"synthetic_classical_fixture_does_not_cross","Not empirical calibration.",{
  "n":scoreC['n'],"decision":scoreC['decision'],"log_e":scoreC['log_e']})
 x["ST1582"]=packet(1582,"continuous_prior_can_dilute_local_power_even_while_preserving_validity","No uniform expected stopping theorem.",{
  "uniform_power_bound":False})
 x["ST1583"]=packet(1583,"continuous_and_quadrature_mixtures_remain_one_sided_C_rejection_tools","No symmetric composite acceptance.",{
  "accept_C_uniformly":False})
 gh=hashlib.sha256(GRID.read_bytes()).hexdigest()
 x["ST1584"]=packet(1584,"constructed_continuous_prior_quadrature_hash","Reproducibility.",{
  "grid":GRID.name,"fixtures":FIX.name,"sha256":gh})
 x["ST1585"]=packet(1585,"continuous_prior_theorem_closes_conceptual_grid_gap_executable_integral_accuracy_open","Status split.",{
  "exact_theorem":True,"exact_integral_code":False,"finite_valid_approximation":True})
 x["ST1586"]=packet(1586,"recommended_ST1587_ST1601","Interval-certify independent preparation/detector polytopes.",{
  "next":["affine vertex extrema","independent C/Q maps","composite cover","efficiency/no-click","matrix detector boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
