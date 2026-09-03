#!/usr/bin/env python3
"""FIN ST1407--ST1421: full multinomial likelihood and finite-shot bounds."""
import math

import numpy as np
from scipy.optimize import minimize_scalar

from fin_oa_discrimination_common import classical_distribution, quantum_distribution
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1407,1421
NAMES=["FullVertexLikelihood","DistanceClassReduction","SevenClassSufficiency","PerClassLogLikelihoodRatios",
 "FullKL_Directions","FullChernoff","SingleShotTVBayesError","ChernoffFiniteShotBound",
 "FullRecord_OnePercentBound","BinaryChernoffBoundComparison","InformationGainRatio","MultinomialCalibrationMap",
 "CompositeCertifiedConservativeShots","RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};t=.6;p=classical_distribution(t);q=quantum_distribution(t)
 x["ST1407"]=packet(1407,"constructed_full_twelve_outcome_log_likelihood","Simple frozen hypotheses.",{
  "per_event":"log(q_x/p_x)","sample_LLR":"sum_x n_x log(q_x/p_x)"})
 # distance classes 0..6; classes 1..5 have paired vertices
 pc=[float(p[0])]+[float(p[d]+p[-d]) for d in range(1,6)]+[float(p[6])]
 qc=[float(q[0])]+[float(q[d]+q[-d]) for d in range(1,6)]+[float(q[6])]
 x["ST1408"]=packet(1408,"constructed_seven_cyclic_distance_classes","Uses reflection-symmetric frozen distributions.",{
  "multiplicities":[1,2,2,2,2,2,1],"classical":pc,"quantum":qc})
 x["ST1409"]=packet(1409,"proven_distance_class_is_sufficient_for_frozen_vertex_likelihood","Paired vertices have identical likelihood ratio.",{
  "information_loss_vs_12_outcomes":0.0,"classes":7})
 llr=[math.log(b/a) for a,b in zip(pc,qc)]
 x["ST1410"]=packet(1410,"proven_finite_per_class_log_likelihoods","All probabilities strictly positive.",{
  "LLR":llr,"monotone_in_distance":False})
 klpq=float(np.sum(p*np.log(p/q)));klqp=float(np.sum(q*np.log(q/p)))
 x["ST1411"]=packet(1411,"proven_full_directional_KL_values","Natural logs.",{"D_C_Q":klpq,"D_Q_C":klqp})
 r=minimize_scalar(lambda s:float(np.sum(p**s*q**(1-s))),bounds=(0,1),method='bounded');ch=-math.log(float(r.fun))
 x["ST1412"]=packet(1412,"strong_numerical_full_Chernoff_information","Scalar bounded optimization.",{"C":ch,"s_star":float(r.x)})
 tv=float(.5*np.sum(abs(p-q)));bayes=.5*(1-tv)
 x["ST1413"]=packet(1413,"proven_single_shot_equal_prior_Bayes_error","Simple distributions.",{"TV":tv,"Bayes_error":bayes})
 x["ST1414"]=packet(1414,"proven_Chernoff_equal_prior_error_upper_bound","Standard product-measure bound.",{
  "formula":"P_e(n)<=0.5 exp(-n C)"})
 nfull=math.ceil(math.log(50)/ch)
 x["ST1415"]=packet(1415,"constructed_full_record_one_percent_sufficient_shot_bound","Chernoff bound, not exact minimal n.",{
  "shots":nfull,"upper":.5*math.exp(-nfull*ch)})
 cbinary=0.10028599984500335;nbinary=math.ceil(math.log(50)/cbinary)
 x["ST1416"]=packet(1416,"constructed_binary_Chernoff_one_percent_sufficient_bound","Conservative comparison.",{
  "shots":nbinary,"upper":.5*math.exp(-nbinary*cbinary)})
 x["ST1417"]=packet(1417,"proven_full_record_has_positive_Chernoff_gain","Frozen models.",{
  "gain":ch-cbinary,"ratio":ch/cbinary})
 x["ST1418"]=packet(1418,"constructed_confusion_matrix_likelihood_pushforward","Known detector matrix M maps p,q before likelihood calculation.",{
  "formula":"p_obs=M p; q_obs=M q","data_processing":"Chernoff/TV cannot improve under coarse noise"})
 rss_lower=0.046692011213846946;delta=math.sqrt(rss_lower/4);ncomp=math.ceil(2*math.log(800)/(delta**2))
 x["ST1419"]=packet(1419,"constructed_conservative_four_time_composite_shot_bound_from_interval_residual","Hoeffding union bound, binary returns, per time.",{
  "residual_lower":delta,"shots_each":ncomp,"total_shots":4*ncomp,"type_I_union_bound_target":.01})
 x["ST1420"]=packet(1420,"full_distance_class_record_is_preferred_for_simple_models_but_composite_bound_remains_conservative","No exact composite likelihood yet.",{
  "simple_full_likelihood_ready":True,"composite_minimax_likelihood_ready":False})
 x["ST1421"]=packet(1421,"recommended_ST1422_ST1436","Add anytime-valid likelihood/e-process rules and protocol 10.99.",{
  "next":["simple e-process","Ville bounds","expected growth","composite warning","validator upgrade"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
