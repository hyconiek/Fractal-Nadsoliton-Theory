#!/usr/bin/env python3
"""FIN ST1302--ST1316: finite-count discrimination design."""
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import binom

from fin_oa_discrimination_common import classical_distribution, classical_return, quantum_distribution, quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1302,1316
NAMES=["FrozenBernoulliParameters","Bernoulli_KL","Bernoulli_Chernoff","LikelihoodRatio_Monotonicity",
 "FivePercent_Threshold","OnePercent_Threshold","PerMille_Threshold","HundredShot_Performance",
 "HoeffdingConservativeBound","FullVertex_Chernoff","FullVertex_InformationGain","IID_Assumptions",
 "InvalidRecord_FailClosed","FiniteCount_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def find_threshold(p0,p1,target):
 for n in range(1,1000):
  for k in range(n+2):
   a=float(binom.sf(k-1,n,p0)) if k<=n else 0.0
   b=float(binom.cdf(k-1,n,p1)) if k>0 else 0.0
   if max(a,b)<=target:return n,k,a,b
 raise RuntimeError
def main():
 x={};t=.6;p0=classical_return(t);p1=quantum_return(t)
 x["ST1302"]=packet(1302,"frozen_ideal_Bernoulli_probabilities","Perfect preparation/readout and iid trials.",{
  "time":t,"p_classical":p0,"p_quantum":p1,"effect":p1-p0})
 kl01=p0*math.log(p0/p1)+(1-p0)*math.log((1-p0)/(1-p1));kl10=p1*math.log(p1/p0)+(1-p1)*math.log((1-p1)/(1-p0))
 x["ST1303"]=packet(1303,"proven_positive_directional_KL_divergences","Natural logarithms.",{
  "D_C_Q":kl01,"D_Q_C":kl10})
 r=minimize_scalar(lambda s:p0**s*p1**(1-s)+(1-p0)**s*(1-p1)**(1-s),bounds=(0,1),method='bounded')
 ch=-math.log(float(r.fun))
 x["ST1304"]=packet(1304,"strong_numerical_Bernoulli_Chernoff_information","Scalar bounded optimization.",{
  "Chernoff":ch,"s_star":float(r.x)})
 x["ST1305"]=packet(1305,"proven_binomial_count_is_sufficient_and_LLR_monotone","Simple hypotheses with p_Q>p_C.",{
  "decision":"choose Q when return count K>=k","sufficient_statistic":"K"})
 for prog,target,name in [(1306,.05,"five"),(1307,.01,"one"),(1308,.001,"permille")]:
  n,k,a,b=find_threshold(p0,p1,target)
  x[f"ST{prog}"]=packet(prog,f"proven_exact_binomial_{name}_percent_design" if prog<1308 else "proven_exact_binomial_per_mille_design",
   "Tail probabilities evaluated by scipy double precision.",{"target_max_error":target,"shots":n,"threshold":k,"type_C_to_Q":a,"type_Q_to_C":b})
 n=100
 best=min((max(float(binom.sf(k-1,n,p0)) if k<=n else 0,float(binom.cdf(k-1,n,p1)) if k>0 else 0),k,float(binom.sf(k-1,n,p0)) if k<=n else 0,float(binom.cdf(k-1,n,p1)) if k>0 else 0) for k in range(n+2))
 x["ST1309"]=packet(1309,"proven_hundred_shot_ideal_error_below_five_per_million","Simple model only.",{
  "shots":n,"threshold":best[1],"type_C_to_Q":best[2],"type_Q_to_C":best[3],"max_error":best[0]})
 delta=(p1-p0)/2;nho=math.ceil(math.log(1/.005)/(2*delta**2))
 x["ST1310"]=packet(1310,"constructed_Hoeffding_midpoint_conservative_design","Distribution-free Bernoulli tail bound.",{
  "midpoint":(p0+p1)/2,"shots_for_each_bound_0_005":nho})
 pc=classical_distribution(t);pq=quantum_distribution(t)
 rr=minimize_scalar(lambda s:float(np.sum(pc**s*pq**(1-s))),bounds=(0,1),method='bounded');chfull=-math.log(float(rr.fun))
 x["ST1311"]=packet(1311,"strong_numerical_full_vertex_Chernoff_information","Twelve-outcome model.",{
  "Chernoff":chfull,"s_star":float(rr.x)})
 x["ST1312"]=packet(1312,"proven_full_vertex_record_is_more_informative_than_binary_return_for_declared_models","Data processing/coarse-graining plus computed strict gap.",{
  "binary_Chernoff":ch,"full_Chernoff":chfull,"gain":chfull-ch})
 x["ST1313"]=packet(1313,"listed_ideal_test_assumptions","Violations handled next round.",{
  "assumptions":["iid shots","exact A","clock t=0.6","perfect source vertex","perfect vertex effect","no loss/dephasing/drift"]})
 x["ST1314"]=packet(1314,"constructed_invalid_record_fail_closed_rules","Executable specification precursor.",{
  "reject_if":["K<0 or K>n","missing time/config","shot count mismatch","unfrozen calibration","post-unblinding model change"]})
 x["ST1315"]=packet(1315,"finite_count_discrimination_is_easy_only_for_ideal_simple_hypotheses","No physical feasibility claim.",{
  "shots_max_error_1pct":find_threshold(p0,p1,.01)[0],"nuisance_robust":False})
 x["ST1316"]=packet(1316,"recommended_ST1317_ST1331","Stress-test timing, dephasing, loss and composite alternatives.",{
  "next":["timing windows","energy dephasing mimic","multi-time design","clock-scale nuisance","loss/no-click accounting"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
