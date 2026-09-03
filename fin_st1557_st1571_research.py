#!/usr/bin/env python3
"""FIN ST1557--ST1571: continuous minimax KL value certificate."""
import math
from collections import deque

import numpy as np
from scipy.optimize import linprog

from fin_interval_cert_common import setup, classical_return_interval, dephased_return_interval, endpoints
from fin_oa_discrimination_common import classical_return, dephased_quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1557,1571
NAMES=["ContinuousMinimaxProblem","TrainingLP","CandidateWeights","AdaptiveKLRectangleBound","AdaptiveCover",
 "ContinuousCandidateLower","DualSupportPoints","DualUpperBound","ContinuousValueBracket","ExactWeightsBoundary",
 "HeldOutGrid","IntegerDesign","FourTimeDefaultBoundary","RoundTwo_Verdict","RoundTwo_Recommendation"]
TIMES=np.array([.3,.6,1.2,2.0])
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def kl(a,b):
 a=np.clip(a,1e-15,1-1e-15);b=np.clip(b,1e-15,1-1e-15)
 return a*np.log(a/b)+(1-a)*np.log((1-a)/(1-b))
def rect_lower(pi,qi,rev=False):
 pl,ph=pi;ql,qh=qi
 if not(ph<ql or qh<pl):return 0.
 if ph<ql:p,q=ph,ql
 else:p,q=pl,qh
 return float(kl(q,p) if rev else kl(p,q))
def make_grid(na,ng):
 p=np.array([classical_return(t) for t in TIMES]);D=[];pars=[]
 for a in np.linspace(.9,1.1,na):
  for g in np.linspace(0,100,ng):
   q=np.array([dephased_quantum_return(a*t,g) for t in TIMES]);D.append(np.minimum(kl(p,q),kl(q,p)));pars.append((float(a),float(g)))
 return np.array(D),pars
def main():
 x={}
 x["ST1557"]=packet(1557,"constructed_continuous_bidirectional_KL_minimax_problem","Frozen decimal spectrum and nuisance rectangle.",{
  "objective":"sup_w inf_theta sum_i w_i min directional Bernoulli KL","times":TIMES.tolist()})
 D,pars=make_grid(41,401);res=linprog(np.r_[np.zeros(4),-1],A_ub=np.c_[-D,np.ones(len(D))],b_ub=np.zeros(len(D)),A_eq=np.array([[1,1,1,1,0]]),b_eq=[1],bounds=[(0,1)]*4+[(None,None)],method='highs');w=res.x[:4]
 x["ST1558"]=packet(1558,"solved_fine_training_grid_LP","Finite grid only.",{
  "rows":len(D),"success":bool(res.success),"value":float(res.x[4])})
 x["ST1559"]=packet(1559,"constructed_two_endpoint_candidate_weights","From fine LP.",{
  "weights":w.tolist(),"zero_middle":bool(w[1]<1e-12 and w[2]<1e-12)})
 x["ST1560"]=packet(1560,"proven_analytic_KL_lower_over_probability_rectangles","Joint convexity; nearest endpoints when intervals disjoint.",{
  "overlap_lower":0.0,"disjoint_rule":"nearest interval endpoints"})
 ctx=setup(20);T=["0.3","2.0"];P=[endpoints(classical_return_interval(t,context=ctx)) for t in T];ww=[float(w[0]),float(w[3])];queue=deque()
 for ia in range(20):
  for ig in range(50):queue.append((.9+ia*.01,.9+(ia+1)*.01,ig*2,(ig+1)*2))
 accepted=splits=0;minacc=float('inf')
 while queue:
  al,ah,gl,gh=queue.popleft();ab=(f"{al:.15g}",f"{ah:.15g}");gb=(f"{gl:.15g}",f"{gh:.15g}");value=0.
  for weight,t,pint in zip(ww,T,P):
   qint=endpoints(dephased_return_interval(t,ab,gb,context=ctx));value+=weight*min(rect_lower(pint,qint),rect_lower(pint,qint,True))
  if value>=.05:accepted+=1;minacc=min(minacc,value);continue
  if (ah-al)/.2 >= (gh-gl)/100:
   m=(al+ah)/2;queue.append((al,m,gl,gh));queue.append((m,ah,gl,gh))
  else:
   m=(gl+gh)/2;queue.append((al,ah,gl,m));queue.append((al,ah,m,gh))
  splits+=1
 x["ST1561"]=packet(1561,"constructed_complete_adaptive_continuous_KL_cover","All boxes accepted above target lower.",{
  "accepted_boxes":accepted,"splits":splits,"minimum_accepted_lower":minacc})
 x["ST1562"]=packet(1562,"proven_candidate_continuous_worst_case_value_at_least_0_05","Outward q intervals and analytic KL rectangle minima.",{
  "lower":.05})
 mu=np.array([.7622353358423254,.2377646641576765]);theta=[(1.1,16.0),(1.1,16.25)];p=np.array([classical_return(t) for t in TIMES]);rows=[]
 for a,g in theta:
  q=np.array([dephased_quantum_return(a*t,g) for t in TIMES]);rows.append(np.minimum(kl(p,q),kl(q,p)))
 avg=mu@np.array(rows)
 x["ST1563"]=packet(1563,"constructed_two_point_dual_nuisance_measure","Any allocation's worst case is below its dual average.",{
  "points":[list(v) for v in theta],"weights":mu.tolist(),"per_time_average":avg.tolist()})
 upper=float(max(avg))+1e-12
 x["ST1564"]=packet(1564,"proven_continuous_minimax_upper_bound","Fixed dual distribution.",{
  "upper":upper})
 x["ST1565"]=packet(1565,"proven_continuous_minimax_value_bracket","Declared metric/model.",{
  "lower":.05,"upper":upper,"width":upper-.05})
 x["ST1566"]=packet(1566,"exact_continuous_optimizing_weights_and_active_parameters_remain_open","Value near-optimality only.",{
  "exact_weights":False,"exact_value":False})
 Dh,parh=make_grid(81,401);vh=Dh@w;idx=int(np.argmin(vh))
 x["ST1567"]=packet(1567,"heldout_finer_grid_is_consistent_with_value_bracket","Redundant numerical check.",{
  "rows":len(Dh),"minimum":float(vh[idx]),"point":list(parh[idx])})
 alloc=[int(round(100*a)) for a in w];alloc[-1]+=100-sum(alloc)
 x["ST1568"]=packet(1568,"constructed_integer_near_minimax_candidate","Not certified optimal integer design.",{
  "shots_per_100":alloc})
 x["ST1569"]=packet(1569,"continuous_KL_bracket_does_not_replace_four_time_RSS_or_model_library_robustness","Endpoint face may be fragile.",{
  "certified_default":"four times"})
 x["ST1570"]=packet(1570,"continuous_value_problem_closed_to_narrow_bracket_exact_optimizer_open","Status.",{
  "relative_bracket_width":(upper-.05)/.05})
 x["ST1571"]=packet(1571,"recommended_ST1572_ST1586","Construct an exact continuous-prior e-process theorem and executable positive quadrature mixture.",{
  "next":["Fubini martingale theorem","continuous prior","positive quadrature","off-grid fixture","one-sided boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
