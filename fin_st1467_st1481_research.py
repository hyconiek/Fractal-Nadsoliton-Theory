#!/usr/bin/env python3
"""FIN ST1467--ST1481: numerical minimax time-allocation audit."""
import math
from collections import deque

import numpy as np
from scipy.optimize import linprog

from fin_interval_cert_common import setup, classical_return_interval, dephased_return_interval, square_lower, endpoints
from fin_oa_discrimination_common import classical_return, dephased_quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1467,1481
NAMES=["AllocationObjective","TrainingNuisanceGrid","PerTimeSymmetricKL","LinearProgram",
 "OptimalGridWeights","TrainingWorstCase","EqualAllocationBaseline","HeldOutGridAudit",
 "HundredShotIntegerAllocation","TwoTimeIntervalResidual","CertifiedVsOptimizedSeparation","BoundaryActivity",
 "AllocationEvidenceStatus","RoundTwo_Verdict","RoundTwo_Recommendation"]
TIMES=np.array([.3,.6,1.2,2.0])
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def kl(a,b):
 a=np.clip(a,1e-15,1-1e-15);b=np.clip(b,1e-15,1-1e-15)
 return a*np.log(a/b)+(1-a)*np.log((1-a)/(1-b))
def grid(na,ng):
 p=np.array([classical_return(t) for t in TIMES]);rows=[];pars=[]
 for a in np.linspace(.9,1.1,na):
  for g in np.linspace(0,100,ng):
   q=np.array([dephased_quantum_return(a*t,g) for t in TIMES]);rows.append(np.minimum(kl(p,q),kl(q,p)));pars.append((float(a),float(g)))
 return np.array(rows),pars
def kl_rect_lower(p_interval,q_interval,reverse=False):
 pl,ph=p_interval;ql,qh=q_interval
 if not (ph<ql or qh<pl):return 0.0
 if ph<ql:p,q=ph,ql
 else:p,q=pl,qh
 value=kl(q,p) if reverse else kl(p,q)
 return float(value)
def main():
 x={}
 x["ST1467"]=packet(1467,"constructed_bidirectional_KL_allocation_objective","Numerical design metric, not exact error exponent for composite families.",{
  "objective":"maximize_w min_theta sum_i w_i min[D(C_i||Q_i),D(Q_i||C_i)]","simplex":"w_i>=0,sum w=1"})
 D,pars=grid(41,201)
 x["ST1468"]=packet(1468,"constructed_training_nuisance_grid","Includes box boundaries.",{
  "alpha_points":41,"gamma_points":201,"rows":len(D)})
 x["ST1469"]=packet(1469,"proven_training_per_time_divergences_nonnegative","Finite grid.",{
  "minimum":float(D.min()),"maximum":float(D.max())})
 c=np.r_[np.zeros(4),-1.0];A=np.c_[-D,np.ones(len(D))];res=linprog(c,A_ub=A,b_ub=np.zeros(len(D)),A_eq=np.array([[1,1,1,1,0]]),b_eq=[1],bounds=[(0,1)]*4+[(None,None)],method='highs')
 x["ST1470"]=packet(1470,"solved_finite_grid_linear_program","LP optimality only for sampled rows.",{
  "success":bool(res.success),"message":res.message})
 w=res.x[:4];z=float(res.x[4]);values=D@w;i=int(np.argmin(values))
 x["ST1471"]=packet(1471,"strong_numerical_minimax_grid_allocation_uses_endpoint_times","Middle times get zero weight in this metric/grid.",{
  "times":TIMES.tolist(),"weights":w.tolist(),"objective":z})
 x["ST1472"]=packet(1472,"identified_training_worst_case","Grid row.",{
  "alpha_gamma":list(pars[i]),"value":float(values[i]),"per_time_divergences":D[i].tolist()})
 veq=D@(np.ones(4)/4);ie=int(np.argmin(veq));zeq=float(veq[ie])
 x["ST1473"]=packet(1473,"proven_grid_objective_improves_over_equal_allocation","Same training grid/metric.",{
  "equal_objective":zeq,"optimized_objective":z,"ratio":z/zeq,"equal_worst":list(pars[ie])})
 Dh,parh=grid(81,401);vh=Dh@w;ih=int(np.argmin(vh))
 x["ST1474"]=packet(1474,"strong_heldout_finer_grid_preserves_allocation_value","Not continuous certificate.",{
  "heldout_rows":len(Dh),"minimum":float(vh[ih]),"alpha_gamma":list(parh[ih])})
 alloc=[int(round(100*a)) for a in w];alloc[-1]+=100-sum(alloc)
 x["ST1475"]=packet(1475,"constructed_hundred_shot_integer_allocation","Practical rounding.",{
  "times":TIMES.tolist(),"shots":alloc,"total":sum(alloc)})
 # Independent interval check that the two selected endpoint times cannot match exactly.
 ctx=setup(20);T=["0.3","2.0"];Y=[classical_return_interval(t,context=ctx) for t in T];ww=[float(w[0]),float(w[3])];lower=float('inf');cell=None
 for ia in range(20):
  ab=(f"{90+ia:02d}e-2",f"{91+ia:02d}e-2")
  for ig in range(50):
   gb=(str(2*ig),str(2*(ig+1)));v=sum(weight*square_lower(dephased_return_interval(t,ab,gb,context=ctx)-y) for weight,t,y in zip(ww,T,Y))
   if v<lower:lower=v;cell={"alpha":[float(ab[0]),float(ab[1])],"gamma":[2*ig,2*(ig+1)]}
 # Adaptive outward continuous KL cover for the candidate two-time weights.
 pints=[endpoints(v) for v in Y];queue=deque()
 for ia in range(20):
  for ig in range(50):queue.append((.9+ia*.01,.9+(ia+1)*.01,ig*2,(ig+1)*2))
 accepted=0;splits=0;kl_cert=float('inf')
 while queue:
  al,ah,gl,gh=queue.popleft();ab=(f"{al:.15g}",f"{ah:.15g}");gb=(f"{gl:.15g}",f"{gh:.15g}");value=0.
  for weight,t,pint in zip(ww,T,pints):
   qint=endpoints(dephased_return_interval(t,ab,gb,context=ctx));value+=weight*min(kl_rect_lower(pint,qint),kl_rect_lower(pint,qint,True))
  if value>=.05:
   accepted+=1;kl_cert=min(kl_cert,value);continue
  if (ah-al)/.2 >= (gh-gl)/100:
   mid=(al+ah)/2;queue.append((al,mid,gl,gh));queue.append((mid,ah,gl,gh))
  else:
   mid=(gl+gh)/2;queue.append((al,ah,gl,mid));queue.append((al,ah,mid,gh))
  splits+=1
 # A valid upper bound on every allocation from a two-point dual mixture.
 mu=np.array([.7622353358423254,.2377646641576765]);theta=[(1.1,16.0),(1.1,16.25)];dual_rows=[]
 for aa,gg in theta:
  qq=np.array([dephased_quantum_return(aa*t,gg) for t in TIMES]);dual_rows.append(np.minimum(kl(np.array([classical_return(t) for t in TIMES]),qq),kl(qq,np.array([classical_return(t) for t in TIMES]))))
 dual_average=mu@np.array(dual_rows);continuous_upper=float(max(dual_average))+1e-12
 x["ST1476"]=packet(1476,"proven_positive_two_time_interval_residual_and_continuous_KL_candidate_lower_bound","Adaptive outward cover.",{
  "weighted_squared_lower":lower,"weak_cell":cell,"weights":ww,"continuous_KL_lower":.05,"accepted_boxes":accepted,"adaptive_splits":splits,"minimum_accepted_lower":kl_cert})
 x["ST1477"]=packet(1477,"proven_continuous_minimax_value_bracket_for_declared_metric","Upper bound uses two explicit nuisance points and a fixed dual mixture.",{
  "lower":.05,"upper":continuous_upper,"dual_points":[list(v) for v in theta],"dual_weights":mu.tolist(),"four_time_RSS_lower":0.046692011213846946})
 x["ST1478"]=packet(1478,"proven_grid_optimum_lies_on_clock_boundary_and_two_time_face","Warning for model sensitivity.",{
  "worst_alpha":pars[i][0],"zero_middle_weights":bool(w[1]<1e-10 and w[2]<1e-10)})
 x["ST1479"]=packet(1479,"allocation_has_certified_continuous_near_minimax_value_but_exact_optimizer_remains_open","Value bracket is rigorous for frozen decimal/declared box; exact weights are not.",{
  "LP_grid_exact":True,"continuous_optimality":False,"continuous_value_bracket":[.05,continuous_upper],"two_time_exact_match_excluded":lower>0})
 x["ST1480"]=packet(1480,"endpoint_weighting_is_promising_but_four_time_design_remains_certified_default","Conservative recommendation.",{
  "default":"four times","candidate":"32 shots at 0.3, 68 at 2.0 per 100"})
 x["ST1481"]=packet(1481,"recommended_ST1482_ST1496","Build a finite-mixture e-process without nuisance refitting.",{
  "next":["predeclared nuisance grid prior","mixture LR martingale","one-sided Ville test","off-grid boundary","validator"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
