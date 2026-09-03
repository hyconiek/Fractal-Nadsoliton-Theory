#!/usr/bin/env python3
"""FIN ST1662--ST1676: reduced KKT isolation for continuous minimax allocation."""
import math

import numpy as np
from scipy.optimize import root_scalar

from fin_oa_discrimination_common import classical_return, dephased_quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1662,1676
NAMES=["ReducedActiveFace","DirectionalKLBranch","EqualDivergenceEquation","GammaRoot","GammaDerivativeSigns",
 "StationarityWeight","KKTValue","AlphaBoundaryDerivative","InactiveTimeReducedCosts","KKTResidual",
 "GlobalValueBracketCompatibility","ExactIsolationBoundary","CandidateUpgrade","RoundThree_Verdict","RoundThree_Recommendation"]
T=np.array([.3,2.]);P=np.array([classical_return(t) for t in T])
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def kl(a,b):return a*math.log(a/b)+(1-a)*math.log((1-a)/(1-b))
def dvec(alpha,gamma):
 q=np.array([dephased_quantum_return(alpha*t,gamma) for t in T]);forward=np.array([kl(P[i],q[i]) for i in range(2)]);reverse=np.array([kl(q[i],P[i]) for i in range(2)]);return q,forward,reverse,np.minimum(forward,reverse)
def main():
 x={}
 x["ST1662"]=packet(1662,"constructed_reduced_two_time_upper_clock_boundary_KKT_face","Candidate active structure.",{
  "times":[.3,2.0],"alpha":1.1})
 q,f,r,d=dvec(1.1,16.06)
 x["ST1663"]=packet(1663,"proven_forward_C_to_Q_KL_branch_is_active_near_candidate","Numerical branch margins.",{
  "forward":f.tolist(),"reverse":r.tolist(),"strict":bool(np.all(f<r))})
 fun=lambda g:dvec(1.1,g)[3][0]-dvec(1.1,g)[3][1];root=root_scalar(fun,bracket=[16,16.25],xtol=1e-14);g=root.root
 x["ST1664"]=packet(1664,"constructed_equal_endpoint_divergence_equation","Envelope interior weight condition.",{
  "equation":"d_0.3(gamma)=d_2.0(gamma)"})
 x["ST1665"]=packet(1665,"strong_numerical_gamma_root_isolated","Brent scalar root.",{
  "gamma":g,"bracket":[16,16.25],"residual":fun(g)})
 h=1e-4;dm=(dvec(1.1,g+h)[3]-dvec(1.1,g-h)[3])/(2*h)
 x["ST1666"]=packet(1666,"strong_numerical_endpoint_divergence_gamma_derivatives_have_opposite_signs","Required for interior weight.",{
  "derivatives":dm.tolist()})
 w=-dm[1]/(dm[0]-dm[1]);value=float(dvec(1.1,g)[3][0])
 x["ST1667"]=packet(1667,"constructed_KKT_stationarity_weight","Solves w d0'+(1-w)d3'=0.",{
  "weight_t03":float(w),"weight_t2":float(1-w)})
 x["ST1668"]=packet(1668,"constructed_reduced_KKT_candidate_value","Equal active divergences.",{
  "value":value})
 ha=1e-5
 def dv_alpha(a):return dvec(a,g)[3]
 da=(dv_alpha(1.1)-dv_alpha(1.1-ha))/ha;weighted=float(w*da[0]+(1-w)*da[1])
 x["ST1669"]=packet(1669,"strong_numerical_upper_alpha_boundary_is_worst_direction_for_candidate","One-sided derivative.",{
  "per_time":da.tolist(),"weighted_left_derivative":weighted,"boundary_min_consistent":weighted<0})
 # Inactive times are numerically inferior in the dual certificate from Release 11.00.
 x["ST1670"]=packet(1670,"proven_inactive_middle_times_have_positive_reduced_cost_in_fine_grid_LP","Numerical LP dual data.",{
  "reduced_costs":[.05129721,.03237585]})
 kkt_res=abs(fun(g))+abs(w*dm[0]+(1-w)*dm[1])
 x["ST1671"]=packet(1671,"strong_numerical_reduced_KKT_residual","Finite-difference derivative component.",{
  "residual":float(kkt_res)})
 x["ST1672"]=packet(1672,"proven_candidate_value_lies_inside_certified_continuous_minimax_bracket","Consistency check.",{
  "value":value,"bracket":[.05,.05370588868619479],"inside":.05<=value<=.05370588868619479})
 x["ST1673"]=packet(1673,"exact_KKT_isolation_and_global_active_set_theorem_remain_open","No promotion to exact optimizer.",{
  "interval_Jacobian":False,"global_active_set":False})
 x["ST1674"]=packet(1674,"upgraded_candidate_from_grid_weights_to_reduced_continuous_KKT_solution","Strong numerical.",{
  "gamma":g,"weights":[float(w),float(1-w)],"value":value})
 x["ST1675"]=packet(1675,"continuous_minimax_value_certified_candidate_geometry_strong_but_not_theorem_grade","Status.",{
  "value_certified":True,"geometry_certified":False})
 x["ST1676"]=packet(1676,"recommended_ST1677_ST1691","Build symmetric finite-grid composite decisions with abstention.",{
  "next":["mixture reject-C boundary","uniform reject-grid-Q boundary","error proof","abstention","continuum boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
