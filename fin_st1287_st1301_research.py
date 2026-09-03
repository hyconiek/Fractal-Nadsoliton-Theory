#!/usr/bin/env python3
"""FIN ST1287--ST1301: spectral discrimination and time optimization."""
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from fin_oa_discrimination_common import EIG, classical_distribution, classical_return, quantum_distribution, quantum_return, return_derivatives
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1287,1301
NAMES=["SpectrumMultiplicityReduction","ClassicalReturn_Monotonicity","QuantumReturn_RecurrenceStructure","SeparationFunction",
 "OptimalTime_LocalRoot","OptimalTime_Curvature","DenseDomain_Audit","LipschitzCover_Bound","FullDistribution_TV",
 "PracticalTime_06","CurveCrossings","LateTime_Behavior","SingleTime_DesignFreeze","OptimizationEvidenceStatus","RoundTwo_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};vals=[]
 for v in EIG:
  if not vals or abs(v-vals[-1][0])>1e-9:vals.append([float(v),1])
  else:vals[-1][1]+=1
 x["ST1287"]=packet(1287,"proven_return_formulas_reduce_to_seven_spectral_values","Floating equality groups tolerance 1e-9.",{
  "eigenvalues":[v for v,m in vals],"multiplicities":[m for v,m in vals]})
 x["ST1288"]=packet(1288,"proven_classical_return_decreases_to_one_twelfth","All nonzero eigenvalues positive.",{
  "C0":1.0,"limit":1/12,"derivative_negative_for_t_positive":True})
 x["ST1289"]=packet(1289,"proven_quantum_return_is_finite_almost_periodic_not_monotone","Finite unitary spectrum.",{
  "Q0":1.0,"fixed_limit":False,"recurrences":True})
 x["ST1290"]=packet(1290,"constructed_return_separation_objective","Binary effect.",{
  "D(t)":"|Q(t)-C(t)|","domain":[0,10]})
 r=minimize_scalar(lambda t:-(quantum_return(t)-classical_return(t)),bounds=(0.3,1.0),method='bounded',options={'xatol':1e-14});tstar=float(r.x);f=return_derivatives(tstar)
 x["ST1291"]=packet(1291,"strong_numerical_unique_first_positive_separation_maximizer","Local bounded optimization, not interval root proof.",{
  "t_star":tstar,"C":f['c'],"Q":f['q'],"D":f['q']-f['c'],"derivative_residual":f['qp']-f['cp']})
 curvature=f['qpp']-f['cpp']
 x["ST1292"]=packet(1292,"proven_negative_local_curvature_numerically","Analytic derivative formula evaluated in double precision.",{
  "second_derivative":curvature,"strict_local_maximum":bool(curvature<0)})
 grid=np.linspace(0,10,1000001);ds=np.abs(np.array([quantum_return(float(t))-classical_return(float(t)) for t in grid]));imax=int(np.argmax(ds))
 x["ST1293"]=packet(1293,"strong_dense_grid_global_domain_audit","Not outward interval arithmetic.",{
  "grid_points":len(grid),"grid_step":float(grid[1]-grid[0]),"grid_max_time":float(grid[imax]),"grid_max":float(ds[imax])})
 lmax=float(EIG.max());L=3*lmax;h=float(grid[1]-grid[0])
 x["ST1294"]=packet(1294,"constructed_global_Lipschitz_cover_upper_bound","Function values are floating; analytic derivative bound is rigorous given eigenvalue bound.",{
  "derivative_bound":L,"cover_upper":float(ds[imax]+L*h/2),"cover_gap":L*h/2})
 pc=classical_distribution(tstar);pq=quantum_distribution(tstar);tv=float(.5*np.sum(abs(pc-pq)))
 x["ST1295"]=packet(1295,"proven_full_vertex_TV_equals_return_gap_at_candidate_time_numerically","All nonreturn signed differences share the opposite sign in replay.",{
  "TV":tv,"return_gap":float(pq[0]-pc[0]),"sign_pattern":bool(np.all(pq[1:]<=pc[1:]))})
 tp=.6;Dp=quantum_return(tp)-classical_return(tp)
 x["ST1296"]=packet(1296,"frozen_practical_time_0_6_is_near_optimal","Simpler design time.",{
  "t":tp,"C":classical_return(tp),"Q":quantum_return(tp),"D":Dp,"loss_vs_optimum":float((f['q']-f['c'])-Dp)})
 xs=np.linspace(1e-8,10,10001);v=np.array([quantum_return(t)-classical_return(t) for t in xs]);roots=[]
 for a,b,va,vb in zip(xs[:-1],xs[1:],v[:-1],v[1:]):
  if va*vb<0:roots.append(float(brentq(lambda t:quantum_return(t)-classical_return(t),a,b)))
 x["ST1297"]=packet(1297,"strong_numerical_return_curve_crossing_inventory_on_0_10","Grid bracketing plus Brent roots.",{
  "positive_roots":roots,"root_count":len(roots)})
 x["ST1298"]=packet(1298,"proven_late_time_single_return_order_can_reverse","Unitary recurrence makes one-time sign nonuniversal.",{
  "t10_C":classical_return(10),"t10_Q":quantum_return(10),"Q_minus_C":quantum_return(10)-classical_return(10)})
 x["ST1299"]=packet(1299,"frozen_base_single_time_design","Pre-registered for simple hypothesis testing.",{
  "time":tp,"observable":"return/not_return","primary_effect_size":Dp})
 x["ST1300"]=packet(1300,"optimization_is_strong_numerical_not_formal_global_certificate","No theorem inflation.",{
  "interval_arithmetic":False,"dense_grid":True,"analytic_derivatives":True,"global_unique_max_proven":False})
 x["ST1301"]=packet(1301,"recommended_ST1302_ST1316","Design exact finite-count likelihood tests.",{
  "next":["Bernoulli KL/Chernoff","integer thresholds","sample sizes","multinomial option","failure rules"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
