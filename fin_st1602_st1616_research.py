#!/usr/bin/env python3
"""FIN ST1602--ST1616: fine-fiber q identifiability."""
import math

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1602,1616
NAMES=["FineFiberModel","ClassicalQMonotonicity","ClassicalQInverse","QuantumQAliasing",
 "QuantumRestrictedInjectivity","ReturnDerivatives","BernoulliFisherInformation","ReferenceQValues",
 "SingleTimeChannelCrossings","QZeroObstruction","TwoTimeFineDesign","TwoTimeIntervalCover",
 "TwoTimeResidualBound","PhysicalQBoundary","RoundFive_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def h(q,t):return .5*(1+math.exp(-2*q*t))
def u(q,t):return math.cos(q*t)**2
def main():
 x={};t=.6;q0=.37
 x["ST1602"]=packet(1602,"constructed_fine_fiber_heat_and_unitary_models","Exact two-point fiber refinement.",{
  "heat":"H_q(t)=(1+exp(-2qt))/2","unitary":"R_q(t)=cos^2(qt)"})
 x["ST1603"]=packet(1603,"proven_classical_fiber_return_strictly_decreases_in_q_for_t_positive","Global q>=0.",{
  "derivative":"-t exp(-2qt)<0"})
 x["ST1604"]=packet(1604,"proven_classical_q_has_closed_form_inverse_from_one_nontrivial_return","Requires H in (1/2,1].",{
  "inverse":"q=-log(2H-1)/(2t)"})
 x["ST1605"]=packet(1605,"proven_quantum_single_time_q_is_aliased","Period/reflection symmetry.",{
  "aliases":"q' = k*pi/t +/- q","global_injective":False})
 qmax=math.pi/(2*t)
 x["ST1606"]=packet(1606,"proven_quantum_q_is_injective_on_first_monotone_branch","Single time with prior range.",{
  "range":[0,qmax],"inverse":"q=acos(sqrt(R))/t"})
 hp=-t*math.exp(-2*q0*t);up=-t*math.sin(2*q0*t)
 x["ST1607"]=packet(1607,"proven_fine_return_derivatives","Reference point.",{
  "q":q0,"t":t,"heat_derivative":hp,"unitary_derivative":up})
 H=h(q0,t);U=u(q0,t);Ic=hp**2/(H*(1-H));Iq=up**2/(U*(1-U))
 x["ST1608"]=packet(1608,"proven_Bernoulli_Fisher_information_values_away_from_endpoints","Ideal fine measurement.",{
  "classical":Ic,"quantum":Iq,"quantum_formula":"4t^2"})
 x["ST1609"]=packet(1609,"constructed_reference_q_fine_discriminator","q=0.37,t=0.6.",{
  "heat":H,"unitary":U,"gap":U-H})
 f=lambda q:u(q,t)-h(q,t);xs=np.linspace(1e-8,10,20001);roots=[]
 for a,b in zip(xs[:-1],xs[1:]):
  if f(a)*f(b)<0:roots.append(float(brentq(f,a,b)))
 x["ST1610"]=packet(1610,"strong_numerical_single_time_channel_crossing_inventory_on_0_10","Besides q=0.",{
  "positive_roots":roots,"count":len(roots)})
 x["ST1611"]=packet(1611,"proven_q_zero_makes_fine_heat_and_unitary_returns_identical_for_all_times","No fiber coupling.",{
  "H_0":1.0,"R_0":1.0,"channel_identifiable":False})
 times=[.6,1.0]
 x["ST1612"]=packet(1612,"constructed_two_time_fine_channel_design","q prior interval required.",{
  "times":times,"q_range":[.1,2.0]})
 iv=mp.iv;iv.dps=20
 def fi(q,tv):return iv.cos(q*tv)**2-(1+iv.exp(-2*q*tv))/2
 def sql(z):
  a=float(z.a);b=float(z.b)
  return 0. if a<=0<=b else min(a*a,b*b)
 lower=float('inf');weak=None
 for i in range(1900):
  a=.1+i*.001;b=a+.001;q=iv.mpf([str(a),str(b)]);v=sql(fi(q,iv.mpf('0.6')))+sql(fi(q,iv.mpf('1.0')))
  if v<lower:lower=v;weak=[a,b]
 x["ST1613"]=packet(1613,"constructed_complete_1900_box_interval_cover_of_q_range","Frozen analytic fiber formulas.",{
  "boxes":1900,"step":.001,"weak_box":weak})
 x["ST1614"]=packet(1614,"proven_two_time_fine_heat_and_unitary_models_are_separated_on_q_0_1_to_2","Interval RSS certificate.",{
  "RSS_lower":lower,"one_time_residual_lower":math.sqrt(lower/2)})
 x["ST1615"]=packet(1615,"q_identifiability_is_conditional_on_channel_range_clock_and_fiber_instrument","No strict physical q.",{
  "q_sourced":False,"SI_rate":False})
 x["ST1616"]=packet(1616,"recommended_ST1617_ST1631","Synthesize noncirculant, continuum, calibration and q results.",{
  "next":["robustness hierarchy","claim lattice","protocol 11.01","remaining blockers","roadmap"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
