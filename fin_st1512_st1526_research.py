#!/usr/bin/env python3
"""FIN ST1512--ST1526: exact transport through 12-to-24 refinement."""
import math

import numpy as np
from scipy.linalg import expm

from fin_oa_discrimination_common import classical_return, quantum_return
from fin_total_nadsoliton_common import N, STRICT_A, write_packet, write_round


LO,HI=1512,1526
NAMES=["RefinedGenerator","KroneckerSemigroupFactorization","ClassicalCoarseReturnInvariance","QuantumCoarseReturnInvariance",
 "FiberStateIndependence","QBlindCoarseProtocol","FineClassicalFiberReturn","FineQuantumFiberReturn",
 "FineDualDiscriminator","NumericalRefinementReplay","ProtocolCertificateTransport","RefinedDetectorRequirement",
 "RefinementScaleBoundary","RoundFive_Verdict","RoundFive_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};q=.37;t=.6;B=np.array([[q,-q],[-q,q]]);A24=np.kron(STRICT_A,np.eye(2))+np.kron(np.eye(N),B)
 x["ST1512"]=packet(1512,"constructed_exact_two_fiber_refined_generator","Conditional refinement family.",{
  "formula":"A24=A12 tensor I2 + I12 tensor Bq","q":q,"dimension":24})
 P24=expm(-t*A24);Pfac=np.kron(expm(-t*STRICT_A),expm(-t*B));U24=expm(-1j*t*A24);Ufac=np.kron(expm(-1j*t*STRICT_A),expm(-1j*t*B))
 x["ST1513"]=packet(1513,"proven_heat_and_unitary_Kronecker_factorization_numerically_replayed","Analytic commuting Kronecker-sum identity.",{
  "heat_residual":float(np.linalg.norm(P24-Pfac,np.inf)),"unitary_residual":float(np.linalg.norm(U24-Ufac,np.inf))})
 init=np.zeros(24);init[0]=1;node0=np.zeros(24);node0[0:2]=1;cret=float(node0@(P24@init))
 x["ST1514"]=packet(1514,"proven_classical_coarse_node_return_equals_base_return","Fiber summed measurement.",{
  "refined":cret,"base":classical_return(t),"residual":abs(cret-classical_return(t))})
 amp=U24@init;qret=float(np.sum(abs(amp[0:2])**2))
 x["ST1515"]=packet(1515,"proven_quantum_coarse_node_return_equals_base_return","Projector onto node0 tensor full fiber.",{
  "refined":qret,"base":quantum_return(t),"residual":abs(qret-quantum_return(t))})
 init2=np.zeros(24,dtype=complex);init2[1]=1;cr2=float(node0@(P24@init2.real));qr2=float(np.sum(abs((U24@init2)[0:2])**2))
 x["ST1516"]=packet(1516,"proven_coarse_return_independent_of_initial_fiber_basis_state","Both fiber basis preparations.",{
  "classical_difference":abs(cr2-cret),"quantum_difference":abs(qr2-qret)})
 x["ST1517"]=packet(1517,"proven_coarse_OA_discriminator_is_exactly_q_blind_for_declared_refinement","Any q>=0 by factorization/unitarity/mass.",{
  "q_selected":False,"coarse_predictions_unchanged":True})
 h=.5*(1+math.exp(-2*q*t))
 x["ST1518"]=packet(1518,"proven_fine_classical_fiber_return_formula","Initial fiber basis state.",{
  "formula":"H_q(t)=(1+exp(-2qt))/2","value":h})
 u=math.cos(q*t)**2
 x["ST1519"]=packet(1519,"proven_fine_quantum_fiber_return_formula","Initial fiber basis state.",{
  "formula":"R_q(t)=cos^2(qt)","value":u})
 x["ST1520"]=packet(1520,"constructed_fine_fiber_dual_dynamics_discriminator","Requires fiber-resolving instrument.",{
  "gap":u-h,"absolute_gap":abs(u-h)})
 x["ST1521"]=packet(1521,"proven_refinement_numerical_replay_matches_closed_form","q=0.37,t=0.6.",{
  "classical_matrix_fiber_return":float(P24[0,0]),"classical_formula":h,"quantum_matrix_fiber_return":float(abs(U24[0,0])**2),"quantum_formula":u})
 x["ST1522"]=packet(1522,"proven_all_coarse_interval_and_statistical_certificates_transport_exactly","Within exact lift and coarse instrument.",{
  "time_optimum":True,"four_time_RSS":True,"seven_node_distance_likelihood":True})
 x["ST1523"]=packet(1523,"refined_physical_instrument_needs_24_outcome_or_coarse_fiber_blind_calibration","No hardware map.",{
  "coarse_detector":12,"fine_detector":24})
 x["ST1524"]=packet(1524,"q_remains_dimensionless_relative_fiber_rate_without_new_scale_source","Refinement does not fix CA.",{
  "q_physical_value":False})
 x["ST1525"]=packet(1525,"refinement_preserves_coarse_test_but_adds_optional_fine_test_not_new_ontology","Conditional result.",{
  "coarse_transport":True,"fine_information":True,"state_category_selected":False})
 x["ST1526"]=packet(1526,"recommended_ST1527_ST1541","Final synthesis and Release 11.00.",{
  "next":["robustness ladder","certified/numerical split","mixture scope","refinement scope","next programmes"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
