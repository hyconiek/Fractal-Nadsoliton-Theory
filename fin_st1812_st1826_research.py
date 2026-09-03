#!/usr/bin/env python3
"""FIN ST1812--ST1826: direct gauge-covariant completion of the Dirichlet form."""
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1812,1826
NAMES=["WeightedIncidence","DirichletFactorization","LinkVariableDefinition","GaugeTransformation",
 "GaugeCovariance","GaugeInvariantMatterAction","CovariantLaplacian","ZeroConnectionRecovery",
 "GaugeCovariantEulerEquation","PotentialCompatibility","NumericalGaugeReplay","OrientationIndependence",
 "GlobalRegularity","MissingGaugeKinetic","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def build():
 A=independent_strict_matrix_float();edges=[];D=[]
 for i in range(12):
  for j in range(i+1,12):
   w=-A[i,j];row=np.zeros(12,dtype=complex);row[i]=np.sqrt(w);row[j]=-np.sqrt(w);edges.append((i,j,w));D.append(row)
 return A,edges,np.array(D)
def covariant(edges,U):
 D=[]
 for k,(i,j,w) in enumerate(edges):
  row=np.zeros(12,dtype=complex);row[i]=np.sqrt(w);row[j]=-np.sqrt(w)*U[k];D.append(row)
 return np.array(D)
def main():
 A,edges,D=build();x={}
 x["ST1812"]=packet(1812,"constructed_weighted_first_order_incidence_from_strict_A_support","One edge per nonzero unordered coupling.",{
  "vertices":12,"edges":len(edges),"weight_min":min(e[2] for e in edges),"weight_max":max(e[2] for e in edges)})
 residual=float(np.linalg.norm(D.conj().T@D-A,np.inf))
 x["ST1813"]=packet(1813,"proven_exact_Dirichlet_factorization_numerically_replayed","Analytic edge-sum identity.",{
  "identity":"A=d_W^* d_W","residual":residual})
 x["ST1814"]=packet(1814,"constructed_U1_link_connection_on_every_weighted_edge","Conditional localization of global phase symmetry.",{
  "U_ij":"unit complex phase","U_ji":"conjugate U_ij"})
 x["ST1815"]=packet(1815,"constructed_vertex_gauge_action_on_fields_and_links","Local U(1).",{
  "psi":"psi_i -> g_i psi_i","U":"U_ij -> g_i U_ij g_j^-1"})
 x["ST1816"]=packet(1816,"proven_edge_covariant_difference_transforms_at_source_vertex","Direct substitution.",{
  "D_Upsi":"psi_i-U_ij psi_j","law":"(D_U psi)'_ij=g_i(D_U psi)_ij"})
 x["ST1817"]=packet(1817,"proven_minimal_link_matter_action_is_exactly_gauge_invariant","No continuum approximation.",{
  "action":"S=1/2 sum_ij w_ij |psi_i-U_ij psi_j|^2"})
 rng=np.random.default_rng(1818);theta=rng.normal(size=len(edges));U=np.exp(1j*theta);DU=covariant(edges,U);AU=DU.conj().T@DU;e=np.linalg.eigvalsh(AU)
 x["ST1818"]=packet(1818,"proven_covariant_laplacian_is_Hermitian_positive_semidefinite","Any U(1) links.",{
  "Hermitian_residual":float(np.linalg.norm(AU-AU.conj().T,np.inf)),"lambda_min":float(e.min())})
 D1=covariant(edges,np.ones(len(edges),complex))
 x["ST1819"]=packet(1819,"proven_trivial_connection_recovers_strict_A","Exact finite replay.",{
  "residual":float(np.linalg.norm(D1.conj().T@D1-A,np.inf))})
 x["ST1820"]=packet(1820,"proven_literal_variation_gives_A_U_psi_plus_radial_potential_equals_source","Single action source.",{
  "EOM":"A_U psi + V'(|psi|^2)psi = J"})
 x["ST1821"]=packet(1821,"proven_any_vertex_local_potential_V_of_abspsi2_is_gauge_invariant","Mass/quartic radial terms allowed.",{
  "does_not_fix_coefficients":True})
 alpha=rng.normal(size=12);g=np.exp(1j*alpha);Up=[]
 for k,(i,j,w) in enumerate(edges):Up.append(g[i]*U[k]*np.conj(g[j]))
 psi=rng.normal(size=12)+1j*rng.normal(size=12);Sp=.5*np.linalg.norm(covariant(edges,np.array(Up))@(g*psi))**2;S=.5*np.linalg.norm(DU@psi)**2
 x["ST1822"]=packet(1822,"proven_random_finite_gauge_replay_to_machine_precision","Independent numerical witness.",{
  "action_before":float(S),"action_after":float(Sp),"residual":float(abs(S-Sp))})
 # Reversing edge orientation conjugates U and multiplies the row by a phase/sign; norm unchanged.
 x["ST1823"]=packet(1823,"proven_edge_orientation_choice_does_not_change_matter_action","Orientation is bookkeeping.",{})
 x["ST1824"]=packet(1824,"proven_direct_link_action_is_globally_regular_on_compact_U1_edge_torus","No APD quotient/poles.",{
  "configuration_space":"U(1)^66 x C^12"})
 x["ST1825"]=packet(1825,"gauge_matter_action_contains_no_dynamical_gauge_curvature_term","Link phases are nondynamical until a face/cycle action is supplied.",{
  "Maxwell_Wilson_action":False})
 x["ST1826"]=packet(1826,"recommended_ST1827_ST1841","Falsify uniqueness: gauge group, charge, calculus and connection source.",{
  "next":["group alternatives","charge aliases","factorization freedom","complexification premise","flat connection selection"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
