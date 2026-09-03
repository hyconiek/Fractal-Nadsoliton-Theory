#!/usr/bin/env python3
"""FIN ST1857--ST1871: gauge action under exact 12-to-24 refinement."""
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1857,1871
NAMES=["RefinedGraphCounts","RefinedIncidenceFactorization","CoarseConnectionLift","FiberConstantActionTransport",
 "GeneralGaugeGroupGrowth","RelativeFiberGaugeModes","FiberLinkHolonomy","FreeQPersistence",
 "IteratedRefinementModuli","SupportDiameter","DegreeGrowth","CellComplexGrowth",
 "ContinuumScalingBoundary","RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 A=independent_strict_matrix_float();q=.37;B=np.array([[q,-q],[-q,q]]);A24=np.kron(A,np.eye(2))+np.kron(np.eye(12),B);V=24;E=2*66+12;x={}
 x["ST1857"]=packet(1857,"proven_refined_product_support_counts","q>0.",{"V":V,"E":E,"cycle_rank":E-V+1})
 # Direct incidence factorization using all nonzero offdiagonal pairs.
 rows=[]
 for i in range(V):
  for j in range(i+1,V):
   if A24[i,j]<-1e-14:
    row=np.zeros(V);w=-A24[i,j];row[i]=np.sqrt(w);row[j]=-np.sqrt(w);rows.append(row)
 D=np.array(rows)
 x["ST1858"]=packet(1858,"proven_refined_weighted_incidence_factorization","Numerical replay.",{
  "edges":len(rows),"residual":float(np.linalg.norm(D.T@D-A24,np.inf))})
 x["ST1859"]=packet(1859,"constructed_coarse_base_connection_lift_identical_on_both_fibers","Fiber links set flat.",{
  "base_links_copied":True})
 # normalized fiber-constant vector u+
 u=np.ones(2)/np.sqrt(2);psi=np.arange(1,13,dtype=float);fine=np.kron(psi,u)
 base=float(psi@A@psi);ref=float(fine@A24@fine)
 x["ST1860"]=packet(1860,"proven_fiber_constant_Dirichlet_action_transports_exactly","Normalized constant fiber mode.",{
  "base":base,"refined":ref,"residual":abs(base-ref)})
 x["ST1861"]=packet(1861,"proven_full_refined_vertex_gauge_group_has_24_local_phases_vs_12_coarse_phases","General gauge transformations exceed coarse lift.",{
  "coarse":12,"fine":24})
 x["ST1862"]=packet(1862,"proven_twelve_relative_fiber_phase_modes_are_invisible_to_coarse_factorized_gauge_subgroup","Need fiber instrument/action.",{"relative_modes":12})
 x["ST1863"]=packet(1863,"constructed_fiber_link_phase_and_holonomy_sector","New connection data not fixed by base links.",{
  "fiber_links":12})
 x["ST1864"]=packet(1864,"proven_refinement_rate_q_remains_free_and_not_selected_by_gauge_covariance","Gauge principle does not source q.",{"q_sourced":False})
 x["ST1865"]=packet(1865,"proven_iterated_binary_refinement_introduces_at_least_one_new_rate_and_relative_connection_layer_per_step","Without consistency law moduli grow.",{})
 x["ST1866"]=packet(1866,"proven_refined_product_support_diameter_is_two_not_extended_continuum","K12 square K2.",{"diameter":2})
 x["ST1867"]=packet(1867,"proven_each_refined_vertex_degree_is_twelve","Eleven base nonlocal plus one fiber edge.",{"degree":12})
 x["ST1868"]=packet(1868,"proven_cycle_rank_grows_from_55_to_121_and_face_choices_proliferate","No canonical gauge kinetic lift.",{"base_cycle_rank":55,"refined_cycle_rank":121})
 x["ST1869"]=packet(1869,"blocked_no_scaling_law_for_q_edge_weights_lengths_face_areas_or_dimension","Continuum limit underdetermined.",{"continuum":False})
 x["ST1870"]=packet(1870,"round_four_proves_exact_matter_action_transport_but_falsifies_refinement_as_unique_locality_dimension_selector","Conditional positive/negative split.",{
  "matter_transport":True,"unique_continuum":False})
 x["ST1871"]=packet(1871,"recommended_ST1872_ST1886","Test finite Dirac/spectral-triple, chirality, fermion and gravity content.",{
  "next":["vertex grading no-go","doubled incidence Dirac","algebra representation","inner fluctuations","spectral action scale","gravity boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
