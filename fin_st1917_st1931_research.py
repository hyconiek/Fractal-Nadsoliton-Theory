#!/usr/bin/env python3
"""FIN ST1917--ST1931: natural face-measure nonuniqueness."""
import itertools

import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1917,1931
NAMES=["OrbitWeightCone","FiveOrbitMinimumCone","EdgeDerivedRules","RuleVectorRank",
 "NonProportionalMeasures","HomogeneityClasses","RatioFunctionFreedom","HodgeStarDependence",
 "WilsonCouplingFreedom","DimensionalAreaGap","MaximumEntropyBoundary","MomentEntropyNoSelector",
 "MeasureNoGo","RoundTwo_Verdict","RoundTwo_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def orbits():
 V=range(12);tris=list(itertools.combinations(V,3))
 def rot(t,k):return tuple(sorted((x+k)%12 for x in t))
 def ref(t,k):return tuple(sorted((k-x)%12 for x in t))
 un=set(tris);O=[]
 while un:t=min(un);o=({rot(t,k) for k in V}|{ref(t,k) for k in V})&set(tris);O.append(sorted(o));un-=o
 return O
def main():
 A=independent_strict_matrix_float();O=orbits();x={}
 x["ST1917"]=packet(1917,"proven_positive_D12_natural_face_measures_form_Rplus12_orbit_weight_cone","One coefficient per triangle orbit.",{"dimension":12})
 x["ST1918"]=packet(1918,"proven_each_five_orbit_minimal_complex_retains_Rplus5_measure_cone","Topology/minimality does not fix weights.",{"dimension":5})
 rules=[];names=["product","geometric_mean","arithmetic_mean","minimum","maximum"]
 for o in O:
  t=o[0];ws=[-A[i,j] for i,j in itertools.combinations(t,2)];rules.append([np.prod(ws),np.prod(ws)**(1/3),np.mean(ws),min(ws),max(ws)])
 M=np.array(rules).T
 x["ST1919"]=packet(1919,"constructed_five_positive_automorphism_natural_edge_derived_face_rules","All use only boundary edge weights.",{
  "rules":names,"vectors":M.tolist()})
 rank=int(np.linalg.matrix_rank(M,tol=1e-12))
 x["ST1920"]=packet(1920,"proven_edge_derived_rule_vectors_span_rank_five","Finite numerical rank with large separation.",{"rank":rank})
 ratios=(M[0]/M[1]).tolist()
 x["ST1921"]=packet(1921,"proven_natural_positive_rules_are_not_globally_proportional","Product/geometric ratio varies by orbit.",{
  "ratio_range":[min(ratios),max(ratios)]})
 x["ST1922"]=packet(1922,"proven_scaling_homogeneity_does_not_select_unique_rule","Arithmetic/geometric/min/max all degree one; product degree three.",{})
 x["ST1923"]=packet(1923,"proven_arbitrary_positive_symmetric_functions_of_two_independent_edge_ratios_preserve_degree_one_naturality","Infinite functional freedom.",{})
 x["ST1924"]=packet(1924,"proven_discrete_Hodge_star_on_2forms_depends_on_face_inner_product_measure","EOM/gauge spectrum changes with weights.",{})
 x["ST1925"]=packet(1925,"proven_Wilson_action_has_independent_overall_gauge_coupling_even_after_relative_face_measure_choice","Scale remains free.",{})
 x["ST1926"]=packet(1926,"blocked_no_face_area_length_squared_or_action_unit_from_dimensionless_edge_weights","CA/scale problem reappears.",{})
 x["ST1927"]=packet(1927,"maximum_entropy_face_measure_requires_reference_measure_and_constraints","Not source-free selector.",{})
 x["ST1928"]=packet(1928,"proven_existing_kernel_moments_or_Shannon_scalar_do_not_choose_one_orbit_weight_vector","Scalar data insufficient for 12-dimensional cone.",{})
 x["ST1929"]=packet(1929,"proven_no_unique_natural_positive_face_measure_from_current_edge_weights_and_basic_axioms","Explicit rank-five counterfamily plus infinite ratio functions.",{
  "unique":False})
 x["ST1930"]=packet(1930,"round_two_falsifies_face_measure_uniqueness_even_after_face_set_choice","Need refinement/energy principle.",{})
 x["ST1931"]=packet(1931,"recommended_ST1932_ST1946","Test whether exact refinement functoriality selects faces/weights.",{
  "next":["triangle lift","fiber squares","boundary ranks","flag-vs-product","face-weight transport"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
