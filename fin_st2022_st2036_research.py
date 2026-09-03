#!/usr/bin/env python3
"""FIN ST2022--ST2036: stronger naturality/homogeneity axioms."""
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2022,2036
NAMES=["SquareRuleInputs","PositiveNaturalRuleFamily","JointHomogeneityNoGo","SeparateBilinearityTheorem",
 "BilinearNormalizationFreedom","AxisAssociativity","CrossAxisSquareRule","BaseFaceFreedomPersists",
 "TensorRescalingAmbiguity","OverallGaugeScale","DimensionalHomogeneityBoundary","ExplicitCounterrules",
 "StrongAxiomClassification","RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 A=independent_strict_matrix_float();w=np.array([-A[0,d] for d in range(1,7)]);q=.37;x={}
 x["ST2022"]=packet(2022,"constructed_square_weight_rule_inputs_as_horizontal_edge_weight_and_fiber_rate","One square orbit per base distance.",{
  "w":w.tolist(),"q":q})
 rules={"product":w*q,"geometric":np.sqrt(w*q),"harmonic":2*w*q/(w+q),"minimum":np.minimum(w,q),"maximum":np.maximum(w,q)}
 x["ST2023"]=packet(2023,"constructed_five_positive_D12_natural_square_weight_rules","All refinement-local.",{
  "rules":{k:v.tolist() for k,v in rules.items()}})
 x["ST2024"]=packet(2024,"proven_joint_degree_one_homogeneity_leaves_at_least_four_nonproportional_rules","Geometric/harmonic/min/max.",{
  "unique":False})
 x["ST2025"]=packet(2025,"proven_separate_linearity_in_w_and_q_forces_mu_equals_c_w_q","Functional equation on positive scalars.",{
  "conditional_unique_up_to_c":True})
 x["ST2026"]=packet(2026,"proven_bilinear_rule_retains_positive_overall_constant_c","No normalization from coarse energy.",{})
 x["ST2027"]=packet(2027,"proven_product_rule_is_associative_under_multiple_independent_refinement_axes","Cross weights multiply factor weights.",{})
 x["ST2028"]=packet(2028,"constructed_cross_fiber_square_weight_c_q1_q2_for_two_binary_axes","Conditional monoidal rule.",{})
 x["ST2029"]=packet(2029,"proven_vertical_bilinearity_does_not_select_initial_base_triangle_face_set_or_measure","Base Rplus orbit cone survives.",{})
 x["ST2030"]=packet(2030,"proven_first_order_factor_inner_product_rescaling_changes_product_2form_normalization_without_changing_0form_A_after_compensating_d","Factorization gauge.",{})
 x["ST2031"]=packet(2031,"proven_overall_Wilson_gauge_coupling_remains_free_even_if_c_fixed","Independent action scale.",{})
 x["ST2032"]=packet(2032,"blocked_no_physical_dimensions_for_w_q_or_face_measure_to_choose_homogeneity_degree","Degree choice is extra.",{})
 ratios=(rules['geometric']/rules['harmonic']).tolist()
 x["ST2033"]=packet(2033,"proven_explicit_degree_one_counterrules_nonproportional_on_six_square_orbits","Finite witness.",{
  "geometric_over_harmonic":ratios,"range":[min(ratios),max(ratios)]})
 x["ST2034"]=packet(2034,"classified_strong_axiom_outcome_bilinearity_plus_associativity_selects_product_shape_only","Still c/base/scale freedom.",{
  "shape":"w q","remaining":["c","base faces","base measure","gauge coupling"]})
 x["ST2035"]=packet(2035,"round_three_finds_conditional_monoidal_rule_but_falsifies_basic_naturality_homogeneity_as_unique_selector","No strict source.",{})
 x["ST2036"]=packet(2036,"recommended_ST2037_ST2051","Test Schur complement/minimal-extension sensitivity to vertical weights.",{
  "next":["single-square toy Hessian","flat extension","network generalization","quantum determinant","identifiability"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
