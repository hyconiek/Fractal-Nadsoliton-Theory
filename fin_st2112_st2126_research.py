#!/usr/bin/env python3
"""FIN ST2112--ST2126: refinement-intertwining d1 classification."""
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2112,2126
NAMES=["CochainRefinementMaps","ChainIntertwiningEquation","LiftedTriangleCoefficient","SquarePullbackZero",
 "LevelOneVerticalFamily","FiberSwapEffect","LevelTwoInheritedCoefficients","LevelTwoNewVerticalFamily",
 "TwoLevelDimensionLowerBound","PrimitiveIntegralRepair","ProductCoboundary","Associativity",
 "RefinementNoUniqueness","RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={}
 x["ST2112"]=packet(2112,"constructed_edge_and_face_cochain_pullbacks_for_binary_product_refinement","Base edges/faces copied to both fibers; vertical/square coarse components zero.",{})
 x["ST2113"]=packet(2113,"constructed_second_differential_intertwining_equation","R2 d1_base=d1_fine R1.",{})
 x["ST2114"]=packet(2114,"proven_intertwining_fixes_each_lifted_triangle_row_scale_equal_to_base_orbit_scale","For duplicate pullback convention.",{})
 x["ST2115"]=packet(2115,"proven_product_square_rows_annihilate_pulled_back_base_edge_cochains","Therefore their coefficients are absent from intertwining equations.",{})
 x["ST2116"]=packet(2116,"proven_level1_boundary_local_equivariant_intertwining_family_retains_six_vertical_square_orbit_scales","Given base d1.",{"free_vertical":6})
 x["ST2117"]=packet(2117,"fiber_swap_symmetry_relates_copies_but_does_not_fix_square_orbit_scales","No additional reduction.",{})
 x["ST2118"]=packet(2118,"proven_level2_intertwining_inherits_all_level1_coefficients_on_lifted_faces","Associative pullback.",{})
 x["ST2119"]=packet(2119,"proven_level2_adds_at_least_seven_new_vertical_edge_type_coefficients","Same orbit lower bound as Hodge weights.",{"new_lower_bound":7})
 x["ST2120"]=packet(2120,"proven_two_level_vertical_d1_moduli_dimension_at_least13_without_integrality","Six plus seven.",{})
 x["ST2121"]=packet(2121,"primitive_integral_cellular_incidence_sets_each_present_face_row_to_standard_plusminus_boundary","Conditional collapse of scales.",{})
 x["ST2122"]=packet(2122,"constructed_canonical_product_cellular_coboundary_given_base_complex_and_refinement_history","All incidence coefficients plusminus one.",{})
 x["ST2123"]=packet(2123,"proven_product_coboundary_is_associative_and_satisfies_d1d0_zero_at_all_levels","Cellular chain theorem.",{})
 x["ST2124"]=packet(2124,"proven_refinement_intertwining_alone_does_not_select_vertical_d1_coefficients_or_base_face_complex","Main no-go.",{})
 x["ST2125"]=packet(2125,"round_three_preserves_standard_cellular_d1_only_conditionally_on_integral_cell_structure","Strict source open.",{})
 x["ST2126"]=packet(2126,"recommended_ST2127_ST2141","Test positivity/Hodge row-rescaling gauge and spectral consequences.",{
  "next":["Q=d1*Md1","rescaling quotient","positive cone","H1 invariance","spectral nonuniqueness"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
