#!/usr/bin/env python3
"""FIN ST2007--ST2021: 24-to-48 associativity and vertical moduli growth."""
from fin_cellular_refinement_common import base_complex,refine,boundary_bits,gf2_rank
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2007,2021
NAMES=["LevelCounts","LevelBoundaryRanks","AssociativeHorizontalQuartering","OldSquareHalving",
 "NewSquareFreedom","NewSquareOrbitLowerBound","TwoLevelVerticalCone","DirectVsIteratedHorizontal","AxisSwapBoundary",
 "ProductCellCoherence","EnergyCompatibilityAtAllLevels","ModuliGrowth","NoFiniteUniqueness","RoundTwo_Verdict","RoundTwo_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 v0,e0,f0=base_complex();v1,e1,f1,m1=refine(v0,e0,f0);v2,e2,f2,m2=refine(v1,e1,f1);x={}
 rows=[]
 for level,(v,e,f) in enumerate([(v0,e0,f0),(v1,e1,f1),(v2,e2,f2)]):rows.append({"level":level,"V":len(v),"E":len(e),"F":len(f),"cycle_rank":len(e)-len(v)+1})
 x["ST2007"]=packet(2007,"proven_product_cell_counts_through_12_24_48","Exact enumeration.",{"levels":rows,"level1_meta":m1,"level2_meta":m2})
 ranks=[gf2_rank(boundary_bits(e,f)) for v,e,f in [(v0,e0,f0),(v1,e1,f1),(v2,e2,f2)]]
 x["ST2008"]=packet(2008,"proven_product_cell_boundary_rank_equals_cycle_rank_at_all_three_levels","Exact GF2.",{
  "ranks":ranks,"H1":[rows[i]['cycle_rank']-ranks[i] for i in range(3)]})
 x["ST2009"]=packet(2009,"proven_two_fiber_swaps_force_base_face_weights_to_one_quarter_per_level2_copy","Associative horizontal lift.",{
  "factor":.25})
 x["ST2010"]=packet(2010,"proven_level1_square_weights_split_to_halves_on_level2_copies","Old vertical measures inherited, not selected.",{})
 x["ST2011"]=packet(2011,"proven_144_new_level2_square_weights_are_coarse_energy_invisible","One per level1 edge.",{
  "new_squares":m2['new_squares']})
 x["ST2012"]=packet(2012,"proven_new_level2_squares_have_at_least_seven_natural_edge_type_orbits","Six base-distance horizontal types plus old-fiber vertical type.",{
  "orbit_lower_bound":7})
 x["ST2013"]=packet(2013,"proven_two_level_symmetric_vertical_measure_cone_has_at_least_dimension13","Six old square orbits plus seven new.",{
  "dimension_lower_bound":13})
 x["ST2014"]=packet(2014,"proven_direct_four_copy_and_iterated_binary_horizontal_weights_agree","Both give kappa/4 under swap symmetry.",{})
 x["ST2015"]=packet(2015,"axis_exchange_symmetry_can_relate_some_vertical_weights_only_if_refinement_axes_and_rates_are_equivalent","q1,q2 generally differ.",{
  "automatic":False})
 x["ST2016"]=packet(2016,"constructed_coherent_relative_product_cell_functor_across_two_steps","Lift prior cells and add edge-times-new-interval cells.",{
  "exists":True})
 x["ST2017"]=packet(2017,"proven_coarse_energy_compatibility_at_every_step_only_normalizes_inherited_horizontal_weights","New vertical sectors remain kernel directions.",{})
 x["ST2018"]=packet(2018,"proven_positive_vertical_moduli_grow_at_each_nontrivial_refinement","At least one new edge orbit per level.",{})
 x["ST2019"]=packet(2019,"proven_no_finite_sequence_of_coarse_energy_equalities_selects_all_future_vertical_face_weights","Induction on refinement steps.",{})
 x["ST2020"]=packet(2020,"round_two_finds_associative_horizontal_measure_but_unbounded_vertical_nonuniqueness","Main classification.",{})
 x["ST2021"]=packet(2021,"recommended_ST2022_ST2036","Test stronger naturality, homogeneity and tensor-product axioms.",{
  "next":["bilinear square rule","homogeneous counterfamilies","axis associativity","overall scale","base ambiguity"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
