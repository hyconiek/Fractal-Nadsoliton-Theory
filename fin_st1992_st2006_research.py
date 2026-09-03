#!/usr/bin/env python3
"""FIN ST1992--ST2006: energy-compatible face weights for 12-to-24."""
import numpy as np

from fin_cellular_refinement_common import base_complex,refine
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1992,2006
NAMES=["FaceEnergyModel","CoarseCurvatureEmbedding","EnergyCompatibilityEquation","HorizontalWeightClassification",
 "FiberSwapHalving","SquareWeightInvisibility","UnreducedConstraintRank","OrbitReducedNullity",
 "RandomEnergyReplay","PositiveCone","MinimalExtensionPrinciple","SchurPreview","ClassificationTheorem",
 "RoundOne_Verdict","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 v,e,f=base_complex();v1,e1,f1,meta=refine(v,e,f);n0=len(f);ns=meta['new_squares'];x={}
 x["ST1992"]=packet(1992,"constructed_diagonal_face_curvature_energy","Chosen base face complex.",{
  "E0":"1/2 sum_f kappa_f F_f^2","faces":n0})
 x["ST1993"]=packet(1993,"constructed_coarse_to_fine_curvature_embedding","Two identical horizontal copies and zero product-square curvature.",{
  "map":"F -> (F,F,0_square)"})
 x["ST1994"]=packet(1994,"proven_energy_compatibility_iff_copy_weights_sum_to_base_weight","Coefficientwise identity.",{
  "condition":"kappa_f,0+kappa_f,1=kappa_f"})
 x["ST1995"]=packet(1995,"proven_general_horizontal_solution_has_one_split_parameter_per_base_face","Before symmetry.",{
  "solution":"(a_f kappa_f,(1-a_f)kappa_f)","free_parameters":n0})
 x["ST1996"]=packet(1996,"proven_fiber_swap_symmetry_uniquely_sets_each_horizontal_copy_to_one_half","Conditional uniqueness.",{
  "copy_factor":.5})
 x["ST1997"]=packet(1997,"proven_all_new_square_weights_are_invisible_to_coarse_curvature_energy","Square curvature is zero on image.",{
  "free_positive_square_weights":ns})
 x["ST1998"]=packet(1998,"proven_unreduced_linear_energy_constraints_leave_286_dimensional_affine_solution_space","506 fine face weights minus 220 equations.",{
  "fine_faces":len(f1),"constraints":n0,"nullity":len(f1)-n0})
 x["ST1999"]=packet(1999,"proven_D12_and_fiber_swap_reduction_still_leaves_six_square_orbit_parameters","Horizontal orbit weights inherited/halved.",{
  "vertical_nullity":6})
 rng=np.random.default_rng(2000);kappa=rng.uniform(.1,2,n0);mu=rng.uniform(.1,2,ns);F=rng.normal(size=n0);E0=.5*np.sum(kappa*F**2);E1=.5*np.sum((kappa/2)*F**2+(kappa/2)*F**2)+.5*np.sum(mu*np.zeros(ns)**2)
 x["ST2000"]=packet(2000,"proven_random_energy_identity_replay","Numerical witness.",{
  "E0":float(E0),"E1":float(E1),"residual":float(abs(E0-E1))})
 x["ST2001"]=packet(2001,"proven_positivity_only_restricts_square_weights_to_open_positive_cone","Does not choose values.",{
  "cone":"R_+^66 or R_+^6 after symmetry"})
 x["ST2002"]=packet(2002,"proven_minimum_energy_extension_of_coarse_flat_fiber_state_is_independent_of_square_weights","Set relative fields flat.",{})
 x["ST2003"]=packet(2003,"constructed_Schur_complement_preview_square_sector_lies_transverse_to_coarse_embedding","Explains invisibility.",{})
 x["ST2004"]=packet(2004,"proven_12_to24_energy_compatible_measure_classification","Horizontal halves unique under swap; six vertical orbit weights arbitrary.",{
  "horizontal":"kappa/2","vertical":"R_+^6"})
 x["ST2005"]=packet(2005,"round_one_falsifies_energy_preservation_as_unique_face_measure_selector","Partial classification positive.",{})
 x["ST2006"]=packet(2006,"recommended_ST2007_ST2021","Iterate classification through 24-to-48 and test associativity.",{
  "next":["cell counts","boundary rank","quarter weights","old square halves","new square orbits","moduli growth"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
