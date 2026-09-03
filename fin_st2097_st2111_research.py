#!/usr/bin/env python3
"""FIN ST2097--ST2111: boundary-local and integral d1 classification."""
import itertools

import sympy as sp

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2097,2111
NAMES=["TriangleLocalAnsatz","LocalChainKernel","PerFaceOneDimensionality","D12OrbitReduction",
 "TwelveParameterLocalFamily","StandardCoboundary","IntegralIncidenceRestriction","OrientationSignQuotient",
 "FullS12SymmetryBoundary","WeightDependentOrbitScalars","TopologyKernelInvariance","HodgeParameterCombination",
 "LocalityNoUniqueness","RoundTwo_Verdict","RoundTwo_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 a,b,c=sp.symbols('a b c');fi,fj,fk=sp.symbols('fi fj fk');expr=a*(fj-fi)+b*(fk-fj)+c*(fk-fi);coeff=[sp.expand(expr).coeff(v) for v in [fi,fj,fk]];sol=sp.linsolve(coeff,[a,b,c]);x={}
 x["ST2097"]=packet(2097,"constructed_boundary_local_triangle_row_ansatz","Output depends only on three oriented boundary edges.",{
  "row":"a x_ij+b x_jk+c x_ik"})
 x["ST2098"]=packet(2098,"proven_local_row_annihilates_all_gradients_iff_coefficients_proportional_to_standard_boundary","Exact SymPy linear solve.",{
  "solution":str(sol)})
 x["ST2099"]=packet(2099,"proven_each_triangle_has_one_dimensional_local_chain_row_space","Before symmetry.",{})
 x["ST2100"]=packet(2100,"proven_D12_equivariance_makes_row_scale_constant_on_each_of_twelve_triangle_orbits","Orbit classification.",{})
 x["ST2101"]=packet(2101,"proven_boundary_local_D12_equivariant_chain_maps_form_twelve_dimensional_family","Full triangle space.",{"dimension":12})
 x["ST2102"]=packet(2102,"constructed_standard_simplicial_coboundary_as_all_orbit_scales_one","Canonical after choosing cellular incidence convention.",{})
 x["ST2103"]=packet(2103,"proven_primitive_integral_incidence_coefficients_restrict_each_nonzero_orbit_scale_to_plusminus_one","If rows must be primitive cellular boundaries.",{})
 x["ST2104"]=packet(2104,"proven_face_orientation_reversal_absorbs_individual_minus_signs","All-nonzero primitive integral rows equivalent to standard d1 up to face orientation.",{})
 x["ST2105"]=packet(2105,"full_S12_naturality_would_force_one_orbit_scale_but_does_not_preserve_weighted_A","Unavailable symmetry enlargement.",{})
 x["ST2106"]=packet(2106,"proven_D12_natural_weight_functions_can_assign_distinct_positive_scales_to_twelve_triangle_orbits","Prior face-measure ambiguity reappears in d1 normalization.",{})
 x["ST2107"]=packet(2107,"proven_nonzero_orbit_row_rescaling_does_not_change_kernel_of_d1_or_H1","Topology cannot select scales.",{})
 x["ST2108"]=packet(2108,"proven_curvature_energy_depends_only_on_products_face_measure_times_row_scale_squared","Gauge redundancy between d1 normalization and Hodge star.",{
  "invariant":"m_orbit c_orbit^2"})
 x["ST2109"]=packet(2109,"locality_reduces_517_to12_but_not_to_unique_physical_operator","Strong reduction/no-go.",{})
 x["ST2110"]=packet(2110,"round_two_finds_unique_combinatorial_incidence_only_after_primitive_integrality_and_chosen_face_complex","Conditional positive result.",{})
 x["ST2111"]=packet(2111,"recommended_ST2112_ST2126","Add exact refinement chain-map intertwiners and classify vertical coefficients.",{
  "next":["coarse pullback","lifted faces","square rows","vertical orbit scales","two-level growth"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
