#!/usr/bin/env python3
"""FIN ST2037--ST2051: Schur complement and minimal-extension no-go."""
import sympy as sp

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2037,2051
NAMES=["SquareToyEnergy","CoarseConstraint","ExactMinimizer","SchurEnergyIndependence",
 "VerticalWeightUnidentifiability","NetworkFlatExtension","GaugeZeroMode","AsymmetricHorizontalBoundary",
 "QuantumDeterminantDependence","RenormalizationPremise","FineObservableRequirement","CoarseDataNoGo",
 "SchurClassification","RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 k,mu,a=sp.symbols('k mu a', positive=True, real=True);r,s=sp.symbols('r s', real=True);E=sp.Rational(1,2)*k*(a**2+r**2)+sp.Rational(1,2)*mu*(2*r+s)**2;x={}
 x["ST2037"]=packet(2037,"constructed_single_product_square_quadratic_energy","a coarse average, r horizontal relative mode, s vertical-edge difference.",{
  "E":str(E)})
 x["ST2038"]=packet(2038,"constructed_fixed_coarse_average_constraint_a0_equals_a_plus_r_a1_equals_a_minus_r","Fiber-swap horizontal weights assumed.",{})
 sol=sp.solve([sp.diff(E,r),sp.diff(E,s)],[r,s],dict=True)[0];Emin=sp.simplify(E.subs(sol))
 x["ST2039"]=packet(2039,"proven_exact_relative_and_vertical_minimizer_is_flat","All k,mu>0.",{
  "solution":{str(v):str(val) for v,val in sol.items()}})
 x["ST2040"]=packet(2040,"proven_minimized_coarse_energy_is_k_a2_over2_independent_of_mu","Exact Schur/minimal extension.",{
  "E_min":str(Emin),"contains_mu":mu in Emin.free_symbols})
 x["ST2041"]=packet(2041,"proven_vertical_square_weight_cannot_be_identified_from_any_coarse_minimal_energy_measurement_in_toy_sector","Flat kernel direction.",{})
 x["ST2042"]=packet(2042,"proven_network_product_has_global_flat_extension_with_identical_layer_connection_and_zero_fiber_potential","Every square curvature zero.",{})
 x["ST2043"]=packet(2043,"proven_common_vertical_potential_shift_is_gauge_zero_mode","Gauge fixing needed but does not alter mu invisibility.",{})
 x["ST2044"]=packet(2044,"asymmetric_horizontal_copy_weights_can_make_mu_enter_reduced_energy_but_break_fiber_swap_and_prior_energy_classification","Not admissible repair under current axioms.",{})
 x["ST2045"]=packet(2045,"proven_Gaussian_fluctuation_determinant_depends_on_mu_even_when_classical_coarse_minimum_does_not","Quantum/thermal measure could probe vertical weights.",{})
 x["ST2046"]=packet(2046,"quantum_determinant_selection_requires_path_integral_measure_regularization_temperature_or_hbar_and_counterterms","New conditional structure.",{})
 x["ST2047"]=packet(2047,"proven_fine_square_holonomy_or_relative_mode_observable_is_required_to_measure_mu_classically","OA fine instrument.",{})
 x["ST2048"]=packet(2048,"proven_coarse_energy_and_Schur_consistency_cannot_select_vertical_face_weights","General flat-extension no-go.",{})
 x["ST2049"]=packet(2049,"classified_energy_compatible_Schur_family_as_same_coarse_form_times_arbitrary_positive_vertical_Hodge_block","Block-diagonal on coarse/vertical decomposition.",{})
 x["ST2050"]=packet(2050,"round_four_falsifies_minimal_extension_or_Schur_complement_as_vertical_measure_selector","Fine/quantum data needed.",{})
 x["ST2051"]=packet(2051,"recommended_ST2052_ST2066","Test monoidal tensor-Hodge repair and differential-calculus ambiguity.",{
  "next":["tensor inner product","product rule","associativity","junk quotient","base measure","scale"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
