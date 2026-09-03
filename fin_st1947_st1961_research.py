#!/usr/bin/env python3
"""FIN ST1947--ST1961: Wilson action, locality and continuum scaling."""
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1947,1961
NAMES=["WilsonActionPositivity","SmallCurvatureQuadraticLimit","FaceScaleFreedom","FiniteHeatTrace",
 "FiniteSpectralDimensionLimits","IteratedRefinementHeatTrace","SpectralDimensionFormula","ConstantQGrowth",
 "GeometricQFreedom","DimensionFourTuningBoundary","SupportLocality","EuclideanSignature","MaxwellNormalizationGap",
 "RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 A=independent_strict_matrix_float();lam=np.linalg.eigvalsh(A);x={}
 x["ST1947"]=packet(1947,"proven_Wilson_face_action_nonnegative_for_positive_face_weights","1-Re U_f>=0.",{
  "action":"sum kappa_f(1-Re U_f)"})
 x["ST1948"]=packet(1948,"proven_small_holonomy_limit_is_positive_quadratic_curvature","For U_f=exp(iF_f).",{
  "limit":"1-cos F = F^2/2+O(F^4)"})
 x["ST1949"]=packet(1949,"proven_face_weights_and_overall_gauge_coupling_set_curvature_scale","Not fixed by matter Dirichlet weights.",{})
 ts=[1e-4,.01,.1,1,10]
 Z=[float(np.sum(np.exp(-t*lam))) for t in ts]
 x["ST1950"]=packet(1950,"constructed_exact_finite_heat_trace_from_A_spectrum","Dimensionless.",{"times":ts,"Z":Z})
 x["ST1951"]=packet(1951,"proven_finite_heat_trace_spectral_dimension_tends_zero_at_both_time_extremes","No asymptotic continuum power law for fixed finite graph.",{})
 x["ST1952"]=packet(1952,"proven_iterated_binary_refinement_heat_trace_factorization","Kronecker sums.",{
  "formula":"Z_k(t)=Z_12(t) product_l(1+exp(-2q_l t))"})
 x["ST1953"]=packet(1953,"proven_effective_spectral_dimension_depends_explicitly_on_all_free_q_l","Log derivative.",{
  "formula":"d_s=2t[sum lambda e^-tlambda/Z12 + sum 2q_l/(exp(2q_l t)+1)]"})
 t=.5;rows=[]
 for k in [1,2,4,8]:
  base=2*t*float(np.sum(lam*np.exp(-t*lam))/np.sum(np.exp(-t*lam)));fiber=2*t*k*(2/(np.exp(2*t)+1));rows.append({"k":k,"d_s":base+fiber})
 x["ST1954"]=packet(1954,"proven_constant_q_refinement_spectral_dimension_grows_with_number_of_layers_at_fixed_time","No stable four-dimensional limit.",{"rows":rows})
 x["ST1955"]=packet(1955,"proven_geometric_or_other_q_sequences_can_shape_intermediate_plateaus_arbitrarily","Requires scaling choice.",{})
 x["ST1956"]=packet(1956,"dimension_four_plateau_can_be_fitted_but_is_not_derived_without_q_l_and_length_scaling_law","Receiver fit not source.",{})
 x["ST1957"]=packet(1957,"proven_base_complete_support_keeps_eleven_nonlocal_neighbors_per_vertex_at_every_refinement","Product adds rather than removes them.",{})
 x["ST1958"]=packet(1958,"proven_constructed_Dirichlet_and_Wilson_actions_are_Euclidean_and_do_not_source_Lorentz_signature","Analytic continuation would be extra.",{})
 x["ST1959"]=packet(1959,"blocked_no_face_area_lattice_spacing_action_unit_or_gauge_coupling_for_Maxwell_normalization","CA remains external.",{})
 x["ST1960"]=packet(1960,"round_four_falsifies_finite_or_freely_refined_heat_trace_as_unique_3plus1_continuum_source","Conditional plateau models remain possible.",{
  "unique_dimension":False})
 x["ST1961"]=packet(1961,"recommended_ST1962_ST1976","Test canonical flag/Hodge repair against refinement functoriality.",{
  "next":["base flag Hodge","refined flag zero modes","product squares","weighted Hodge ambiguity","spectral action"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
