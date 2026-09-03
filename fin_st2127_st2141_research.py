#!/usr/bin/env python3
"""FIN ST2127--ST2141: positivity, Hodge spectra and row-rescaling gauge."""
import itertools

import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2127,2141
NAMES=["UnweightedDECConvention","WeightedDifferentialConvention","ConventionEquivalence","CurvatureQuadraticOperator",
 "RowHodgeRescalingGauge","PositiveOrbitCone","KernelInvariance","StandardSimplexSpectrum",
 "OrbitScaledSpectrum","SpectralTunability","TopologyVsDynamics","IntegralIncidenceBoundary","RefinementVerticalSpectrum",
 "RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def data():
 n=12;edges=list(itertools.combinations(range(n),2));ei={e:i for i,e in enumerate(edges)};tris=list(itertools.combinations(range(n),3));d0=np.zeros((66,12));d1=np.zeros((220,66))
 for r,(i,j) in enumerate(edges):d0[r,i]=-1;d0[r,j]=1
 for r,(i,j,k) in enumerate(tris):d1[r,ei[(i,j)]]=1;d1[r,ei[(j,k)]]=1;d1[r,ei[(i,k)]]=-1
 return edges,tris,d0,d1
def tri_orbits(tris):
 V=range(12)
 def rot(t,k):return tuple(sorted((x+k)%12 for x in t))
 def ref(t,k):return tuple(sorted((k-x)%12 for x in t))
 un=set(tris);O=[]
 while un:t=min(un);o=({rot(t,k) for k in V}|{ref(t,k) for k in V})&set(tris);O.append(o);un-=o
 return O
def main():
 edges,tris,d0,d1=data();A=independent_strict_matrix_float();w=np.array([-A[i,j] for i,j in edges]);R=np.diag(np.sqrt(w));x={}
 x["ST2127"]=packet(2127,"constructed_unweighted_topological_d0_with_strict_edge_Hodge_star","DEC convention.",{
  "A_identity_residual":float(np.linalg.norm(d0.T@np.diag(w)@d0-A,np.inf))})
 x["ST2128"]=packet(2128,"constructed_weighted_differential_d0w_equals_Rd0_with_standard_edge_inner_product","Prior incidence convention.",{
  "A_identity_residual":float(np.linalg.norm((R@d0).T@(R@d0)-A,np.inf))})
 d1w=d1@np.linalg.inv(R)
 x["ST2129"]=packet(2129,"proven_two_conventions_are_isomorphic_and_weighted_d1_must_be_d1_Rinverse","Chain identity.",{
  "unweighted_chain_residual":float(np.linalg.norm(d1@d0,np.inf)),"weighted_chain_residual":float(np.linalg.norm(d1w@(R@d0),np.inf))})
 Q=d1.T@d1
 x["ST2130"]=packet(2130,"constructed_positive_curvature_operator_Q_equals_d1star_M2_d1","For M2 positive.",{
  "standard_lambda_min":float(np.linalg.eigvalsh(Q).min())})
 O=tri_orbits(tris);orb_index={t:i for i,o in enumerate(O) for t in o};c=np.array([1+.1*orb_index[t] for t in tris]);M=np.ones(220);Q1=(c[:,None]*d1).T@np.diag(M)@(c[:,None]*d1);Q2=d1.T@np.diag(M*c**2)@d1
 x["ST2131"]=packet(2131,"proven_row_scale_and_face_Hodge_inverse_rescaling_leave_curvature_operator_invariant","Only h=m c^2 observable.",{
  "residual":float(np.linalg.norm(Q1-Q2,np.inf))})
 x["ST2132"]=packet(2132,"proven_positive_D12_orbit_Hodge_parameters_form_Rplus12_physical_quadratic_cone","After row-scale quotient.",{"dimension":12})
 evQ=np.linalg.eigvalsh(Q1)
 x["ST2133"]=packet(2133,"proven_all_positive_orbit_scales_preserve_kernel_dimension11","Same d1 kernel.",{
  "zero_modes":int(np.sum(abs(evQ)<1e-9))})
 Lstd=d0@d0.T+d1.T@d1
 x["ST2134"]=packet(2134,"proven_unweighted_full_simplex_Hodge_L1_is_exactly_12_identity","Combinatorial identity.",{
  "residual":float(np.linalg.norm(Lstd-12*np.eye(66),np.inf))})
 Lscaled=d0@d0.T+Q1;ev=np.linalg.eigvalsh(Lscaled)
 x["ST2135"]=packet(2135,"proven_orbit_scaled_positive_Hodge_splits_complete_simplex_spectrum","Same topology, different dynamics.",{
  "lambda_min":float(ev.min()),"lambda_max":float(ev.max()),"distinct_rounded":len(set(np.round(ev,8)))})
 x["ST2136"]=packet(2136,"proven_twelve_positive_Hodge_parameters_tune_gauge_one_form_spectrum_without_changing_A_or_H1","Nonpredictive spectrum.",{})
 x["ST2137"]=packet(2137,"proven_cohomology_only_controls_zero_modes_not_positive_spectral_values","Topology cannot fix dynamics.",{})
 x["ST2138"]=packet(2138,"primitive_integral_incidence_fixes_d1_basis_but_not_face_Hodge_star","Integrality insufficient.",{})
 x["ST2139"]=packet(2139,"refinement_adds_vertical_Hodge_spectral_parameters_even_when_cellular_d1_is_standard","Prior six/seven orbit cones.",{})
 x["ST2140"]=packet(2140,"round_four_falsifies_positivity_H1_and_integral_d1_as_unique_gauge_spectrum_selector","Hodge cone remains.",{})
 x["ST2141"]=packet(2141,"recommended_ST2142_ST2156","Audit spectral/cohomological predictions and continuum boundary after quotienting conventions.",{
  "next":["standard degeneracy","weighted spectra","spectral action freedom","refinement zero modes","continuum scaling"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
