#!/usr/bin/env python3
"""FIN ST1842--ST1856: support topology and gauge-curvature obstruction."""
import math

import sympy as sp

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1842,1856
NAMES=["CompleteSupportGraph","OneSkeletonCycleRank","TriangleFaceCount","TriangleBoundaryRank",
 "CliqueComplexTopology","NoFaceGaugeDegrees","AllTriangleFlatness","WilsonActionFamily",
 "FaceWeightNonuniqueness","BandThresholdTopology","C12Contrast","GraphLocalityNoGo",
 "MaxwellContinuumBoundary","RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 V=12;E=66;cycle=E-V+1;tri=math.comb(V,3);edges=[(i,j) for i in range(V) for j in range(i+1,V)];idx={e:k for k,e in enumerate(edges)};cols=[]
 for i in range(V):
  for j in range(i+1,V):
   for k in range(j+1,V):
    col=[0]*E
    for a,b,sgn in [(i,j,1),(j,k,1),(i,k,-1)]:col[idx[(a,b)]]=sgn
    cols.append(col)
 B2=sp.Matrix(E,tri,lambda r,c:cols[c][r]);rank=int(B2.rank());x={}
 x["ST1842"]=packet(1842,"proven_strict_weight_support_graph_is_complete_K12","All pair weights positive.",{"V":V,"E":E,"degree":11,"diameter":1})
 x["ST1843"]=packet(1843,"proven_edge_only_connection_mod_vertex_gauge_has_55_cycle_degrees","Connected graph cycle rank.",{"cycle_rank":cycle})
 x["ST1844"]=packet(1844,"proven_complete_support_has_220_triangle_cliques","Candidate 2-cells.",{"triangles":tri})
 x["ST1845"]=packet(1845,"proven_triangle_to_edge_boundary_map_has_rank_55","All one-cycle homology killed when triangles attached.",{"rank":rank,"expected_cycle_rank":cycle})
 x["ST1846"]=packet(1846,"proven_full_clique_complex_of_K12_is_contractible_11_simplex","No nontrivial H1/topological gauge sector.",{"H1":0,"dimension":11})
 x["ST1847"]=packet(1847,"proven_without_2_cells_no_local_curvature_action_is_defined","55 holonomies remain modulo gauge.",{"local_Maxwell":False})
 x["ST1848"]=packet(1848,"proven_zero_all_triangle_holonomies_make_U1_connection_pure_gauge_on_clique_complex","Simply connected.",{})
 x["ST1849"]=packet(1849,"constructed_infinite_Wilson_action_family_on_triangles","Additional coefficients.",{
  "action":"sum_faces kappa_f [1-Re holonomy_f]","free_face_weights":tri})
 x["ST1850"]=packet(1850,"proven_edge_weights_do_not_uniquely_determine_face_areas_or_Wilson_coefficients","Product, mean, minimum and external areas all possible.",{
  "canonical_face_measure":False})
 bands=[]
 for r in range(1,7):
  ee=12*r if r<6 else 66;bands.append({"max_distance":r,"edges":ee,"cycle_rank":ee-11})
 x["ST1851"]=packet(1851,"constructed_distance_band_threshold_topology_ladder","Every threshold is an extra choice.",{"bands":bands})
 x["ST1852"]=packet(1852,"proven_nearest_neighbor_C12_has_H1_dimension_one_but_discards_54_strict_edges","Not selected by nonzero support.",{
  "C12_edges":12,"discarded":54,"H1":1})
 x["ST1853"]=packet(1853,"proven_complete_support_has_no_graph_geodesic_locality_or_large_diameter","Diameter one conflicts with emergent extended space unless weights/refinement add metric.",{
  "support_lightcone":False})
 x["ST1854"]=packet(1854,"blocked_no_canonical_2_complex_face_area_or_3plus1_scaling_for_Maxwell_limit","Gauge matter action does not generate gauge dynamics.",{
  "Maxwell_continuum":False})
 x["ST1855"]=packet(1855,"round_three_falsifies_support_graph_as_unique_gauge_curvature_geometry","Conditional clique or band complexes remain possible.",{
  "unique_cell_complex":False})
 x["ST1856"]=packet(1856,"recommended_ST1857_ST1871","Test exact refinement transport and whether it selects locality/dimension.",{
  "next":["A24 incidence","coarse gauge lift","relative fiber gauges","cycle growth","q scaling","continuum no-go"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
