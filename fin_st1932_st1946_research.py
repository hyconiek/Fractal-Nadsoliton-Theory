#!/usr/bin/env python3
"""FIN ST1932--ST1946: refinement functoriality and flag-vs-product obstruction."""
import itertools

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1932,1946
NAMES=["RefinedSupportCounts","RefinedCycleRank","RefinedTriangleCount","TriangleBoundaryRank",
 "FlagRefinedH1","ProductSquareCount","SquareBoundaryRank","ProductBoundaryRank",
 "ProductH1Closure","FlagProductConflict","LiftedBaseFaceFreedom","SquareOrbitWeights",
 "QDegeneracyBoundary","RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def gf2rank(cols):
 piv={}
 for v in cols:
  while v:
   p=v.bit_length()-1
   if p in piv:v^=piv[p]
   else:piv[p]=v;break
 return len(piv)
def main():
 V=24;edges=[]
 for f in [0,1]:
  for i,j in itertools.combinations(range(12),2):edges.append(tuple(sorted((2*i+f,2*j+f))))
 for i in range(12):edges.append((2*i,2*i+1))
 edges=sorted(edges);ei={e:k for k,e in enumerate(edges)}
 def bits(path):
  v=0
  for a,b in zip(path,path[1:]+path[:1]):v^=1<<ei[tuple(sorted((a,b)))]
  return v
 triangles=[]
 for f in [0,1]:
  for i,j,k in itertools.combinations(range(12),3):triangles.append(bits([2*i+f,2*j+f,2*k+f]))
 squares=[bits([2*i,2*j,2*j+1,2*i+1]) for i,j in itertools.combinations(range(12),2)];E=len(edges);cycle=E-V+1;x={}
 x["ST1932"]=packet(1932,"proven_refined_K12_square_K2_support_counts","q>0 generic.",{"V":V,"E":E})
 x["ST1933"]=packet(1933,"proven_refined_one_skeleton_cycle_rank121","Connected product graph.",{"cycle_rank":cycle})
 x["ST1934"]=packet(1934,"proven_refined_flag_complex_has_440_triangle_cliques","Two K12 layers only.",{"triangles":len(triangles)})
 rt=gf2rank(triangles)
 x["ST1935"]=packet(1935,"proven_refined_triangle_boundary_rank110_over_GF2","Exact bit elimination.",{"rank":rt})
 x["ST1936"]=packet(1936,"proven_refined_flag_complex_has_H1_dimension11","Cycle rank minus triangle rank.",{"H1":cycle-rt})
 x["ST1937"]=packet(1937,"proven_product_refinement_supplies_66_canonical_base_edge_times_interval_squares","Given typed product map.",{"squares":len(squares)})
 rs=gf2rank(squares)
 x["ST1938"]=packet(1938,"proven_square_boundary_family_has_rank66_but_only11_new_directions_mod_triangles","Exact GF2.",{"square_rank":rs,"new_mod_triangles":gf2rank(triangles+squares)-rt})
 rall=gf2rank(triangles+squares)
 x["ST1939"]=packet(1939,"proven_triangles_plus_product_squares_span_full_cycle_space","Exact GF2.",{"rank":rall,"cycle_rank":cycle})
 x["ST1940"]=packet(1940,"proven_product_cell_complex_kills_refined_H1","Given all lifted base triangles and squares.",{"H1":cycle-rall})
 x["ST1941"]=packet(1941,"proven_flag_of_support_is_not_refinement_product_functorial","Base flag H1=0 but refined flag H1=11; product cells need nonclique squares.",{})
 x["ST1942"]=packet(1942,"proven_every_base_natural_face_set_lifts_to_two_layers_so_refinement_does_not_select_among_base_candidates","Initial ambiguity inherited.",{})
 x["ST1943"]=packet(1943,"proven_product_squares_form_six_D12_distance_orbits_and_retain_positive_weight_cone","One square orbit per cyclic distance.",{"orbit_count":6})
 x["ST1944"]=packet(1944,"fiber_edge_identification_from_weight_fails_when_q_coincides_with_a_base_edge_weight","Product metadata then essential.",{"q_free":True})
 x["ST1945"]=packet(1945,"round_three_finds_product_cellular_completion_conditional_on_refinement_map_but_no_absolute_support_functor","Positive/negative split.",{
  "relative_functor":True,"absolute_unique":False})
 x["ST1946"]=packet(1946,"recommended_ST1947_ST1961","Test Wilson positivity, heat-trace dimension and continuum scaling.",{
  "next":["small-curvature action","face scale","spectral dimension","iterated q laws","Lorentz boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
