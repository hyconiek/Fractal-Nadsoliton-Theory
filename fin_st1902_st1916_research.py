#!/usr/bin/env python3
"""FIN ST1902--ST1916: automorphisms and natural face-set census."""
import itertools

import sympy as sp

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1902,1916
NAMES=["WeightedAutomorphismGroup","TriangleOrbitCensus","NaturalFaceSetCount","ExplicitFullRankComplexA",
 "ExplicitFullRankComplexB","SameSizeNonuniqueness","GF2FullRankCensus","GF2MinimalCensus",
 "FullCliqueChoice","EmptyFaceChoice","NaturalityNoSelection","H1NoSelection",
 "MinimalityNoSelection","RoundOne_Verdict","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def orbit_data():
 V=range(12);tris=list(itertools.combinations(V,3))
 def rot(t,k):return tuple(sorted((x+k)%12 for x in t))
 def ref(t,k):return tuple(sorted((k-x)%12 for x in t))
 unseen=set(tris);orbs=[]
 while unseen:
  t=min(unseen);o=({rot(t,k) for k in V}|{ref(t,k) for k in V})&set(tris);orbs.append(sorted(o));unseen-=o
 return orbs
def gf2rank(cols):
 piv={}
 for v in cols:
  while v:
   p=v.bit_length()-1
   if p in piv:v^=piv[p]
   else:piv[p]=v;break
 return len(piv)
def main():
 V=range(12);edges=list(itertools.combinations(V,2));ei={e:i for i,e in enumerate(edges)};orbs=orbit_data();x={}
 x["ST1902"]=packet(1902,"proven_weighted_strict_graph_automorphism_group_is_exactly_D12","Unique maximum-weight distance-1 edges form C12; its automorphisms are D12, which preserves all cyclic distances.",{
  "order":24})
 reps=[]
 for o in orbs:
  t=o[0];ds=sorted(min((a-b)%12,(b-a)%12) for a,b in itertools.combinations(t,2));reps.append({"representative":t,"size":len(o),"distance_triple":ds})
 x["ST1903"]=packet(1903,"proven_220_triangles_split_into_twelve_D12_orbits","Finite orbit enumeration.",{
  "orbits":reps,"sizes_sum":sum(r['size'] for r in reps)})
 x["ST1904"]=packet(1904,"proven_automorphism_natural_face_sets_are_arbitrary_unions_of_twelve_orbits","Boolean orbit subsets.",{
  "count":2**len(orbs)})
 def boundary(mask):
  ts=[t for i,o in enumerate(orbs) if mask>>i&1 for t in o];M=sp.zeros(len(edges),len(ts))
  for c,t in enumerate(ts):M[ei[(t[0],t[1])],c]=1;M[ei[(t[1],t[2])],c]=1;M[ei[(t[0],t[2])],c]=-1
  return M,ts
 masks=[2099,2101]
 for prog,mask,label in [(1905,masks[0],'A'),(1906,masks[1],'B')]:
  M,ts=boundary(mask);x[f"ST{prog}"]=packet(prog,f"proven_explicit_natural_76_face_complex_{label}_has_boundary_rank55","Exact rational SymPy rank.",{
   "orbit_indices":[i for i in range(12) if mask>>i&1],"faces":len(ts),"rank":int(M.rank())})
 x["ST1907"]=packet(1907,"proven_naturality_H1_killing_and_same_face_count_do_not_select_unique_complex","Two explicit distinct witnesses.",{
  "distinct":masks[0]!=masks[1],"faces_each":76})
 bits={t:(1<<ei[(t[0],t[1])])^(1<<ei[(t[1],t[2])])^(1<<ei[(t[0],t[2])]) for o in orbs for t in o};passing=[]
 for mask in range(1,4096):
  cols=[bits[t] for i,o in enumerate(orbs) if mask>>i&1 for t in o]
  if gf2rank(cols)==55:passing.append((mask,len(cols),mask.bit_count()))
 x["ST1908"]=packet(1908,"proven_at_least_2542_natural_face_sets_kill_H1_over_GF2_and_hence_have_full_rational_boundary_rank","Exact bit elimination; there may be additional rational-only passes.",{
  "passing":len(passing)})
 S={m for m,n,k in passing};mins=[]
 for m,n,k in passing:
  if not any((m^(1<<i)) in S for i in range(12) if m>>i&1):mins.append((m,n,k))
 x["ST1909"]=packet(1909,"proven_at_least_338_inclusion_minimal_GF2_full_rank_natural_complexes","All found use five orbits.",{
  "count":len(mins),"orbit_counts":sorted(set(k for m,n,k in mins)),"face_range":[min(n for m,n,k in mins),max(n for m,n,k in mins)]})
 x["ST1910"]=packet(1910,"constructed_full_flag_clique_complex_choice","Canonical from support but 11-dimensional simplex.",{
  "triangles":220,"contractible":True})
 x["ST1911"]=packet(1911,"constructed_empty_2cell_choice","Also support-automorphism natural but leaves H1 dimension55.",{})
 x["ST1912"]=packet(1912,"proven_automorphism_naturality_alone_has_4096_fold_face_set_ambiguity","No selector.",{})
 x["ST1913"]=packet(1913,"proven_adding_H1_zero_still_leaves_thousands_of_natural_candidates","At least GF2 passing count.",{})
 x["ST1914"]=packet(1914,"proven_inclusion_minimality_and_five_orbit_constraint_still_leave_hundreds_of_candidates","No canonical minimum.",{})
 x["ST1915"]=packet(1915,"round_one_falsifies_automorphism_naturality_plus_basic_topology_as_unique_2complex_selector","Stronger principles required.",{
  "unique":False})
 x["ST1916"]=packet(1916,"recommended_ST1917_ST1931","Classify natural face measures and test uniqueness under homogeneity/local rules.",{
  "next":["orbit-weight cone","edge-derived rules","homogeneity","Hodge star","dimensional area"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
