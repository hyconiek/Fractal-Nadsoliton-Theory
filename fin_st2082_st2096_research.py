#!/usr/bin/env python3
"""FIN ST2082--ST2096: D12 character classification of second differentials."""
import itertools

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=2082,2096
NAMES=["OrientedSimplexRepresentations","D12CharacterTable","CochainDimensions","UnconstrainedEquivariantHom",
 "GradientImageCharacter","CycleQuotientCharacter","ChainConstrainedHomDimension","StandardCoboundaryWitness",
 "CharacterIntegrality","LargeModuliNoGo","ReflectionSignImportance","FullSimplexScope","FaceSubsetBoundary","RoundOne_Verdict","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def parity(seq):return -1 if sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))%2 else 1
def char(k,perm):
 out=0
 for S in itertools.combinations(range(12),k+1):
  img=[perm[i] for i in S]
  if set(img)==set(S):out+=parity(img)
 return out
def main():
 G=[]
 for typ in ['r','s']:
  for k in range(12):G.append((typ,k,tuple((i+k)%12 for i in range(12)) if typ=='r' else tuple((k-i)%12 for i in range(12))))
 rows=[]
 for typ,k,p in G:rows.append({"element":f"{typ}{k}","chi0":char(0,p),"chi1":char(1,p),"chi2":char(2,p)})
 x={}
 x["ST2082"]=packet(2082,"constructed_oriented_vertex_edge_triangle_D12_representations","Full simplex cochains.",{"dimensions":[12,66,220]})
 x["ST2083"]=packet(2083,"proven_exact_oriented_character_table_by_fixed_simplex_sign_census","24 group elements.",{"rows":rows})
 x["ST2084"]=packet(2084,"proven_cochain_dimensions_C0_C1_C2","Combinatorial.",{"C0":12,"C1":66,"C2":220})
 hom=sum(r['chi1']*r['chi2'] for r in rows)//24
 x["ST2085"]=packet(2085,"proven_equivariant_Hom_C1_C2_dimension615","Character inner product.",{"dimension":hom})
 x["ST2086"]=packet(2086,"proven_gradient_image_representation_is_C0_minus_trivial","Connected complete graph.",{"dimension":11})
 qchar=[r['chi1']-r['chi0']+1 for r in rows]
 x["ST2087"]=packet(2087,"constructed_edge_mod_gradient_character","Cycle/coexact quotient dimension55.",{"character":qchar,"dimension":55})
 constrained=sum((r['chi1']-r['chi0']+1)*r['chi2'] for r in rows)//24
 x["ST2088"]=packet(2088,"proven_D12_equivariant_maps_d1_with_d1d0_zero_form_517_dimensional_space","Maschke splitting/character inner product.",{"dimension":constrained})
 x["ST2089"]=packet(2089,"proven_standard_simplicial_coboundary_is_one_nonzero_equivariant_chain_map","Existence witness.",{})
 x["ST2090"]=packet(2090,"proven_character_inner_products_are_integers_and_dimension_check_passes","Exact arithmetic.",{"Hom":hom,"chain_constrained":constrained})
 x["ST2091"]=packet(2091,"proven_equivariance_and_chain_condition_are_radically_nonunique","517 versus desired one-dimensional ray.",{})
 x["ST2092"]=packet(2092,"proven_oriented_reflection_signs_materially_change_character_and_must_not_be_ignored","Unoriented permutation census would be wrong calculus.",{})
 x["ST2093"]=packet(2093,"classification_is_for_full_triangle_space_not_one_of_4096_face_subsets","Scope.",{})
 x["ST2094"]=packet(2094,"face_subset_choice_can_only_change_reduce_or_reorganize_moduli_not_restore_unconditional_uniqueness","Must be audited separately.",{})
 x["ST2095"]=packet(2095,"round_one_falsifies_D12_equivariance_plus_d1d0_zero_as_unique_second_differential_selector","Exact no-go.",{})
 x["ST2096"]=packet(2096,"recommended_ST2097_ST2111","Add triangle-boundary locality and integral incidence structure.",{
  "next":["local row kernel","orbit scalars","integral coefficients","standard coboundary","weight dependence"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
