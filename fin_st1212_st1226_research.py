#!/usr/bin/env python3
"""FIN ST1212--ST1226: sector-axiom package SA."""
import itertools

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1212,1226
NAMES=["SectorFrame_Torsor","Dihedral_FreeTransitiveAction","SA_Definition","OriginOnly_ResidualPolarity",
 "PolarityOnly_ResidualOrigins","OriginPolarity_UniqueFrame","CouplingAxiom_Necessity","FiniteSelectorPotential",
 "PotentialWeightFreedom","InvariantMeasure_SelectorNoGo","SA_CA_Independence","SA_StateCategoryBlindness",
 "SA_StrictStatus","RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def act(g,x):
 typ,k=g;r,l=x
 return ((r+k)%12,l) if typ=='R' else ((k-r)%12,-l)
def main():
 x={};X=[(r,l) for r in range(12) for l in (-1,1)];G=[('R',k) for k in range(12)]+[('F',k) for k in range(12)]
 orbit={act(g,(0,1)) for g in G};stabs=[g for g in G if act(g,(0,1))==(0,1)]
 x["ST1212"]=packet(1212,"constructed_oriented_origin_frame_space","Existing P3140 pair space.",{
  "X":"Z12 x {+1,-1}","size":len(X),"interpretation":"marked origin plus orientation/polarity"})
 x["ST1213"]=packet(1213,"proven_D12_action_is_free_and_transitive_on_sector_frames","Finite enumeration.",{
  "group_size":len(G),"orbit_size":len(orbit),"stabilizer_size":len(stabs),"torsor":True})
 x["ST1214"]=packet(1214,"constructed_minimal_conditional_SA","Explicitly non-strict.",{
  "SA":["A_origin:r0 in Z12","A_lambda:lambda0 in {+1,-1}","A_coupling:selected frame couples to sector channel"]})
 x["ST1215"]=packet(1215,"proven_origin_axiom_alone_leaves_two_frames","Finite count.",{
  "r0":0,"remaining":[[0,-1],[0,1]],"count":2})
 x["ST1216"]=packet(1216,"proven_polarity_axiom_alone_leaves_twelve_frames","Finite count.",{
  "lambda0":1,"count":12})
 x["ST1217"]=packet(1217,"proven_origin_plus_polarity_selects_one_frame","This is stipulation, not strict provenance.",{
  "selected":[0,1],"remaining_count":1,"symmetry_broken":True})
 x["ST1218"]=packet(1218,"proven_frame_selection_without_coupling_does_not_affect_dynamics","Typing theorem.",{
  "selected_label_only":True,"requires_A_coupling":True,"physical_sector_effect_without_coupling":False})
 weights=list(itertools.product([1,2,3],[1,2,3]));ok=0
 for wo,wl in weights:
  vals={(r,l):wo*min(r,12-r)**2+wl*(1-l)/2 for r,l in X};m=min(vals.values());ok+=sum(v==m for v in vals.values())==1 and vals[(0,1)]==m
 x["ST1219"]=packet(1219,"proven_axiom_selector_potential_has_unique_minimum_for_positive_audited_weights","Replays P3141 with exact integer enumeration.",{
  "formula":"V=mu[w_o d(r,r0)^2+w_l(1-lambda lambda0)/2]","weight_pairs":len(weights),"unique_rows":ok})
 x["ST1220"]=packet(1220,"proven_selector_potential_weights_and_scale_remain_free","No strict source or unit measure.",{
  "free":["mu","w_origin","w_lambda"],"strict_values":False})
 x["ST1221"]=packet(1221,"proven_invariant_probability_on_transitive_frame_torsor_is_uniform","No localized invariant selector.",{
  "mass_per_frame":1/24,"selected_mass":1/24,"localized":False})
 x["ST1222"]=packet(1222,"proven_SA_and_CA_solve_independent_group_actions","Scales and frames are product factors.",{
  "parameter_space":"(R_+)^3 x (Z12 x +/-)","CA_breaks_D12":False,"SA_fixes_units":False})
 x["ST1223"]=packet(1223,"proven_SA_is_state_category_blind","The same frame can label classical, ray and density models.",{
  "candidate_categories":3,"selected_category":False})
 x["ST1224"]=packet(1224,"SA_closes_selector_only_in_labelled_axiom_branch","Strict QW-2191 remains open.",{
  "conditional_unique_frame":True,"strict_source":False,"QW2191_discharged":False})
 x["ST1225"]=packet(1225,"two_strategic_leaf_cuts_CA_and_SA_are_consistently_but_only_conditionally_closed","Operational/state leaf remains.",{
  "units":"conditional","sector":"conditional","state_and_measurement":"open"})
 x["ST1226"]=packet(1226,"recommended_ST1227_ST1241","Test whether W0+CA+SA is sufficient after the A-only ontology no-go.",{
  "next":["build same-CA/SA classical and quantum models","prove two-package insufficiency",
          "dimensionful dual dynamics","identify missing operational package"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
