#!/usr/bin/env python3
"""FIN ST1197--ST1211: conversion-axiom package CA."""
import itertools
import numpy as np

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1197,1211
NAMES=["CA_Definition","CA_RankThree_Minimality","PhysicalHamiltonian_UnitLift","PhysicalHeatTime_UnitLift",
 "DerivedSpeedUnit","DerivedEnergyMomentumMassUnits","ActionLift","CA_PositiveScaleTorsor","DimensionlessPredictions_Invariant",
 "Thermodynamic_kB_Gap","CA_SelectorBlindness","CA_StateCategoryBlindness","CA_AnnihilationNeutrality",
 "RoundTwo_Verdict","RoundTwo_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={}
 # columns are ell, tau, hbar in (H,L,T) exponents
 basis=np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=float) # rows H,L,T; columns ell,tau,hbar after reordering below
 # Explicit vectors: ell=(0,1,0), tau=(0,0,1), hbar=(1,0,0)
 vecs={"ell":[0,1,0],"tau":[0,0,1],"hbar":[1,0,0]}
 x["ST1197"]=packet(1197,"constructed_CA_positive_scale_triple","Explicit non-strict package.",{
  "CA":["ell_*>0","tau_*>0","hbar_*>0"],"dimensions":{"ell":"L","tau":"T","hbar":"H"}})
 rows=[]
 for r in range(4):
  for sub in itertools.combinations(vecs,r):
   if sub:
    M=np.array([vecs[k] for k in sub],float).T
    rank=int(np.linalg.matrix_rank(M))
   else:
    rank=0
   rows.append((sub,rank))
 x["ST1198"]=packet(1198,"proven_full_CA_triple_is_unique_rank_three_subset","Replays and independently checks P3146 dimensional algebra.",{
  "subsets":len(rows),"rank3_subsets":[list(s) for s,r in rows if r==3],"proper_subset_max_rank":max(r for s,r in rows if len(s)<3)})
 x["ST1199"]=packet(1199,"constructed_dimensionful_quantum_generator","State/channel choice remains conditional.",{
  "formula":"H_phys=(hbar_*/tau_*) A","unit":"energy","unitary":"exp[-i H_phys t_phys/hbar_*]=exp[-i(t_phys/tau_*)A]"})
 x["ST1200"]=packet(1200,"constructed_dimensionful_heat_generator","Classical/open channel choice remains conditional.",{
  "formula":"L_heat=-(1/tau_*)A","semigroup":"exp[-(t_phys/tau_*)A]","unit":"inverse time"})
 x["ST1201"]=packet(1201,"constructed_derived_speed_scale","Not a light-cone theorem.",{
  "formula":"c_*=ell_*/tau_*","unit":"L/T","equals_observed_c":False})
 x["ST1202"]=packet(1202,"proven_CA_derives_consistent_mechanical_unit_family","Pure dimensional consequences.",{
  "energy":"E_*=hbar_*/tau_*","momentum":"p_*=hbar_*/ell_*","mass":"m_*=hbar_* tau_*/ell_*^2",
  "identity":"E_*=m_* c_*^2"})
 x["ST1203"]=packet(1203,"constructed_action_unit_lift","Does not supply a unique Lagrangian or EOM.",{
  "formula":"S_phys=hbar_* S_dimensionless","path_weight":"exp(i S_phys/hbar_*)=exp(i S_dimensionless)"})
 x["ST1204"]=packet(1204,"proven_CA_values_form_free_positive_three_torsor","W0 contains no accepted scale-charged selector.",{
  "group":"(R_+)^3 acts freely by independent rescaling of ell,tau,hbar","canonical_numeric_triple":False})
 x["ST1205"]=packet(1205,"proven_dimensionless_W0_predictions_are_CA_rescaling_invariant","Units label physical magnitudes without changing A's ratios.",{
  "invariants":["eigenvalue ratios","dimensionless phase","Markov transition at fixed t/tau","entropy in nats"],
  "absolute_SI_values":False})
 x["ST1206"]=packet(1206,"proven_CA_does_not_convert_Shannon_to_thermodynamic_entropy_without_kB_or_temperature_structure",
  "hbar, ell and tau are insufficient for temperature dimension.",{
  "Shannon_H":"dimensionless","physical_entropy":"S=k_B H","k_B_in_CA":False,"Landauer_bridge_complete":False})
 x["ST1207"]=packet(1207,"proven_CA_is_selector_blind","Positive scales are unchanged under Z12 origin/polarity action.",{
  "selects_r0":False,"selects_lambda0":False,"QW2191_discharged":False})
 x["ST1208"]=packet(1208,"proven_CA_is_state_category_blind","All Release 10.96 candidates accept the same unit triple.",{
  "classical":True,"pure_quantum":True,"density_quantum":True,"selects_one":False})
 x["ST1209"]=packet(1209,"proven_CA_only_dimensionalizes_existing_blockers","It does not create conservation or completeness.",{
  "unitary_blocker":"norm with dimensional time","heat_blocker":"mass with dimensional time",
  "new_ontological_persistence":False})
 x["ST1210"]=packet(1210,"CA_is_minimal_for_HLT_units_but_insufficient_for_physics","Conditional success and strict-source failure separated.",{
  "dimension_rank":3,"strict_source":False,"selector":False,"state_category":False,"instrument":False})
 x["ST1211"]=packet(1211,"recommended_ST1212_ST1226","Audit the sector-axiom package as a finite frame torsor.",{
  "next":["D12 action on 24 branches","origin/polarity minimality","selector potential",
          "invariant-measure no-go","units/category independence"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
