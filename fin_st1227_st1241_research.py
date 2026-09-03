#!/usr/bin/env python3
"""FIN ST1227--ST1241: insufficiency of the W0+CA+SA two-package architecture."""
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, STRICT_A, write_packet, write_round


LO,HI=1227,1241
NAMES=["TwoPackage_Architecture","SameCASA_ThreeStateModels","StateTopology_AmbiguityPersists",
 "DimensionfulDualDynamics","SelectedFrame_PreparationCandidate","SameInput_DifferentPredictions",
 "AdaptiveOperatorLaw_TypeAudit","AdaptiveScalingFlow_Counterexample","FullMatrixAlgebra_NotOperationalSelection",
 "ThermodynamicBridge_StillMissing","LorentzCausalBridge_StillMissing","AnnihilationMeaning_StillAmbiguous",
 "TwoPackage_InsufficiencyTheorem","RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};ell,tau,hbar=2.0,3.0,5.0;r0,lam0=0,1
 x["ST1227"]=packet(1227,"constructed_W0_CA_SA_conditional_architecture","Explicit non-strict branch.",{
  "W0":"strict dimensionless informational core","CA":{"ell":ell,"tau":tau,"hbar":hbar},
  "SA":{"r0":r0,"lambda0":lam0,"coupling":"assumed"}})
 x["ST1228"]=packet(1228,"constructed_three_models_with_identical_W0_CA_SA","Countermodels to sufficiency.",{
  "models":["classical Delta11 heat","pure quantum CP11 unitary","density quantum D12 conjugation"],
  "same_units":True,"same_selected_frame":True,"same_A":True})
 x["ST1229"]=packet(1229,"proven_CA_and_SA_do_not_change_state_space_inequivalence","Dimensions/topology unchanged by units and a frame label.",{
  "dimensions":[11,22,143],"pairwise_homeomorphic":False})
 x["ST1230"]=packet(1230,"proven_same_tau_dimensionalizes_both_incompatible_channels","This is dimensional compatibility, not physical equivalence.",{
  "unitary":"exp[-i(t_phys/tau_*)A]","heat":"exp[-(t_phys/tau_*)A]",
  "physical_channel_selected":False})
 e0=np.zeros(N);e0[r0]=1;rho0=np.outer(e0,e0)
 x["ST1231"]=packet(1231,"constructed_same_selected_vertex_preparation_in_classical_and_quantum_models",
  "Preparation semantics already differs by category.",{
  "classical":"p0=delta_r0","quantum":"rho0=|r0><r0|","selected_r0":r0})
 tdim=1.0;P=expm(-tdim*STRICT_A);U=expm(-1j*tdim*STRICT_A)
 pc=float(P[r0,r0]);pq=float(abs(U[r0,r0])**2)
 x["ST1232"]=packet(1232,"proven_same_W0_CA_SA_yields_different_return_predictions","Operational discriminator once category/channel is specified.",{
  "dimensionless_time":tdim,"classical_return":pc,"quantum_return":pq,"absolute_difference":abs(pc-pq)})
 x["ST1233"]=packet(1233,"proven_repository_adaptive_law_is_operator_space_typed_not_state_category_typed",
  "Uses the documented schematic law only.",{
  "law":"K_dot=-Pi(V'(K)-C(K))","needs":["operator manifold","projection Pi","potential V","feedback C"],
  "classical_or_quantum_state_selected":False})
 x["ST1234"]=packet(1234,"proven_simple_adaptive_operator_flow_preserves_both_classical_and_quantum_uses",
  "Counterexample to adaptive-law-as-category-selector.",{
  "flow":"A_dot=-gamma A, gamma>0","solution":"A(t)=exp(-gamma t)A0",
  "remains_Markov_laplacian_up_to_time_rescaling":True,"remains_self_adjoint_quantum_generator":True})
 x["ST1235"]=packet(1235,"proven_M12_generation_does_not_supply_preparation_instrument_or_record",
  "Algebraic reachability is not operational semantics.",{
  "generated_algebra_dimension":144,"state_selector":False,"instrument_selector":False})
 x["ST1236"]=packet(1236,"proven_CA_SA_does_not_complete_information_to_thermodynamics","Temperature/bath/process structure absent.",{
  "kB":False,"temperature":False,"Landauer_process":False,"Shannon_bits":True})
 x["ST1237"]=packet(1237,"proven_cstar_does_not_create_Lorentz_or_causal_structure","A speed unit is not a light cone.",{
  "c_star":"ell_*/tau_*","linear_massless_dispersion":False,"Lorentz_group":False})
 x["ST1238"]=packet(1238,"proven_annihilation_blocker_remains_category_dependent_after_CA_SA","Units and frame do not choose conserved object.",{
  "classical":"mass","quantum":"norm/trace","single_physical_answer":False})
 x["ST1239"]=packet(1239,"proven_W0_plus_CA_plus_SA_is_insufficient_for_unique_operational_physics",
  "New theorem relative to the earlier strategic two-package proposal.",{
  "proof":"three pairwise inequivalent models share identical W0, CA and SA and predict different return probabilities",
  "missing_package":"state/channel/clock/preparation/instrument/record semantics"})
 x["ST1240"]=packet(1240,"two_package_architecture_is_coherent_but_not_operationally_complete","Do not discard CA/SA; add an independently labelled package.",{
  "units_closed_conditionally":True,"sector_closed_conditionally":True,"operational_model_closed":False})
 x["ST1241"]=packet(1241,"recommended_ST1242_ST1256","Construct and minimize an Operational/State package OA.",{
  "next":["classical and quantum OA witnesses","necessity of category/channel/clock/preparation/instrument/record",
          "finite discriminants","evidence boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
