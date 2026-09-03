#!/usr/bin/env python3
"""FIN ST1122--ST1136: internal currents, entropy and dilation."""
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, STRICT_A, UNIFORM, write_packet, write_round


LO,HI=1122,1136
NAMES=["EdgeCurrent_Definition","ContinuityEquation_Certificate","ShannonEntropyProduction",
 "UniformDetailedBalance","StationaryCurrent_NoGo","HeatChannel_KrausConstruction",
 "HeatChannel_ChoiRank","KrausDilation_Nonuniqueness","FiniteEnvironment_SemigroupNoGo",
 "InfiniteOrResetEnvironment_Requirement","GlobalVsLocalInformation","ExternalPump_Rejection",
 "StrictA_GradientNotCirculation","Internalization_CurrentVerdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def entropy(p):
 p=np.asarray(p);q=p[p>0];return float(-np.sum(q*np.log(q)))
def main():
 x={};W=-STRICT_A.copy();np.fill_diagonal(W,0)
 rng=np.random.default_rng(1122);p=rng.random(N);p/=p.sum()
 J=W*(p[None,:]-p[:,None]) # J_ij: into i from j
 dp=-STRICT_A@p
 x["ST1122"]=packet(1122,"constructed_strict_edge_current","Classical Markov-state candidate required.",{
  "definition":"J_ij=W_ij(p_j-p_i)","antisymmetry_residual":float(np.linalg.norm(J+J.T,np.inf)),
  "edge_weight_min":float(W[W>0].min())})
 x["ST1123"]=packet(1123,"proven_vertex_continuity_equation","Finite exact formula, floating replay.",{
  "equation":"p_dot_i=sum_j J_ij=-(Ap)_i","residual":float(np.linalg.norm(J.sum(1)-dp,np.inf)),
  "global_rate":float(dp.sum())})
 dH=float(-dp@np.log(p)); edge=float(0.5*np.sum(W*(p[:,None]-p[None,:])*(np.log(p)[:,None]-np.log(p)[None,:])))
 x["ST1124"]=packet(1124,"proven_nonnegative_Shannon_entropy_production","Only for the classical symmetric heat candidate.",{
  "dH_dt":dH,"edge_formula":edge,"residual":abs(dH-edge),"nonnegative":bool(dH>=-1e-13)})
 P=expm(-STRICT_A)
 x["ST1125"]=packet(1125,"proven_uniform_detailed_balance","Uniform state follows from symmetric positive rates.",{
  "stationary_residual":float(np.linalg.norm(STRICT_A@UNIFORM,np.inf)),
  "detailed_balance_residual":float(np.max(abs(UNIFORM[:,None]*W-UNIFORM[None,:]*W.T)))})
 x["ST1126"]=packet(1126,"proven_no_nonzero_stationary_edge_circulation_for_strict_heat",
  "Irreducibility uses all positive offdiagonal W entries.",{
  "stationary_distribution":"uniform","stationary_current_norm":0.0,
  "sustained_internal_pump_from_A":False,"reason":"gradient detailed-balance flow has zero equilibrium currents"})
 colres=float(np.linalg.norm(P.sum(0)-1,np.inf)); positive=int(np.sum(P>0))
 x["ST1127"]=packet(1127,"constructed_exact_fixed_time_CPTP_heat_embedding",
  "Embedding discards input coherences; it is one channel choice.",{
  "kraus":"K_ij=sqrt(P_ij)|i><j|","P":"exp(-A) at t=1",
  "kraus_completeness_residual":colres,"positive_kraus_terms":positive,
  "action":"rho -> diag(P diag(rho))"})
 x["ST1128"]=packet(1128,"proven_declared_measure_prepare_channel_has_full_Choi_rank",
  "For t=1 every P_ij is strictly positive.",{
  "choi_rank":N*N,"explicit_environment_dimension":N*N,"dimension":N,
  "boundary":"rank is for this coherence-destroying embedding, not every extension of classical P"})
 x["ST1129"]=packet(1129,"proven_Stinespring_representation_is_nonunique",
  "Minimal dilations are unique only up to environment isometry, not as physical apparatus.",{
  "unitary_mixing_of_kraus_operators":True,"A_selects_environment_basis":False,
  "A_selects_preparation":False})
 x["ST1130"]=packet(1130,"proven_fixed_finite_closed_environment_cannot_realize_mixing_heat_semigroup_for_all_times",
  "Assumes one finite environment, fixed initial environment state and time-independent total Hamiltonian.",{
  "proof":"finite closed unitary matrix elements are almost periodic; nonstationary strict heat trajectories converge to uniform and are not almost periodic",
  "exact_all_time_finite_dilation":False,"single_time_dilation":True})
 x["ST1131"]=packet(1131,"conditional_infinite_or_reset_internal_environment_needed_for_exact_irreversible_semigroup",
  "This is an added operational structure, not an external ontological layer.",{
  "options":["infinite internal bath","collision model with fresh internal ancillas","singular/continuum limit"],
  "strict_source":False})
 x["ST1132"]=packet(1132,"proven_local_entropy_growth_compatible_with_global_unitarity",
  "Requires a system-environment partition.",{
  "classical_heat_entropy_t0":entropy(p),"classical_heat_entropy_t1":entropy(P@p),
  "global_information_loss_required":False})
 x["ST1133"]=packet(1133,"ontology_rejects_external_resource_but_allows_internal_partition",
  "No actual FIN partition or bath is derived.",{
  "external_pump_fundamental":False,"internal_environment_conditional":True,
  "global_closed_state_required":True})
 x["ST1134"]=packet(1134,"proven_strict_A_supplies_relaxational_gradient_flow_not_sustained_circulation",
  "A different antisymmetric/non-detailed-balance generator would be extra.",{
  "equilibrium_current":0.0,"entropy_production_until_equilibrium":True,
  "self_sustaining_pump_from_heat":False})
 x["ST1135"]=packet(1135,"internal_heat_interpretation_requires_noncanonical_environment_structure",
  "The mathematical heat semigroup itself remains exact.",{
  "strict_heat_valid":True,"closed_finite_microdynamics_from_A_alone":False,
  "state_ontology_selection":False})
 x["ST1136"]=packet(1136,"recommended_ST1137_ST1151","Next round tests refinement/fiber nonuniqueness and dual compatibility.",{
  "next":["12-to-24 exact lifts","coarse intertwining","free fiber rate","hidden state entropy",
          "unitary/heat refinement compatibility","scale-source verdict"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
