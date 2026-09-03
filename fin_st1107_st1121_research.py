#!/usr/bin/env python3
"""FIN ST1107--ST1121: completeness and no-killing classification."""
import math
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, ONES, STRICT_A, write_packet, write_round


LO,HI=1107,1121
NAMES=[
 "FiniteLinear_CompleteFlows","MarkovGenerator_IffCriterion","StrictHeat_ConservativeCertificate",
 "UniformKilling_Counterfamily","HamiltonianDensity_CompleteFlow","A_Dephasing_LindbladFamily",
 "RelaxationJump_SelectorObstruction","TracePreservation_HeisenbergCriterion","InfiniteChannelFamily_Underdetermination",
 "NoKilling_LinearClassification","NonlinearBoundary_NotControlled","CompactState_NotEnoughForVectorField",
 "StrictNoKilling_CurrentScope","CompletenessBridge_Verdict","RoundTwo_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={}
 x["ST1107"]=packet(1107,"proven_finite_matrix_exponentials_are_global","Only the selected linear laws.",{
  "unitary":"exp(-itA) exists for all real t","heat":"exp(-tA) exists for all t>=0","finite_time_blowup":False})
 x["ST1108"]=packet(1108,"proven_finite_CTMC_generator_criterion_paid","Standard finite-state theorem applied to Q=-A.",{
  "criterion":"Q_ij>=0 for i!=j and Q1=0","offdiag_min":float((-STRICT_A)[~np.eye(N,dtype=bool)].min()),
  "row_residual":float(np.linalg.norm((-STRICT_A)@ONES,np.inf))})
 x["ST1109"]=packet(1109,"proven_strict_heat_is_positive_and_mass_preserving_for_all_t",
  "Follows from the Metzler-generator theorem, stronger than sampling times.",{
  "semigroup":"P_t=exp(-tA)","entrywise_nonnegative":True,"stochastic":True,"killing_rate":0.0})
 rates=[0.01,0.1,1.0]
 x["ST1110"]=packet(1110,"proven_one_parameter_killing_extensions_exist",
  "Each gamma is additional physics/mathematics.",{
  "family":"Q_gamma=-A-gamma I","gammas":rates,"survival_t1":[math.exp(-g) for g in rates],
  "A_selects_gamma":False})
 x["ST1111"]=packet(1111,"proven_hamiltonian_density_flow_is_complete_CPTP_group",
  "Quantum density-state category supplied.",{
  "generator":"L_H(rho)=-i[A,rho]","trace_preserving":True,"positive":True,"invertible":True})
 x["ST1112"]=packet(1112,"constructed_A_sourced_dephasing_family",
  "Rate gamma remains free; this is not selected by A alone.",{
  "generator":"L_gamma(rho)=-i[A,rho]+gamma(A rho A-1/2{A^2,rho})",
  "gamma_domain":"gamma>=0","trace_preserving":True,"stationary_commutant":True})
 x["ST1113"]=packet(1113,"proven_canonical_relaxation_jump_requires_extra_selector_inside_degenerate_shells",
  "A has five two-dimensional eigenspaces.",{
  "spectral_multiplicities":[1,2,2,2,2,2,1],
  "issue":"mapping each excited shell to the ground ray needs a covector/basis or summed noise convention",
  "A_unique_relaxation_law":False,"QW2191_discharged":False})
 x["ST1114"]=packet(1114,"proven_trace_preservation_equals_Lstar_identity_zero",
  "Finite quantum dynamical semigroup.",{
  "criterion":"L*(I)=0","killing_or_leakage":"L*(I)<0 in a subunital survival model",
  "A_commutator_pays_criterion":True})
 x["ST1115"]=packet(1115,"proven_infinitely_many_complete_channels_share_A",
  "Hamiltonian, dephasing rates, Markov heat and killed laws are inequivalent.",{
  "families":["unitary","classical heat","A-dephasing gamma>=0","uniform killing gamma>0"],
  "unique_channel":False})
 x["ST1116"]=packet(1116,"proven_linear_no_killing_is_a_generator_constraint",
  "Not an ontological consequence of self-reference.",{
  "classical":"Q1=0","quantum":"L*(I)=0","survival_extension":"adjoin cemetery state for missing mass/trace"})
 x["ST1117"]=packet(1117,"blocked_no_full_nonlinear_completeness_from_A",
  "A controls only the linear part; nonlinear growth and boundary behavior require a specified law.",{
  "full_vector_field":False,"global_Lipschitz_or_coercive_bound":False,"finite_time_collapse_excluded":False})
 x["ST1118"]=packet(1118,"proven_compact_state_space_needs_tangent_invariant_vector_field",
  "Compactness alone does not stop a proposed ambient vector field from leaving the set.",{
  "required":["tangency/inward condition","regularity","global solution"],"supplied_by_A_for_unknown_nonlinearity":False})
 x["ST1119"]=packet(1119,"strict_no_killing_proven_only_for_declared_linear_channels",
  "No universal theorem over future nonlinear/open extensions.",{
  "closed":["unitary norm","classical heat mass","density unitary trace","A-dephasing trace"],
  "open":["full nadsoliton nonlinearity","state-dependent environment","boundary sectors"]})
 x["ST1120"]=packet(1120,"kernel_to_complete_full_flow_bridge_open_and_nonunique",
  "The linear channels are complete but mutually incompatible as state ontologies.",{
  "linear_completeness":True,"full_FIN_completeness":False,"selection_from_A":False})
 x["ST1121"]=packet(1121,"recommended_ST1122_ST1136","Next round tests whether heat can be internalized without an external pump.",{
  "next":["edge-current continuity","entropy production","stationary-current no-go","Kraus dilation",
          "finite-environment recurrence obstruction","internal resource verdict"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
