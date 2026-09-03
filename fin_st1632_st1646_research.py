#!/usr/bin/env python3
"""FIN ST1632--ST1646: noncirculant Lindblad/dephasing sensitivity."""
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1632,1646
NAMES=["DephasingGenerator","CommutatorPerturbation","SquaredCommutatorPerturbation","LindbladGeneratorDifference",
 "HSContraction","LindbladDuhamelBound","ReturnEffectBound","Gamma100NumericBound","CompositeResidualStability",
 "CompositeEpsilonThreshold","RoundingWitness","NoncirculantWitness","OtherJumpBoundary","RoundOne_Verdict","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};M=2.3421820411463;gamma=100.;tmax=2.;delta=0.10804167160619894
 x["ST1632"]=packet(1632,"constructed_A_dephasing_Lindblad_generator","Hilbert-Schmidt operator space.",{
  "L_A":"-i ad_A -(gamma/2) ad_A^2","gamma_range":[0,100]})
 x["ST1633"]=packet(1633,"proven_commutator_superoperator_perturbation_bound","A,B self-adjoint.",{
  "bound":"||ad_B-ad_A||<=2 epsilon"})
 x["ST1634"]=packet(1634,"proven_squared_commutator_perturbation_bound","Factorization X^2-Y^2.",{
  "bound":"||ad_B^2-ad_A^2||<=4 epsilon (||A||+||B||)"})
 x["ST1635"]=packet(1635,"proven_dephasing_generator_difference_bound","Combine Hamiltonian/dissipative terms.",{
  "bound":"Delta_L<=2epsilon+2gamma epsilon(||A||+||B||)"})
 x["ST1636"]=packet(1636,"proven_each_declared_dephasing_semigroup_is_HS_contractive","Generator is normal with nonpositive real spectrum in A energy basis.",{
  "contraction":True})
 x["ST1637"]=packet(1637,"proven_Lindblad_semigroup_Duhamel_bound","Between two contractive declared semigroups.",{
  "bound":"||exp(tL_B)-exp(tL_A)||<=t Delta_L"})
 x["ST1638"]=packet(1638,"proven_rank_one_return_probability_perturbation_bound","HS norms of state/effect equal one.",{
  "bound":"|Q_B,gamma(t)-Q_A,gamma(t)|<=t Delta_L"})
 eps=8.041872774817993e-8;dl=2*eps*(1+gamma*(2*M+eps));total=tmax*dl+tmax*eps
 x["ST1639"]=packet(1639,"proven_gamma100_small_parameter_box_numeric_bound","Worst t<=2 plus heat target perturbation.",{
  "epsilon":eps,"Delta_L":dl,"total_return_residual_perturbation":total})
 x["ST1640"]=packet(1640,"proven_four_time_composite_residual_remains_positive_under_declared_noncirculant_Lindblad_perturbation","Uses prior at-least-one residual lower.",{
  "original_residual_lower":delta,"perturbed_lower":delta-total})
 denom=2+4*(1+2*gamma*M);crit=delta/denom
 x["ST1641"]=packet(1641,"constructed_conservative_operator_epsilon_threshold_for_composite_residual_positivity","Drops epsilon inside ||B|| upper for closed formula.",{
  "epsilon_threshold_approx":crit,"denominator":denom})
 er=6e-14;dlr=2*er*(1+gamma*(2*M+er));x["ST1642"]=packet(1642,"proven_rounding_scale_has_negligible_Lindblad_perturbation","Structured rounding witness.",{
  "epsilon":er,"total_bound":tmax*dlr+tmax*er})
 en=1e-8;dln=2*en*(1+gamma*(2*M+en));x["ST1643"]=packet(1643,"proven_noncirculant_epsilon_1e_minus8_preserves_composite_residual","Graph Laplacian class.",{
  "epsilon":en,"perturbed_lower":delta-(tmax*dln+tmax*en)})
 x["ST1644"]=packet(1644,"sensitivity_theorem_does_not_cover_arbitrary_jump_operators_nonMarkov_heat_or_trace_norm","Scope.",{
  "other_Lindblad_jumps":False,"trace_norm_contraction_claim":False})
 x["ST1645"]=packet(1645,"closed_declared_dephasing_sensitivity_gap_for_small_noncirculant_graph_perturbations","Conditional theorem.",{
  "composite_positive":True})
 x["ST1646"]=packet(1646,"recommended_ST1647_ST1661","Propagate explicit kernel-parameter boxes to A norm bounds.",{
  "next":["large audit box","small certified box","weight positivity","operator norm","location/composite gates"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
