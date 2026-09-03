#!/usr/bin/env python3
"""FIN ST1542--ST1556: noncirculant graph-Laplacian perturbation bounds."""
from fin_oa_discrimination_common import classical_return, quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1542,1556
NAMES=["PerturbationClass","HeatDuhamelBound","UnitaryDuhamelBound","ClassicalReturnBound","QuantumReturnBound",
 "GapPerturbationBound","BaseTimeRobustGap","GlobalLocationMarginFormula","OperatorNormThreshold","EntrywiseThreshold",
 "OneEminus8Witness","RoundingScaleWitness","UniquenessBoundary","CompositeDephasingBoundary","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};tc=.5945308;T=10.;margin=3.089572412462438e-6;base=quantum_return(.6)-classical_return(.6)
 x["ST1542"]=packet(1542,"constructed_noncirculant_graph_laplacian_operator_norm_class","Conditional class.",{
  "A_B":"self-adjoint positive semidefinite graph Laplacians","row_sum_zero":True,"B_minus_A_norm_le_epsilon":True,"circulant_required":False})
 x["ST1543"]=packet(1543,"proven_heat_semigroup_Duhamel_contraction_bound","Requires A,B positive semidefinite.",{
  "bound":"||exp(-tB)-exp(-tA)||<=t epsilon"})
 x["ST1544"]=packet(1544,"proven_unitary_group_Duhamel_bound","Self-adjoint A,B.",{
  "bound":"||exp(-itB)-exp(-itA)||<=t epsilon"})
 x["ST1545"]=packet(1545,"proven_classical_return_probability_perturbation_bound","Same vertex effect.",{
  "bound":"|C_B(t)-C_A(t)|<=t epsilon"})
 x["ST1546"]=packet(1546,"proven_quantum_return_probability_perturbation_bound","Uses |u|+|v|<=2.",{
  "bound":"|Q_B(t)-Q_A(t)|<=2t epsilon"})
 x["ST1547"]=packet(1547,"proven_return_gap_perturbation_bound","Triangle inequality.",{
  "bound":"|D_B(t)-D_A(t)|<=3t epsilon"})
 eps=1e-8
 x["ST1548"]=packet(1548,"proven_base_t06_gap_stays_positive_under_operator_perturbation","Noncirculant class.",{
  "epsilon":eps,"lower_gap":base-1.8*eps})
 x["ST1549"]=packet(1549,"proven_global_target_location_margin_formula","Compares candidate time with every outside time <=10.",{
  "lower_margin":"m-3(tc+10)epsilon","m":margin,"tc":tc})
 crit=margin/(3*(tc+T))
 x["ST1550"]=packet(1550,"proven_operator_norm_threshold_for_all_global_maximizers_to_remain_in_target_interval","Does not prove uniqueness inside.",{
  "epsilon_strict_upper":crit,"target":[.59,.60]})
 x["ST1551"]=packet(1551,"constructed_entrywise_sufficient_threshold_via_norm_le_12delta","Arbitrary symmetric entrywise perturbation bound.",{
  "delta_strict_upper":crit/12})
 x["ST1552"]=packet(1552,"proven_epsilon_1e_minus8_retains_large_exclusion_margin","Explicit witness.",{
  "epsilon":eps,"margin":margin-3*(tc+T)*eps})
 er=6e-14
 x["ST1553"]=packet(1553,"proven_previous_rounding_scale_is_far_inside_noncirculant_location_threshold","Operator norm witness.",{
  "epsilon":er,"margin":margin-3*(tc+T)*er,"threshold_ratio":crit/er})
 x["ST1554"]=packet(1554,"noncirculant_bound_localizes_global_maxima_but_does_not_preserve_unique_stationary_root","Derivative/eigenvector perturbation not certified.",{
  "global_location":True,"unique_maximum":False})
 x["ST1555"]=packet(1555,"Duhamel_lift_does_not_cover_dephasing_generator_sensitivity_or_nonPSD_heat","Scope boundary.",{
  "composite_RSS_lift":False,"arbitrary_open_quantum":False})
 x["ST1556"]=packet(1556,"recommended_ST1557_ST1571","Export the newly obtained continuous minimax value certificate under new numbering.",{
  "next":["adaptive interval KL cover","dual upper bound","value bracket","exact-weight boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
