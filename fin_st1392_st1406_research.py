#!/usr/bin/env python3
"""FIN ST1392--ST1406: preparation, detector and loss robustness."""
import math

from fin_oa_discrimination_common import classical_return, quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1392,1406
NAMES=["BaseGap","CommonUniformPreparationContamination","IndependentPreparationTVBound","PreparationIdentifiabilityThreshold",
 "BinaryDetectorConfusionMap","KnownConfusionGap","UnknownDetectorErrorBound","CombinedAdversarialErrorBound",
 "TenPercentCommonContamination","KnownUniformLoss","UnknownLossNoClick","VertexDependentEfficiency",
 "RobustnessResponsibilityMatrix","RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};t=.6;c=classical_return(t);q=quantum_return(t);gap=q-c
 x["ST1392"]=packet(1392,"proven_frozen_base_return_gap","Float replay of exact model formulas.",{"time":t,"C":c,"Q":q,"gap":gap})
 x["ST1393"]=packet(1393,"proven_common_uniform_source_contamination_shrinks_gap_linearly","Uniform state is stationary for both declared channels.",{
  "formula":"p'_H=(1-epsilon)p_H+epsilon/12","gap_factor":"1-epsilon"})
 x["ST1394"]=packet(1394,"proven_independent_return_probability_errors_reduce_gap_by_at_most_two_epsilon","Triangle inequality.",{
  "assumption":"each model return prediction perturbed by at most epsilon","lower_gap":"D-2epsilon"})
 threshold=gap/2
 x["ST1395"]=packet(1395,"proven_adversarial_preparation_identifiability_threshold","Binary return only.",{
  "epsilon_must_be_below":threshold,"at_threshold_lower_gap":0.0})
 x["ST1396"]=packet(1396,"constructed_binary_detector_confusion_channel","Common calibrated rates.",{
  "formula":"p_obs=a+(1-a-b)p","a":"false positive","b":"false negative"})
 x["ST1397"]=packet(1397,"proven_known_common_detector_confusion_scales_gap","Invertible when a+b!=1.",{
  "gap_factor":"abs(1-a-b)","nonidentifying_surface":"a+b=1"})
 x["ST1398"]=packet(1398,"proven_independent_unknown_detector_bias_has_D_minus_2delta_bound","Worst-case scalar probability errors.",{
  "lower_gap":"D-2delta","delta_threshold":threshold})
 prep=.05;det=.05;lower=gap-2*(prep+det)
 x["ST1399"]=packet(1399,"constructed_combined_adversarial_error_budget","Simple additive worst-case budget.",{
  "preparation_error_each":prep,"detector_error_each":det,"certified_lower_gap":lower,"positive":lower>0})
 eps=.1
 x["ST1400"]=packet(1400,"proven_ten_percent_common_uniform_preparation_retains_large_gap","Common contamination model.",{
  "epsilon":eps,"gap_after":(1-eps)*gap})
 eta=.8
 x["ST1401"]=packet(1401,"proven_known_uniform_loss_preserves_conditional_vertex_distribution","No-click probability remains informative.",{
  "efficiency":eta,"click_conditioned_gap":gap,"attempt_level_click_probability":eta})
 x["ST1402"]=packet(1402,"proven_unknown_uniform_loss_is_unidentifiable_if_no_clicks_discarded","Postselection obstruction.",{
  "requires_attempt_count":True,"requires_no_click_record":True})
 x["ST1403"]=packet(1403,"constructed_vertex_dependent_efficiency_as_higher_risk_nuisance","Not covered by scalar loss model.",{
  "formula":"p_obs(x) proportional eta_x p(x)","requires_calibration_matrix":True})
 x["ST1404"]=packet(1404,"constructed_OA_calibration_responsibility_matrix","Separates assumptions from predictions.",{
  "source":"preparation contamination bound","detector":"confusion/efficiency matrix","clock":"alpha band","environment":"dephasing family","record":"attempts+no-clicks"})
 x["ST1405"]=packet(1405,"base_discriminator_has_explicit_positive_error_budgets_but_not_complete_hardware_robustness","No apparatus characterization.",{
  "common_uniform_models":True,"adversarial_scalar_bound":True,"full_detector_drift":False})
 x["ST1406"]=packet(1406,"recommended_ST1407_ST1421","Use the full vertex likelihood and symmetry-reduced sufficient statistic.",{
  "next":["seven distance classes","likelihood ratios","Chernoff bounds","sample design","coarse-return information loss"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
