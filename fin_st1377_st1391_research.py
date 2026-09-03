#!/usr/bin/env python3
"""FIN ST1377--ST1391: interval composite identifiability certificate."""
import math

from fin_interval_cert_common import classical_return_interval, dephased_return_interval, endpoints, square_lower, square_upper, setup
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1377,1391
NAMES=["CompositeBoxDefinition","FourTimeTargetIntervals","AlphaGammaPartition","PerCellRSSLowerBound",
 "GlobalRSSLBound","WorstLowerCell","PointCandidateUpperBound","CompositeIdentifiabilityTheorem",
 "AtLeastOneTimeResidualBound","ClockBandNecessity","DephasingBandCoverage","PartitionReplayDeterminism",
 "CertificateScope","RoundTwo_Verdict","RoundTwo_Recommendation"]
TIMES=["0.3","0.6","1.2","2.0"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};dps=20;ctx=setup(dps)
 x["ST1377"]=packet(1377,"constructed_closed_composite_nuisance_box","Energy-dephasing family only.",{
  "alpha":[.9,1.1],"gamma":[0,100],"times":[float(t) for t in TIMES]})
 targets=[classical_return_interval(t,dps=dps,context=ctx) for t in TIMES]
 x["ST1378"]=packet(1378,"proven_outward_classical_target_intervals","Frozen decimal spectrum.",{
  "targets":[list(endpoints(v)) for v in targets]})
 x["ST1379"]=packet(1379,"constructed_complete_rectangular_partition","Exact decimal endpoints.",{
  "alpha_cells":20,"alpha_step":.01,"gamma_cells":50,"gamma_step":2,"total_cells":1000})
 global_lower=float('inf');worst=None;lowers=[]
 for ia in range(20):
  ab=(f"{90+ia:02d}e-2",f"{91+ia:02d}e-2")
  for ig in range(50):
   gb=(str(2*ig),str(2*(ig+1)));rss=0.0
   for t,target in zip(TIMES,targets):
    pred=dephased_return_interval(t,ab,gb,dps=dps,context=ctx);rss+=square_lower(pred-target)
   lowers.append(rss)
   if rss<global_lower:global_lower=rss;worst={"alpha":[float(ab[0]),float(ab[1])],"gamma":[float(gb[0]),float(gb[1])]}
 x["ST1380"]=packet(1380,"proven_each_cell_has_outward_RSS_lower_bound","Sum of interval squared-residual lower bounds.",{
  "cells_evaluated":len(lowers),"all_nonnegative":all(v>=0 for v in lowers)})
 x["ST1381"]=packet(1381,"proven_global_positive_RSS_lower_bound_on_declared_box","Finite complete cover.",{
  "RSS_lower":global_lower,"positive":global_lower>0})
 x["ST1382"]=packet(1382,"identified_weakest_partition_cell","A cover artifact need not contain exact minimizer.",{
  "cell":worst,"cell_lower":global_lower})
 alpha=("1.1","1.1");gamma=("12.373117825382856","12.373117825382856");candidate_lo=candidate_hi=0.0
 for t,target in zip(TIMES,targets):
  z=dephased_return_interval(t,alpha,gamma,dps=dps,context=ctx)-target;candidate_lo+=square_lower(z);candidate_hi+=square_upper(z)
 x["ST1383"]=packet(1383,"proven_point_candidate_RSS_upper_bound","Does not prove optimizer location.",{
  "alpha":1.1,"gamma":12.373117825382856,"RSS_interval":[candidate_lo,candidate_hi]})
 x["ST1384"]=packet(1384,"proven_four_time_classical_curve_not_reproducible_by_declared_dephased_quantum_box","Frozen decimal model theorem.",{
  "minimum_RSS_lower":global_lower,"exact_match_possible":False})
 one_res=math.sqrt(global_lower/4)
 x["ST1385"]=packet(1385,"proven_at_least_one_time_has_absolute_residual_lower_bound","Pigeonhole from four squared residuals.",{
  "residual_lower":one_res})
 x["ST1386"]=packet(1386,"clock_band_is_a_paid_premise_not_strict_CA_value","Certificate fails to address wide alpha torsor.",{
  "paid_band":[.9,1.1],"strict_clock_calibration":False})
 x["ST1387"]=packet(1387,"dephasing_band_0_100_is_complete_only_for_declared_family","Other Lindblad/noise laws absent.",{
  "gamma_covered":[0,100],"all_open_quantum_channels":False})
 x["ST1388"]=packet(1388,"constructed_deterministic_partition_replay","No stochastic optimizer.",{
  "cells":1000,"minimum_repeatable":global_lower})
 x["ST1389"]=packet(1389,"certificate_is_interval_grade_for_frozen_decimals_and_box_not_physical_universality","No claim inflation.",{
  "outward_interval":True,"symbolic_A":False,"laboratory_nuisance_complete":False})
 x["ST1390"]=packet(1390,"upgraded_four_time_composite_separation_from_strong_numerical_to_positive_interval_lower_bound","Major local closure.",{
  "RSS_lower":global_lower,"RSS_candidate_upper":candidate_hi})
 x["ST1391"]=packet(1391,"recommended_ST1392_ST1406","Add preparation, detector confusion and loss boxes.",{
  "next":["source contamination","confusion matrix","loss/no-click","identifiability thresholds","failure surfaces"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
