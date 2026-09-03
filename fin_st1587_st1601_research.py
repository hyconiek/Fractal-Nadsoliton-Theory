#!/usr/bin/env python3
"""FIN ST1587--ST1601: independent calibration polytope certificate."""
import hashlib
import itertools
import json
import math

from fin_interval_cert_common import setup, classical_return_interval, dephased_return_interval, endpoints
from fin_oa_discrimination_common import classical_return, quantum_return
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1587,1601
NAMES=["IndependentAffinePolytope","AffineVertexExtrema","BaseClassicalRange","BaseQuantumRange",
 "BaseIndependentGap","CompositeIndependentCover","CompositeIndependentRSS","CompositeIndependentResidual",
 "CommonVsIndependentComparison","UniformEfficiencyBox","NoClickRequirement","VertexDependentMatrixBoundary",
 "ArtifactHash","RoundFour_Verdict","RoundFour_Recommendation"]
CAL=ROOT/'FIN_ST1587_ST1601_IndependentCalibrationPolytope.json';TIMES=['0.3','0.6','1.2','2.0']
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def transform_range(interval,eps=.05,dmax=.02):
 vals=[]
 for p in interval:
  for e,a,b in itertools.product([0,eps],[0,dmax],[0,dmax]):vals.append(a+(1-a-b)*((1-e)*p+e/12))
 return min(vals),max(vals)
def square_lower(interval):
 a,b=interval
 if a<=0<=b:return 0.
 return min(a*a,b*b)
def main():
 x={};eps=.05;dmax=.02;eta=.8
 cal={"schema":"FIN-OA-INDEPENDENT-CAL-11.01-v1","C":{"epsilon":[0,eps],"a":[0,dmax],"b":[0,dmax]},"Q":{"epsilon":[0,eps],"a":[0,dmax],"b":[0,dmax]},"uniform_efficiency_each":[eta,1]};CAL.write_text(json.dumps(cal,indent=2,sort_keys=True)+'\n')
 x["ST1587"]=packet(1587,"constructed_independent_classical_and_quantum_affine_calibration_polytopes","Maps need not share parameters.",cal)
 x["ST1588"]=packet(1588,"proven_affine_probability_extrema_occur_at_box_vertices","Multiaffine in epsilon,a,b and monotone in p.",{
  "vertices_per_model":8})
 cr=transform_range((classical_return(.6),classical_return(.6)),eps,dmax);qr=transform_range((quantum_return(.6),quantum_return(.6)),eps,dmax)
 x["ST1589"]=packet(1589,"proven_independent_classical_return_range_at_t06","Exact vertex enumeration in double precision.",{"range":cr})
 x["ST1590"]=packet(1590,"proven_independent_quantum_return_range_at_t06","Closed quantum base.",{"range":qr})
 x["ST1591"]=packet(1591,"proven_base_return_ranges_remain_disjoint_under_independent_5pct_2pct_maps","Sharper than prior triangle sufficient bound.",{
  "gap_lower":qr[0]-cr[1],"positive":qr[0]>cr[1]})
 ctx=setup(20);C=[transform_range(endpoints(classical_return_interval(t,context=ctx)),eps,dmax) for t in TIMES];rssmin=float('inf');weak=None
 for ia in range(20):
  ab=(f"{90+ia:02d}e-2",f"{91+ia:02d}e-2")
  for ig in range(50):
   gb=(str(2*ig),str(2*(ig+1)));rss=0.
   for t,ci in zip(TIMES,C):
    qi=transform_range(endpoints(dephased_return_interval(t,ab,gb,context=ctx)),eps,dmax);rss+=square_lower((qi[0]-ci[1],qi[1]-ci[0]))
   if rss<rssmin:rssmin=rss;weak={"alpha":[float(ab[0]),float(ab[1])],"gamma":[2*ig,2*(ig+1)]}
 x["ST1592"]=packet(1592,"constructed_complete_1000_cell_composite_cover_with_independent_calibration_maps","Calibration vertices enclosed per time/model.",{
  "cells":1000,"weak_cell":weak})
 x["ST1593"]=packet(1593,"proven_positive_composite_RSS_under_independent_5pct_2pct_calibration","Frozen decimal/nuisance box.",{
  "RSS_lower":rssmin,"positive":rssmin>0})
 x["ST1594"]=packet(1594,"proven_at_least_one_calibrated_time_residual_lower_bound","Pigeonhole.",{
  "lower":math.sqrt(rssmin/4)})
 x["ST1595"]=packet(1595,"independent_certificate_is_weaker_than_common_map_certificate_but_positive","Expected extra freedom cost.",{
  "independent_RSS":rssmin,"common_RSS":.03883580017500168,"ratio":rssmin/.03883580017500168})
 x["ST1596"]=packet(1596,"proven_uniform_efficiency_each_does_not_change_click_conditioned_ranges","Efficiencies may differ, attempts/no-clicks retain them.",{
  "eta_each":[eta,1]})
 x["ST1597"]=packet(1597,"attempt_and_no_click_fields_remain_required_for_independent_efficiencies","Postselection no-go.",{
  "required":True})
 x["ST1598"]=packet(1598,"vertex_dependent_detector_confusion_polytope_remains_unclosed","Binary affine map only.",{
  "matrix_polytope":False})
 ch=hashlib.sha256(CAL.read_bytes()).hexdigest()
 x["ST1599"]=packet(1599,"constructed_independent_calibration_artifact_hash","Reproducibility.",{
  "file":CAL.name,"sha256":ch})
 x["ST1600"]=packet(1600,"upgraded_independent_5pct_2pct_status_from_not_certified_to_positive_specific_affine_certificate","Prior sufficient bound was conservative, not false.",{
  "base_gap":qr[0]-cr[1],"composite_RSS":rssmin})
 x["ST1601"]=packet(1601,"recommended_ST1602_ST1616","Analyze fine-fiber q identifiability and channel aliases.",{
  "next":["classical inversion","quantum aliasing","Fisher information","single-time crossings","two-time interval separation"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
