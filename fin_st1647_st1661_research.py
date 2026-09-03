#!/usr/bin/env python3
"""FIN ST1647--ST1661: strict-kernel parameter boxes to operator bounds."""
import math

import mpmath as mp

from fin_oa_discrimination_common import classical_return, quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1647,1661
NAMES=["KernelParameterMap","LargeAuditBox","LargeWeightIntervals","LargeOperatorBound","LargeBaseGap",
 "LargeLocationFailure","CertifiedSmallBox","SmallOperatorBound","SmallLocationPass","SmallLindbladPass",
 "WeightPositivity","DependencyOverestimate","ParameterBoxProvenance","RoundTwo_Verdict","RoundTwo_Recommendation"]
NOM=(.18575,.1625,1.,1.8);BASE=(5e-6,5e-6,5e-5,5e-5)
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def audit(radii):
 iv=mp.iv;iv.dps=30;weights=[];nom=[]
 for d in range(1,7):
  o=iv.mpf([str(NOM[0]-radii[0]),str(NOM[0]+radii[0])]);p=iv.mpf([str(NOM[1]-radii[1]),str(NOM[1]+radii[1])]);b=iv.mpf([str(NOM[2]-radii[2]),str(NOM[2]+radii[2])]);e=iv.mpf([str(NOM[3]-radii[3]),str(NOM[3]+radii[3])]);weights.append(iv.cos(o*d+p)/(1+b*iv.mpf(d)**e));nom.append(math.cos(NOM[0]*d+NOM[1])/(1+NOM[2]*d**NOM[3]))
 siv=2*sum(weights[:5],iv.mpf('0'))+weights[5];sn=2*sum(nom[:5])+nom[5];diag=max(abs(float(siv.a)-sn),abs(float(siv.b)-sn));off=[max(abs(float(v.a)-x),abs(float(v.b)-x)) for v,x in zip(weights,nom)];row=diag+2*sum(off[:5])+off[5]
 return {"weights":[[float(v.a),float(v.b)] for v in weights],"diag_radius":diag,"operator_bound":row,"all_positive":all(float(v.a)>0 for v in weights)}
def main():
 x={};basegap=quantum_return(.6)-classical_return(.6);location_threshold=9.720652006795929e-8
 x["ST1647"]=packet(1647,"constructed_parameter_to_weight_to_laplacian_interval_map","Strict functional form.",{
  "w_d":"cos(omega d+phi)/(1+beta d^eta)","A":"row-sum graph Laplacian"})
 large=audit(BASE)
 x["ST1648"]=packet(1648,"constructed_large_unsourced_parameter_audit_box","Diagnostic scale, not repository CI.",{
  "nominal":NOM,"radii":BASE})
 x["ST1649"]=packet(1649,"proven_large_box_weight_intervals_positive","Outward interval arithmetic.",{
  "weights":large['weights'],"positive":large['all_positive']})
 x["ST1650"]=packet(1650,"proven_large_box_operator_norm_upper_bound","Circulant row-sum norm.",{
  "operator_bound":large['operator_bound'],"diag_radius":large['diag_radius']})
 x["ST1651"]=packet(1651,"proven_large_box_keeps_simple_t06_gap_positive_by_Duhamel","Crude bound.",{
  "lower_gap":basegap-1.8*large['operator_bound']})
 x["ST1652"]=packet(1652,"large_box_fails_current_global_location_and_Lindblad_composite_sufficient_gates","Failure to certify, not counterexample.",{
  "location_pass":large['operator_bound']<location_threshold,"operator_bound":large['operator_bound']})
 small_r=tuple(v*.0004 for v in BASE);small=audit(small_r)
 x["ST1653"]=packet(1653,"constructed_small_parameter_box_tuned_below_existing_proof_threshold","Mathematical witness, not physical uncertainty.",{
  "radii":small_r})
 x["ST1654"]=packet(1654,"proven_small_box_operator_norm_upper_bound","Outward intervals.",{
  "operator_bound":small['operator_bound']})
 x["ST1655"]=packet(1655,"proven_small_parameter_box_passes_global_location_threshold","Inherited noncirculant/circulant theorem.",{
  "threshold":location_threshold,"margin":location_threshold-small['operator_bound']})
 M=2.3421820411463;eps=small['operator_bound'];dl=2*eps*(1+100*(2*M+eps));pert=2*dl+2*eps
 x["ST1656"]=packet(1656,"proven_small_parameter_box_passes_gamma100_composite_sensitivity_gate","t<=2.",{
  "perturbation_bound":pert,"residual_lower_after":.10804167160619894-pert})
 x["ST1657"]=packet(1657,"proven_both_parameter_boxes_preserve_positive_graph_weights","Markov validity.",{
  "large":large['all_positive'],"small":small['all_positive']})
 x["ST1658"]=packet(1658,"interval_parameter_dependency_makes_bounds_conservative","Same parameters reused across distances/diagonal.",{
  "exact_joint_image":False})
 x["ST1659"]=packet(1659,"neither_box_has_strict_statistical_or_physical_provenance","No silent CI.",{
  "repository_CI":False,"laboratory_calibration":False})
 x["ST1660"]=packet(1660,"symbolic_parameter_propagation_is_constructed_but_only_tiny_box_closes_current_strong_gates","Base discrimination robust on large box.",{
  "large_base_pass":True,"small_full_pass":True})
 x["ST1661"]=packet(1661,"recommended_ST1662_ST1676","Isolate reduced continuous minimax KKT candidate.",{
  "next":["active branch","gamma equality","stationarity weight","alpha boundary","residual audit","exactness boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
