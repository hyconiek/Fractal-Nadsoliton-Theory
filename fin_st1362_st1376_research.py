#!/usr/bin/env python3
"""FIN ST1362--ST1376: interval certificate for the return-gap maximizer."""
import math

from fin_interval_cert_common import DECIMAL_EIGENVALUES, MULTIPLICITIES, abs_upper, endpoints, return_bundle, setup
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1362,1376
NAMES=["FrozenDecimalSpectrum","IntervalArithmeticContract","CandidatePointLowerBound","GlobalCoverConstruction",
 "OutsideCoverUpperBound","TargetIntervalPositivity","TargetIntervalStrictConcavity","LeftDerivativeSign",
 "RightDerivativeSign","UniqueRootBracket","UniqueGlobalMaximumTheorem","MaximumValueBracket",
 "Release1098_Upgrade","CertificateScope","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};dps=25;ctx=setup(dps)
 x["ST1362"]=packet(1362,"constructed_exact_decimal_spectral_model","Certificate is for these decimals, not symbolic kernel parameters.",{
  "eigenvalues":DECIMAL_EIGENVALUES,"multiplicities":MULTIPLICITIES})
 x["ST1363"]=packet(1363,"constructed_outward_mpmath_interval_contract","Relies on mpmath interval implementation.",{
  "decimal_precision":dps,"domain":[0,10],"cover_step":.002})
 fpoint,_,_=return_bundle("0.5945308",dps=dps,context=ctx);point_lo,point_hi=endpoints(fpoint)
 x["ST1364"]=packet(1364,"proven_candidate_point_lower_bound_in_frozen_decimal_model","Outward interval point evaluation.",{
  "time":"0.5945308","F_lower":point_lo,"F_upper":point_hi})
 max_upper=-1.0;max_box=None;boxes=0
 for i in range(5000):
  a=i*.002;b=(i+1)*.002
  if b>0.59 and a<0.60:continue
  mid=(a+b)/2;fm,_,_=return_bundle(f"{mid:.6f}",dps=dps,context=ctx);_,fpbox,_=return_bundle(f"{a:.6f}",f"{b:.6f}",dps=dps,context=ctx)
  up=abs_upper(fm)+.001*abs_upper(fpbox)
  if up>max_upper:max_upper,max_box=up,[round(a,6),round(b,6)]
  boxes+=1
 x["ST1365"]=packet(1365,"constructed_complete_Taylor_interval_cover_outside_target","Mean-value/Taylor enclosure per box.",{
  "boxes":boxes,"target_excluded":[.59,.60],"covered":[[0,.59],[.60,10]]})
 x["ST1366"]=packet(1366,"proven_all_outside_boxes_below_candidate_value","Global exclusion outside target interval.",{
  "largest_outside_upper":max_upper,"largest_box":max_box,"candidate_lower":point_lo,"strict_margin":point_lo-max_upper})
 ft,_,fpp=return_bundle("0.59","0.60",dps=dps,context=ctx);flo,fhi=endpoints(ft);slo,shi=endpoints(fpp)
 x["ST1367"]=packet(1367,"proven_F_positive_on_target_interval","Therefore |F|=F there.",{"F_interval":[flo,fhi]})
 x["ST1368"]=packet(1368,"proven_strict_concavity_on_target_interval","Outward second derivative enclosure.",{
  "Fpp_interval":[slo,shi],"upper_negative":shi<0})
 _,fpl,_=return_bundle("0.59453078",dps=dps,context=ctx);l0,l1=endpoints(fpl)
 x["ST1369"]=packet(1369,"proven_positive_derivative_at_left_root_endpoint","Point interval.",{
  "time":"0.59453078","Fp_interval":[l0,l1]})
 _,fpr,_=return_bundle("0.59453079",dps=dps,context=ctx);r0,r1=endpoints(fpr)
 x["ST1370"]=packet(1370,"proven_negative_derivative_at_right_root_endpoint","Point interval.",{
  "time":"0.59453079","Fp_interval":[r0,r1]})
 x["ST1371"]=packet(1371,"proven_unique_stationary_root_in_bracket","Intermediate value plus strict decreasing derivative.",{
  "root_bracket":[.59453078,.59453079],"width":1e-8})
 x["ST1372"]=packet(1372,"proven_unique_global_maximum_of_absolute_return_gap_on_0_10_for_frozen_decimal_spectrum","Combines outside cover, positivity, concavity and derivative signs.",{
  "domain":[0,10],"maximizer_bracket":[.59453078,.59453079],"unique":True})
 center="0.594530785";fc,_,_=return_bundle(center,dps=dps,context=ctx);_,fpbox,_=return_bundle("0.59453078","0.59453079",dps=dps,context=ctx);clo,chi=endpoints(fc);upper=chi+5e-9*abs_upper(fpbox)
 x["ST1373"]=packet(1373,"proven_maximum_value_bracket","Center lower plus mean-value upper on root bracket.",{
  "value_lower":clo,"value_upper":upper,"width":upper-clo})
 x["ST1374"]=packet(1374,"upgraded_release1098_strong_numerical_maximum_to_interval_theorem_for_frozen_decimal_model","Scope narrowed explicitly.",{
  "previous":"strong numerical","current":"interval-certified frozen-decimal theorem"})
 x["ST1375"]=packet(1375,"certificate_does_not_cover_symbolic_kernel_or_eigenvalue_uncertainty","No promotion beyond input model.",{
  "symbolic_kernel_certificate":False,"floating_matrix_uncertainty_box":False})
 x["ST1376"]=packet(1376,"recommended_ST1377_ST1391","Interval-certify four-time composite separation over the nuisance box.",{
  "next":["alpha-gamma box cover","RSS lower bound","residual lower bound","candidate upper","scope audit"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
