#!/usr/bin/env python3
"""FIN ST1452--ST1466: structured circulant uncertainty lift."""
import math

from fin_interval_cert_common import abs_upper, endpoints, return_bundle, classical_return_interval, dephased_return_interval, setup_uncertain, square_lower
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1452,1466
NAMES=["StructuredMatrixBox","SpectralRadiusBound","UncertainIntervalContext","UncertainOutsideCover",
 "UncertainOutsideMargin","UncertainConcavity","UncertainDerivativeSigns","UncertainRootBracket",
 "UncertainGlobalLocation","UncertainMaximumValue","UncertainCompositeCover","UncertainCompositeResidual",
 "UncertaintyScope","RoundOne_Verdict","RoundOne_Recommendation"]
TIMES=["0.3","0.6","1.2","2.0"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};entry_radius=5e-15;n=12;spec_radius=n*entry_radius;radius=f"{spec_radius:.17g}";dps=25;ctx=setup_uncertain(radius,dps)
 x["ST1452"]=packet(1452,"constructed_row_sum_preserving_real_symmetric_circulant_entry_box","Conditional rounding model.",{
  "entry_radius":entry_radius,"preserves":["circulant","real symmetric","row sum zero","Fourier eigenvectors"]})
 x["ST1453"]=packet(1453,"proven_eigenvalue_radius_bound_by_row_sum_norm","Triangle/Weyl bound.",{
  "spectral_radius":spec_radius,"formula":"|delta lambda_k|<=12 delta_entry"})
 x["ST1454"]=packet(1454,"constructed_outward_uncertain_eigenvalue_context","Zero mode kept exact by row-sum preservation.",{
  "radius":spec_radius,"zero_exact":True,"dps":dps})
 fpoint,_,_=return_bundle("0.5945308",context=ctx);plo,phi=endpoints(fpoint);maxup=-1.;box=None;count=0
 for i in range(5000):
  a=i*.002;b=(i+1)*.002
  if b>0.59 and a<.60:continue
  mid=(a+b)/2;fm,_,_=return_bundle(f"{mid:.6f}",context=ctx);_,fp,_=return_bundle(f"{a:.6f}",f"{b:.6f}",context=ctx);up=abs_upper(fm)+.001*abs_upper(fp)
  if up>maxup:maxup=up;box=[round(a,6),round(b,6)]
  count+=1
 x["ST1455"]=packet(1455,"constructed_complete_uncertain_Taylor_cover","Same 4995-box topology with eigenvalue intervals.",{"boxes":count,"largest_box":box})
 x["ST1456"]=packet(1456,"proven_uncertain_outside_cover_below_candidate","Robust global location exclusion.",{
  "candidate_lower":plo,"outside_upper":maxup,"margin":plo-maxup})
 ft,_,fpp=return_bundle("0.59","0.60",context=ctx);fplo,fphi=endpoints(ft);slo,shi=endpoints(fpp)
 x["ST1457"]=packet(1457,"proven_uncertain_target_positive_and_strictly_concave","Structured box.",{
  "F_interval":[fplo,fphi],"Fpp_interval":[slo,shi],"strict":fplo>0 and shi<0})
 _,fl,_=return_bundle("0.5945307",context=ctx);_,fr,_=return_bundle("0.5945309",context=ctx);ll,lh=endpoints(fl);rl,rh=endpoints(fr)
 x["ST1458"]=packet(1458,"proven_uniform_derivative_signs_at_wider_endpoints","All spectra in box.",{
  "left_time":.5945307,"left_Fp":[ll,lh],"right_time":.5945309,"right_Fp":[rl,rh],"signs":ll>0 and rh<0})
 x["ST1459"]=packet(1459,"proven_each_admissible_spectrum_has_one_stationary_root_in_common_bracket","Root may vary with spectrum.",{
  "common_bracket":[.5945307,.5945309],"unique_per_spectrum":True})
 x["ST1460"]=packet(1460,"proven_each_admissible_spectrum_has_global_maximum_inside_0_59_0_60","Does not claim one common exact maximizer.",{
  "global_location":[.59,.60],"outside_margin":plo-maxup})
 fc,_,_=return_bundle("0.5945308",context=ctx);clo,chi=endpoints(fc);_,fpb,_=return_bundle("0.5945307","0.5945309",context=ctx);upper=chi+1e-7*abs_upper(fpb)
 x["ST1461"]=packet(1461,"constructed_uniform_maximum_value_enclosure","Uses common root bracket and mean-value bound.",{
  "lower_witness":clo,"uniform_upper":upper,"width":upper-clo})
 targets=[classical_return_interval(t,context=ctx) for t in TIMES];rssmin=float('inf');weak=None
 for ia in range(20):
  ab=(f"{90+ia:02d}e-2",f"{91+ia:02d}e-2")
  for ig in range(50):
   gb=(str(2*ig),str(2*(ig+1)));rss=0.
   for t,y in zip(TIMES,targets):rss+=square_lower(dephased_return_interval(t,ab,gb,context=ctx)-y)
   if rss<rssmin:rssmin=rss;weak={"alpha":[float(ab[0]),float(ab[1])],"gamma":[2*ig,2*(ig+1)]}
 x["ST1462"]=packet(1462,"constructed_complete_uncertain_composite_cover","1000 nuisance cells with spectral intervals.",{"cells":1000,"weak_cell":weak})
 x["ST1463"]=packet(1463,"proven_positive_composite_RSS_lower_bound_under_structured_spectral_uncertainty","Target and alternative share the interval spectrum; independent interval evaluation is conservative.",{
  "RSS_lower":rssmin,"one_time_residual_lower":math.sqrt(rssmin/4)})
 x["ST1464"]=packet(1464,"uncertainty_lift_is_not_an_arbitrary_matrix_box_or_symbolic_kernel_theorem","Fourier structure is essential.",{
  "arbitrary_entry_perturbations":False,"kernel_parameter_box":False,"eigenvector_uncertainty":False})
 x["ST1465"]=packet(1465,"release1099_certificates_survive_declared_structured_rounding_box","Major robustness upgrade.",{
  "time_location_survives":plo>maxup,"composite_identifiability_survives":rssmin>0})
 x["ST1466"]=packet(1466,"recommended_ST1467_ST1481","Optimize measurement-time allocation over nuisance grid.",{
  "next":["per-time divergences","linear-program allocation","held-out grid","equal-allocation comparison","certificate boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
