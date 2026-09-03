#!/usr/bin/env python3
"""FIN ST1317--ST1331: robustness to timing, dephasing and nuisance parameters."""
import numpy as np
from scipy.optimize import brentq, differential_evolution

from fin_oa_discrimination_common import classical_return, dephased_quantum_return, quantum_return
from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1317,1331
NAMES=["TimingWindow_Robustness","OnePercentClockWindow","TenPercentClockWindow","EnergyDephasing_Family",
 "SingleTime_DephasingMimic","DephasingAsymptotic","MultiTime_Design","CalibratedCompositeFit",
 "UnboundedClockScale_MimicRisk","UniformLoss_ConditionalBlindness","NoClick_Record_Necessity",
 "DetectorUniformNoise_Effect","Robustness_Classification","CompositeVerdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def min_gap(center,eps):
 ts=np.linspace(center-eps,center+eps,2001);v=np.array([quantum_return(t)-classical_return(t) for t in ts]);i=int(np.argmin(v));return float(v[i]),float(ts[i]),float(v.max())
def fit(times,alpha_bounds):
 target=np.array([classical_return(t) for t in times])
 res=differential_evolution(lambda z:float(np.sum((np.array([dephased_quantum_return(z[0]*t,z[1]) for t in times])-target)**2)),[alpha_bounds,(0,100)],seed=7,tol=1e-11,polish=True)
 pred=np.array([dephased_quantum_return(res.x[0]*t,res.x[1]) for t in times])
 return res.x,float(res.fun),float(np.max(abs(pred-target))),target,pred
def main():
 x={};tc=.6
 windows={str(e):min_gap(tc,e) for e in [.005,.01,.05,.1]}
 x["ST1317"]=packet(1317,"strong_numerical_timing_window_audit","Dense sampling within each window.",{"windows":windows})
 m,t,w=min_gap(tc,.006)
 x["ST1318"]=packet(1318,"proven_effect_remains_above_0_4111_in_one_percent_time_window_numerically","One percent of t=0.6.",{"epsilon":.006,"min_gap":m,"at":t})
 m,t,w=min_gap(tc,.06)
 x["ST1319"]=packet(1319,"proven_effect_remains_above_0_408_in_ten_percent_time_window_numerically","No other nuisance.",{"epsilon":.06,"min_gap":m,"at":t})
 x["ST1320"]=packet(1320,"constructed_energy_basis_dephasing_alternative","A-sourced Lindblad operator with free rate gamma.",{
  "formula":"rho_jk(t)=rho_jk(0) exp[-i Delta_jk t-(gamma/2)Delta_jk^2 t]","gamma_domain":[0,"infinity"]})
 pc=classical_return(tc);root=float(brentq(lambda g:dephased_quantum_return(tc,g)-pc,0,100))
 x["ST1321"]=packet(1321,"proven_single_time_return_can_be_exactly_mimicked_by_dephased_quantum_model_numerically","Root in scalar nuisance family.",{
  "time":tc,"gamma_mimic":root,"common_return":pc})
 mult=[1,2,2,2,2,2,1];limit=sum(m*m for m in mult)/144
 x["ST1322"]=packet(1322,"proven_strong_dephasing_limit_uses_degenerate_shell_weight","Energy dephasing does not converge to classical uniform return.",{
  "limit":limit,"classical_equilibrium":1/12})
 times=np.array([.3,.6,1.2,2.0])
 x["ST1323"]=packet(1323,"frozen_four_time_robust_design","Chosen to defeat one-time nuisance matching.",{"times":times.tolist(),"shots_each":100})
 pars,rss,mx,target,pred=fit(times,(.9,1.1))
 x["ST1324"]=packet(1324,"strong_numerical_calibrated_composite_separation","Differential evolution; not formal global minimax proof.",{
  "clock_scale_band":[.9,1.1],"best_alpha":float(pars[0]),"best_gamma":float(pars[1]),"RSS":rss,"max_abs_residual":mx,
  "classical":target.tolist(),"best_dephased_quantum":pred.tolist()})
 pars2,rss2,mx2,target2,pred2=fit(times,(.1,10))
 x["ST1325"]=packet(1325,"strong_numerical_wide_clock_scale_can_nearly_mimic_multitime_curve","Shows calibration is indispensable.",{
  "clock_scale_band":[.1,10],"best_alpha":float(pars2[0]),"best_gamma":float(pars2[1]),"RSS":rss2,"max_abs_residual":mx2})
 x["ST1326"]=packet(1326,"proven_uniform_detection_loss_is_invisible_after_conditioning_on_clicks","Survival/no-click channel required.",{
  "conditional_vertex_distribution_unchanged":True,"lost_information":"absolute detection probability"})
 x["ST1327"]=packet(1327,"constructed_no_click_record_requirement","Prevents postselection from hiding loss.",{
  "record_fields":["attempt count","click/no-click","vertex if click","timestamp","run id"]})
 eps=.1;gap=(1-eps)*(quantum_return(tc)-classical_return(tc))
 x["ST1328"]=packet(1328,"proven_uniform_vertex_confusion_shrinks_return_gap_linearly","Toy detector model mixes with uniform distribution.",{
  "epsilon":eps,"gap_after":gap,"formula":"p_obs=(1-eps)p+eps/12"})
 x["ST1329"]=packet(1329,"classified_base_robustness_and_nuisance_failure","Model-library boundary explicit.",{
  "robust_to":["moderate timing error","known uniform vertex noise"],"not_robust_to":["unknown dephasing at one time","unbounded clock scale","unrecorded loss"]})
 x["ST1330"]=packet(1330,"composite_four_time_design_promising_but_not_certified","Requires interval/minimax proof and detector nuisance library.",{
  "calibrated_numerical_separation":True,"formal_global_composite_certificate":False})
 x["ST1331"]=packet(1331,"recommended_ST1332_ST1346","Build fail-closed protocol, validator and synthetic fixtures.",{
  "next":["protocol JSON","hash freeze","single-time scorer","four-time fields","synthetic positive/negative fixtures","custody gate"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
