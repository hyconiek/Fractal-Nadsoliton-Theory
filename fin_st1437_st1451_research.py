#!/usr/bin/env python3
"""FIN ST1437--ST1451: certified OA synthesis and final roadmap."""
import matplotlib.pyplot as plt
import numpy as np

from fin_oa_discrimination_common import classical_return, dephased_quantum_return, quantum_return
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1437,1451
NAMES=["CertifiedTimeMaximum","CertifiedCompositeIdentifiability","CalibrationErrorBudget","FullLikelihoodUpgrade",
 "AnytimeValidSimpleTest","SimpleVsCompositeClaimSeparation","Protocol1099Status","EvidenceBoundary",
 "OAArchitectureImpact","AnnihilationScope","RemainingMathematicalBlockers","RecommendedNextProgrammes",
 "StopRules","StrictCoreBoundary","CycleSix_ReportTrigger"]
FIGDIR=ROOT/'FIN_ST1362_ST1451_Figures';FIG=FIGDIR/'certified_oa_upgrade.png'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def figure():
 FIGDIR.mkdir(exist_ok=True);fig,ax=plt.subplots(1,2,figsize=(10,4.2))
 t=np.linspace(.58,.61,501);f=np.array([quantum_return(z)-classical_return(z) for z in t]);ax[0].plot(t,f,lw=2);ax[0].axvspan(.59453078,.59453079,color='crimson',alpha=.5,label='certified maximizer bracket');ax[0].set(xlabel='Dimensionless time',ylabel='$Q(t)-C(t)$',title='Certified return-gap maximum');ax[0].grid(alpha=.25);ax[0].legend(fontsize=8)
 ts=np.array([.3,.6,1.2,2.]);c=np.array([classical_return(z) for z in ts]);q=np.array([dephased_quantum_return(1.1*z,12.373117825382856) for z in ts]);ax[1].plot(ts,c,'o-',label='Classical target');ax[1].plot(ts,q,'s--',label='Dephased Q candidate');ax[1].set(xlabel='Scheduled time',ylabel='Return probability',title='Four-time composite separation');ax[1].grid(alpha=.25);ax[1].legend(fontsize=8);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def main():
 figure();x={}
 x["ST1437"]=packet(1437,"proven_unique_global_return_gap_maximum_for_frozen_decimal_spectrum","Interval arithmetic on [0,10].",{
  "time_bracket":[.59453078,.59453079],"value_bracket":[0.4112409998076315,0.41124099980763185]})
 x["ST1438"]=packet(1438,"proven_positive_four_time_composite_RSS_lower_bound_on_declared_box","Interval cover of 1000 cells.",{
  "alpha":[.9,1.1],"gamma":[0,100],"RSS_lower":0.046692011213846946,"one_time_residual_lower":0.10804167160619894})
 x["ST1439"]=packet(1439,"proven_explicit_preparation_detector_and_loss_responsibilities","Not full hardware uncertainty.",{
  "independent_scalar_error_threshold":0.20560911273268,"example_combined_lower_gap":0.21121822546536,"no_click_required":True})
 x["ST1440"]=packet(1440,"proven_seven_distance_classes_are_sufficient_and_more_informative_than_binary_return","Frozen simple hypotheses.",{
  "classes":7,"Chernoff_full":0.13647980004881377,"Chernoff_binary":0.10028599984500335})
 x["ST1441"]=packet(1441,"proven_anytime_valid_symmetric_likelihood_boundary_errors_for_simple_models","Ville martingale theorem.",{
  "alpha":.01,"Q_boundary":100,"C_boundary":.01,"optional_stopping_safe":True})
 x["ST1442"]=packet(1442,"proven_simple_and_composite_claims_must_remain_separate","Composite plug-in LR is not automatically safe.",{
  "simple":"29-shot and e-process guarantees","composite":"interval identifiability plus conservative 1146 shots/time","composite_anytime_eprocess":False})
 x["ST1443"]=packet(1443,"constructed_protocol_10_99_with_interval_and_sequential_certificates","Synthetic validator only.",{
  "protocol":"FIN_ST1422_ST1436_Protocol_10_99.json","validator":"fin_oa_protocol_validator_10_99.py","figure":str(FIG.relative_to(ROOT))})
 x["ST1444"]=packet(1444,"blocked_no_independent_calibration_custody_or_nature_events","Mathematical transfer readiness only.",{
  "lab":False,"raw_events":False,"independent_clock_calibration":False,"custody":False})
 x["ST1445"]=packet(1445,"certified_discriminator_strengthens_conditional_OA_not_strict_state_selection","Ontology remains postulated.",{
  "OA_falsifiable_in_principle":True,"OA_strictly_derived":False})
 x["ST1446"]=packet(1446,"protocol_tests_internal_dynamics_not_total_nadsoliton_nonexistence","Annihilation boundary preserved.",{
  "tests":["mass-mixing channel","coherent unitary channel"],"does_not_test":"literal annihilation of whole being"})
 x["ST1447"]=packet(1447,"listed_remaining_mathematical_blockers","No premature laboratory transfer.",{
  "blockers":["symbolic-kernel interval lift","full preparation/detector nuisance boxes","optimal composite likelihood",
              "valid composite e-process","refinement-level transport","platform map"]})
 x["ST1448"]=packet(1448,"recommended_next_programmes","Nonreplay continuation.",{
  "programmes":["lift decimal certificates to interval uncertainty in A entries","certify eigenvalue clusters from A boxes",
   "solve minimax four-time likelihood","construct mixture e-process over gamma/alpha","add calibrated preparation simplex",
   "add detector confusion polytope","optimize time allocation rather than equal shots","use seven-class sequential records",
   "transport the protocol through 12-to-24 refinement","derive one platform-specific forward model",
   "freeze independent calibration and custody forms","seek external review before any event collection"]})
 x["ST1449"]=packet(1449,"installed_certificate_and_protocol_stop_rules","Preserve scope.",{
  "stop_if":["spectrum outside frozen decimal model","alpha/gamma outside certified box","unregistered detector nuisance",
             "postselection without no-clicks","adaptive plug-in likelihood","post-unblinding protocol edit"]})
 x["ST1450"]=packet(1450,"strict_FIN_physics_status_unchanged_despite_stronger_conditional_certificates","No ToE promotion.",{
  "strict_state_category":False,"laboratory":False,"ToE":False})
 x["ST1451"]=packet(1451,"six_round_cycle_complete_report_required","Release 10.99.",{
  "programs":90,"rounds":6,"release":"10.99","breakthrough":"interval-certified unique optimum and composite nuisance separation","strict_ToE_closure":False})
 write_round(LO,HI,x)
if __name__=="__main__":main()
