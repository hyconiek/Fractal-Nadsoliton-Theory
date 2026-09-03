#!/usr/bin/env python3
"""FIN ST1347--ST1361: final OA discrimination verdict and roadmap."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fin_oa_discrimination_common import classical_return, dephased_quantum_return, quantum_return
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1347,1361
NAMES=["ClaimHierarchy","IdealSimpleHypothesisTheorem","FullVertexAdvantage","TimingRobustnessVerdict",
 "SingleTimeCompositeNoGo","FourTimeCompositeEvidence","ClockCalibrationNecessity","LossAccountingNecessity",
 "ExecutableSpec_Status","LaboratoryBoundary","ConditionalPhysicsMeaning","AnnihilationRelevance",
 "RecommendedNextProgrammes","ResearchStopRules","CycleSix_ReportTrigger"]
FIGDIR=ROOT/'FIN_ST1272_ST1361_Figures';FIG=FIGDIR/'oa_return_discrimination.png'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def make_figure():
 FIGDIR.mkdir(exist_ok=True);t=np.linspace(0,3,601);c=np.array([classical_return(z) for z in t]);q=np.array([quantum_return(z) for z in t]);qd=np.array([dephased_quantum_return(z,12.3) for z in t])
 fig,ax=plt.subplots(figsize=(8.2,4.8));ax.plot(t,c,label='Classical heat return $C(t)$',lw=2);ax.plot(t,q,label='Closed quantum return $Q(t)$',lw=2);ax.plot(t,qd,label='Dephased quantum, $\\gamma=12.3$',lw=1.6,ls='--')
 for z in [.3,.6,1.2,2.0]:ax.axvline(z,color='0.75',lw=.8)
 ax.scatter([.6],[classical_return(.6)],s=35);ax.scatter([.6],[quantum_return(.6)],s=35);ax.set(xlabel='Dimensionless time $t$',ylabel='Return probability',ylim=(0,1.03));ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def main():
 make_figure();x={}
 x["ST1347"]=packet(1347,"constructed_epistemic_claim_hierarchy","Prevents ideal predictions from becoming evidence.",{
  "levels":["exact finite model","strong numerical optimization","synthetic protocol validation","laboratory evidence absent"]})
 x["ST1348"]=packet(1348,"proven_ideal_simple_hypotheses_have_sub_one_percent_29_shot_test","Exact binomial-tail numerical evaluation.",{
  "time":.6,"shots":29,"threshold":19,"errors":[0.0071823135380980115,0.007797504846906937]})
 x["ST1349"]=packet(1349,"proven_full_vertex_readout_has_higher_Chernoff_information","Declared frozen distributions.",{
  "binary":0.10028599984500335,"full_vertex":0.13647980004881377})
 x["ST1350"]=packet(1350,"strong_timing_robustness_for_base_models","No dephasing/loss.",{
  "ten_percent_window_min_gap":0.408,"base_gap":quantum_return(.6)-classical_return(.6)})
 x["ST1351"]=packet(1351,"proven_single_time_return_is_nonidentifying_with_free_dephasing","Scalar nuisance root exists.",{
  "gamma_mimic_near":12.6,"consequence":"primary rule tests only frozen C versus closed Q, not C versus all open quantum models"})
 x["ST1352"]=packet(1352,"strong_numerical_four_time_composite_separation_with_calibrated_clock","Not interval-certified minimax.",{
  "times":[.3,.6,1.2,2.0],"clock_band":[.9,1.1],"best_RSS":0.05239469837328462,"best_max_residual":0.15060643505744295})
 x["ST1353"]=packet(1353,"proven_clock_calibration_is_structurally_necessary_for_claimed_composite_test","Wide scale freedom nearly mimics curve.",{
  "wide_band":[.1,10],"best_max_residual":0.01950198,"CA_value_must_be_external_or_calibrated":True})
 x["ST1354"]=packet(1354,"proven_no_click_accounting_is_necessary_to_detect_loss","Conditional clicks erase uniform survival loss.",{
  "attempts_and_clicks_required":True,"postselected_only_record_invalid":True})
 x["ST1355"]=packet(1355,"constructed_executable_spec_and_fail_closed_validator","Synthetic fixtures only.",{
  "protocol":"FIN_ST1332_ST1346_Protocol.json","validator":"fin_oa_protocol_validator.py","figure":str(FIG.relative_to(ROOT))})
 x["ST1356"]=packet(1356,"blocked_no_apparatus_raw_nature_record_or_independent_custody","Physical boundary unchanged.",{
  "laboratory":False,"provider_registrar_analyst_split":False,"holdout":False})
 x["ST1357"]=packet(1357,"conditional_OA_makes_dual_dynamics_empirically_distinguishable_in_principle","Only after preparation/clock/instrument assumptions.",{
  "mathematical_prediction":True,"physical_realization":False})
 x["ST1358"]=packet(1358,"proven_protocol_discriminates_internal_channel_models_not_literal_total_annihilation","Scope of relevance.",{
  "tests":"heat mass redistribution versus unitary norm dynamics","does_not_test":"nonexistence of the whole nadsoliton"})
 x["ST1359"]=packet(1359,"recommended_next_programmes","Ranked nonreplay continuation.",{
  "programmes":["interval-certify the time optimum","solve exact composite minimax discrimination","add state-preparation contamination",
   "add detector confusion matrix and drift","use full multinomial likelihood","derive confidence sequences for sequential runs",
   "certify four-time identifiability over nuisance boxes","map OA fields to one candidate platform without claiming realization",
   "freeze custody/hash/holdout transfer package","obtain independent laboratory review before data",
   "only then collect raw event-level records","compare failure with predeclared alternative library"]})
 x["ST1360"]=packet(1360,"installed_claim_and_research_stop_rules","No p-hacking or channel repair.",{
  "stop_if":["matrix hash mismatch","clock outside band","unrecorded no-clicks","invalid counts","post-unblinding model edit","nuisance class not predeclared"]})
 x["ST1361"]=packet(1361,"six_round_cycle_complete_report_required","No external evidence.",{
  "programs":90,"rounds":6,"release":"10.98","main_result":"large ideal OA discriminator but single-time composite no-go; calibrated four-time route remains strong numerical",
  "strict_ToE_closure":False})
 write_round(LO,HI,x)
if __name__=="__main__":main()
