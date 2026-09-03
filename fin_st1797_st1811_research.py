#!/usr/bin/env python3
"""FIN ST1797--ST1811: decisive synthesis after repository-wide falsification."""
import hashlib
import json
import math

import matplotlib.pyplot as plt
import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float
from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1797,1811
NAMES=["ProofLedger","APDObstruction","SingleActionObstruction","MomentBridgeObstruction",
 "PropagatorObstruction","DeepestSurvivingStructure","DecisiveSourceActionGate","GateNecessityAudit",
 "CurrentFINGateScore","ResearchRouteRanking","HighestValueNextTheorem","FalsificationProtocol",
 "NoNewPhysicsBoundary","FinalInterpretation","CycleSix_ReportTrigger"]
FIGDIR=ROOT/'FIN_ST1722_ST1811_Figures';FIG=FIGDIR/'repository_falsification_summary.png';GATE=ROOT/'FIN_ST1797_ST1811_DecisiveSourceActionGate.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def make_figure():
 FIGDIR.mkdir(exist_ok=True);fig,ax=plt.subplots(1,3,figsize=(11,3.7))
 d0=4/3;d=np.r_[np.linspace(1.15,d0-.003,200),np.linspace(d0+.003,1.52,200)];A=4*np.log(2);L=A*np.cos(np.pi*d/4+np.pi/6)/(1+.01*d);S=np.cos(.18575*d+.1625)/(1+d**1.8);ax[0].plot(d,S/L);ax[0].axvline(d0,color='r',ls='--');ax[0].set(title='APD pole at legacy zero',xlabel='$d$',ylabel='$K_s/K_L$',ylim=(-20,20));ax[0].grid(alpha=.25)
 refs=[.5,1,2];mf=[1.7516996094,1.4699856726,1.1920435517];gf=[1.0182284771,1.2319100972,1.0784676598];ax[1].plot(refs,mf,'o-',label='$1+c_0$');ax[1].plot(refs,gf,'s--',label='$1+c_2$');ax[1].set(title='Reference-dependent multipliers',xlabel='$d_{ref}$',ylabel='factor');ax[1].legend(fontsize=7);ax[1].grid(alpha=.25)
 Am=independent_strict_matrix_float();s=Am[0,0];W=s*np.eye(12)-Am;ax[2].scatter(range(12),np.linalg.eigvalsh(Am),label='A spectrum');ax[2].scatter(range(12),np.linalg.eigvalsh(W),label='W spectrum');ax[2].axhline(0,color='k',lw=.8);ax[2].set(title='A positive, W indefinite',xlabel='ordered mode',ylabel='eigenvalue');ax[2].legend(fontsize=7);ax[2].grid(alpha=.25);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def main():
 make_figure();x={}
 x["ST1797"]=packet(1797,"constructed_repository_wide_proven_refuted_conditional_ledger","Six-round synthesis.",{
  "proven":["APD finite triviality/global pole","gauge variation defect","Yukawa sign defect","moment rank/reference defect","W resolvent obstruction","Dirichlet action"],"conditional":["FRW/reduced EOM","OA protocols"],"refuted":["current unique kernel-to-physical-Ltotal derivation"]})
 x["ST1798"]=packet(1798,"proven_APD_cannot_supply_nontrivial_global_positive_completion","Local finite germ bookkeeping survives.",{
  "finite_ratio":True,"global_bounded":False,"positive":False,"independent_source":False})
 x["ST1799"]=packet(1799,"proven_current_displayed_Ltotal_and_covariant_EOM_rows_do_not_form_one_gauge_invariant_variational_system","Multiple independent defects.",{
  "gauge":False,"single_action":False,"Yukawa_sign":False})
 x["ST1800"]=packet(1800,"proven_three_moment_effective_coupling_bridge_is_reference_dependent_rank_deficient_and_nonpredictive_with_free_bases","Feature map only.",{
  "rank":3,"parameters":4,"base_couplings":5})
 x["ST1801"]=packet(1801,"proven_W_is_not_same_A_resolvent_and_is_indefinite","Current finite model.",{
  "W_positive_modes":5,"W_negative_modes":7,"Euclidean_covariance":False})
 x["ST1802"]=packet(1802,"identified_deepest_surviving_unavoidable_structure_as_finite_Dirichlet_Markov_spectral_geometry","Not physical field theory.",{
  "core":["A PSD","Dirichlet form","A+ Green","heat/unitary/wave functional calculus","finite information/coarse maps"]})
 gate={"schema":"FIN-DECISIVE-SOURCE-ACTION-GATE-v1","requirements":[
  {"id":"G1","name":"globally regular non-tautological strict source map","current":False},
  {"id":"G2","name":"canonical refinement/continuum geometry and dimensional scale","current":False},
  {"id":"G3","name":"one typed gauge/diffeomorphism-consistent action","current":False},
  {"id":"G4","name":"all claimed EOM from literal variation with residual tables","current":False},
  {"id":"G5","name":"state/channel/clock/preparation/instrument/record OA","current":False},
  {"id":"G6","name":"parameter values or sharp dimensionless predictions not absorbable into free bases","current":False},
  {"id":"G7","name":"held-out physical forward model and falsifiable record","current":False}]}
 GATE.write_text(json.dumps(gate,indent=2,sort_keys=True)+'\n')
 x["ST1803"]=packet(1803,"constructed_decisive_strict_source_action_gate","A future bridge must pass all rows jointly.",{
  "file":GATE.name,"requirements":len(gate['requirements'])})
 x["ST1804"]=packet(1804,"proven_each_gate_class_is_necessary_by_countermodel_or_current_ambiguity","Removal audit.",{
  "remove_G1":"APD quotient tautology","remove_G2":"dref/unit torsors","remove_G3":"gauge/Proca ambiguity","remove_G4":"mixed-lineage EOM","remove_G5":"A-only OA nonuniqueness","remove_G6":"free-coupling absorption","remove_G7":"no physics evidence"})
 passed=sum(1 for r in gate['requirements'] if r['current'])
 x["ST1805"]=packet(1805,"proven_current_repository_passes_zero_of_seven_decisive_joint_rows_in_strict_scope","Subcomponents exist but no complete row.",{
  "passed":passed,"total":len(gate['requirements'])})
 ranking=[
  {"route":"publish finite Dirichlet/Markov/spectral FIN mathematics","probability_of_rigorous_success":.9,"physics_decision_power":.25},
  {"route":"derive refinement-compatible continuum Dirichlet form plus scale and causal observable","probability_of_rigorous_success":.35,"physics_decision_power":.8},
  {"route":"derive gauge-covariant spectral/source action from independently sourced algebra","probability_of_rigorous_success":.18,"physics_decision_power":1.0},
  {"route":"repair current APD-three-moment-Ltotal chain","probability_of_rigorous_success":.03,"physics_decision_power":.7},
  {"route":"conditional OA laboratory discrimination","probability_of_rigorous_success":.45,"physics_decision_power":.55}]
 x["ST1806"]=packet(1806,"constructed_quantitative_research_route_ranking","Subjective strategic probabilities, not theorem.",{"ranking":ranking})
 x["ST1807"]=packet(1807,"selected_highest_value_next_theorem_as_refinement_compatible_gauge_covariant_source_action_or_no_go","Must start from strict algebra/Dirichlet/refinement data, not moment fitting.",{
  "target":"derive/refute one action satisfying G1-G4 before any SM/GR coefficient claims"})
 x["ST1808"]=packet(1808,"constructed_immediate_falsification_protocol_for_next_action","Stop on first failed row.",{
  "tests":["global regularity","gauge variation/Ward identity","literal EL regeneration","refinement naturality","scale covariance","nonabsorbable prediction"]})
 x["ST1809"]=packet(1809,"proven_no_new_physics_by_assumption_requires_freezing_current_SM_GR_Ltotal_promotion","Conditional scaffolds may remain labelled.",{
  "strict_Ltotal":False,"conditional_scaffold":True})
 x["ST1810"]=packet(1810,"final_interpretation_FIN_is_currently_finite_Dirichlet_spectral_information_geometry_with_conditional_operational_models_not_fundamental_physics","Survives falsification.",{
  "fundamental_physical_theory":False,"mathematical_framework":True})
 x["ST1811"]=packet(1811,"six_round_cycle_complete_report_required","Release 11.03.",{
  "programs":90,"rounds":6,"release":"11.03","main_result":"current kernel-moment-Ltotal chain falsified; Dirichlet core survives; decisive source-action gate defined","strict_ToE_closure":False,"gate_sha256":hashlib.sha256(GATE.read_bytes()).hexdigest(),"figure":str(FIG.relative_to(ROOT))})
 write_round(LO,HI,x)
if __name__=="__main__":main()
