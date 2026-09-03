#!/usr/bin/env python3
"""FIN ST1707--ST1721: Release 11.02 synthesis."""
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np

from fin_total_nadsoliton_common import ROOT, write_packet, write_round


LO,HI=1707,1721
NAMES=["ExtendedRobustnessLattice","LindbladSensitivityVerdict","KernelParameterBoxVerdict","KKTGeometryVerdict",
 "SymmetricGridDecisionVerdict","DetectorMatrixVerdict","JointQClockVerdict","Protocol1102Freeze",
 "EvidenceBoundary","AnnihilationScope","RemainingBlockers","RecommendedNextProgrammes","StopRules","FinalScientificVerdict","CycleSix_ReportTrigger"]
FIGDIR=ROOT/'FIN_ST1632_ST1721_Figures';FIG=FIGDIR/'lindblad_kernel_qclock_summary.png';PROTO=ROOT/'FIN_ST1707_ST1721_Protocol_11_02.json'
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def h(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def figure():
 FIGDIR.mkdir(exist_ok=True);fig,ax=plt.subplots(1,3,figsize=(11,3.7))
 eps=np.logspace(-14,-4,250);M=2.3421820411463;g=100;pert=2*(2*eps*(1+g*(2*M+eps)))+2*eps;ax[0].loglog(eps,pert);ax[0].axhline(.1080416716,color='r',ls='--');ax[0].set(title='Dephasing perturbation bound',xlabel='$\\epsilon$',ylabel='worst return perturbation');ax[0].grid(alpha=.25)
 vals=[.0002010564907,8.041872775e-8];ax[1].bar(['large audit','small certified'],vals);ax[1].axhline(9.720652007e-8,color='r',ls='--',label='location threshold');ax[1].set_yscale('log');ax[1].set(title='Kernel parameter boxes',ylabel='$\\|\\Delta A\\|$ bound');ax[1].legend(fontsize=7);ax[1].grid(axis='y',alpha=.25)
 tau=np.linspace(.5,2,200);kappa=.37;q=kappa*tau;ax[2].plot(tau,q);ax[2].set(title='$q/\\tau_*$ scaling orbit',xlabel='$\\tau_*$',ylabel='$q$');ax[2].grid(alpha=.25);fig.tight_layout();fig.savefig(FIG,dpi=180);plt.close(fig)
def main():
 figure();x={}
 x["ST1707"]=packet(1707,"constructed_generator_parameter_inference_detector_clock_robustness_lattice","Conditional OA hierarchy.",{
  "layers":["propagator","Lindblad generator","kernel parameters","minimax geometry","composite decision","detector matrix","q-clock"]})
 x["ST1708"]=packet(1708,"proven_declared_A_dephasing_composite_separation_survives_small_noncirculant_graph_perturbations","Other jumps open.",{
  "epsilon_example":8.041872774817993e-8,"residual_lower":.10789050485132269})
 x["ST1709"]=packet(1709,"kernel_parameter_propagation_has_large_base_robust_box_and_tiny_full_certificate_box","Neither box physical CI.",{
  "large_operator_bound":.00020105649070913632,"small_operator_bound":8.041872774817993e-8})
 x["ST1710"]=packet(1710,"reduced_continuous_KKT_candidate_isolated_strong_numerically_value_already_certified","Exact active-set theorem open.",{
  "gamma":16.05917256231611,"weights":[.31859567830451857,.6814043216954815],"value":.053705210722682,"certified_value_bracket":[.05,.05370588868619479]})
 x["ST1711"]=packet(1711,"finite_grid_symmetric_decision_has_two_valid_error_boundaries_and_abstention","Off-grid Q not covered.",{
  "Q_boundary":"mixture LR>=100","C_boundary":"max component LR<=.01","alpha":.01})
 x["ST1712"]=packet(1712,"arbitrary_vertex_misclassification_ball_preserves_simple_full_distribution_separation","Independent 5% detector matrices.",{
  "ideal_TV":.4112182254653598,"lower_TV":.3112182254653598})
 x["ST1713"]=packet(1713,"proven_joint_q_clock_scaling_torsor_blocks_separate_identification","Any number of physical times.",{
  "invariant":"kappa=q/tau_*","transformation":"(q,tau)->(cq,c tau)","requires_anchor":True})
 protocol={"schema":"FIN-OA-11.02-v1","parents":{"11.01":h('FIN_ST1617_ST1631_Protocol_11_01.json'),"continuous_prior":h('FIN_ST1572_ST1586_ContinuousPriorQuadrature.json'),"independent_calibration":h('FIN_ST1587_ST1601_IndependentCalibrationPolytope.json'),"detector_qclock":h('FIN_ST1692_ST1706_DetectorMatrixAndQClock.json')},"lindblad":{"gamma_max":100,"small_epsilon":8.041872774817993e-8,"residual_lower":.10789050485132269},"kernel_boxes":{"large_radii":[5e-6,5e-6,5e-5,5e-5],"small_radii":[2e-9,2e-9,2e-8,2e-8]},"minimax":{"value_bracket":[.05,.05370588868619479],"KKT_candidate":{"gamma":16.05917256231611,"weights":[.31859567830451857,.6814043216954815]}},"symmetric_grid":{"alpha":.01,"abstention":True},"qclock":{"identifiable":"q/tau_*","q_and_tau_separate":False}}
 PROTO.write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n')
 x["ST1714"]=packet(1714,"constructed_protocol_11_02_hash_freeze","Composed parent chain.",{
  "file":PROTO.name,"sha256":h(PROTO.name),"figure":str(FIG.relative_to(ROOT))})
 x["ST1715"]=packet(1715,"blocked_no_physical_parameter_CI_clock_anchor_platform_or_events","Evidence boundary.",{
  "parameter_CI":False,"clock_anchor":False,"platform":False,"events":False})
 x["ST1716"]=packet(1716,"all_new_tests_remain_internal_channel_parameter_inference_not_total_annihilation","Ontology.",{
  "literal_nonexistence":False})
 x["ST1717"]=packet(1717,"listed_remaining_mathematical_frontier","Nonreplay.",{
  "blockers":["arbitrary Lindblad jumps","physical kernel CI","exact KKT active set","continuous symmetric decision",
              "detector plus loss matrix box","independent q/tau anchor","platform law"]})
 x["ST1718"]=packet(1718,"recommended_next_programmes","Next six-round candidates.",{
  "programmes":["general GKSL jump-operator ball","interval parameter CI from kernel provenance","Krawczyk KKT isolation",
   "continuous sup-LR certification","detector-loss augmented matrix polytope","joint q/tau with second scaling observable",
   "multi-fiber q consistency","information-optimal continuous prior","platform forward model","external clock calibration","independent audit","raw events after freeze"]})
 x["ST1719"]=packet(1719,"installed_release1102_stop_rules","No scope inflation.",{
  "stop_if":["A-dephasing bound used for arbitrary jumps","small box called measured CI","KKT candidate called exact optimizer",
             "finite grid called continuum decision","TV ball called calibrated apparatus","kappa called separate q and tau"]})
 x["ST1720"]=packet(1720,"final_verdict_conditional_OA_now_has_generator_parameter_inference_and_detector_level_theorems_plus_qclock_no_go","Strict physics unchanged.",{
  "conditional_strength":True,"strict_state_category":False,"laboratory":False})
 x["ST1721"]=packet(1721,"six_round_cycle_complete_report_required","Release 11.02.",{
  "programs":90,"rounds":6,"release":"11.02","breakthrough":"Lindblad sensitivity, parameter boxes, symmetric grid decision, detector ball, q-clock no-go","strict_ToE_closure":False})
 write_round(LO,HI,x)
if __name__=="__main__":main()
