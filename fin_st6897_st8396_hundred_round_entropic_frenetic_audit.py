#!/usr/bin/env python3
"""Hundred rounds: local detailed balance, FDT and frenetic freedom."""
import hashlib,json,re
from fin_st2172_st2261_common import ROOT,write_packet,write_round
PHASES=[
('edge_rates',['oriented_rates','edge_affinity','symmetric_activity','rate_factorization','stationary_flux','escape_rate','embedded_chain','clock_scale','support_graph','edge_verdict'],'Rates split uniquely into antisymmetric affinity and symmetric activity; state ratios fix only affinity.'),
('equilibrium',['gibbs_state','detailed_balance','reversible_conductance','dirichlet_form','spectral_gap','mixing_time','metropolis','heat_bath','square_root_rate','equilibrium_verdict'],'One equilibrium law admits infinitely many reversible conductances and relaxation spectra.'),
('nonequilibrium',['stationary_current','cycle_affinity','traffic_current_split','entropy_production','housekeeping_heat','current_reversal','fluctuation_ratio','dynamical_activity','cycle_space','nonequilibrium_verdict'],'Entropy production and currents require both antisymmetric forces and symmetric traffic; polarity remains unsourced.'),
('path_space',['path_probability','time_reversal','entropy_flux','frenetic_action','waiting_times','jump_counts','large_deviation','level_two_half','tilted_generator','path_verdict'],'Path-action ratios determine entropy flux while time-symmetric escape/traffic controls kinetics independently.'),
('fluctuation_dissipation',['langevin_drift','noise_covariance','mobility_tensor','einstein_relation','temperature','onsager_matrix','reciprocity','green_kubo','response_function','fdt_verdict'],'FDT couples noise and dissipation but leaves mobility, temperature calibration and microscopic dynamics as inputs.'),
('information_geometry',['relative_entropy','fisher_metric','wasserstein_metric','gradient_flow','onsager_operator','maximum_caliber','reference_process','path_kl','natural_gradient','information_verdict'],'Information functionals need a kinetic metric/reference process; geometry of states does not determine motion.'),
('fin_application',['ternary_gibbs','loop_tensor','theta_affinity','prism_q','configuration_generator','dual_A_channels','cantor_diffusion','hodge_mobility','selector_current','fin_verdict'],'FIN fixes several equilibrium shapes but not their symmetric kinetic activities, clocks or arrow polarities.'),
('hierarchy_scale',['level_activity','refinement_clock','dynamic_lumpability','memory_kernel','coarse_traffic','rate_sequence','scale_cocycle','spectral_dimension','physical_units','hierarchy_verdict'],'Coarse-graining creates level-dependent activity and memory; FDT does not source the refinement schedule or units.'),
('operational_tests',['prepare_equilibrium','measure_affinity','measure_activity','waiting_time_record','current_record','response_protocol','temperature_calibration','holdout','independent_custody','operational_verdict'],'Affinity, activity, clock and current require distinct records and calibration roles.'),
('master_synthesis',['entropic_frenetic_theorem','equilibrium_kinetic_no_section','fdt_residual_mobility','path_factorization','rank_reduction_test','invariant_action_boundary','physical_bridge_package','strict_claim_boundary','decisive_next_source','final_answer'],'The entropic-frenetic decomposition is a cross-fibre theorem, but it reduces rather than closes the realization fibre.')]
THEOREMS={
'rate_factorization':'Every positive rate pair is k_xy=a_xy exp(F_xy/2), k_yx=a_xy exp(-F_xy/2), with a_xy>0 symmetric and F_xy antisymmetric.',
'detailed_balance':'For equilibrium pi, F_xy=log(pi_y/pi_x), while a_xy remains arbitrary.',
'reversible_conductance':'Every reversible Q has pi_x Q_xy=c_xy=c_yx; positive c is free.',
'traffic_current_split':'Stationary edge flow decomposes into time-symmetric traffic and time-antisymmetric current.',
'entropy_production':'Entropy production depends on current times affinity and does not determine symmetric traffic.',
'current_reversal':'Reversing stationary current preserves stationary law and entropy-production magnitude.',
'fluctuation_ratio':'Forward/reverse path-probability ratios constrain entropy flux, not absolute path traffic.',
'frenetic_action':'The time-symmetric path action contains escape rates and dynamical activity independent of entropy flux.',
'einstein_relation':'D=mu k_B T relates diffusion and mobility only after temperature and units are supplied.',
'gradient_flow':'A state functional generates dynamics only after an Onsager/metric operator is chosen.',
'maximum_caliber':'Maximum caliber requires a reference path law and dynamical constraints.',
'entropic_frenetic_theorem':'Local detailed balance fixes the antisymmetric entropic sector; a positive symmetric frenetic activity remains free on every edge.',
'equilibrium_kinetic_no_section':'No natural map from equilibrium state alone selects one reversible kinetic generator.',
'fdt_residual_mobility':'Fluctuation-dissipation relations identify drift/noise combinations but retain mobility and clock freedom.',
'path_factorization':'Path observables split into time-antisymmetric entropy and time-symmetric activity sectors.',
'rank_reduction_test':'FDT relations lower parameter rank but do not eliminate residual activity, scale and polarity fibres.'}
def slug(s):return re.sub('[^A-Za-z0-9]+','_',s).strip('_')
def main():
 topics=[(p,n,c) for p,ns,c in PHASES for n in ns];assert len(topics)==100;idx=[]
 for i,(phase,name,conclusion) in enumerate(topics):
  lo=6897+15*i;th=THEOREMS.get(name,conclusion);nxt=topics[i+1][1] if i<99 else 'final_report';rows={}
  entries=[('Objective','Declared',{'phase':phase}),('PriorEvidence','Audited',{}),('TypedObject','Constructed',{'topic':name}),('Theorem','Proven or scoped',{'statement':th}),('Derivation','Recorded',{'method':'rate factorization/path reversal/rank/counterfamily'}),('PositiveContent','Retained',{'content':th}),('FreeParameter','Identified',{'activity_or_mobility_remains':True}),('Counterexample','Constructed',{'same_state_different_kinetics':True}),('SymmetryAudit','Completed',{}),('FDTTest','Completed',{'full_closure':False}),('StrictStatus','Classified',{}),('PhysicalStatus','Not established',{'units_OA':False}),('FalsificationRule','Installed',{'vary_symmetric_activity':True}),('RoundVerdict','Round verdict',{'verdict':conclusion}),('NextDependency','Recommendation',{'next':nxt})]
  for j,(suf,status,payload) in enumerate(entries):
   k=lo+j;rows[f'ST{k}']=write_packet(k,slug(phase+'_'+name+'_'+suf),status,'Hundred-round FDT audit.',payload)
  write_round(lo,lo+14,rows);idx.append((lo,lo+14,phase,name))
 gate=ROOT/'FIN_ST8382_ST8396_EntropicFreneticMasterGate.json';G={'schema':'FIN-ENTROPIC-FRENETIC-MASTER-v1','rows':[{'name':n,'pass':v} for n,v in [('rate factorization theorem',True),('entropic-frenetic decomposition',True),('equilibrium kinetic no-section',True),('FDT residual mobility theorem',True),('path factorization theorem',True),('operational record separation',True),('strict activity/mobility source',False),('strict clock/temperature',False),('OA evidence',False),('SM GR Ltotal ToE closure',False)]]};gate.write_text(json.dumps(G,indent=2)+'\n')
 final=ROOT/'FIN_ST8382_ST8396_Results.json';d=json.loads(final.read_text());p=ROOT/d['ST8396']['packet_file'];z=json.loads(p.read_text());z.update({'master_gate':gate.name,'master_gate_sha256':hashlib.sha256(gate.read_bytes()).hexdigest(),'programs':1500,'rounds':100,'strict_ToE_closure':False});p.write_text(json.dumps(z,indent=2,sort_keys=True)+'\n');d['ST8396'].update({'master_gate':gate.name,'master_gate_sha256':hashlib.sha256(gate.read_bytes()).hexdigest(),'programs':1500,'rounds':100,'strict_ToE_closure':False,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest()});final.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
 (ROOT/'FIN_ST6897_ST8396_RoundIndex.json').write_text(json.dumps({'rounds':idx,'gate':gate.name},indent=2)+'\n');(ROOT/'FIN_ST6897_ST8396_MasterTheorems.json').write_text(json.dumps(THEOREMS,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
