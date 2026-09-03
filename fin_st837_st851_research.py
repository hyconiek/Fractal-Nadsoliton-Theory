#!/usr/bin/env python3
"""FIN ST837--ST851: conditional nonlinear nadsoliton candidate audit."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST837_ST851_Results.json';SUM=R/'FIN_ST837_ST851_Summary.csv';names=['DNLS_StateSpace','Defocusing_EnergyBound','Focusing_NormConstrainedBound','Uniform_StationaryBranch','LocalizedBranch_ExistenceRoute','Linearized_StabilityMatrix','PumpLoss_AmplitudeBalance','SeedRequirement','GainSourceRequirement','SelectorRequirement','SolitonTerminologyGate','ContinuumScaling_Obligation','HeldOutCarrier_Transfer','StrictSource_Gate','Evidence_Gate'];N={837+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 x={}
 x['ST837']=f(837,'constructed_finite_DNLS_model','Conditional.',{'equation':'i psi_dot=A psi+g |psi|^2 psi','conserved':['norm','Hamiltonian']})
 x['ST838']=f(838,'proven_defocusing_energy_bounded_below','g supplied.',{'H':'psi*A*psi+g/2 sum|psi|^4','condition':'g>=0'})
 x['ST839']=f(839,'proven_focusing_bound_at_fixed_norm','No unconstrained bound.',{'condition':'g<0, fixed norm N0','lower_bound':'g N0^2/2 plus 0'})
 x['ST840']=f(840,'proven_uniform_stationary_branch','Not localized.',{'psi':'sqrt(N0/12) exp(-i mu t) 1','mu':'g N0/12'})
 x['ST841']=f(841,'conditional_localized_branch_route','Proof not executed.',{'method':['variational minimizer at fixed norm','implicit continuation','concentration compactness finite']})
 x['ST842']=f(842,'constructed_Bogoliubov_stability_matrix','No interval spectrum.',{'blocks':'A+nonlinear diagonal and pairing'})
 x['ST843']=f(843,'proven_conditional_pump_loss_balance','pump/loss supplied.',{'equation':'N_dot=2(P-gamma N)N','nonzero_equilibrium':'N=P/gamma'})
 x['ST844']=f(844,'proven_exact_symmetric_seed_no_go_persists','Noise/seed needed.',{'reason':'equivariant deterministic flow preserves uniform orbit'})
 x['ST845']=f(845,'blocked_no_strict_gain_coefficient_source','g/P/gamma unsourced.',{'source':False})
 x['ST846']=f(846,'blocked_no_branch_selector','D12 orbit not member.',{'QW_2191':'open'})
 x['ST847']=f(847,'soliton_name_not_yet_earned','Need traveling/localized stability.',{'requirements':['localization','orbital stability','collision/scattering or coherent persistence']})
 x['ST848']=f(848,'continuum_scaling_open','Refinement and nonlinear scaling missing.',{'needed':'g_N scaling with a_N and norm density'})
 x['ST849']=f(849,'heldout_carrier_transfer_required','No transfer.',{'carriers':['C16','C24','C48']})
 x['ST850']=f(850,'blocked_no_strict_nonlinear_source','Conditional model only.',{'strict':False})
 x['ST851']=f(851,'blocked_no_evidence','Round two local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
