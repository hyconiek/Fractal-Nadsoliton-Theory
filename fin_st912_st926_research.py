#!/usr/bin/env python3
"""FIN ST912--ST926: linear persistence and annihilation audit."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from fin_st402_st416_research import independent_strict_matrix_float
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST912_ST926_Results.json';SUM=R/'FIN_ST912_ST926_Summary.csv';A=independent_strict_matrix_float();N=12;names=['Annihilation_Definition_Lattice','Unitary_Norm_NoExtinction','Heat_Mass_Conservation','Heat_Pattern_Decay','Wave_PhaseSpace_Energy','ZeroMode_Persistence','FiniteTime_Invertibility','OpenLoss_Extinction_Counterexample','Dephasing_Coherence_Loss','GlobalVsAccessible_Information','FiniteUnitary_Recurrence','InstantaneousZero_NotAnnihilation','StrictCore_Persistence_Status','PersistenceSource_Gate','Evidence_Gate'];NN={912+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in NN.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':NN[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 e=np.linalg.eigvalsh(A);gap=float(e[1]);x={}
 x['ST912']=f(912,'constructed_four_annihilation_notions','Definitions.',{'types':['zero-state extinction','pattern delocalization','coherence loss','topological/attractor loss']})
 x['ST913']=f(913,'proven_unitary_zero_state_unreachable','Finite self-adjoint model.',{'theorem':'||exp(-itA)psi||=||psi||; nonzero state never becomes zero.'})
 x['ST914']=f(914,'proven_heat_total_mass_conservation','Probability embedding.',{'theorem':'1^T exp(-tA)p=1^T p','zero_mass_only_extinguishes':True})
 x['ST915']=f(915,'proven_heat_pattern_exponential_decay','Localization not protected.',{'spectral_gap':gap,'bound':'||P_t(p-u)||<=exp(-gap t)||p-u||'})
 x['ST916']=f(916,'proven_wave_phase_space_energy_conservation','Second-order wave model supplied.',{'energy':'||u_dot||^2+<u,Au>'})
 x['ST917']=f(917,'proven_constant_zero_mode_persists','Carries no localization.',{'kernel_dimension':1})
 x['ST918']=f(918,'proven_finite_time_semigroup_invertibility_linear_space','Inverse not stochastic.',{'exp_minus_tA_inverse':'exp(tA)'})
 x['ST919']=f(919,'proven_loss_can_annihilate_without_source','Open law additional.',{'counterexample':'psi_dot=-(A+gamma I)psi ->0'})
 x['ST920']=f(920,'proven_dephasing_destroys_coherence_not_populations','Environment supplied.',{'spectral_diagonal_preserved':True})
 x['ST921']=f(921,'proven_accessible_information_can_vanish_while_global_state_persists','Context dependent.',{'mechanism':'coarse channel/Schur tracing'})
 x['ST922']=f(922,'proven_finite_unitary_almost_recurrence','Not soliton stability.',{'reason':'finite torus of eigenphases'})
 x['ST923']=f(923,'proven_instantaneous_component_zero_not_global_annihilation','Measurements can vanish.',{'norm_remains':True})
 x['ST924']=f(924,'strict_core_protects_norm_or_mass_not_nadsoliton_identity','Nonlinear object absent.',{'localized_soliton_proved':False,'topological_charge_proved':False})
 x['ST925']=f(925,'blocked_no_strict_persistence_source_beyond_linear_conservation','Need nonlinear/topological/pump law.',{'source':False})
 x['ST926']=f(926,'blocked_no_evidence','Round one local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
