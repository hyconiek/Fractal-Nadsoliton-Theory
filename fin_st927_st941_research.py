#!/usr/bin/env python3
"""FIN ST927--ST941: charge, winding and nonlinear persistence audit."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST927_ST941_Results.json';SUM=R/'FIN_ST927_ST941_Summary.csv';names=['U1_NoetherCharge','ChargeBlocksZeroState','DiscreteWinding_Definition','WindingHomotopyProtection','AmplitudeZero_UnwindingGate','LinearUnitary_WindingNotProtected','DNLS_ChargeEnergy','DiscreteBreather_Route','OrbitalStability_Obligation','TopologicalCharge_NotExported','PumpedAttractor_Persistence','AttractorNotConservation','CollisionAnnihilation_Question','StrictProtection_Gate','Evidence_Gate'];N={927+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 x={}
 x['ST927']=f(927,'proven_U1_charge_for_phase_invariant_Hamiltonian','Model conditional.',{'Q':'sum |psi_x|^2'})
 x['ST928']=f(928,'proven_positive_charge_blocks_global_zero','Does not protect shape.',{'Q>0':'psi cannot be identically zero'})
 x['ST929']=f(929,'constructed_discrete_winding','Branch convention needed.',{'W':'(1/2pi) sum principal phase increments'})
 x['ST930']=f(930,'proven_winding_constant_under_nonzero_continuous_homotopy','Requires gap from zero/cut.',{'integer':True})
 x['ST931']=f(931,'proven_unwinding_requires_amplitude_zero_or_phasecut_event','Finite lattice.',{'gate':'min |psi_x|=0 or increment hits branch cut'})
 x['ST932']=f(932,'linear_unitary_does_not_guarantee_winding_protection','Components may cross zero.',{'protected':False})
 x['ST933']=f(933,'proven_DNLS_norm_energy_conservation','g supplied.',{'invariants':['Q','H']})
 x['ST934']=f(934,'conditional_discrete_breather_existence_route','No proof for FIN parameters.',{'methods':['anti-continuum continuation','variational constrained minimizer']})
 x['ST935']=f(935,'blocked_no_orbital_stability_theorem','Need spectral/GSS analysis.',{'pass':False})
 x['ST936']=f(936,'blocked_no_strict_topological_charge','No sourced phase field/nonzero constraint.',{'charge':False})
 x['ST937']=f(937,'proven_pump_loss_nonzero_attractor_conditionally','Pump supplied.',{'Qstar':'P/gamma'})
 x['ST938']=f(938,'proven_attractor_persistence_not_conservation','Pump removal permits decay.',{'external_resource':True})
 x['ST939']=f(939,'collision_annihilation_not_defined_current_model','Need multi-soliton/scattering law.',{'open':True})
 x['ST940']=f(940,'blocked_no_strict_identity_protection','Only norm/mass conditional mechanisms.',{'nadsoliton_identity':False})
 x['ST941']=f(941,'blocked_no_evidence','Round two local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
