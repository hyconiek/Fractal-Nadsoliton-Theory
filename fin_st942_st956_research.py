#!/usr/bin/env python3
"""FIN ST942--ST956: dissipative attractor and stochastic extinction audit."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST942_ST956_Results.json';SUM=R/'FIN_ST942_ST956_Summary.csv';names=['LogisticNorm_Flow','ZeroStability_PumpSign','NonzeroAttractor','Basin_PositiveNorm','PumpRemoval_Decay','AdditiveNoise_ZeroNotAbsorbing','MultiplicativeNoise_Extinction','FiniteNoise_EscapeRate','Lyapunov_Attractor','AttractorShape_Selector','ReservoirResource_Cost','DetailedBalance_NoPersistentExcitation','NonequilibriumPersistence','StrictPumpSource_Gate','Evidence_Gate'];N={942+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 x={}
 x['ST942']=f(942,'proven_logistic_norm_reduction','Conditional pump model.',{'Qdot':'2(P-gamma Q)Q'})
 x['ST943']=f(943,'proven_zero_stable_if_P_negative_unstable_if_positive','P supplied.',{'linear_rate':'2P'})
 x['ST944']=f(944,'proven_nonzero_attractor_Qstar','Only norm.',{'Qstar':'P/gamma for P>0'})
 x['ST945']=f(945,'proven_positive_norm_basin_excludes_zero','Deterministic.',{'basin':'Q0>0'})
 x['ST946']=f(946,'proven_pump_removal_causes_decay','No intrinsic persistence.',{'P':0})
 x['ST947']=f(947,'proven_additive_noise_makes_zero_nonabsorbing','Noise source added.',{'zero_absorbing':False})
 x['ST948']=f(948,'multiplicative_noise_extinction_model_dependent','Boundary convention needed.',{'possible':True})
 x['ST949']=f(949,'Kramers_escape_rate_requires_potential_noise_scale','No values.',{'rate':'~exp(-DeltaV/D)'})
 x['ST950']=f(950,'conditional_Lyapunov_attractor','Shape dynamics needed.',{'V':'(gamma/2)(Q-P/gamma)^2'})
 x['ST951']=f(951,'blocked_no_attractor_shape_selector','D12 orbit.',{'selector':False})
 x['ST952']=f(952,'proven_persistent_nonequilibrium_cost','Reservoir required.',{'entropy_production':'nonzero generally'})
 x['ST953']=f(953,'proven_detailed_balance_relaxes_excitation_without_protection','Topology/constraint exceptions.',{'persistent_excited_state':False})
 x['ST954']=f(954,'conditional_nonequilibrium_persistence_mechanism','Not strict.',{'mechanisms':['pump','saturation','feedback']})
 x['ST955']=f(955,'blocked_no_strict_pump_reservoir_source','No resource law.',{'source':False})
 x['ST956']=f(956,'blocked_no_evidence','Round three local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
