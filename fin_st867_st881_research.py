#!/usr/bin/env python3
"""FIN ST867--ST881: nonproxy metric-response gravity intake."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST867_ST881_Results.json';SUM=R/'FIN_ST867_ST881_Summary.csv';names=['MetricFamily_Type','LocalWeightMetricLaw','DirichletMetricVariation','StressTensor_Candidate','DiscreteConservation_Test','BianchiI_Anisotropic_Test','FRW_Baseline_Test','EinsteinHilbert_Comparison','NewtonAction_UnitGate','EquivalencePrinciple_Gate','HeldOutGeometry_Test','LegacyGravityRole_Stop','MetricResponse_AcceptanceMatrix','StrictGravity_Gate','Evidence_Gate'];N={867+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 x={}
 x['ST867']=f(867,'metric_family_not_sourced','Need g carrier.',{'required':'g->A(g)'})
 x['ST868']=f(868,'constructed_conditional_local_weight_law','Not strict.',{'example':'w_xy(g)=w0 exp(-distance_g/ell)'})
 x['ST869']=f(869,'proven_variation_exists_for_supplied_law','Law imported.',{'deltaE':'1/2 sum delta w |fx-fy|^2'})
 x['ST870']=f(870,'constructed_formal_stress_candidate','No spacetime tensor.',{'T':'metric derivative of action'})
 x['ST871']=f(871,'conservation_not_exported','Ward/diffeomorphism structure missing.',{'pass':False})
 x['ST872']=f(872,'BianchiI_test_blocked_no_anisotropic_source','Matches P2687.',{'source':False})
 x['ST873']=f(873,'Minkowski_zero_only_FRW_baseline','Nonflat needs source.',{'closure':False})
 x['ST874']=f(874,'EH_comparison_receiver_only','No nonproxy equality.',{'variation_bundle':False})
 x['ST875']=f(875,'blocked_no_Newton_action_unit','Scale torsor.',{'G_source':False,'hbar_action':False})
 x['ST876']=f(876,'blocked_no_equivalence_principle_derivation','Matter coupling absent.',{'pass':False})
 x['ST877']=f(877,'heldout_geometry_test_not_reached','No model.',{'geometry':'anisotropic unseen'})
 x['ST878']=f(878,'legacy_gravity_role_transfer_forbidden','No bridge.',{'beta_power':False})
 x['ST879']=f(879,'acceptance_matrix_zero_rows','No complete package.',{'accepted':0})
 x['ST880']=f(880,'blocked_no_strict_gravity_bridge','A8 remains partial/frozen.',{'closure':False})
 x['ST881']=f(881,'blocked_no_evidence','Round four local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
