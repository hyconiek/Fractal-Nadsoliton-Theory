#!/usr/bin/env python3
"""FIN ST1242--ST1256: minimal operational/state package OA."""
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, STRICT_A, write_packet, write_round


LO,HI=1242,1256
NAMES=["OA_Definition","ClassicalOA_Witness","QuantumOA_Witness","ReturnProbability_Discriminator",
 "StateCategory_Necessity","Channel_Necessity","ClockMap_Necessity","Preparation_Necessity",
 "Instrument_Necessity","Record_Necessity","CompositionEnvironment_Boundary","OperationalDiscriminant_Library",
 "Evidence_Status","OA_StrictStatus","RoundFive_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};r0=0;t=1.0;P=expm(-t*STRICT_A);U=expm(-1j*t*STRICT_A)
 x["ST1242"]=packet(1242,"constructed_minimal_OA_interface","Package is explicit and non-strict.",{
  "OA_fields":["StateCategory","PhysicalChannel","ClockMap","Preparation","Instrument","RecordRule"]})
 x["ST1243"]=packet(1243,"constructed_classical_OA_witness","One complete finite operational candidate.",{
  "StateCategory":"Delta11","PhysicalChannel":"P_t=exp(-tA)","ClockMap":"t=t_phys/tau_*",
  "Preparation":"delta_r0","Instrument":"vertex readout","RecordRule":"timestamped outcome counts"})
 x["ST1244"]=packet(1244,"constructed_quantum_OA_witness","One complete finite operational candidate.",{
  "StateCategory":"D12","PhysicalChannel":"rho->U_t rho U_t*","ClockMap":"t=t_phys/tau_*",
  "Preparation":"|r0><r0|","Instrument":"vertex PVM with declared update","RecordRule":"timestamped outcomes"})
 pc=float(P[r0,r0]);pq=float(abs(U[r0,r0])**2)
 x["ST1245"]=packet(1245,"proven_OA_candidates_are_operationally_distinguishable_in_model","No laboratory data supplied.",{
  "test":"prepare r0, evolve dimensionless time 1, measure r0 return",
  "classical":pc,"quantum":pq,"difference":abs(pc-pq)})
 x["ST1246"]=packet(1246,"proven_StateCategory_field_is_necessary","Without it simplex/ray/density domains remain.",{
  "models_remaining_without_field":3})
 x["ST1247"]=packet(1247,"proven_PhysicalChannel_field_is_necessary","A supports incompatible heat, unitary and dephasing channels.",{
  "channels_remaining_without_field":"at least continuum","same_state_carrier_can_have_multiple_channels":True})
 x["ST1248"]=packet(1248,"proven_ClockMap_field_is_necessary_for_operational_time","tau supplies a unit but not a readout protocol.",{
  "missing_without_field":["start event","elapsed-time observable","calibration record"]})
 x["ST1249"]=packet(1249,"proven_Preparation_field_is_necessary","A channel without an input ensemble has no unique experimental prediction.",{
  "possible_initial_states":"continuum","selected_by_A":False})
 x["ST1250"]=packet(1250,"proven_Instrument_field_is_necessary","Same POVM admits different postmeasurement maps.",{
  "example":["Luders","measure-and-prepare"],"same_probabilities_single_shot":True,"different_sequential_statistics":True})
 x["ST1251"]=packet(1251,"proven_RecordRule_field_is_necessary_for_falsification","Outcomes not registered under a frozen schema cannot test the model.",{
  "requires":["event fields","run/config id","timestamps","immutable analysis rule"]})
 x["ST1252"]=packet(1252,"proven_composition_environment_is_needed_beyond_single_use_but_not_minimal_for_one_finite_test","Boundary of minimal OA.",{
  "single_test_OA_complete":True,"long_history_or_open_system_requires_extension":True})
 x["ST1253"]=packet(1253,"constructed_local_discriminant_library","Predictions remain conditional.",{
  "tests":["return curve","reversibility/echo","coherence interference","heat entropy monotonicity","sequential instrument disturbance"]})
 x["ST1254"]=packet(1254,"blocked_no_external_evidence_or_independent_custody","OA witnesses are executable specifications only.",{
  "apparatus":False,"raw_events":False,"independent_holdout":False})
 x["ST1255"]=packet(1255,"OA_is_unavoidable_but_not_strict_sourced","Release 10.96 no-go prevents silent selection.",{
  "conditional_package":True,"strict_source":False,"physical_category_selected":False})
 x["ST1256"]=packet(1256,"recommended_ST1257_ST1271","Prove the three-package architecture and audit necessity/independence.",{
  "next":["W0+CA+SA+OA theorem","remove-one-package countermodels","annihilation interpretation",
          "strategic stop rules","post-verdict research programme"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
