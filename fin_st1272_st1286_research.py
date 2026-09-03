#!/usr/bin/env python3
"""FIN ST1272--ST1286: freeze classical and quantum OA hypotheses."""
import hashlib
import json

import numpy as np

from fin_oa_discrimination_common import EIG, classical_distribution, classical_return, quantum_distribution, quantum_return
from fin_total_nadsoliton_common import N, STRICT_A, write_packet, write_round


LO,HI=1272,1286
NAMES=["OA_ComparisonScope","FrozenClassicalOA","FrozenQuantumOA","CommonPreparation","CommonVertexInstrument",
 "ReturnBernoulliObservable","FullVertexObservable","SpectralReturnFormulas","TranslationEquivariance",
 "ProbabilityNormalization","MatrixFingerprint","TimeDomainFreeze","NoNuisanceBaseHypotheses","BlindRoleBoundary","RoundOne_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 x={};source=0
 x["ST1272"]=packet(1272,"constructed_two_hypothesis_OA_scope","Conditional model comparison, not nature data.",{
  "H_C":"classical heat on Delta11","H_Q":"closed unitary on D12","shared_A":True,"source_vertex":source})
 x["ST1273"]=packet(1273,"frozen_classical_OA_hypothesis","CA and SA remain conditional.",{
  "state":"p in Delta11","channel":"p(t)=exp(-tA)p0","clock":"t=t_phys/tau_*","measurement":"vertex return"})
 x["ST1274"]=packet(1274,"frozen_quantum_OA_hypothesis","No dephasing in the base alternative.",{
  "state":"rho in D12","channel":"rho(t)=U_t rho0 U_t*","U_t":"exp(-itA)","measurement":"vertex PVM"})
 x["ST1275"]=packet(1275,"constructed_common_selected_vertex_preparation","Preparation implementation is not laboratory-sourced.",{
  "classical":"delta_0","quantum":"|0><0|","SA_frame":[0,1]})
 x["ST1276"]=packet(1276,"constructed_common_vertex_effect","Same outcome label, category-specific state update.",{
  "effect":"E0=|0><0|","outcomes":["return","not_return"],"instrument_update_frozen":False})
 x["ST1277"]=packet(1277,"proven_return_measurement_is_Bernoulli_under_each_frozen_model","Finite categorical result.",{
  "p_C(t)":"(exp(-tA))_00","p_Q(t)":"|(exp(-itA))_00|^2"})
 t=1.0;pc=classical_distribution(t);pq=quantum_distribution(t)
 x["ST1278"]=packet(1278,"constructed_full_vertex_multinomial_observable","Return statistic is a coarse binary reduction.",{
  "classical_probabilities_t1":pc.tolist(),"quantum_probabilities_t1":pq.tolist(),
  "TV_t1":float(0.5*np.sum(abs(pc-pq)))})
 x["ST1279"]=packet(1279,"proven_circulant_spectral_return_formulas","All Fourier modes have equal source weight 1/12.",{
  "classical":"C(t)=12^-1 sum_k exp(-lambda_k t)",
  "quantum":"Q(t)=|12^-1 sum_k exp(-i lambda_k t)|^2"})
 x["ST1280"]=packet(1280,"proven_return_curves_are_source_vertex_independent","Circulant translation symmetry only.",{
  "vertices":N,"same_diagonal_all_vertices":True,"selects_absolute_origin":False})
 x["ST1281"]=packet(1281,"proven_both_frozen_output_distributions_are_normalized","Numerical replay at t=1.",{
  "classical_sum":float(pc.sum()),"quantum_sum":float(pq.sum()),"classical_min":float(pc.min()),"quantum_min":float(pq.min())})
 matrix_text=json.dumps(np.round(STRICT_A,15).tolist(),separators=(',',':'))
 x["ST1282"]=packet(1282,"constructed_frozen_matrix_fingerprint","Fingerprint is reproducibility metadata, not exact symbolic provenance.",{
  "rounding_decimals":15,"sha256":hashlib.sha256(matrix_text.encode()).hexdigest(),"dimension":N})
 x["ST1283"]=packet(1283,"frozen_design_time_domain","Chosen before synthetic or future event records.",{
  "domain":[0.0,10.0],"optimization_target":"absolute return-probability separation"})
 x["ST1284"]=packet(1284,"base_hypotheses_exclude_unfrozen_nuisance_parameters","Robust alternatives are deferred to round four.",{
  "dephasing":0.0,"clock_scale":1.0,"loss":0.0})
 x["ST1285"]=packet(1285,"blocked_independent_provider_registrar_analyst_roles","Local generation cannot satisfy external custody.",{
  "provider_ne_registrar_ne_analyst":False,"raw_nature_events":False})
 x["ST1286"]=packet(1286,"recommended_ST1287_ST1301","Derive curves and optimize discrimination times.",{
  "next":["spectral multiplicity reduction","stationary derivative","global grid/Lipschitz audit","crossings","practical time freeze"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
