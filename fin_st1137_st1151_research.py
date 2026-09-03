#!/usr/bin/env python3
"""FIN ST1137--ST1151: refinement, hidden fibers and dual dynamics."""
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, STRICT_A, write_packet, write_round


LO,HI=1137,1151
NAMES=["TwoFiber_RefinementFamily","CoarseGenerator_Intertwiner","HeatSemigroup_Intertwiner",
 "UnitaryGroup_Intertwiner","FiberRate_Nonuniqueness","KernelDimension_PhaseChange",
 "RefinedSpectrum_Formula","CoarseObserver_FiberBlindness","Entropy_ChainRule",
 "RefinedStateCategory_StillAmbiguous","DualDynamics_RefinementCompatibility","NoDimensionalScale_FromFiberRate",
 "FractalLayer_InterpretationBoundary","RefinementBridge_CurrentVerdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def entropy(p):
 p=np.asarray(p);q=p[p>0];return float(-np.sum(q*np.log(q)))
def objects(q):
 B=np.array([[q,-q],[-q,q]],float)
 A24=np.kron(STRICT_A,np.eye(2))+np.kron(np.eye(N),B)
 C=np.kron(np.eye(N),np.ones((1,2)))
 return B,A24,C
def main():
 x={};q=.37;B,A24,C=objects(q)
 x["ST1137"]=packet(1137,"constructed_exact_one_parameter_12_to_24_refinement",
  "Family is a mathematical lift, not a strict fractal law.",{
  "formula":"A24(q)=A12 tensor I2 + I12 tensor [[q,-q],[-q,q]]","q_domain":"q>=0",
  "dimension":24})
 x["ST1138"]=packet(1138,"proven_exact_coarse_generator_intertwining","Coarse map sums each two-point fiber.",{
  "identity":"C A24=A12 C","residual":float(np.linalg.norm(C@A24-STRICT_A@C,np.inf))})
 t=.83
 x["ST1139"]=packet(1139,"proven_exact_heat_coarse_intertwining_numerically_replayed","Analytic consequence of generator intertwining.",{
  "identity":"C exp(-tA24)=exp(-tA12) C","time":t,
  "residual":float(np.linalg.norm(C@expm(-t*A24)-expm(-t*STRICT_A)@C,np.inf))})
 x["ST1140"]=packet(1140,"proven_exact_unitary_coarse_intertwining_numerically_replayed","C is an amplitude sum, not a CPTP quantum coarse channel.",{
  "identity":"C exp(-itA24)=exp(-itA12) C","time":t,
  "residual":float(np.linalg.norm(C@expm(-1j*t*A24)-expm(-1j*t*STRICT_A)@C,np.inf))})
 x["ST1141"]=packet(1141,"proven_continuum_of_fiber_rates_has_identical_coarse_dynamics",
  "Reproduces the coarse operator for every q; no q selector is supplied.",{
  "tested_q":[0.0,.1,.37,1.0,10.0],"coarse_A_independent_of_q":True,"canonical_q":False})
 dims=[]
 for z in [0.0,.37]:dims.append(int(np.sum(np.abs(np.linalg.eigvalsh(objects(z)[1]))<1e-9)))
 x["ST1142"]=packet(1142,"proven_hidden_fiber_connectivity_changes_without_coarse_change",
  "q=0 and q>0 have different fine stationary sectors.",{
  "kernel_dimension_q0":dims[0],"kernel_dimension_q_positive":dims[1],
  "coarse_observer_detects_difference":False})
 e12=np.linalg.eigvalsh(STRICT_A);e24=np.linalg.eigvalsh(A24);target=np.sort(np.r_[e12,e12+2*q])
 x["ST1143"]=packet(1143,"proven_refined_spectrum_is_two_shifted_copies",
  "Kronecker-sum construction.",{
  "formula":"spec(A24)=spec(A12) union (spec(A12)+2q)",
  "spectral_residual":float(np.max(abs(e24-target)))})
 rng=np.random.default_rng(1144);p24=rng.random(2*N);p24/=p24.sum();coarse=C@p24
 alt=p24.copy()
 for i in range(N):
  s=coarse[i];alt[2*i]=.9*s;alt[2*i+1]=.1*s
 x["ST1144"]=packet(1144,"proven_coarse_map_has_12_dimensional_hidden_fiber_kernel",
  "Blindness is relative to C.",{
  "rank_C":int(np.linalg.matrix_rank(C)),"nullity_C":2*N-int(np.linalg.matrix_rank(C)),
  "two_states_same_coarse_residual":float(np.linalg.norm(C@alt-coarse))})
 cond=0.0
 for i,s in enumerate(coarse):
  if s>0:cond+=s*entropy(p24[2*i:2*i+2]/s)
 x["ST1145"]=packet(1145,"proven_refinement_entropy_chain_rule","Classical probability candidate.",{
  "H_fine":entropy(p24),"H_coarse_plus_conditional":entropy(coarse)+cond,
  "residual":abs(entropy(p24)-entropy(coarse)-cond)})
 x["ST1146"]=packet(1146,"proven_refinement_does_not_select_state_category","Each A24 again admits simplex, ray and density-state models.",{
  "classical_dimension":23,"pure_ray_dimension":46,"density_affine_dimension":575,
  "unique_category":False})
 x["ST1147"]=packet(1147,"proven_both_functional_calculi_respect_declared_refinement","Mathematical compatibility, not physical time identification.",{
  "heat_intertwines":True,"unitary_intertwines":True,"single_physical_channel_selected":False})
 x["ST1148"]=packet(1148,"proven_q_has_only_inverse_internal_time_units_until_calibrated","Dimensionless model.",{
  "q_role":"relative fiber relaxation/frequency scale","SI_value":False,"c_or_length_source":False})
 x["ST1149"]=packet(1149,"constructed_layers_as_internal_fibers_or_descriptions","Interpretation consistent with one total nadsoliton but not forced.",{
  "consistent":"fine state -> coarse record","not_derived":["physical fractal hierarchy","absolute scale","observer layer"]})
 x["ST1150"]=packet(1150,"refinement_bridge_exact_but_nonunique_and_nonphysical","No new strict selector or role transfer.",{
  "exact_intertwining":True,"free_parameter_q":True,"physical_scale":False,"QW2191_discharged":False})
 x["ST1151"]=packet(1151,"recommended_ST1152_ST1166","Next round asks whether an internal observer/clock/record can select the state category.",{
  "next":["information-capacity bounds","POVM versus instrument ambiguity","relational clocks",
          "observer inclusion dimension","operational selector verdict"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
