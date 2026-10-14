#!/usr/bin/env python3
"""RG-E01..05: bounded finite computations, produces only designs/model evidence."""
import numpy as np, json, math, hashlib
from pathlib import Path
from scipy.sparse.linalg import expm_multiply
from source_model import *
HERE=Path(__file__).resolve().parent
MODELS=all_models()
KINETICS=('heat_bath','metropolis','barker')
def out(rel,obj):
    p=HERE/rel;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
def stat(A,g):
    z=pair_difference(A,g)
    return [float(x) for x in z]

def theoretical_two_time(A,g=G_N2,t1=.4,t2=1.):
    """Prepared at independent uniform g0 equilibrium then switched to g, ideal nondisturbing readout.
    13x13 joint with explicit invalid outcome; no clipping and no fitted rate.
    """
    outcome=[];checks=[]
    d=np.fromfunction(lambda i,j: (j-i)%Q,(Q,Q),dtype=int).astype(int).ravel()
    p0=np.ones(Q*Q)/(Q*Q)
    for kinetic in KINETICS:
        L=build_n2_generator(A,g,kinetic)
        p1=expm_multiply(L.T*t1,p0)
        inp=np.zeros((Q*Q,Q))
        for j in range(Q): inp[:,j]=p1*(d==j)
        p2=expm_multiply(L.T*(t2-t1),inp)
        joint=np.zeros((Q+1,Q+1))
        for k in range(Q):
          for j in range(Q):joint[k,j]=p2[d==j,k].sum()
        pstat=pi_n2(A,g)
        res=np.asarray(pstat@L).ravel()
        rev=(pstat[:,None]*L.toarray());db=np.max(np.abs(rev-rev.T))
        g0=build_n2_generator(A,0,kinetic)
        checks.append({'kinetic':kinetic,'max_rowsum_abs':float(np.max(np.abs(np.asarray(L.sum(axis=1)).ravel()))),
                       'max_stationary_residual':float(np.max(np.abs(res))),
                       'max_detailed_balance_residual':float(db),'g0_calibration_same_heatbath':bool(np.max(np.abs((g0-build_n2_generator(A,0,'heat_bath')).toarray()))<1e-14),
                       'min_joint_probability':float(joint.min()),'total_joint_mass':float(joint.sum()),
                       'final_invalid_probability':float(joint[-1,:].sum()+joint[:,-1].sum())})
        outcome.append({'kinetic':kinetic,'joint_13x13':joint.tolist(), 'first_marginal_12':[float(x) for x in joint[:12,:].sum(axis=1)],
                        'final_marginal_12':[float(x) for x in joint[:,:12].sum(axis=0)]})
    by={d['kinetic']:np.array(d['joint_13x13']) for d in outcome}
    sep={f'{x}_vs_{y}':tv(by[x],by[y]) for i,x in enumerate(KINETICS) for y in KINETICS[i+1:]}
    return outcome,checks,sep

def main():
    models=MODELS
    ids=['FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE']
    for name,A in models.items():
      assert A.shape==(12,12) and np.max(np.abs(A-A.T))<1e-12
    source={'date':'2026-10-14','status':'NO_NEW_PHYSICAL_INPUT',
      'model_supplied':{'q':12,'M3_g':G_M3,'N2_g':G_N2,'A7':'rank-seven Fourier projection from supplied exact-decimal strict kernel',
        'kinetics_M3':'leave-one-out heat bath rate 1 per copy, supplied law'},
      'capabilities':{'M3_ordered_first_two':'mathematical simulator only; no demonstrated physical preparation or addressable hidden copy',
                      'N2_difference':'PHYS-014 design-only, two q12 nodes, g3, difference-label readout; physical apparatus/records not provided'},
      'physical_records_received':False,'independent_A7_source_received':False,
      'clock_units':'not established by numerical g0 control',
      'readout':'ideal diagonal model maps; real confusion, latency, disturbance, invalid and drift unmeasured'}
    out('RG-E01/OPERATIONAL_CONTRACT.json',source)
    # calculate exact formula predictions and obstructions
    shape={}; probs={}; contrast={}
    for name,A in models.items():
      probs[name]={'N2_diff_g3':stat(A,G_N2), 'M3_ordered_piY_g5':m3_visible(A,G_M3).tolist()}
      shape[name]=static_mode_shape(A)
      contrast[name]={}
      for j,l,b,bp,c in [(0,1,0,1,0),(2,5,0,3,7),(0,4,3,5,0)]:
        r0=conditional_R(A,b,c);r1=conditional_R(A,bp,c)
        lhs=math.log(r0[j]/r0[l])-math.log(r1[j]/r1[l]);rhs=(G_M3/3)*(A[j,b]-A[l,b]-A[j,bp]+A[l,bp])
        contrast[name][f'{j}_{l}_{b}_{bp}_{c}']={'difference_log_odds':float(lhs),'matrix_contrast':float(rhs),'identity_error':float(lhs-rhs)}
    out('RG-E02/SHAPE_PREDICTIONS.json',{'shape':shape,'conditional_contrasts':contrast,'observable_laws':probs,'input_scope':'nominal binary64 supplied kernel'} )
    # static N2 TV and certified conditional minimax lower bounds in a calibrated *class* via TV triangle
    n2={x:np.array(probs[x]['N2_diff_g3']) for x in ids}
    m3={x:np.array(probs[x]['M3_ordered_piY_g5']) for x in ids}
    dist_N2={f'{x}_vs_{y}':tv(n2[x],n2[y]) for i,x in enumerate(ids) for y in ids[i+1:]}
    dist_M3={f'{x}_vs_{y}':tv(m3[x],m3[y]) for i,x in enumerate(ids) for y in ids[i+1:]}
    # numeric input perturbation already declared trace preserving, not an alternative independent source
    pert_tv=tv(n2['FIN_A7'],np.array(probs['PERT_34_10P']['N2_diff_g3']))
    gain_fraction=.01; instrument_error=.005
    uncertainty_by_model={}
    for x in ids:
        A=models[x];w=float(np.ptp(A[0]))
        bg=math.tanh((G_N2*gain_fraction*w)/8)
        uncertainty_by_model[x]={'N2_logweight_oscillation_gain_1pct':float(G_N2*gain_fraction*w/2),
         'gain_tv_cap_from_tanh':float(bg),'hypothetical_detector_tv_cap':instrument_error,
         'total_hypothetical_model_radius':float(bg+instrument_error)}
    pair_bounds={}
    for i,x in enumerate(ids):
      for y in ids[i+1:]:
        key=f'{x}_vs_{y}'; d=dist_N2[key]
        r=uncertainty_by_model[x]['total_hypothetical_model_radius']+uncertainty_by_model[y]['total_hypothetical_model_radius']
        pair_bounds[key]={'nominal_tv':d,'lower_bound_if_calibrated':max(0.,d-r-1e-12),
          'two_sided_model_uncertainty_budget':r,'max_common_tv_model_error_for_separation':d/2,
          'hypothetical_gate':'positive' if d>r+1e-12 else 'NOT_GUARANTEED',
          'actual_calibration_gate':'OPEN_INPUT'}
    discr={'N2_pairwise_nominal_tv':dist_N2,'M3_ordered_pairwise_nominal_tv':dist_M3,
        'FIN_vs_trace_perturbation_N2_tv':pert_tv,'hypothetical_design_only_uncertainties':{'fractional_gain':gain_fraction,'detector_TV':instrument_error},
        'hypothetical_by_model':uncertainty_by_model,'N2_nominal_plus_hypothetical_robust_separation':pair_bounds,
        'unconditional_minimax_if_arbitrary_unknown_readout':0,
        'explanation':'A common detector channel K maps all input states to an identical output, causing TV=0. A numerical positive minimax lower bound is prohibited without channel calibration.',
        'open_input':'instrument confusion/invalid probabilities, true gain uncertainty, source independence, sampling precision'}
    out('RG-E03/SEPARATION_BOUNDS.json',discr)
    out('RG-E03/DISCRIMINATION_CONTRACT.json',{'selected_candidate':'N2_difference_distribution_g3',
       'models':ids,'same_gain_trace_and_readout':True,'allowable_calibration_set':'HYPOTHETICAL_ONLY',
       'actual_minimax_result':'NO_SEPARATION_WITH_UNRESTRICTED_READOUT',
       'guaranteed_lower_bound_given_explicit_TV_error':'max(0,TV_nominal - eps_model1 - eps_model2)',
       'unrestricted_rank7_alternatives':'no positive uniform separation from arbitrarily close alternatives'})
    # dynamic predictions: N2 PHYS014 topology, distinct from earlier M3 ROM
    multi={};diag={};kinetic_sep={}
    for x in ids:
        outs,checks,seps=theoretical_two_time(models[x]);multi[x]=outs;diag[x]=checks;kinetic_sep[x]=seps
    out('RG-E04/KINETIC_CONTRACTS.json',{'N2_scope':'two labelled q12 nodes, distinct from M3',
      'kinetic_definitions':{'heat_bath':'rate copy update 1; categorical Gibbs refresh with self draws',
           'metropolis':'proposal among 11 distinct labels; attempt 11/12 per copy; Metropolis accept',
           'barker':'proposal among 11 distinct labels; attempt 22/12 per copy; Barker accept'},
      'g0_rate_calibration':'all g0 generators equal and nonzero single-copy Fourier autocorr exp(-t); t dimensionless',
      'g0_physical_clock':'OPEN_INPUT: no independent measured calibration of time units',
      'checks':diag})
    out('RG-E04/MULTITIME_PREDICTION_DESIGN.json',{
      'protocol':'HYPOTHETICAL_FORWARD_COMPUTATIONAL_DESIGN (no physical validation)',
      'N2_two_time':{'coupling':3.,'t1':.4,'t2':1.,'prep':'uniform joint over 144 pair states, g0 before switch',
       'measurement':'noninvasive ideal difference d=(b-a) mod12, 13th invalid outcome retained',
       'time_unit':'generator attempt clock calibrated in dimensionless g0 model only',
       'joint_13x13_by_model_kinetic':multi,'pairwise_kinetic_joint_TV':kinetic_sep},
      'empirical_g0_clock_unavailable':True,'physical_confirmation_status':'DESIGN_ONLY'})
    out('RG-E04/INSTRUMENT_MAPS.json',{'abstract_ideal_K':'For each state (a,b), output d=(b-a)%12 with probability 1; invalid=0',
       'finite_alphabet':list(range(12))+['INVALID'],
       'finite_sequence':'p -> D_d exp(tL^T)p; state-dependent instrument noninvasiveness only a premise',
       'imperfect_K':'13x144 stochastic matrix, entries must be independently calibrated; UNKNOWN',
       'detector_latency_drift_disturbance':'UNKNOWN; do not equate classical noninvasive filter with physical instrument',
       'unknown_all_other_outcomes':'retain invalid or unreadable record as explicit 13th outcome'})
    out('RG-E05/PREREGISTRATION.json',{'status':'DESIGN_ONLY_NOT_EVIDENCE','system':'proposed PHYS014 N2 q12 pair, not validated apparatus',
       'baseline':'independent baseline g0 uniform at equilibrium, record raw pair states and timestamps',
       'static_observable':'12 outcomes of pair label difference at g3',
       'kinetic_observable':'13x13 two-time joint difference/INVALID at .4 and 1.0 model time',
       'clock':'use independent physical g0 relaxation calibration and document its uncertainty before interpreting model time',
       'kinetic_candidates':list(KINETICS),'interaction_candidates':ids,
       'calibration_and_validation_split':'freeze physical channels / compilation/temperature/g0 clock/confusion/invalid before independent validation; keep disjoint records',
       'unavailable_data':'no selected hardware channel logs, no measured pair distribution, no g0 relaxation clock, no detector error or drift',
       'primary_decision_rule':'calibrated nuisance-model observed-law sets must have nonzero separated TV; otherwise no discrimination',
       'replication':'serial correlation and trajectory-level effective sample size required; all outcomes incl invalid',
       'source_independence':'programmed A7 does not explain physical source',
       'reuses_old_D_windows':False,'no_new_exp_performed':True})
    out('RG-E05/STOP_CERTIFICATE.json',{'status':'DESIGN_ONLY / NO_NEW_SOURCE / NO_NEW_PHYSICAL_INPUT',
       'cause':'Neither independent A7 source nor real q12 validation records/clock/readout calibrations are supplied',
       'requires':'one independently calibrated physical pair interface with frozen channel IDs/weights and drift/invalid statistics',
       'do_not_claim':['physical A7 selection','physical dynamics identification','independent experimental validation','N2 implies M3 memory', 'uniqueness over all rank-seven interactions']})
    out('RG-E05/PHYSICAL_CLAIM_PATH.json',{'M':'q12 generalized Curie-Weiss supplied A7, exact finite model identities',
       'PM':'real programmable q12 CTMC pair with measured stability, NOT ADMITTED',
       'PB':'independent mapping physical 24 channels to two q12 states, detector error/clock, NOT ADMITTED',
       'PP':'N2 static difference and N2 two-time joint model predictions; numerical designs only',
       'EV':'no independent validation records; OPEN'})
    print(json.dumps({'N2_TV':dist_N2,'M3_TV':dist_M3,'kinetic_TV':kinetic_sep,
     'hypothetical_robust':{k:v['lower_bound_if_calibrated'] for k,v in pair_bounds.items()},
     'perturbation_tv':pert_tv},indent=2))
if __name__=='__main__':main()
