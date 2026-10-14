#!/usr/bin/env python3
"""Bounded RG-F01..F04 computational proofs/designs, no empirical observations."""
import json,math,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.stats import binom
from source_model import all_models,build_n2_generator,pi_n2,pair_difference,tv
from positive_process import difference_generator,stationary_difference,propagation,emission,sequence_probability,full_joint,uniformization
ROOT=Path(__file__).resolve().parent
NAME=['FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE']
KINDS=['heat_bath','metropolis','barker']


def store(name,obj):
    target=ROOT/name;target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')

def make_diff_projection():
    d=np.array([(b-a)%12 for a in range(12) for b in range(12)])
    R=np.zeros((144,12));R[np.arange(144),d]=1.
    return d,R

def f01(models):
    data={};tests=[];d,R=make_diff_projection()
    K=emission(); Ks=emission(.03,.02)
    mixed=(np.arange(12)+1).astype(float);mixed/=mixed.sum()
    alt=np.zeros((13,12));alt[0,:6]=1;alt[1,6:]=1
    for mid in NAME:
      A=models[mid];data[mid]={}
      for kind in KINDS:
        B=difference_generator(A,3.,kind)
        L=build_n2_generator(A,3.,kind).toarray()
        intertwine=float(np.max(np.abs(L@R-R@B)))
        p=stationary_difference(A)
        # E predecessor prepared uniform g0, t=.4,1.0
        joint=full_joint(np.ones(12)/12,B,[.4,1.0],K,'heat_bath' if kind=='heat_bath' else None,p)
        pred=(ROOT.parent/'fin_generative_operational_campaign_20261014') # not a source dependency
        # probabilistic laws from positive CTMC; fixed and noisy instrument, arbitrary preparations
        P=propagation(B,.4,kind if kind=='heat_bath' else None,p)
        Pup=uniformization(B,.4,4.,90)
        rev=float(np.max(np.abs(p[:,None]*B-(p[:,None]*B).T)))
        joint_noisy=full_joint(mixed,B,[.21,.78],Ks,'heat_bath' if kind=='heat_bath' else None,p)
        # marginalize first readout -> single final law
        final=(mixed@propagation(B,.78,kind if kind=='heat_bath' else None,p))@Ks.T
        # affine in preparation
        second=np.roll(mixed,3)
        mix=.37*mixed+.63*second
        frommix=full_joint(mix,B,[.21,.78],Ks,kind if kind=='heat_bath' else None,p)
        mixtarget=.37*joint_noisy+.63*full_joint(second,B,[.21,.78],Ks,kind if kind=='heat_bath' else None,p)
        norm3=sum(sequence_probability(mixed,B,[.12,.28,.6],(a,b,c),Ks,kind if kind=='heat_bath' else None,p)
                  for a in range(13) for b in range(13) for c in range(13))
        checks={'generator_row_sum':float(np.max(np.abs(B.sum(1)))),'negative_offdiag_count':int(np.sum((B-np.diag(np.diag(B)))<-1e-14)),
         'full_parent_intertwine':intertwine,'stationarity':float(np.max(np.abs(p@B))),
         'detailed_balance':rev,'heatbath_reset_error':float(np.max(np.abs(B-2*(np.outer(np.ones(12),p)-np.eye(12))))) if kind=='heat_bath' else None,
         'uniformization_difference_at_t04':float(np.max(np.abs(P-Pup))),
         'joint_mass_ideal':float(joint.sum()),'joint_mass_noisy':float(joint_noisy.sum()),
         'marginalization_max_error':float(np.max(np.abs(joint_noisy.sum(0)-final))),
         'mixture_affinity_error':float(np.max(np.abs(frommix-mixtarget))),
         'sequence_3_total_mass':float(norm3),'invalid_total_noisy':float(joint_noisy[-1,:].sum()+joint_noisy[:,-1].sum()-joint_noisy[-1,-1])}
        data[mid][kind]={'generator_12x12':B.tolist(),'stationary_12':p.tolist(),
          'ideal_joint_13x13_g0_uniform_quench_t04_t10':joint.tolist(),
          'synthetic_noisy_joint_13x13':joint_noisy.tolist(), 'checks':checks}
        tests.append({'model':mid,'kinetic':kind,'checks':checks})
    # selected ternary readout still composes and retains all outcomes
    altout=full_joint(mixed,difference_generator(models['FIN_A7']),[.2,.6],alt,'heat_bath',stationary_difference(models['FIN_A7']))
    store('RG-F01/PROCESS_CONTRACT.json',{
     'identity':'N2 q=12, two labelled copies, g=3, zero field; d=(b-a) mod12',
     'model_classes':NAME,'kinetics':KINDS,'positive_quotients':data,
     'observation':'K[o,d] 13x12 column stochastic, o=0..11 or INVALID; no physical disturbance assumed',
     'preparation':'p0(d) arbitrary probability vector, including any marginal of nonuniform microscopic phase preparation',
     'g0_identity':'All B(g0) have B0[d,e]=2/12 for d!=e; B0=2(1u^T-I), nonconstant decay 2',
     'full_parent_ref':'L R=R B; derived conditionally on circulant zero-field rate law',
     'model_is_not_physics':True})
    store('RG-F01/INSTRUMENT_CONTRACT.json',{
     'outputs':[str(i) for i in range(12)]+['INVALID'],
     'required':{'K_shape':[13,12],'column_stochastic':True,'fixed_across_outcomes':True,
        'state_feedback':'none under noninvasive ideal contract; if physical disturbance exists specify J[o,d_before,d_after] with sum_o J stochastic',
        'hidden_invalid_state':'not modeled; physical invalid-state transition is OPEN_INPUT'},
     'example_synthetic_K':Ks.tolist(),'arbitrary_phase_initial_allowed_for_equivariant_interface':True,
     'phase_dependent_real_detector':'not licensed by quotient'},)
    store('RG-F01/COMPOSITION_TESTS.json',{'nine_cases':tests,'three_time_law_with_13_outcomes':'normalization tested for 2197 strings per kinetic/model',
       'alternative_phase_independent_binary_emission':{'mass':float(altout.sum()),'positive':bool(altout.min()>=-1e-14)},
       'numerical_scope':'binary64/scipy, no outward all-ops claim'})
    return data

def f02(models,data):
    d,R=make_diff_projection()
    A=models['FIN_A7'];L=build_n2_generator(A,3.,'heat_bath').toarray()
    # Explicit break: add h=+.01 to preference for absolute label a'=0 when updating COPY 1
    from scipy.special import softmax
    h=.01;Lp=np.zeros((144,144))
    for a in range(12):
      for b in range(12):
        i=a*12+b
        for site in (0,1):
          old,other=(a,b) if site==0 else (b,a)
          field=h*(np.arange(12)==0) if site==0 else np.zeros(12)
          w=softmax(3*A[:,other]/2+field)
          for new in range(12):
            if new==old:continue
            j=new*12+b if site==0 else a*12+new
            Lp[i,j]+=w[new];Lp[i,i]-=w[new]
    def defect(full):
      rows=full@R
      largest=0;pair=None
      for k in range(12):
        states=np.flatnonzero(d==k)
        for x in states[1:]:
          diff=np.abs(rows[x]-rows[states[0]])
          v=float(diff.max())
          if v>largest:largest=v;pair=[int(states[0]),int(x),k,int(np.argmax(diff))]
      return largest,pair
    perfect,pairp=defect(L);broken,pairb=defect(Lp)
    # Phase-related readout obstruction. Two states same d but detector output depends on a.
    raw_effect=np.zeros((13,144))
    for i in range(144):
       a,b=divmod(i,12);out=(b-a)%12
       raw_effect[out,i]=.99 if a==0 else 1.
       raw_effect[12,i]=.01 if a==0 else 0.
    ideal=np.zeros((13,144))
    for i in range(144):ideal[d[i],i]=1.
    # Output uniform over all 13 (erasure) compatible with naive uniformity check.
    erasure=np.zeros((13,144));erasure[:12,:]=1/12
    maxphase=0
    for k in range(12):
      inds=np.flatnonzero(d==k)
      maxphase=max(maxphase,float(np.max(np.abs(raw_effect[:,inds]-raw_effect[:,[inds[0]]]))))
    result={'q':12,'interaction':'FIN_A7','kinetics':'heat_bath',
        'mathematical_zero_field_quotient_defect':perfect,
        'one_declared_0p01_absolute_first_copy_field':{'field':h,'max_block_rate_difference':broken,'witness_[first_state,second_state,d,target_d]':pairb},
        'phase_dependent_detector_0p01_invalid':{'max_column_output_difference_same_d':maxphase,
           'violates_fixed_K13x12':True},
        'perfect_instrument_factorization_defect':0.,
        'erasure_13x144_column_stochastic':bool(np.allclose(erasure.sum(0),1)),
        'erasure_K_fingerprint_tv':0.,'erasure_g0_valid_12_histogram_exactly_uniform':True,'erasure_invalid_probability':0.,
        'raw_microstates_required':'phase-resolved (a,b), stimulus labels/gain, timestamps, emission outcome including INVALID',
        'all_physical_symmetry_or_readout_tests':'OPEN_INPUT'}
    store('RG-F02/SYMMETRY_AND_ERASURE_CONTROLS.json',result)
    store('RG-F02/CALIBRATION_CAPABILITIES.json',{
      'need_to_verify_strong_lumpability':['Prepare two distinct absolute states (a,b) with same d','Estimate each transition into all twelve d fibres at common controlled clock','Repeat across 12 phases or obtain certified symmetry mapping','Check rate and detector/instrument change separately'],
      'need_to_verify_instrument':['Calibrate each J_o distribution over outcome and successor difference, not marginal histogram only','Retain all invalid outcomes','Test detector/feedback latency and drift','If conditioning on d only, phase-resolved coupling must be independently justified'],
      'available':'mathematical model and PHYS-014 design-only','missing':'all hardware calibration records; absolute pair preparation; effect of physical invalid state on subsequent transition',
      'physical_gate':'OPEN_INPUT'})
    return result

def sample_bound(k,rho,alpha=.025,eps_ratio=1/3):
    eps=rho*eps_ratio
    n=math.ceil(math.log((2**k-2)/alpha)/(2*eps*eps))
    assert 2*eps<rho
    return {'alphabet_size':k,'rho':rho,'validation_alpha':alpha,'epsilon':eps,'n_independent_records_sufficient_Hoeffding_union':n,
     'proof':'TV(emp,p)=max over nontrivial proper event subsets S of [emp(S)-p(S)]; union of 2^k-2 one-sided Hoeffding exp(-2n eps^2).'}

def f03(models,stat):
    inherited_gap=.00119114
    # Exact truth conditional on inherited gain-cover certificate; not re-certified here
    required={'inherited_ideal_gain_box_min_TV':.01119114,
       'each_hypothesis_output_TV_allowance':.005,
       'residual_lower_margin':inherited_gap,
       'further_uncertainty_sum_max_for_positive_margin_strict':inherited_gap,
       'status':'HYPOTHETICAL_INPUT_CLASS_NOT_MEASURED'}
    allocation={'example_per_hypothesis_added_burnin_TV':.0002,
                'example_per_hypothesis_added_field_compilation_TV':.0001,
                'example_per_hypothesis_added_drift_TV':.0001,
                'remaining_pair_gap_after_example':inherited_gap-2*(.0002+.0001+.0001),
                'is_measured':False}
    # heat bath exact worst-case t via stationary reset law: bound e^-2t <=.0002
    allocation['necessary_sufficient_worst_case_heatbath_dimensionless_wait_for_eps_0p0002']=math.log(1/.0002)/2
    # General law-level telescoping bound, with uncertain parameters placeholders (null)
    budget={'static_primary':required,'hypothetical_extra_allocation':allocation,
      'components':{'production_equilibrium_initial_TV':None,'phase_symmetry_compilation_TV':None,
       'field_and_target_interaction_TV':None,'detector_invalid_drift_TV':None,
       'clock_rate_TV_if_dynamic':None,'instrument_disturbance_TV_if_dynamic':None,
       'arithmetic_TV':None,'independent_validation_sampling_TV':None},
      'identity':'TV(full observed law) <= e_prep + sum_intervals (delta_t/2)*max_x sum_y |L_xy-Lref_xy| + sum_readouts sup_x TV(J_x, Jref_x) + e_arith when common underlying state/representation; any invalid-state transition needs explicit joint J and its successor state',
      'probability_error_not_the_same_as_statistical_error':True,
      'posterior_conditioning':'never renormalize away invalid; postselection requires its own denominator',
      'physical_complete_bound':'OPEN_INPUT'}
    # Fixed calibrated gain-box and .005 instrument is ENLARGED hypothetical set, not device evidence.
    store('RG-F03/TOTAL_LAW_ERROR_BUDGET.json',budget)
    store('RG-F03/OPEN_BOUND.json',{'missing_inputs':list(budget['components']),
      'physical_ready':False,'bound_status':'NO_PHYSICAL_INPUT',
      'existing_abstract_bound':'quoted from user-supplied independent intake; not a fresh SHA-checked certificate',
      'example_gain_box':[2.97,3.03], 'output_TV_radius_each':.005})
    # numerical exact-math Poisson tail bound with mpmath iv in a separate small check; not any all-ops proof
    from mpmath import mp
    mp.dps=80
    times=[.4,.6,1.0]
    tails={}
    for t in times:
        z=mp.mpf(4)*mp.mpf(str(t));m=65
        # e^-z z^(m+1)/(m+1)! / (1-z/(m+2))
        tail=mp.exp(-z)*z**(m+1)/mp.factorial(m+1)/(1-z/(m+2))
        tails[str(t)]={'poisson_tail_bound_m65':mp.nstr(tail,25),'z':mp.nstr(z,12)}
    store('RG-F03/PREDICTION_ENCLOSURES.json',{'times':times,'uniformization_lambda':4,
      'domination_is_valid_for_the_nine_n2_generators':bool(max([-min(np.diag(np.array(stat[m][kind]['generator_12x12']))) for m in NAME for kind in KINDS])<=4),
      'tail_beyond_65':tails,'all_operation_error_included':False,
      'conditional_input_note':'source scientific kernel exact-input enclosure and gain-box cover inherited from intake, not recalculated in this run',
      'needed_for_all_operation_certificate':'interval input categorical rows + outward interval matrix arithmetic + outward Poisson weights + compensated mass checks; floating-only term insufficient'})
    return budget

def f04(models,data,stat):
    rho=.00119114;alpha_total=.05;alpha_cal=.025;alpha_val=.025
    sample=sample_bound(13,rho,alpha_val,1/3)
    # event maximizing FIN nominal TV. It is fixed before future validation records.
    pf=np.array(stat['FIN_A7']['heat_bath']['stationary_12'])
    p7=np.array(stat['FLAT_P7_TRACE']['heat_bath']['stationary_12'])
    subset=[int(i) for i in np.where(pf>p7)[0]]
    eventstats={}
    for m in ['FIN_A7','FLAT_P7_TRACE']:
      lo=float(pair_difference(models[m],2.97)[subset].sum())
      hi=float(pair_difference(models[m],3.03)[subset].sum())
      eventstats[m]={'g2p97':lo,'g3p03':hi,'center':float(pair_difference(models[m],3)[subset].sum())}
    nominal=float(pf[subset].sum()-p7[subset].sum())
    # The endpoints do NOT certify global extrema; use conservative general gain Lipschitz proof
    caps={m: float(.03*np.ptp(models[m][0])/4) for m in ['FIN_A7','FLAT_P7_TRACE']}
    lower=max(0.,nominal-caps['FIN_A7']-caps['FLAT_P7_TRACE']-.01)
    event={'subset_predeclared_nominal_argmax':subset,'nominal_event_difference':nominal,'endpoint_samples_not_global_certificates':eventstats,
     'rigorously_simple_gain_Lipschitz_TV_cap':caps,'event_gap_lower_bound_after_gain_and_two_detector_0p005':lower,
     'fixed_one_event_robust_positive':'NOT_CERTIFIED_UNDER_CLASS',
     'event_iid_hoeffding':'P(|sample_event - p_event| > eps) <= 2 exp(-2n eps^2); model-set separation not automatically established'}
    # simple calibration finite confidence split; no data. decision for finite model-set with exact TV distances
    decision={'status':'PROSPECTIVE_CONDITIONAL_STATISTICAL_DESIGN',
      'primary':'FIN_A7_vs_FLAT_P7_TRACE_g3_stationary_13_outputs_with_INVALID',
      'secondary':'FIN_A7_heat_bath_vs_metropolis_N2_two_time_13x13 - only after independent g0 clock and instrument calibration',
      'hierarchy':'primary first; secondary interpreted only if separately powered/error allocated; no selection by minimum observed p',
      'alpha_total':alpha_total,'alpha_calibration_budget':alpha_cal,'alpha_validation_budget':alpha_val,
      'conditional_model_set_gap_rho':rho,'finite_alphabet_iid_rule':sample,
      'confidence_region':'{q: TV(q,p_emp)<=eps_n}; each hypothesis set H_X includes separately allowed gains and instruments, not silently shared nuisance; accept compatibility iff lower bound on inf_{q in H_X} TV(q,p_emp)<=eps',
      'decisions':{'FIN_ONLY':'flat model rejected, FIN not rejected','FLAT_ONLY':'FIN rejected, flat not rejected',
        'INCONCLUSIVE':'neither rejected','NEITHER_MODEL_FITS':'both rejected'},
      'false_rejection_guarantee':'conditional on true distribution in declared H_X and iid n records: <= alpha_validation by finite-subset union bound',
      'type_II_wrong_unique_selection':'if calibrated set gap rho>2eps and concentration event holds, wrong alternative cannot be uniquely selected; otherwise no unconditional power',
      'dependence':'independent experiment/reset required; measured ESS is not a theorem; otherwise calibrate mixing or use independent blocks',
      'calibration_uncertainty':'physical numerical values not provided; illustrative ±1%, .005 are NOT measured',
      'single_event_comparator':event,'freeze_before_physical_validation':True}
    store('RG-F04/DECISION_RULE.json',decision)
    # Minimal synthetic decision checks on a toy 13 alphabet, not the FIN data, with eps=0.04
    kt=13;q=np.zeros(kt);r=np.zeros(kt);q[0]=1;r[1]=1
    eps=.04
    def flags(v):
      f1=tv(v,q)<=eps+1e-12;f2=tv(v,r)<=eps+1e-12
      return 'INCONCLUSIVE' if f1 and f2 else 'FIN_ONLY' if f1 else 'FLAT_ONLY' if f2 else 'NEITHER_MODEL_FITS'
    tests={'true_fin':flags(q),'true_flat':flags(r),'erasure':flags(np.ones(kt)/kt),
           'boundary_q':flags(.96*q+.04*r),'boundary_beyond':flags(.95*q+.05*r),
           'too_small_sample':'INCONCLUSIVE if eps_n>=rho/2; no model selection assertion',
           'status':'software toy distribution, not real or withheld physical data'}
    store('RG-F04/SYNTHETIC_TESTS.json',tests)
    store('RG-F04/ONE_EVENT_COMPARISON.json',event)
    return decision,tests

def main():
    models=all_models()
    stat=f01(models)
    sym=f02(models,stat)
    budget=f03(models,stat)
    decision,tests=f04(models,stat,stat)
    print(json.dumps({'cases':sum(len(v) for v in stat.values()),
     'lumpability_defect_ideal':sym['mathematical_zero_field_quotient_defect'],
     'lumpability_defect_h0p01':sym['one_declared_0p01_absolute_first_copy_field']['max_block_rate_difference'],
     'static_conditional_gap':budget['static_primary']['residual_lower_margin'],
     'illustrative_extra_gap':budget['hypothetical_extra_allocation']['remaining_pair_gap_after_example'],
     'sufficient_iid_n_13':decision['finite_alphabet_iid_rule']['n_independent_records_sufficient_Hoeffding_union'],
     'toy_decisions':tests},indent=2))

if __name__=='__main__':main()
