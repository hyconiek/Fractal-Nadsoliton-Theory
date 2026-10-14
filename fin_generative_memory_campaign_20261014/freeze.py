"""RG-D04: deterministic freeze BEFORE any reference at t=.35/.8 or delta_2.
This script is the pretest producer. No withheld parent propagation appears herein.
"""
import hashlib,json,time,platform
from model import *
from bounds import *
from train import diagnostic_twotime
ORDER=6; FIRST=.35;SECOND=.8;PREPS=prep_names()+['delta_2']
REFERENCE_ALLOWANCE_TV=1e-8

def main():
    # Only training result can influence choice; verify it exists
    training=json.loads((HERE/'RG-D03/TRAINING_PREDICTIONS.json').read_text())
    train_choice={'selection':'channel_conditioned_Krylov_Galerkin','depth':ORDER,'cap':32,'readout':'first-stage event Y=(0,0) vs explicit complement; second-stage full ordered Y','basis':'Q S^ell D_E S^k z_p and Q S^k z_p, fixed graded ordering; rank selected by linear algebra, not trajectory fitting','training_origin':'training only t=.15,.5,1,1.5 and two-time .15->.5, preparations conditional_00/delta0/delta1','scope':'preparation-conditioned response ROM; coefficients derived on demand by the SAME FROZEN algorithm from parent S and stated p0, not a universal state-independent 27-pole memory kernel','reference_tolerance_TV':REFERENCE_ALLOWANCE_TV,'primary_TV_target':.01,'primary_preparation_class':PREPS,'frozen_single_time_tests':[FIRST,SECOND], 'frozen_joint_test':{'t1':FIRST,'t2':SECOND,'first_outcomes':['Y00','not_Y00'],'second_outcomes':'full 144 ordered Y','postselection':'none; first complement retained','parent':'q12 M3 exact-decimal FIN or same trace-matched countermodel g=5.145228719489142','comparison_baseline':'stationary-conditional averaged 144-state Markov generator of the SAME three-copy parent','withheld_reference':'not accessed until verify_frozen.py runs'}}
    # independently of any source trajectory at reserved time, compute predicted values and a priori analytic envelopes
    pred={};coeff={};modelbounds={}
    for name,A in models().items():
        m=build(A);pred[name]={};coeff[name]={};modelbounds[name]={}
        for prep in PREPS:
            B,labels=channel_conditioned_basis(m,prep,ORDER,32)
            V=basis_matrix(m,B);T,wi,wr,value=reduced(m,B)
            z=z_from_p(m,make_prep(m,prep))
            single={};a_priori={}
            for t in (FIRST,SECOND):
                y=value(prep,t)
                single[str(t)]=[float(a) for a in y]
                bd=bound_onetime(m,V,T,z,t)
                a_priori[str(t)]={'TV_upper_formula':.5*bd['l2_error_upper_formula'],'moment_tail_l2':bd['tail_formula'],'numerical_prediction_min':float(y.min()), 'prediction_mass':float(y.sum())}
            joint=diagnostic_twotime(m,T,B,make_prep(m,prep),FIRST,SECOND)
            jb=bound_twotime(m,V,T,z,FIRST,SECOND)
            pred[name][prep]={'single_time_full144':single,'joint_time_first_binary_second144':joint.tolist()}
            modelbounds[name][prep]={'single_time':a_priori,'two_time':jb,'joint_probability_min':float(joint.min()),'joint_probability_mass':float(joint.sum())}
            coeff[name][prep]={'hidden_rank':int(B.shape[1]),'source_labels':labels,'resolved_rank':144,'symmetric_generator_shape':list(T.shape),'T_trace':float(np.trace(T)),'T_eig_min':float(np.linalg.eigvalsh(T)[0]),'T_eig_max':float(np.linalg.eigvalsh(T)[-1]),'memory_K0_channel00':float(memory_operator_stats(m,B)['k0_00']),'projected_K0_channel00':float(memory_operator_stats(m,B)['reduced_k0_00']),'preparation_map':'w0 = [U^T (p0/sqrt pi), B(p0)^T (p0/sqrt pi)]','observation_map':'sqrt(piY) * first_144_components(w)'}
        print('frozen predicted',name,flush=True)
    inputs=['model.py','train.py','bounds.py','interval_inputs.py','source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md','source_context/INTAKE_20261014.md','source_context/FIN_RG_C01_C05_FULL_HANDOFF_20261014.zip','RG-D02/INPUT_ENCLOSURES.json','RG-D03/TRAINING_PREDICTIONS.json']
    hashes={p:sha(HERE/p) for p in inputs}
    d04={'execution_status':'FROZEN_BEFORE_REFERENCE','scientific_verdict':'CONDITIONAL_BOUNDED_PREPARATION_SPECIFIC_RESPONSE_PROPOSAL','protocol':train_choice,'models':list(models()),'coefficients':coeff,'input_hashes':hashes,'fitted_parameters':[],'arithmetic_status':'nominal binary64 S,T and envelope; exact-decimal parent input separately mpmath.iv enclosed; reference allowance 1e-8 conditional on IEEE arithmetic and exp/softmax implementations','bound_formula':'reversible generator spectrum [-3,0], Taylor shifted by 1.5, Poisson tail, triangle/Cauchy bound for two-time filtering','disclaimer':'B depends deterministically on initial preparation, so it is not a single universal kinetic generator for all preparations. Positivity is checked on frozen full output arrays, not asserted globally.'}
    write('RG-D04/FROZEN_PREDICTOR.json',d04)
    write('RG-D04/PREDICTIONS_PREOPEN.json',pred)
    write('RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json',modelbounds)
    pass_models={}
    for name in pred:
        maxbound=max([max(v['single_time'][t]['TV_upper_formula'] for t in ('0.35','0.8')) for v in modelbounds[name].values()]+[v['two_time']['TV_upper_formula'] for v in modelbounds[name].values()])
        minprob=min([min(v['single_time'][t]['numerical_prediction_min'] for t in ('0.35','0.8')) for v in modelbounds[name].values()]+[v['joint_probability_min'] for v in modelbounds[name].values()])
        pass_models[name]={'max_model_TV_upper_formula':maxbound,'input_reference_plus_numerical_allowance_TV':REFERENCE_ALLOWANCE_TV,'max_total_nominal_TV_upper':maxbound+REFERENCE_ALLOWANCE_TV+1e-7,'meets_0.01_prior_to_validation':maxbound+REFERENCE_ALLOWANCE_TV+1e-7<=.01,'all_frozen_probabilities_nonnegative_numerically':minprob>=-1e-12,'min_frozen_probability':minprob}
    prereg={'status':'FROZEN_BEFORE_HELDOUT_REFS','protocol':train_choice,'frozen_bound_pretest':pass_models,'selection_history':'baseline global 29 directions worked on training but did not bound unseen initial state; opted for one algorithmic preparation-conditioned channel-aware ROM, selected solely from training data / structure and declared before withheld reference','new_fitted_parameters':[],'noninference':'finite-model computational validation only','validation_abort_rule':'if any parent/countermodel selected test fails 0.01 or a frozen predicted probability is negative, retain failure; do not refit under this identifier'}
    write('RG-D04/PREREGISTRATION.json',prereg)
    locked=['RG-D04/FROZEN_PREDICTOR.json','RG-D04/PREDICTIONS_PREOPEN.json','RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json','RG-D04/PREREGISTRATION.json']
    hashes_locked={p:sha(HERE/p) for p in locked}
    write('RG-D04/LOCK.json',{'pretest_files_sha256':hashes_locked,'status':'LOCKED_BEFORE_RUNNING_VALIDATE','validator_must_check_sha256':True})
    print('PRE-FREEZE COMPLETE\n',json.dumps(pass_models,indent=2))
if __name__=='__main__':main()
