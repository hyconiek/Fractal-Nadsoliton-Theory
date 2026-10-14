"""RG-D05. Opens heldout parent references ONLY AFTER verifying RG-D04 lock.
Uses independent tensor uniformization for reference and sparse transition as crosscheck.
Does not rebuild/select/modify/retune prediction or generator coefficients.
"""
import json,hashlib,math,time
from model import *
from scipy.linalg import expm

def baseline_joint(m,t1,t2):
    s=np.sqrt(m['piY']);a=np.zeros(Y);a[0]=1/s[0]
    av1=expm(m['A']*t1)@a
    ve=np.zeros(Y);ve[0]=av1[0]
    p2=expm(m['A']*(t2-t1))
    return np.r_[s*(p2@ve),s*(p2@(av1-ve))]

def check():
    locked=json.loads((HERE/'RG-D04/LOCK.json').read_text())
    for p,expected in locked['pretest_files_sha256'].items():
        seen=sha(HERE/p)
        if seen!=expected:raise RuntimeError('FREEZE CHECK FAILED '+p)
    frozen_contract=json.loads((HERE/'RG-D04/FROZEN_PREDICTOR.json').read_text())
    for rel,expected in frozen_contract['input_hashes'].items():
        if sha(HERE/rel)!=expected:raise RuntimeError('FROZEN INPUT CHANGED '+rel)
    pred=json.loads((HERE/'RG-D04/PREDICTIONS_PREOPEN.json').read_text())
    env=json.loads((HERE/'RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json').read_text())
    protocol=json.loads((HERE/'RG-D04/PREREGISTRATION.json').read_text())
    assert protocol['status']=='FROZEN_BEFORE_HELDOUT_REFS'
    assert protocol['protocol']['depth']==6
    t1,t2=.35,.8;allowed=1e-8;model_tests={};raw={};summary=[]
    for name,A in models().items():
        m=build(A);model_tests[name]={};raw[name]={}
        mj=baseline_joint(m,t1,t2)
        for prep in (prep_names()+['delta_2']):
            p0=make_prep(m,prep)
            one={};refs={}
            for t in (t1,t2):
                p=uniform_tensor(m,p0,t,65)
                p_ind=uniform_p(m,p0,t,65)
                full=observed(m,p);ind=observed(m,p_ind)
                assert max(abs(full-ind))<1e-12
                q=np.array(pred[name][prep]['single_time_full144'][str(t)])
                baseline=markov(m,t)
                bound=env[name][prep]['single_time'][str(t)]['TV_upper_formula']+allowed+1e-7
                gap=tv(full,q)
                one[str(t)]={'TV_ROM_vs_parent':gap,'TV_averaged_Markov_vs_parent':tv(full,baseline),'a_priori_upper_envelope_with_arithmetic_allowance':bound,'within_envelope':gap<=bound,'within_predeclared_0.01':gap<=.01,'mass_ROM':float(q.sum()),'min_ROM':float(q.min()),'mass_ref':float(full.sum()),'tensor_sparse_max_abs':float(max(abs(full-ind)))}
                refs[str(t)]=full.tolist()
            p1=uniform_tensor(m,p0,t1,65)
            mask=np.zeros(N);mask[:12]=1
            ev=uniform_tensor(m,p1*mask,t2-t1,65)
            nev=uniform_tensor(m,p1*(1-mask),t2-t1,65)
            j_exact=np.r_[observed(m,ev),observed(m,nev)]
            j_pred=np.array(pred[name][prep]['joint_time_first_binary_second144'])
            jb=env[name][prep]['two_time']['TV_upper_formula']+allowed+1e-7
            jj=tv(j_exact,j_pred)
            assert abs(j_exact.sum()-1)<1e-12
            assert np.max(np.abs(j_exact[:Y]+j_exact[Y:]-refs[str(t2)]))<1e-12
            joint={'TV_ROM_vs_parent':jj,'TV_averaged_Markov_vs_parent':tv(j_exact,mj),'a_priori_upper_envelope_with_arithmetic_allowance':jb,'within_envelope':jj<=jb,'within_predeclared_0.01':jj<=.01,'first_event_probability_ref':float(j_exact[:Y].sum()),'first_event_probability_ROM':float(j_pred[:Y].sum()),'min_ROM':float(j_pred.min()),'mass_ROM':float(j_pred.sum()),'mass_ref':float(j_exact.sum()),'full_second_time_marginal_compatibility_inf':float(np.max(np.abs(j_exact[:Y]+j_exact[Y:]-refs[str(t2)])))}
            raw[name][prep]={'one_time_reference_full144':refs,'joint_reference_288':j_exact.tolist()}
            model_tests[name][prep]={'one_time':one,'two_time':joint}
            for t,d in one.items():summary.append({'model':name,'preparation':prep,'scope':'one_time_'+t,'TV_ROM':d['TV_ROM_vs_parent'],'TV_markov':d['TV_averaged_Markov_vs_parent'],'bound':d['a_priori_upper_envelope_with_arithmetic_allowance'],'PASS':d['within_predeclared_0.01'] and d['within_envelope']})
            summary.append({'model':name,'preparation':prep,'scope':'two_time_0.35_0.8','TV_ROM':joint['TV_ROM_vs_parent'],'TV_markov':joint['TV_averaged_Markov_vs_parent'],'bound':jb,'PASS':joint['within_predeclared_0.01'] and joint['within_envelope']})
            print('VALIDATED',name,prep,'one',one['0.8']['TV_ROM_vs_parent'],'two',joint['TV_ROM_vs_parent'],flush=True)
    verdict={'execution_status':'EXECUTED_AFTER_HASHED_FREEZE','scientific_verdict':'COMPUTATIONAL_FINITE_MODEL_TEST_ONLY','validation_kind':'reserved same-g synthetic parent process, no physical measurement','frozen_hashes':locked['pretest_files_sha256'],'model_tests':model_tests,'n_test_laws':len(summary),'n_tests_meeting_bound_and_0.01':sum(x['PASS'] for x in summary),'aggregate_max_TV_rom':max(x['TV_ROM'] for x in summary),'aggregate_max_TV_markov':max(x['TV_markov'] for x in summary),'FIN_max_TV_rom':max(x['TV_ROM'] for x in summary if x['model']=='FIN_A7'),'FIN_max_TV_markov':max(x['TV_markov'] for x in summary if x['model']=='FIN_A7'),'all_probabilities_valid_on_selected_laws':all(z['two_time']['min_ROM']>=-1e-12 and all(v['min_ROM']>=-1e-12 for v in z['one_time'].values()) for model in model_tests.values() for z in model.values()),'reference_analytic_allowance_TV':allowed,'numerical_safety_margin_TV':1e-7,'has_universal_preparation_independent_generator':False,'fitted_after_open':False}
    write('RG-D05/VALIDATION.json',verdict)
    write('RG-D05/RAW_MODEL_OUTPUTS.json',raw)
    write('RG-D05/TEST_SUMMARY.json',summary)
    return verdict
if __name__=='__main__':
    v=check();print('OVERALL',v['n_tests_meeting_bound_and_0.01'],'/',v['n_test_laws'],'worst ROM',v['aggregate_max_TV_rom'],'worst baseline',v['aggregate_max_TV_markov'])
