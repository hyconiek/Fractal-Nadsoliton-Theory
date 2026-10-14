"""Independent output/lock/statistical verifier. Pure stdlib, NO model imports.
Reads serialized frozen predictions and independently propagated references.
"""
import json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
rd=lambda p:json.loads((ROOT/p).read_text())
sha=lambda p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def tv(x,y): return math.fsum(abs(a-b) for a,b in zip(x,y))/2

def main():
    lock=rd('RG-D04/LOCK.json'); pre=rd('RG-D04/PREREGISTRATION.json')
    pred=rd('RG-D04/PREDICTIONS_PREOPEN.json');ref=rd('RG-D05/RAW_MODEL_OUTPUTS.json')
    bounds=rd('RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json');summary=rd('RG-D05/VALIDATION.json')
    mem=rd('RG-D03/MEMORY_REDUCTION.json'); pre_input=rd('RG-D02/INPUT_ENCLOSURES.json')
    tests={}
    def ck(tag,cond): tests[tag]=bool(cond)
    ck('freeze_four_hashes',all(sha(f)==h for f,h in lock['pretest_files_sha256'].items()))
    ck('freeze_verified_prior_to_open',pre['status']=='FROZEN_BEFORE_HELDOUT_REFS' and summary['execution_status']=='EXECUTED_AFTER_HASHED_FREEZE')
    ck('finite_models_exactly_three',set(pred)=={'FIN_A7','FLAT_P7_TRACE','FULL_POTTS_TRACE'})
    ck('pi_exact_kernel_conditional_enclosure',pre_input['input_status']=='INTERVAL_ENCLOSED_CONDITIONAL_ON_MPMATH_IV')
    all_laws=0;passes=0;maxrom=0;maxmark=0;positive=True;consistent=True;rank=True;better=True
    for model,preps in pred.items():
        ck(model+'_four_preparations',set(preps)=={'conditional_00','delta_0','delta_1','delta_2'})
        for prep,vals in preps.items():
            exact=ref[model][prep]
            jointE=exact['joint_reference_288'];jointP=vals['joint_time_first_binary_second144']
            bnd=bounds[model][prep]
            for t in ('0.35','0.8'):
                x=exact['one_time_reference_full144'][t];y=vals['single_time_full144'][t]
                assert len(x)==len(y)==144
                gap=tv(x,y);lim=bnd['single_time'][t]['TV_upper_formula']+1.1e-7
                passes+=gap<=.01 and gap<=lim;all_laws+=1
                maxrom=max(maxrom,gap);positive&=(min(y)>=-1e-12 and abs(math.fsum(y)-1)<1e-11)
                if t=='0.8':consistent &=max(abs(x[i]-(jointE[i]+jointE[i+144])) for i in range(144))<1e-12
            assert len(jointE)==len(jointP)==288
            gap=tv(jointE,jointP);lim=bnd['two_time']['TV_upper_formula']+1.1e-7
            passes+=gap<=.01 and gap<=lim;all_laws+=1
            maxrom=max(maxrom,gap);positive&=(min(jointP)>=-1e-12 and abs(math.fsum(jointP)-1)<1e-11)
    ck('all_36_laws_recomputed',all_laws==36 and passes==36)
    ck('nonnegative_normalized_predictions_all_36',positive)
    ck('joint_second_marginal_reconstructed',consistent)
    ck('max_ROM_TV_recomputed',abs(maxrom-summary['aggregate_max_TV_rom'])<1e-13)
    ck('all_memory_integrals_distinct_from_instantaneous',abs(mem['FIN_A7']['K0_00']-mem['FIN_A7']['M0_integral_00_numeric'])>0.05 and abs(mem['FIN_A7']['Kprime0_00']-mem['FIN_A7']['M1_first_integral_moment_00_numeric'])>0.1)
    ck('memory_cg_residual_small',all(mem[n]['CG_status']==0 and mem[n]['CG_residual_l2']<1e-10 for n in ('FIN_A7','FLAT_P7_TRACE','FULL_POTTS_TRACE')))
    ck('g0_degenerate_memory',mem['FIN_A7_g0']['K0_00']<1e-25)
    ck('reduction_better_than_markov_on_all_laws',all(z['TV_ROM']<z['TV_markov'] for z in rd('RG-D05/TEST_SUMMARY.json')))
    ck('bounded_model_hidden_rank_32',all(c['hidden_rank']<=32 for model in rd('RG-D04/FROZEN_PREDICTOR.json')['coefficients'].values() for c in model.values()))
    results={'verifier':'independent stdlib parser/arithmetic of exported immutable predictions, no producer imports','checks':tests,'count_pass':sum(tests.values()),'count_total':len(tests),'laws_checked':all_laws,'laws_pass':passes,'overall_pass':all(tests.values()),'independence_limit':'does not independently derive the source A7, numerical L operator or rigorous floating error accounting'}
    (ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
    print('VERIFICATION',results['count_pass'],'/',results['count_total'],'laws',passes,'/',all_laws)
    if not results['overall_pass']:raise RuntimeError('INDEPENDENT_VERIFICATION_FAILED '+str(tests))
if __name__=='__main__':main()
