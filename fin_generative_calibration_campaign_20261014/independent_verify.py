#!/usr/bin/env python3
"""Independent RG-F verifier. Does NOT import run_calibration or positive_process.
It constructs a 144-state generator and propagates separate labelled branches.
Source_model is the frozen E finite-model input (not independent A7 physics).
"""
import json, math, zipfile
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.special import softmax
from mpmath import mp
from source_model import all_models,build_n2_generator,tv
R=Path(__file__).resolve().parent

def load(p):return json.loads((R/p).read_text())
def check(cond,name,results,value=None):
    results.append({'test':name,'passed':bool(cond),'value':None if value is None else float(value)})
    if not cond:print('FAIL',name,value)

def independent_B(L,ds):
    B=np.zeros((12,12)); diffs=[]
    for d in range(12):
        ph=np.flatnonzero(ds==d)
        rates=np.array([np.bincount(ds,weights=L[row],minlength=12) for row in ph])
        B[d]=rates[0]
        diffs.append(float(np.max(np.abs(rates-B[d]))))
    return B,max(diffs)

def reference_joint(L,ds):
    p1=expm_multiply(.4*L.T,np.ones(144)/144)
    j=np.zeros((13,13))
    for d1 in range(12):
        branch=p1*(ds==d1)
        end=expm_multiply(.6*L.T,branch)
        j[d1,:12]=np.bincount(ds,weights=end,minlength=12)
    return j

def main():
    src=load('RG-F01/PROCESS_CONTRACT.json')
    checks=load('RG-F01/COMPOSITION_TESTS.json')
    sy=load('RG-F02/SYMMETRY_AND_ERASURE_CONTROLS.json')
    dec=load('RG-F04/DECISION_RULE.json')
    bud=load('RG-F03/TOTAL_LAW_ERROR_BUDGET.json')
    tails=load('RG-F03/PREDICTION_ENCLOSURES.json')
    a=all_models();ds=np.array([(j-i)%12 for i in range(12) for j in range(12)])
    results=[];max_joint=0;max_B=0
    with zipfile.ZipFile(R/'source_context/FIN_RG_E01_E05_FULL_HANDOFF_20261014.zip') as archive:
        predecessor=json.loads(archive.read('fin_generative_operational_campaign_20261014/RG-E04/MULTITIME_PREDICTION_DESIGN.json'))
    prior=predecessor['N2_two_time']['joint_13x13_by_model_kinetic']
    for name in ('FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE'):
      for kind in ('heat_bath','metropolis','barker'):
        L=build_n2_generator(a[name],3.,kind).toarray()
        B,phase=independent_B(L,ds)
        saved=np.array(src['positive_quotients'][name][kind]['generator_12x12'])
        err=float(np.max(np.abs(B-saved)));max_B=max(err,max_B)
        j=reference_joint(L,ds)
        stored=np.array(src['positive_quotients'][name][kind]['ideal_joint_13x13_g0_uniform_quench_t04_t10'])
        e=float(np.max(np.abs(j-stored)));max_joint=max(e,max_joint)
        check(err<3e-14,f'{name}.{kind}.independent_lumped_matrix',results,err)
        check(phase<3e-14,f'{name}.{kind}.phase_lumpability',results,phase)
        check(e<5e-13,f'{name}.{kind}.full_parent_two_time_law',results,e)
        original=np.array(next(row['joint_13x13'] for row in prior[name] if row['kinetic']==kind))
        check(np.max(np.abs(j-original))<5e-13,f'{name}.{kind}.RG_E_regression',results,float(np.max(np.abs(j-original))))
        check(float(np.max(np.abs(B.sum(axis=1))))<1e-13,f'{name}.{kind}.row_sum',results)
        off=B-np.diag(np.diag(B))
        check(off.min()>-1e-14,f'{name}.{kind}.offdiag_nonnegative',results)
        check(j.min()>-2e-14 and abs(j.sum()-1)<2e-13,f'{name}.{kind}.positive_mass',results)
        p=softmax(1.5*a[name][0])
        check(np.max(np.abs(p@B))<2e-13,f'{name}.{kind}.gibbs_stationary',results)
        check(np.max(np.abs(p[:,None]*B-(p[:,None]*B).T))<2e-13,f'{name}.{kind}.detailed_balance',results)
        L0=build_n2_generator(a[name],0,kind).toarray();B0,_=independent_B(L0,ds)
        target=2*(np.ones((12,12))/12-np.eye(12))
        check(np.max(np.abs(B0-target))<3e-14,f'{name}.{kind}.g0_clock_identity',results)
    # gamma = inverse tau in model units; no physical calibration implied
    check(abs(sy['one_declared_0p01_absolute_first_copy_field']['max_block_rate_difference']-.002711545661289394)<1e-12,
          'asymmetric_control_saved_witness',results)
    check(sy['erasure_K_fingerprint_tv']==0 and sy['erasure_g0_valid_12_histogram_exactly_uniform'],
          'erasure_obstruction',results)
    k=dec['finite_alphabet_iid_rule']['alphabet_size'];rho=mp.mpf('0.00119114');eps=rho/3
    n=int(dec['finite_alphabet_iid_rule']['n_independent_records_sufficient_Hoeffding_union'])
    mp.dps=90; tailprob=(2**k-2)*mp.exp(-2*n*eps**2)
    check(tailprob<=mp.mpf('.025') and (2**k-2)*mp.exp(-2*(n-1)*eps**2)>mp.mpf('.025'),
          'finite_sample_40m_integer_threshold',results,float(tailprob))
    check(2*eps<rho,'gap_vs_two_confidence_radii',results,float(2*eps))
    check(abs(bud['hypothetical_extra_allocation']['remaining_pair_gap_after_example']-.00039114)<1e-12,
          'positive_example_budget',results)
    for t in ['0.4','0.6','1.0']:
        z=4*mp.mpf(t);m=65
        true_tail=mp.gammainc(m+1,0,z)/mp.factorial(m) # Poisson P(N>m) = lower gamma(m+1,z)/m!
        conservative=mp.mpf(tails['tail_beyond_65'][t]['poisson_tail_bound_m65'])
        # published bound decimal serialized nearest -> compare with relative tolerance
        check(float(true_tail)<float(conservative)*(1+1e-19),f'poisson_tail_{t}',results,float(true_tail))
    summary={'date':'2026-10-14','independent_verifier':'144 labelled generator + branch-wise sparse uniformization/exp propagation; does not import positive_process',
        'checks_total':len(results),'checks_passed':sum(x['passed'] for x in results),
        'max_144_vs_12_joint_coordinate_error':max_joint,'max_full_parent_quotient_entry_error':max_B,
        'scope':'numerical finite-model checks, not independent physical evidence or all-operation outward interval arithmetic',
        'checks':results}
    (R/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(f"{summary['checks_passed']}/{summary['checks_total']} independent checks PASS; max joint error {max_joint:.3g}")
    if summary['checks_passed']!=summary['checks_total']:raise SystemExit(1)
if __name__=='__main__':main()
