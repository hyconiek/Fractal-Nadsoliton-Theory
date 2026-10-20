#!/usr/bin/env python3
"""Exact finite controls for the reset history/source probe; no source law adopted."""
from collections import Counter, defaultdict
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def block_law(p,h,r,eps):
    out={}
    for seq in itertools.product(range(len(p)),repeat=r):
        base=F(1); tilt=F(1)
        for i in seq: base*=p[i]; tilt*=h[i]
        out[seq]=base*(1+eps*tilt)
    return out


def moment(p,h,power): return sum((p[i]*h[i]**power for i in range(len(p))),F(0))


def mark_moment(indices,p,h,r,eps):
    """Integrate grouped block factors, retaining repeated observations."""
    total=F(0)
    for offset in range(r):
        groups=defaultdict(Counter)
        for i in indices:
            b,position=divmod(i+offset,r)
            groups[b][position]+=1
        ans=F(1)
        for powers in groups.values():
            value=F(1)
            for count in powers.values(): value*=moment(p,h,count)
            if len(powers)==r:
                correction=eps
                for count in powers.values(): correction*=moment(p,h,count+1)
                value+=correction
            ans*=value
        total+=ans/r
    return total


def iid_moment(indices,p,h):
    value=F(1)
    for count in Counter(indices).values(): value*=moment(p,h,count)
    return value


def short_window(p,h,r,eps,length):
    """Brute integration of complete blocks, separate from grouped moments."""
    law=block_law(p,h,r,eps); out=defaultdict(F)
    for offset in range(r):
        nb=(offset+length+r-1)//r
        for blocks in itertools.product(law,repeat=nb):
            weight=F(1,r)
            seq=()
            for b in blocks: weight*=law[b]; seq+=b
            out[seq[offset:offset+length]]+=weight
    return dict(out)


def controls():
    checks=[]
    def ck(name,condition):
        assert condition,name
        checks.append(name)
    p=[F(1,2),F(1,3),F(1,6)]
    h=[F(1,2),-F(1,2),-F(1,2)]
    eps=F(1,2)
    ck('nonuniform_fixture_centered',moment(p,h,1)==0 and sum(p)==1)
    ck('nonuniform_fixture_variance',moment(p,h,2)==F(1,4))
    for r in [3,4,5]:
        law=block_law(p,h,r,eps)
        ck(f'block_{r}_positive_normalized',min(law.values())>0 and sum(law.values())==1)
        for omit in range(r):
            marginal=defaultdict(F)
            for seq,value in law.items(): marginal[seq[:omit]+seq[omit+1:]]+=value
            desired={seq:product(p[i] for i in seq) for seq in itertools.product(range(3),repeat=r-1)}
            ck(f'block_{r}_omit_{omit}_iid',dict(marginal)==desired)
        ck(f'block_{r}_time_reversal',all(law[x]==law[x[::-1]] for x in law))
        full=sum((value*product(h[i] for i in seq) for seq,value in law.items()),F(0))
        ck(f'block_{r}_full_moment',full==eps*moment(p,h,2)**r)
        ck(f'phase_average_{r}_moment',mark_moment(list(range(r)),p,h,r,eps)==full/r)
    window=short_window(p,h,3,eps,4)
    ck('brute_stationary_window_normalization',sum(window.values())==1)
    left=defaultdict(F); right=defaultdict(F)
    for seq,value in window.items(): left[seq[:3]]+=value; right[seq[1:]]+=value
    ck('brute_stationarity_shift',dict(left)==dict(right))
    ck('brute_time_reversal',all(window[x]==window[x[::-1]] for x in window))
    # Sparse finite analysis lift, independent of the window enumeration.
    law3=block_law(p,h,3,eps)
    stationary={(w,j):prob/3 for w,prob in law3.items() for j in range(3)}
    def step(dist):
        out=defaultdict(F)
        for (w,j),prob in dist.items():
            if j<2: out[w,j+1]+=prob
            else:
                for v,pv in law3.items(): out[v,0]+=prob*pv
        return dict(out)
    ck('hidden_lift_stationary_inflow_exact',step(stationary)==stationary)
    start={(next(iter(law3)),1):F(1)}
    for _ in range(3): start=step(start)
    ck('hidden_lift_block_reset_after_r_events',all(start.get((w,1),0)==prob for w,prob in law3.items()) and sum(start.values())==1)
    ck('every_phase_has_one_eligible_refresh_count',all(len([n for n in range(3,6) if (i+n)%3==j])==1 for i in range(3) for j in range(3)))
    for i in range(4):
        for j in range(i+1,4):
            pairs=defaultdict(F)
            for seq,value in window.items(): pairs[(seq[i],seq[j])]+=value
            ck(f'brute_pair_{i}_{j}_iid',all(pairs[a,b]==p[a]*p[b] for a in range(3) for b in range(3)))
    v=moment(p,h,2)
    for a,b in itertools.product(range(4),repeat=2):
        indices=[0,a,a+b]
        diff=mark_moment(indices,p,h,3,eps)-iid_moment(indices,p,h)
        ck(f'Poisson_count_triple_{a}_{b}',diff==(eps*v**3/3 if (a,b)==(1,1) else 0))
    for counts in itertools.product(range(3),repeat=3):
        indices=[0]
        for n in counts: indices.append(indices[-1]+n)
        diff=mark_moment(indices,p,h,4,eps)-iid_moment(indices,p,h)
        ck('Poisson_count_order4_'+'_'.join(map(str,counts)),diff==(eps*v**4/4 if counts==(1,1,1) else 0))
    ck('triple_contrast_coefficient',eps*v**3/3==F(1,384))
    ck('triple_contrast_positive_rational_lower',F(1,384)/9==F(1,3456)>0)
    for i,j,k in itertools.product(range(3),repeat=3):
        ck(f'conditional_mark_prediction_{i}_{j}_{k}',left[(i,j,k)]/(p[i]*p[j])==p[k]*(1+eps*h[i]*h[j]*h[k]/3))
    ck('genuine_non_Markov_prediction',left[(0,0,0)]/(p[0]*p[0])!=p[0])
    # Full instantaneous state-conditioned tensor, from a pair law, not a
    # claim about predictable rates in the whole observed filtration.
    lam=F(2)
    B=[[lam*(p[j]-int(i==j)) for j in range(3)] for i in range(3)]
    gamma=[]
    for i in range(3):
        gi=[[F(0) for _ in range(3)] for _ in range(3)]
        for j in range(3):
            if i==j: continue
            d=[int(a==j)-int(a==i) for a in range(3)]
            for a in range(3):
                for b in range(3): gi[a][b]+=lam*p[j]*d[a]*d[b]
        gamma.append(gi)
    G=[[sum(p[i]*gamma[i][a][b] for i in range(3)) for b in range(3)] for a in range(3)]
    ck('state_conditioned_gamma_matches_reset',all(gamma[i][j][j]==B[i][j] for i in range(3) for j in range(3) if i!=j))
    ck('full_mean_gamma_reset_formula',all(G[i][j]==2*lam*(int(i==j)*p[i]-p[i]*p[j]) for i in range(3) for j in range(3)))
    ck('mean_gamma_inverse_gives_Markov_candidate',all(-G[i][j]/(2*p[i])==B[i][j] for i in range(3) for j in range(3)))
    # Rational proxies for exp(-lambda t) validate kernel semigroup algebra.
    def kernel(z): return [[z*int(i==j)+(1-z)*p[j] for j in range(3)] for i in range(3)]
    A,C=kernel(F(2,5)),kernel(F(3,7)); target=kernel(F(6,35))
    ck('pair_kernel_semigroup_exact',all(sum(A[i][k]*C[k][j] for k in range(3))==target[i][j] for i in range(3) for j in range(3)))
    z=F(2,5); P=kernel(z)
    ck('backward_martingale_drift_cancellation_exact',all(
       -lam*z*(p[j]-int(i==j))+sum(B[i][k]*P[k][j] for k in range(3))==0
       for i in range(3) for j in range(3)))
    # Twelve-state positive rational control; canonical contrast preserves
    # every permutation stabilizing p.
    p12=[F(i,78) for i in range(1,13)]; center=sum(x*x for x in p12)
    h12=[x-center for x in p12]
    t12=block_law(p12,h12,3,eps)
    ck('q12_canonical_contrast_centered',moment(p12,h12,1)==0 and moment(p12,h12,2)>0)
    ck('q12_full_block_positive_normalized',min(t12.values())>0 and sum(t12.values())==1)
    marginal=defaultdict(F)
    for seq,value in t12.items(): marginal[seq[:2]]+=value
    ck('q12_pair_marginal_exact',all(marginal[i,j]==p12[i]*p12[j] for i in range(12) for j in range(12)))
    # Scope countercontrols: a noncentered h does not preserve pair marginals;
    # a constant uniform-p canonical contrast is zero and cannot be a witness.
    bad=block_law(p,[F(1)]*3,3,eps)
    ck('noncentered_contrast_rejected_as_probability',sum(bad.values())!=1)
    uniform=[F(1,12)]*12
    hc=[x-sum(y*y for y in uniform) for x in uniform]
    ck('uniform_canonical_contrast_is_zero_negative_control',all(x==0 for x in hc))
    # Native input binding is nominal numerical diagnostic only.
    import numpy as np
    spec=importlib.util.spec_from_file_location('seed_source_probe',ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py')
    seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,_,_,A7=seed.build_rank7()
    native=seed.pair_difference_distribution(A7,seed.G_FROZEN)
    hn=native-float(native@native); vn=float(native@(hn*hn))
    ck('nominal_native_p_nonuniform_positive',float(min(native))>0 and vn>0 and np.max(native)-np.min(native)>0)
    return checks,{'rational_triple_coefficient':'1/384','rational_triple_lower_at_half_time':'1/3456',
                   'native_nominal_p':native.tolist(),'native_nominal_variance_v':vn,
                   'native_nominal_epsilon_half_third_coefficient':float(vn**3/6),
                   'numeric_scope':'binary64 diagnostic, not an outward kernel/correlation certificate'}


def product(values):
    out=F(1)
    for v in values: out*=v
    return out


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=HERE/'RESULTS.json'); a=parser.parse_args()
    checks,details=controls()
    inputs=['AGENTS.md','FIN_Ten_Research_Programs_Monograph.md',
            'fin_intrinsic_entropy_source_review_20261019/verify_intake.py',
            'fin_intrinsic_entropy_source_review_20261019/PROOF_AUDIT_AND_SCOPE.md',
            'fin_post_rp004_history_review_20261020/MEAN_BRACKET_CORRECTION.md',
            'fin_complete_anonymous_history_review_20261020/PROOF_COMPLETIONS_AND_SCOPE.md',
            'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
            'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py']
    result={'date':'2026-10-20','status':'PASS_BOUNDED_RESET_HISTORY_SOURCE_PROBE',
            'bounded_control_count':len(checks),'bounded_controls':checks,'details':details,
            'input_hashes':{n:sha(ROOT/n) for n in inputs},'F01':'OPEN_SOURCE_NOT_DERIVED',
            'model_status':'DIAGNOSTIC_COUNTERFAMILY_NOT_NATIVE_FIN','formal_proof':False,'physical_validation':False}
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':len(checks),'status':result['status'],'F01':result['F01'],'output':str(a.output)},indent=2))


if __name__=='__main__': main()
