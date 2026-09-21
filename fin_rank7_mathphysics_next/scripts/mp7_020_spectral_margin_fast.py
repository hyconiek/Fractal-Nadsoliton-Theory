#!/usr/bin/env python3
import os
from fractions import Fraction as F
from pathlib import Path
from math import isqrt, nextafter, inf
import json, collections, time
R7N=Path(os.environ.get('R7N_ROOT','/mnt/data/r7n_source/extracted/FIN_R7N_HANDOFF_20260920'))
R7O3=Path(os.environ.get('R7O3_ROOT','/mnt/data/r7o3_source/extracted/FIN_R7O3_TARGETP_HANDOFF_20260920'))
OUT=Path(os.environ.get('MP7_WORK_ROOT','/mnt/data/fin_rank7_mathphysics_next'))/'results/MP7-020_spectral_margin.json'
TAU=F(67,250); SCALE=10**30

def sqrt_upper_rational(x:F, scale=SCALE):
    A=x.numerator*scale*scale; q=isqrt(A//x.denominator)
    while F(q*q,scale*scale)<x:q+=1
    return F(q,scale)

def prev_float_fraction(x):
    f=float(x)
    return F.from_float(nextafter(f,-inf))

def trG_from_record(x):
    den=int(x['basis_den']); B=x['basis_num']
    return sum(F(int(B[r][k])**2,den*den) for r in range(4) for k in range(3))

def witness_margin_saved(x, exact_d3=False, exact_gersh=False, exact_G=None):
    trG=sum(exact_G[i][i] for i in range(3)) if exact_G is not None else trG_from_record(x)
    candidates=[]
    if exact_gersh:
        gl=min(F(v) for v in x['pd_bounds_exact']['gersh_lower_bounds'])
    else:
        gl=prev_float_fraction(x.get('gersh_lower',-1.0))
    if gl>0:candidates.append(('gersh',gl/trG))
    if exact_d3:
        d3lo=F(x['pd_bounds_exact']['d3'][0])
    else:
        d3=x.get('d3')
        d3lo=prev_float_fraction(d3[0]) if d3 and d3[0]>0 else F(0)
    if d3lo>0:
        # K is SPD by the saved accepted certificate. Since K=tau G-E and E>=0,
        # tr K <= tau tr G. For SPD 3x3: lambda_min >= 4 det(K)/tr(K)^2.
        lmin=4*d3lo/(TAU*TAU*trG*trG)
        candidates.append(('det-trace',lmin/trG))
    if not candidates:return None,None
    kind,m=max(candidates,key=lambda kv:kv[1])
    return m,{'bound_kind':kind,'trG':float(trG),'d3_lower':float(d3lo),'gersh_lower':float(gl)}

def newb():return {'count':0,'min':None,'meta':None}
def upd(b,m,meta):
    if m is None:return False
    b['count']+=1
    if b['min'] is None or m<b['min']:b['min']=m;b['meta']=meta
    return True

def proc_reason(bucket,x,label):
    r=x['reason']
    if r=='SAFE_BY_TRACE':m=TAU-F(x['trace_upper'])/2;meta={'path':x['path'],'reason':r}
    elif r=='SAFE_BY_E2':
        q=sqrt_upper_rational(F(x['e2_upper']));m=TAU-q;meta={'path':x['path'],'reason':r,'sqrt_e2_upper':float(q)}
    elif r in ('SAFE_BY_COMPRESSION','SAFE_BY_CENTERED_MOMENT'):
        m,detail=witness_margin_saved(x);meta={'path':x['path'],'reason':r,**(detail or {})}
    else:return None,{'label':label,'path':x.get('path'),'reason':r}
    return (m,meta) if m and m>0 else (None,{'label':label,'path':x.get('path'),'reason':r,'margin':float(m) if m else None})

def main():
    st=time.time();B=collections.defaultdict(newb);fails=[]
    base=json.load(open(R7N/'checkpoints/R7N-020_trace_e2_compression_v1.json'))
    for x in base['trace_terminals']:
        m=TAU-F(x['trace_upper'])/2;upd(B['R7N_TRACE'],m,{'path':x['path'],'reason':'SAFE_BY_TRACE'})
    for x in base['e2_terminals']:
        q=sqrt_upper_rational(F(x['e2_upper']));m=TAU-q;upd(B['R7N_E2'],m,{'path':x['path'],'reason':'SAFE_BY_E2','sqrt_e2_upper':float(q)})
    for x in base['compression_terminals']:
        m,d=witness_margin_saved(x)
        if not m or m<=0:fails.append({'label':'R7N_COMPRESSION','path':x['path'],'detail':d})
        else:upd(B['R7N_COMPRESSION'],m,{'path':x['path'],'reason':'SAFE_BY_COMPRESSION',**d})
    for fn,key in [('R7N-021_refine_once_v1.json','R7N_REFINE_ONCE'),('R7N-021_t_refine_once_v2.json','R7N_T_REFINE_ONCE'),('R7N-021_t_refine_second_cheap_v1.json','R7N_T_REFINE_SECOND')]:
        d=json.load(open(R7N/'checkpoints'/fn))
        for x in d['refined_safe']:
            m,meta=proc_reason(B[key],x,key)
            if m is None:fails.append(meta)
            else:upd(B[key],m,meta)
    # R7O3: exact determinant/Gershgorin endpoints and exact Gram matrix are saved.
    for line in open(R7O3/'certificates/active_leaf_certificates.jsonl'):
        x=json.loads(line); G=[[F(v) for v in row] for row in x['gram_matrix_exact']]
        m,d=witness_margin_saved(x,exact_d3=True,exact_gersh=True,exact_G=G)
        if not m or m<=0:fails.append({'label':'R7O3','leaf_id':x['leaf_id'],'detail':d})
        else:upd(B['R7O3_REPAIRS'],m,{'leaf_id':x['leaf_id'],'original_index':x['original_index'],'pd_method':x['pd_method'],**d})
    total=sum(b['count'] for b in B.values())
    key,bmin=min(((k,v) for k,v in B.items() if v['min'] is not None),key=lambda kv:kv[1]['min'])
    def sb(b):return {'count':b['count'],'min_margin_float':float(b['min']) if b['min'] else None,'min_meta':b['meta']}
    # Accepted unbounded-tail theorems use the sharper sigma_* threshold.  MP7-020
    # does not re-prove those theorems; it audits their exact separation from tau0.
    r17=json.load(open(R7N/'results/R7N-017_target_p_implication.json'))
    tail_sep=F(r17['sigma_separation_lower'])
    fr1=json.load(open(R7N/'inputs/intake_review_20260919/FR1_replay.json'))
    fr42=json.load(open(R7N/'inputs/intake_review_20260919/FR42_replay.json'))
    tail_gate=(fr1['large_J3_tail']['status']=='INTERVAL_CERTIFIED' and
               fr1['large_J4_tail']['status']=='INTERVAL_CERTIFIED' and
               fr1['large_J5_tail']['status']=='INTERVAL_CERTIFIED' and
               fr42.get('status')=='PASS' and r17.get('sigma_lt_tau0') is True)
    global_margin=min(bmin['min'],tail_sep) if tail_gate else None
    out={'task':'MP7-020','tau':'67/250','compact_terminal_count':total,'expected_compact_terminal_count':25656,
         'buckets':{k:sb(v) for k,v in sorted(B.items())},
         'compact_min_margin_exact':str(bmin['min']), 'compact_min_margin_float':float(bmin['min']),
         'compact_lambda2_upper_float':float(TAU-bmin['min']),'compact_min_source':key,'compact_min_meta':bmin['meta'],
         'tail_gate_pass':tail_gate,'tail_sigma_separation_exact':str(tail_sep),'tail_sigma_separation_float':float(tail_sep),
         'tail_sources':['accepted FR1 large-J3/J4/J5 tail theorems','accepted FR42 large-J6 tail theorem'],
         'global_margin_exact':str(global_margin) if global_margin is not None else None,
         'global_margin_float':float(global_margin) if global_margin is not None else None,
         'global_lambda2_upper_float':float(TAU-global_margin) if global_margin is not None else None,
         'global_min_source':key if global_margin==bmin['min'] else ('TAILS' if global_margin is not None else None),
         'failures':fails,'elapsed_seconds':time.time()-st,
         'rigor_note':'For R7N saved float lower endpoints, nextafter(value,-inf) is converted exactly to Fraction, hence is <= the exact rational endpoint that was rounded to that float. For 3D SPD K, tr K<=tau tr G and lambda_min(K)>=4 det(K)/tr(K)^2. lambda_max(G)<=tr G. Tail separation reuses accepted FR1/FR42 theorems and audits only their exact sigma_*-to-tau0 gap.',
         'scope':'complete declared shared-field domain: 25,656-cell compact refined cover plus accepted unbounded FR1/FR42 tails.'}
    OUT.write_text(json.dumps(out,indent=2)); print(json.dumps(out,indent=2))
if __name__=='__main__':main()
