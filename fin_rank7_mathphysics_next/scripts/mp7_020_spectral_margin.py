#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json, sys, itertools, collections, time
sys.set_int_max_str_digits(0)

R7N=Path('/mnt/data/r7n_source/extracted/FIN_R7N_HANDOFF_20260920')
R7O3=Path('/mnt/data/r7o3_source/extracted/FIN_R7O3_TARGETP_HANDOFF_20260920')
OUT=Path('/mnt/data/fin_rank7_mathphysics_next/results/MP7-020_spectral_margin.json')
TAU=F(67,250)

# Source package modules have their historical ROOT hard-coded. A temporary symlink is
# created by the invoking shell, never inside either immutable source archive.
sys.path[:0]=[str(R7N/'src'), str(R7N/'inputs/FR223_20260916/src')]
from intervals import QI
import target_p_trace_cover as tc
import target_p_trace_tight_rounded as tr
import compression_interval_probe as old
import compression_centered_moment as cm

SCALE=10**30
def sqrt_upper_rational(x:F, scale=SCALE):
    # q=ceil(sqrt(x)*scale)/scale exactly, using integer arithmetic.
    A=x.numerator*scale*scale
    q=isqrt(A//x.denominator)
    while F(q*q, scale*scale) < x: q+=1
    while q>0 and F((q-1)*(q-1),scale*scale) >= x: q-=1
    return F(q,scale)

def interval_det3(K):
    return (K[0][0]*K[1][1]*K[2][2]+2*K[0][1]*K[0][2]*K[1][2]
            -K[0][0]*K[1][2]*K[1][2]-K[1][1]*K[0][2]*K[0][2]-K[2][2]*K[0][1]*K[0][1])

def exact_margin_from_K(K,G):
    d1=K[0][0]
    d2=K[0][0]*K[1][1]-K[0][1]*K[0][1]
    d3=interval_det3(K)
    sylv=d1.lo>0 and d2.lo>0 and d3.lo>0
    gl=[]
    for a in range(3):
        offsum=sum(max(abs(K[a][b].lo),abs(K[a][b].hi)) for b in range(3) if b!=a)
        gl.append(K[a][a].lo-offsum)
    gersh=min(gl)
    if not sylv and gersh<=0:
        return None, {'sylv':False,'gersh':str(gersh),'d3lo':str(d3.lo)}
    tr_hi=sum(K[i][i].hi for i in range(3))
    lb_det=F(0)
    if d3.lo>0 and tr_hi>0:
        lb_det=4*d3.lo/(tr_hi*tr_hi)
    lb=max(F(0),gersh,lb_det)
    trG=sum(G[i][i] for i in range(3))
    if lb<=0 or trG<=0: return None, {'sylv':sylv,'gersh':str(gersh),'lb_det':str(lb_det),'trG':str(trG)}
    return lb/trG, {'sylv':sylv,'gersh':str(gersh),'lb_det':str(lb_det),'trK_hi':str(tr_hi),'trG':str(trG)}

def G_from_basis(B):
    return [[sum(B[r][i]*B[r][j] for r in range(4)) for j in range(3)] for i in range(3)]

def pairwise_K(cell, basis_num, basis_den):
    B=[[F(int(basis_num[r][k]),basis_den) for k in range(3)] for r in range(4)]
    G=G_from_basis(B)
    K=[[QI(TAU*G[i][j]) for j in range(3)] for i in range(3)]
    wb=tc.weight_bounds(cell)
    for i,j in itertools.combinations(range(7),2):
        Dhi=sum(b for a,b in wb)
        lo=wb[i][0]*wb[j][0]/(Dhi*Dhi)
        hi=tc.pair_prob_product_upper(wb,i,j)
        q=QI(lo,hi)
        dz=[]
        for k in range(3):
            z=QI(0)
            for r in range(4): z += B[r][k]*(old.OBS[i][r]-old.OBS[j][r])
            dz.append(z)
        for a in range(3):
            for b in range(a,3):
                K[a][b]=K[a][b]-q*dz[a]*dz[b]
                if a!=b:K[b][a]=K[a][b]
    return K,G

def centered_K(cell,basis_num,basis_den,center_c):
    B=[[F(int(basis_num[r][k]),basis_den) for k in range(3)] for r in range(4)]
    G=G_from_basis(B)
    c=[F(x) for x in center_c]
    wb=tr.weights(cell); D=QI(sum(a for a,b in wb),sum(b for a,b in wb))
    Z=[]
    for i in range(7):
        row=[]
        for k in range(3):
            z=QI(0)
            for r in range(4): z += B[r][k]*old.OBS[i][r]
            row.append(z-QI(c[k]))
        Z.append(row)
    E=[[QI(0) for _ in range(3)] for __ in range(3)]
    for a in range(3):
        for b in range(a,3):
            N=QI(0)
            for i in range(7): N += QI(wb[i][0],wb[i][1])*Z[i][a]*Z[i][b]
            E[a][b]=N/D; E[b][a]=E[a][b]
    K=[[QI(TAU*G[a][b])-E[a][b] for b in range(3)] for a in range(3)]
    return K,G

def R7O3_margin(d):
    G=[[F(x) for x in row] for row in d['gram_matrix_exact']]
    E=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            lo,hi=d['moment_entry_enclosures'][i][j]; E[i][j]=QI(F(lo),F(hi))
    K=[[QI(TAU*G[i][j])-E[i][j] for j in range(3)] for i in range(3)]
    return exact_margin_from_K(K,G)

def rec_cell(x): return tuple(tuple(map(F,p)) for p in x['cell'])

def compression_margin_record(x, centered=False):
    cell=rec_cell(x)
    Bn=x.get('basis_num'); den=int(x.get('basis_den',100000))
    c=x.get('center_c')
    if centered:
        # Older records save center_c under centered_c in some lanes; if absent,
        # recompute the historical fixed center exactly via the original routine.
        if c is None: c=x.get('centered_c')
        if Bn is None or c is None:
            B,Bn_np,den,eigs,rows,det,cF=cm.center_data(cell,den)
            Bn=Bn_np.tolist(); c=[str(v) for v in cF]
        K,G=centered_K(cell,Bn,den,c)
    else:
        if Bn is None:
            Bn_np,den,eigs,rows,det=old.center_basis(cell,den); Bn=Bn_np.tolist()
        K,G=pairwise_K(cell,Bn,den)
    return exact_margin_from_K(K,G)

def load(p): return json.load(open(p))

def update_best(bucket, m, meta):
    if m is None: return
    bucket['count']+=1
    if bucket['min_margin'] is None or m < bucket['min_margin']:
        bucket['min_margin']=m; bucket['min_meta']=meta

def newbucket(): return {'count':0,'min_margin':None,'min_meta':None}

def serialize_bucket(b):
    return {'count':b['count'], 'min_margin_exact':str(b['min_margin']) if b['min_margin'] is not None else None,
            'min_margin_float':float(b['min_margin']) if b['min_margin'] is not None else None, 'min_meta':b['min_meta']}

def main():
    st=time.time(); buckets=collections.defaultdict(newbucket); failures=[]
    # R7N base lanes
    d=load(R7N/'checkpoints/R7N-020_trace_e2_compression_v1.json')
    for x in d['trace_terminals']:
        m=TAU-F(x['trace_upper'])/2; update_best(buckets['R7N_TRACE'],m,{'path':x['path'],'trace_upper':x['trace_upper']})
    for x in d['e2_terminals']:
        q=sqrt_upper_rational(F(x['e2_upper'])); m=TAU-q
        update_best(buckets['R7N_E2'],m,{'path':x['path'],'e2_upper':x['e2_upper'],'sqrt_upper':str(q)})
    for idx,x in enumerate(d['compression_terminals']):
        m,det=compression_margin_record(x,False)
        if m is None: failures.append(['R7N_COMPRESSION',x['path'],det])
        else:update_best(buckets['R7N_COMPRESSION'],m,{'path':x['path'],**det})
    # first generic refine
    d1=load(R7N/'checkpoints/R7N-021_refine_once_v1.json')
    for x in d1['refined_safe']:
        r=x['reason']
        if r=='SAFE_BY_TRACE': m=TAU-F(x['trace_upper'])/2
        elif r=='SAFE_BY_E2': m=TAU-sqrt_upper_rational(F(x['e2_upper']))
        elif r=='SAFE_BY_COMPRESSION':
            m,det=compression_margin_record(x,False)
            if m is None: failures.append(['R7N_REFINE_ONCE',x['path'],det]); continue
        else: failures.append(['R7N_REFINE_ONCE_UNKNOWN',x['path'],r]); continue
        update_best(buckets['R7N_REFINE_ONCE'],m,{'path':x['path'],'reason':r})
    # t-refine once, centered or pairwise
    dt=load(R7N/'checkpoints/R7N-021_t_refine_once_v2.json')
    for x in dt['refined_safe']:
        r=x['reason']
        if r=='SAFE_BY_TRACE':m=TAU-F(x['trace_upper'])/2
        elif r=='SAFE_BY_E2':m=TAU-sqrt_upper_rational(F(x['e2_upper']))
        elif r=='SAFE_BY_CENTERED_MOMENT':
            m,det=compression_margin_record(x,True)
            if m is None:failures.append(['R7N_T_REFINE_ONCE_CENTERED',x['path'],det]);continue
        elif r=='SAFE_BY_COMPRESSION':
            m,det=compression_margin_record(x,False)
            if m is None:failures.append(['R7N_T_REFINE_ONCE_COMP',x['path'],det]);continue
        else:failures.append(['R7N_T_REFINE_ONCE_UNKNOWN',x['path'],r]);continue
        update_best(buckets['R7N_T_REFINE_ONCE'],m,{'path':x['path'],'reason':r})
    # second t refine is centered moment only in final safe set
    d2=load(R7N/'checkpoints/R7N-021_t_refine_second_cheap_v1.json')
    for x in d2['refined_safe']:
        r=x['reason']
        if r=='SAFE_BY_TRACE':m=TAU-F(x['trace_upper'])/2
        elif r=='SAFE_BY_E2':m=TAU-sqrt_upper_rational(F(x['e2_upper']))
        elif r=='SAFE_BY_CENTERED_MOMENT':
            m,det=compression_margin_record(x,True)
            if m is None:failures.append(['R7N_T_REFINE_SECOND',x['path'],det]);continue
        else:failures.append(['R7N_T_REFINE_SECOND_UNKNOWN',x['path'],r]);continue
        update_best(buckets['R7N_T_REFINE_SECOND'],m,{'path':x['path'],'reason':r})
    # R7O3 repairs: reconstruct K directly from saved exact G and moment intervals.
    with open(R7O3/'certificates/active_leaf_certificates.jsonl') as f:
        for line in f:
            x=json.loads(line); m,det=R7O3_margin(x)
            if m is None:failures.append(['R7O3',x['leaf_id'],det]);continue
            update_best(buckets['R7O3_REPAIRS'],m,{'leaf_id':x['leaf_id'],'parent':x['original_index'],'pd_method':x['pd_method'],**det})
    total=sum(v['count'] for v in buckets.values())
    global_min=None; global_key=None
    for k,v in buckets.items():
        if v['min_margin'] is not None and (global_min is None or v['min_margin']<global_min): global_min=v['min_margin'];global_key=k
    out={'task':'MP7-020','tau':'67/250','compact_terminal_count':total,
         'expected_compact_terminal_count':25656,
         'buckets':{k:serialize_bucket(v) for k,v in sorted(buckets.items())},
         'compact_min_margin_exact':str(global_min) if global_min is not None else None,
         'compact_min_margin_float':float(global_min) if global_min is not None else None,
         'compact_min_source':global_key,
         'compact_lambda2_upper_float':float(TAU-global_min) if global_min is not None else None,
         'failures':failures,'elapsed_seconds':time.time()-st,
         'scope':'compact refined cover only; tails analyzed separately',
         'method':'TRACE: tau-tr/2; E2: tau-ceil_sqrt(e2); 3D witness: lambda_min(K)/lambda_max(G), with lambda_min(K)>=max(Gershgorin,4 det(K)/tr(K)^2) and lambda_max(G)<=tr(G).'}
    OUT.write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
