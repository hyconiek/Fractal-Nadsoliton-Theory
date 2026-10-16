#!/usr/bin/env python3
"""Exact finite activity/noise research fixtures; no Monte Carlo or fitting.

Standard-library rational checks are separate from the optional native-A7
binary64 diagnostic. Analytical theorems are not proved by flag counts.
"""
import argparse
import ast
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def counts(N,q):
    if q==1:return [(N,)]
    return [(i,)+tail for i in range(N+1) for tail in counts(N-i,q-1)]
def multinomial(n):return math.factorial(sum(n))//math.prod(math.factorial(i) for i in n)
def move(n,i,j):
    m=list(n); m[i]-=1; m[j]+=1; return tuple(m)
def edge(x,y):return (x,y) if x<y else (y,x)
def occupation_drift(Q,n):
    N=sum(n)
    return [sum(v*F(y[i]-n[i],N) for y,v in Q[n].items()) for i in range(len(n))]
def gamma(Q,n):
    N=sum(n); q=len(n)
    return [[sum(v*F((y[i]-n[i])*(y[j]-n[j]),N*N) for y,v in Q[n].items()) for j in range(q)] for i in range(q)]

def exact_hb(N,q,t):
    states=counts(N,q)
    W={n:F(multinomial(n))*t**sum(i*i for i in n) for n in states}
    Q={n:{} for n in states}
    for n in states:
        for i in range(q):
            if n[i]==0:continue
            r=list(n); r[i]-=1
            weights=[t**(2*a) for a in r]; z=sum(weights)
            for j in range(q):
                if i!=j:Q[n][move(n,i,j)]=n[i]*weights[j]/z
    return W,Q

def wheel(centre,labels):
    q=len(centre); a,b,c=labels
    roots=[]
    for i,j in ((a,b),(a,c),(b,c),(b,a),(c,a),(c,b)):
        v=[0]*q; v[i]+=1; v[j]-=1; roots.append(tuple(v))
    ring=[tuple(centre[i]+v[i] for i in range(q)) for v in roots]
    assert all(min(n)>=0 for n in ring)
    stress=defaultdict(F)
    for j,n in enumerate(ring):
        stress[edge(centre,n)]+=1; stress[edge(n,ring[(j+1)%6])]-=1
    return dict(stress)

def d12_average(stress,q):
    out=defaultdict(F)
    for sign,shift in itertools.product((-1,1),range(q)):
        def transform(n):
            target=[0]*q
            for i,v in enumerate(n):target[(sign*i+shift)%q]=v
            return tuple(target)
        for (x,y),v in stress.items():out[edge(transform(x),transform(y))]+=v/F(2*q)
    return {e:v for e,v in out.items() if v}

def adjacency(stress):
    out=defaultdict(dict)
    for (x,y),v in stress.items():out[x][y]=v; out[y][x]=v
    return out

def safe_amplitude(W,Q,stress,factor):
    B=[]; adj=adjacency(stress)
    for (x,y),v in stress.items():
        if v<0:B += [W[x]*Q[x][y]/(2*abs(v)),W[y]*Q[y][x]/(2*abs(v))]
    for x,row in adj.items():
        for i,n_i in enumerate(x):
            if not n_i:continue
            outgoing=sum(v for y,v in Q[x].items() if y[i]==x[i]-1)
            positive=sum(v for y,v in row.items() if v>0 and y[i]==x[i]-1)
            assert outgoing<n_i
            if positive:B.append(W[x]*(n_i-outgoing)/(2*positive))
    assert B and min(B)>0
    return factor*min(B)

def perturb(W,Q,stress,s):
    out={x:dict(row) for x,row in Q.items()}
    for (x,y),v in stress.items():
        out[x][y]+=s*v/W[x]; out[y][x]+=s*v/W[y]
    return out

def reconstruct(W,Gamma,N,allowed_edges=None):
    """Endpoint inversion, rejecting incompatible fibres instead of repairing."""
    if type(N) is not int or N<1 or not W:raise ValueError('INVALID_POPULATION_OR_EMPTY_DOMAIN')
    D=set(W); first=next(iter(D))
    if not isinstance(first,tuple):raise ValueError('INVALID_COUNT_COORDINATES')
    q=len(first); Q={n:{} for n in D}; terminal=0
    if q<2 or set(Gamma)!=D:raise ValueError('INVALID_ALPHABET_OR_TENSOR_DOMAIN')
    for n in D:
        if not isinstance(n,tuple) or len(n)!=q or any(type(v) is not int for v in n):raise ValueError('INVALID_COUNT_COORDINATES')
        if type(W[n]) not in (int,F):raise ValueError('NONEXACT_WEIGHT')
        if W[n]<=0 or min(n)<0 or sum(n)!=N:raise ValueError('INVALID_STATE_OR_WEIGHT')
        tensor=Gamma[n]
        if len(tensor)!=q or any(len(row)!=q for row in tensor):raise ValueError('TENSOR_SHAPE')
        if any(type(v) not in (int,F) for row in tensor for v in row):raise ValueError('NONEXACT_TENSOR')
        if any(sum(row)!=0 for row in tensor):raise ValueError('NONZERO_TENSOR_ROW_SUM')
        if any(tensor[i][j]!=tensor[j][i] or (i!=j and tensor[i][j]>0) for i,j in itertools.product(range(q),repeat=2)):raise ValueError('TENSOR_SIGN_OR_SYMMETRY')
    for i in range(q):
        for j in range(i+1,q):
            fibres=defaultdict(list)
            for n in D:fibres[tuple(n[k] for k in range(q) if k not in (i,j))].append(n)
            for nodes in fibres.values():
                nodes.sort(key=lambda n:n[i]); segments=[]
                for n in nodes:
                    if not segments or n[i]!=segments[-1][-1][i]+1:segments.append([n])
                    else:segments[-1].append(n)
                for segment in segments:
                    prev=F(0)
                    for k,n in enumerate(segment):
                        S=-N*N*W[n]*Gamma[n][i][j]
                        if S<0:raise ValueError('NEGATIVE_NODE_SUM')
                        current=S-prev
                        if k==len(segment)-1:
                            if current!=0:raise ValueError('TERMINAL_RESIDUAL')
                            terminal+=1
                        else:
                            if current<0:raise ValueError('NEGATIVE_CONDUCTANCE')
                            y=segment[k+1]
                            if current and allowed_edges is not None and edge(n,y) not in allowed_edges:raise ValueError('KNOWN_FORBIDDEN_EDGE')
                            if current:Q[n][y]=current/W[n]; Q[y][n]=current/W[y]
                        prev=current
    return Q,terminal

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True); args=p.parse_args(); checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail}); assert ok,name
    W,Q=exact_hb(3,3,F(2)); centre=(1,1,1); S=wheel(centre,(0,1,2))
    for x,row in adjacency(S).items():
        assert all(sum(v*(y[i]-x[i]) for y,v in row.items())==0 for i in range(3))
    ck('exact wheel equilibrium stress on all seven affected nodes',len(adjacency(S))==7 and len(S)==12)
    ck('base exact Gibbs detailed balance',all(W[x]*v==W[y]*Q[y][x] for x,row in Q.items() for y,v in row.items()))
    s=safe_amplitude(W,Q,S,F(1,2)); R=perturb(W,Q,S,s)
    ck('positive drift-preserving reversible perturbation',s>0 and all(v>0 and W[x]*v==W[y]*R[y][x] for x,row in R.items() for y,v in row.items()))
    ck('all pointwise coordinate drifts preserved exactly',all(occupation_drift(Q,n)==occupation_drift(R,n) for n in W))
    ck('all per-slot attempt caps paid',all(sum(v for y,v in R[x].items() if y[i]==x[i]-1)<n_i for x in W for i,n_i in enumerate(x) if n_i))
    G={n:gamma(Q,n) for n in W}; H={n:gamma(R,n) for n in W}
    delta=H[centre][0][0]-G[centre][0][0]
    ck('pointwise covariance and centre holding rate genuinely change',delta==4*s/(9*W[centre]) and delta>0 and sum(R[centre].values())>sum(Q[centre].values()))
    for i,j in itertools.product(range(3),repeat=2):assert sum(W[n]*(H[n][i][j]-G[n][i][j]) for n in W)==0
    ck('all nine invariant averaged covariance identities exact',True)
    inv,terminal=reconstruct(W,G,3); invr,_=reconstruct(W,H,3)
    ck('endpoint inverse recovers both distinct full generators exactly',inv==Q and invr==R,terminal_checks=terminal)
    subset={n:w for n,w in W.items() if n not in ((2,1,0),(1,2,0),(1,1,1))}
    reflected={x:{y:v for y,v in Q[x].items() if y in subset} for x in subset}
    inverse,terminal2=reconstruct(subset,{n:gamma(reflected,n) for n in subset},3)
    ck('reflected disconnected fibre segments reconstructed without omitted exits',inverse==reflected,terminal_checks=terminal2)
    W0,Q0=exact_hb(3,3,F(1)); s0=safe_amplitude(W0,Q0,S,F(0))
    ck('exact g0 family reduction under vanishing amplitude',s0==0 and perturb(W0,Q0,S,s0)==Q0)
    # Nonnegative node sums do not imply a positive or compatible edge reconstruction.
    def scalar_nodes(data):
        prev=F(0); edges=[]
        for value in data[:-1]:prev=value-prev; edges.append(prev)
        return edges,data[-1]-prev
    e,bad=scalar_nodes(list(map(F,(1,1,1))))
    ck('incompatible terminal dataset retained as rejection',bad!=0)
    e,bad=scalar_nodes(list(map(F,(2,1,1,2))))
    ck('terminal-compatible but negative conductance dataset rejected',bad==0 and e==[F(2),F(-1),F(2)])
    for data,reason in (((1,1,1),'TERMINAL_RESIDUAL'),((2,1,1,2),'NEGATIVE_CONDUCTANCE')):
        nbad=len(data)-1; states=[(k,nbad-k) for k in range(nbad+1)]
        weights={n:F(1) for n in states}
        tensors={n:[[F(data[k],nbad*nbad),-F(data[k],nbad*nbad)],[-F(data[k],nbad*nbad),F(data[k],nbad*nbad)]] for k,n in enumerate(states)}
        try:reconstruct(weights,tensors,nbad)
        except ValueError as error:assert str(error)==reason
        else:raise AssertionError('negative compatibility fixture accepted')
    ck('actual inverse API fails closed on both incompatible tensors',True)
    bad_calls=[({}, {}, 1),({(0,0):F(1)}, {(0,0):[[F(0),F(0)],[F(0),F(0)]]},0),
        ({(1,):F(1)}, {(1,):[[F(0)]]},1),
        ({(F(1,2),F(1,2)):F(1)}, {(F(1,2),F(1,2)):[[F(0),F(0)],[F(0),F(0)]]},1),
        ({(1,0):float('nan')}, {(1,0):[[F(0),F(0)],[F(0),F(0)]]},1)]
    for weights,tensors,pop in bad_calls:
        try:reconstruct(weights,tensors,pop)
        except ValueError:pass
        else:raise AssertionError('malformed exact inverse fixture accepted')
    try:reconstruct(W,G,3,allowed_edges=set())
    except ValueError as error:assert str(error)=='KNOWN_FORBIDDEN_EDGE'
    else:raise AssertionError('forbidden support silently changed')
    ck('exact inverse rejects malformed domains and forbidden edges',True)
    # Q12 D12 exact fixture: nonzero stress on uniform HB tests algebra, not g0 preservation.
    W12,Q12=exact_hb(3,12,F(1)); n12=(1,1,1)+(0,)*9
    S12=d12_average(wheel(n12,(0,1,2)),12); adj=adjacency(S12)
    ck('D12 averaged stress counts and noncancellation',len(adj)==60 and len(S12)==132 and sum(v>0 for v in S12.values())==72)
    ck('D12 averaged stress annihilates every coordinate exactly',all(all(sum(v*(y[i]-x[i]) for y,v in row.items())==0 for i in range(12)) for x,row in adj.items()))
    R12=perturb(W12,Q12,S12,F(1,16))
    ck('exact q12 positive rates and labelled per-slot cap',all(v>0 for row in R12.values() for v in row.values()) and all(sum(v for y,v in R12[x].items() if y[i]==x[i]-1)<n_i for x in W12 for i,n_i in enumerate(x) if n_i))
    ck('q12 centre activity and covariance values',sum(R12[n12].values())-sum(Q12[n12].values())==F(1,192) and gamma(R12,n12)[0][0]-gamma(Q12,n12)[0][0]==F(1,2592))
    # Beyond single-copy support, covariance no longer identifies the generator.
    perturbations={(0,1):F(1,10),(1,2):F(-1,10),(2,3):F(1,10),(0,3):F(-1,90)}
    ck('multiple-count-jump covariance inverse negative control',all(sum(v*(b-a)**2 for (a,b),v in perturbations.items() if k in (a,b))==0 for k in range(4)) and any(v!=0 for v in perturbations.values()))
    # Parity-preserving edge thinning: a nontrivial regular same-parity-drift class member.
    sigma=[(-1)**i for i in range(12)]
    Rp={x:{y:(F(1,10)*v if sum((y[i]-x[i])*sigma[i] for i in range(12))==0 else v) for y,v in row.items()} for x,row in Q12.items()}
    def mchar(A,n,power):return sum(v*F(sum(sigma[i]*(y[i]-n[i]) for i in range(12)),3)**power for y,v in A[n].items())
    ck('same-parity activity changes without changing parity drift or bracket',Rp!=Q12 and all(mchar(Rp,n,1)==mchar(Q12,n,1) and mchar(Rp,n,2)==mchar(Q12,n,2) for n in W12))
    # Native A7 finite check, reconstructed from versioned constants, floating only.
    seed=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
    constants={node.targets[0].id:ast.literal_eval(node.value) for node in ast.parse(seed.read_text()).body if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ('OMEGA','PHI','ETA')}
    w={d:math.cos(constants['OMEGA']*d+constants['PHI'])/(1+d**constants['ETA']) for d in range(1,7)}
    lam={r:sum((2 if d<6 else 1)*w[d]*(1-math.cos(2*math.pi*r*d/12)) for d in range(1,7)) for r in (3,4,5,6)}
    ar=[(2*sum(lam[r]*math.cos(2*math.pi*r*d/12) for r in (3,4,5))+lam[6]*(-1)**d)/12 for d in range(12)]
    A=[[ar[(i-j)%12] for j in range(12)] for i in range(12)]; gain=3.; Wn={}; Qn={}
    for n in W12:
        Wn[n]=multinomial(n)*math.exp(gain/6*sum(n[i]*A[i][j]*n[j] for i,j in itertools.product(range(12),repeat=2)))
        row={}
        for i,ni in enumerate(n):
            if not ni:continue
            cavity=list(n); cavity[i]-=1
            field=[gain/3*sum(A[j][k]*cavity[k] for k in range(12)) for j in range(12)]; shift=max(field)
            target=[math.exp(v-shift) for v in field]; norm=sum(target)
            for j in range(12):
                if j!=i:row[move(n,i,j)]=ni*target[j]/norm
        Qn[n]=row
    sn=safe_amplitude(Wn,Qn,S12,gain*gain/(1+gain*gain)); Rn=perturb(Wn,Qn,S12,sn)
    drift_error=max(abs(a-b) for n in Wn for a,b in zip(occupation_drift(Qn,n),occupation_drift(Rn,n)))
    db_error=max(abs(Wn[x]*v-Wn[y]*Rn[y][x]) for x,row in Rn.items() for y,v in row.items())
    ck('native A7 q12 N3 nonzero stress diagnostic',sn>0 and drift_error<2e-14 and db_error<2e-10 and all(v>0 for row in Rn.values() for v in row.values()),drift_residual=drift_error,db_residual=db_error)
    ck('native diagnostic per-slot margins remain positive',all(sum(v for y,v in Rn[x].items() if y[i]==x[i]-1)<ni for x in Wn for i,ni in enumerate(x) if ni))
    inputs=['AGENTS.md','FIN_RATE_ACTIVITY_SOURCE_NEXT_CHAIN_20261016.md',
        'fin_rate_coupling_review_20261016/PROOF_AUDIT_AND_ADDITIONAL_CONTROLS.md',
        'fin_rate_coupling_review_20261016/MANIFEST.sha256',
        'fin_critical_metastability_review_20261015/INTAKE_20261015.md',
        'fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md',
        'fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json',
        'fin_uniform_critical_review_20261014/INTAKE_20261014.md',
        'fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md',
        'fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md',
        'fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md',
        'fin_rank7_followup/inputs/fin_handoff_audit/results.json',
        'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py']
    result={'date':'2026-10-16','status':'PASS_BOUNDED_ACTIVITY_COVARIANCE_RESEARCH_CHECKS','implementation_sha256':sha(Path(__file__)),
        'checks':checks,'count':len(checks),'exact_q3':{'N':3,'states':len(W),'stress_edges':len(S),'amplitude':str(s),'centre_variance_change':str(delta)},
        'exact_q12_uniform_fixture':{'N':3,'states':len(W12),'stress_states':len(adj),'stress_edges':len(S12),'amplitude':'1/16',
            'scope':'nonzero algebraic perturbation at uniform law, NOT the g0-preserving family at g0 or the asymptotic MT domain'},
        'native_A7_diagnostic':{'N':3,'gain':gain,'arithmetic':'binary64 diagnostic, NOT interval proof','amplitude':sn,'drift_residual':drift_error,'db_residual':db_error},
        'inputs':{rel:sha(ROOT/rel) for rel in inputs},
        'source_status':'No independent activity provider derived; AS-002/003 NOT_ACTIVATED; new intrinsic theorems only'}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'exact_q3':result['exact_q3'],'exact_q12':result['exact_q12_uniform_fixture'],'native_A7_diagnostic':result['native_A7_diagnostic']},indent=2))

if __name__=='__main__':main()
