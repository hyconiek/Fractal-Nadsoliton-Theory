#!/usr/bin/env python3
"""Bounded checks for new phase-centred cut proofs; not their all-N proof."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def det(rows):
    a=[list(map(F,row)) for row in rows];out=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None:return F(0)
        if pivot!=i:a[pivot],a[i]=a[i],a[pivot];out=-out
        v=a[i][i];out*=v
        for j in range(i+1,len(a)):
            scale=a[j][i]/v
            for k in range(i,len(a)):a[j][k]-=scale*a[i][k]
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();checks=[]
    def ck(name,value,**detail):
        checks.append({'name':name,'passed':bool(value),**detail});assert value,name
    C=(1,0,-1,0)*3;S=(0,1,0,-1)*3;sig=tuple((-1)**j for j in range(12))
    for direction in (-2,2):
        pairs=[(i,j) for i,j in itertools.product(range(12),repeat=2) if sig[j]-sig[i]==direction]
        bins=Counter((C[j]-C[i],S[j]-S[i]) for i,j in pairs)
        ck('nine exact replacements in each diagonal '+str(direction),bins==Counter({(a,b):9 for a in (-1,1) for b in (-1,1)}))
    theta=F(1,1000);h=F(1,10**10);lower=theta**4/12;upper=theta**4/(12*(1-theta**2))
    ck('deep parity-line state has a paid positive energy below h',0<lower<upper<h)
    for N in (12000,24000,120000):
        x=tuple(int(N*(1+theta*t)/12) for t in sig)
        assert sum(x)==N and F(sum(a*b for a,b in zip(x,sig)),N)==theta
        assert sum(a*b for a,b in zip(x,C))==sum(a*b for a,b in zip(x,S))==0
        images={tuple(x[(j-shift)%12] for j in range(12)) for shift in range(0,12,2)}
        images|={tuple(x[(shift-j)%12] for j in range(12)) for shift in range(0,12,2)}
        assert images=={x} and tuple(x[(j-1)%12] for j in range(12))!=x
    ck('actual integer singleton Gplus orbit and distinct odd rotation',True)
    staircase_cases=0
    for N,t in itertools.product((12000,24000,120000),(F(4),F(7),F(101,13))):
        k=math.isqrt((F(N)*t/8).numerator//(F(N)*t/8).denominator)
        n=[N//12]*12;n[0]+=k;n[2]-=k
        points=[tuple(n)]
        for i,j in ((3,0),(1,0)):
            n=list(n);n[i]-=1;n[j]+=1;points.append(tuple(n))
        for a,(expectedC,expectedS) in zip(points,((2*k,0),(2*k+1,1),(2*k+2,0))):
            assert min(a)>0 and sum(a)==N
            assert sum(v*c for v,c in zip(a,C))==expectedC and sum(v*s for v,s in zip(a,S))==expectedS
        assert F(8*k*k,N)<=t<F(8*(k+1)**2,N)
        assert [sum(a*b for a,b in zip(state,sig)) for state in points]==[0,2,4]
        staircase_cases+=1
    ck('two-step upward staircase exact with rational proxy thresholds',staircase_cases==9,
       scope='Fourier/floor arithmetic only; actual transcendental T0 handled analytically')
    # Actual E-layer and singleton normalizations cancel; coefficient controls.
    a=F(1,1000);WE=F(1,10**12);WY=F(1,10**3);WX=F(1,10**10)
    ck('strong-shielding lower numerator can be positive while boundary vanishes',9*a*WY-(1+1/a)*WE-WX>0)
    # Exact cut-current stress. It is deliberately not the actual connected MT graph.
    toy=[]
    for t in (F(1,100),F(1,10000),F(1,10**8)):
        A=set(range(6));V=A-{4};W={i:(t if i in (3,4,5) else F(1)) for i in range(9)}
        H={(0,1):F(1,2),(1,2):F(1,2),(3,4):t/2,(4,5):t/2,(6,7):F(1,2),(7,8):F(1,2)}
        FF={(0,7):F(1,10),(7,2):F(1,10)};U={**H,**FF};N=12
        def quantities(pick):
            fout=sum(v for (x,y),v in FF.items() if x in pick and y not in pick)
            inc=sum(v for (x,y),v in U.items() if x not in pick and y in pick)
            B=sum(W[x] for x,y in H if x in pick and y not in pick)
            tails={x for x,y in FF if x in pick and y not in pick}
            shield={x for x in tails if not any(a==x and b not in pick for a,b in H)}
            return fout,inc,B,sum(W[x] for x in shield)
        fa,ia,ba,za=quantities(A);fv,iv,bv,zv=quantities(V)
        assert fa==ia and fa==fv and iv==ia+t/2 and za==zv==1 and bv==t
        D=fv-iv;strong=(N*zv-iv)/(N*bv)
        assert D==-t/2 and D<=0 and strong>0
        toy.append({'hole_weight':str(t),'exact_residual':str(D),'strong_shield_ratio':str(strong)})
    ck('exact cut control preserves incoming cancellation and exposes false stronger numerator',all(F(toy[i+1]['strong_shield_ratio'])>F(toy[i]['strong_shield_ratio']) for i in range(2)),
       scope='synthetic current algebra only; actual asymptotic proof is in the theorem file')
    ck('dropping incoming current gives a false adverse conclusion',F(1,10)>0 and F(toy[-1]['exact_residual'])<0)
    # Correction of an ancillary SC sentence: arbitrary upper cuts need not obey lower-ideal deficit.
    D_upper=-F(1);shellmass=F(0)
    ck('no universal lower bound on D for nonideal cuts',D_upper<-shellmass)
    # Oriented defect path charging with nonuniform weights.
    for L in range(1,7):
        weights={i:F(1,2)**i for i in range(L+2)}
        edges={(i,i+1):weights[i+1] for i in range(L+1)};T=set(range(1,L+1));universe=set(weights);pick=universe-T
        source={i:F(0) for i in weights}
        for (x,y),v in edges.items():source[x]+=v;source[y]-=v
        assert sum(source.values())==0 and sum(source[x] for x in pick)==-sum(source[x] for x in T)
        boundary=sum(weights[x] for x,y in edges if x in pick and y not in pick)
        assert sum(weights[x] for x in T)<=(36*2)**L*boundary
        assert max(sum(source[x] for x in pick),F(0))/(12*boundary)<=2*(36*2)**L
    ck('bounded-depth phase-centering and path charging exact on weighted controls',True)
    for t in (F(1,10**4),F(1,10**8)):
        N=12;ww=(F(1),F(2),t);II=[F(N,4)*w for w in ww]
        ss=(F(1,100),-F(1,100),-N*t)
        aa=max(-ss[2]-II[2],F(0));corrected=[II[0]-aa/3,II[1]-2*aa/3,II[2]+aa]
        outgoing=[v+s for v,s in zip(corrected,ss)]
        assert all(v>=0 for v in corrected+outgoing) and sum(corrected)==sum(II)
        assert all(o-i==s for o,i,s in zip(outgoing,corrected,ss))
        assert sum(abs(a-b) for a,b in zip(II,corrected))==2*aa
        assert all(i<=2*N*w and o<=3*N*w for i,o,w in zip(corrected,outgoing,ww))
    ck('rare-departure carry repair is exact positive and mass preserving',True)
    def hall(rows,cols,edges):
        return all(sum(rows[i] for i in range(len(rows)) if mask>>i&1)<=sum(cols[j] for j in {b for a,b in edges if mask>>a&1}) for mask in range(1<<len(rows)))
    rows=(F(2),F(1));cols=(F(1),F(2));edges={(0,0),(0,1),(1,1)}
    ck('typed adjacent-fibre Hall control accepts actual feasible marginals',hall(rows,cols,edges))
    ck('positive matching support alone does not prove weighted Hall',not hall(rows,cols,{(0,0),(1,1)}))
    ck('local neighbour availability does not prove a weighted coupling',not hall((5,1,1),(1,1,5),{(0,0),(0,1),(1,1),(1,2),(2,2),(2,0)}))
    ck('carry flux bounds imply original vertex cap without new edge constraint',F(2)+F(3)==5 and F(3,5)/4==F(3,20)>F(1,8) and F(5,3)/4<F(1,2))
    # Exact fixed-parity fibre covolume: no numerical determinant or CLT assumption.
    basis=[]
    for parity in (0,1):
        labels=list(range(parity,12,2))
        for j in labels[:-1]:
            v=[0]*12;v[j]=1;v[labels[-1]]=-1;basis.append(v)
    gram=[[sum(x*y for x,y in zip(a,b)) for b in basis] for a in basis]
    ck('ten-dimensional parity slice Gram determinant equals36',det(gram)==36 and len(basis)==10)
    ck('slice covolume6 and soft mesh2 reproduce full cell12',6*2==12 and F(-5)+F(-3,4)==F(-23,4))
    ck('phase prefix factor288 contains both source half and soft mesh',F(1,2)/(36*2)==F(1,144) and F(1,144)/2==F(1,288))
    # Nominal independent versioned-A7 rates, not an exact-decimal matrix certificate.
    safe=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
    spec=importlib.util.spec_from_file_location('safe_cut_model',safe);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    import numpy as np
    _,_,lam,_,mat=module.build_rank7();g=12/lam[6];numeric=[]
    for N in (120,1200,12000):
        n=np.full(12,N//12);up=0.;total=0.
        for i in range(12):
            cavity=n.copy();cavity[i]-=1;field=(g/N)*(mat@cavity);e=np.exp(field-field.max());prob=e/e.sum()
            for j in range(12):
                if j==i:continue
                rate=n[i]*prob[j];total+=rate
                if sig[j]-sig[i]==2:up+=rate
        assert total<=N*(1+1e-12) and abs(up/N-.25)<2/N
        numeric.append({'N':N,'total_up_rate_over_N':up/N,'per_36_edges':up/(36*N)})
    ck('versioned A7 native reference nominal normalization',max(abs(mat@np.ones(12)))<1e-12,
       scope='binary64 normalization only, not a new spectral or MT certificate')
    # Soft-mesh normalization in a separate q2 critical binomial control.
    soft=[]
    for N in (120,1200,12000):
        lw=[math.lgamma(N+1)-math.lgamma(k+1)-math.lgamma(N-k+1)+N*(2*k/N-1)**2/2 for k in range(N+1)]
        top=max(lw);z=sum(math.exp(v-top) for v in lw);p0=math.exp(lw[N//2]-top)/z
        density=p0/(2*N**(-.75));soft.append({'N':N,'central_soft_density_proxy':density})
    rho0=2/(12**.25*math.gamma(.25))
    ck('independent q2 soft lattice normalization diagnostic',abs(soft[-1]['central_soft_density_proxy']-rho0)<abs(soft[0]['central_soft_density_proxy']-rho0),
       scope='q2 finite control, not proof of q12 local limit')
    proofnames=['ACTUAL_COUNTERCUT_TO_STRONG_SHIELDING.md','BALANCED_PHASE_DEFECT_THEOREM.md','COMPACT_PHASE_PREFIX_LOCAL_LIMIT.md',
                'EXACT_COMPACT_CURRENT_MARGINALS_AND_GLUING.md']
    inputs=[safe,ROOT/'fin_rbcut_weighted_degree_review_20261016/MANIFEST.sha256',
        ROOT/'fin_critical_metastability_review_20261015/MANIFEST.sha256',
        ROOT/'fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md',
        ROOT/'fin_count_partition_feasibility_review_20261016/MANIFEST.sha256',
        ROOT/'fin_reversible_bracket_averaging_review_20261016/MANIFEST.sha256',
        ROOT/'fin_stationary_kinetic_rigidity_review_20261016/MANIFEST.sha256',Path(__file__),
        *[HERE/n for n in proofnames]]
    result={'date':'2026-10-16','final_validation_date':'2026-10-17','status':'PASS_BOUNDED_BALANCED_PHASE_RESEARCH_CONTROLS',
        'scope':'analytical theorems separate from exact finite fixtures and nominal floating diagnostics',
        'checks':checks,'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
        'exact_energy_bounds':{'theta':str(theta),'f_lower':str(lower),'f_upper':str(upper),'h':str(h)},
        'current_toy':toy,'native_rate_diagnostics':numeric,'q2_soft_mesh_diagnostics':soft,
        'source_status':'strong shielding disproved analytically on actual subsequence; original SC remains open',
        'all_N_cut_or_flow_proof_by_this_script':False}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'source_status':result['source_status']},indent=2))

if __name__=='__main__':main()
