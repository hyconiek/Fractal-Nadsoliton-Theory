"""Independent replay of the Target-P partial tree and its accepted leaves.

No numerical basis generation: the saved rational bases are checked directly.
Weight powers are enclosed anew, retaining the true sqrt(3) exponent.
"""
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json
import math
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parent))
from review import HERE,ROOT,SOURCE,sr,load,save,interval_pair
sr.iv.dps=45
TAU=F(67,250)

def cell(row):return tuple(tuple(map(F,p)) for p in row['cell'])
def volume(box):return math.prod(b-a for a,b in box)
def assembled():
    base=load('checkpoints/R7N-020_trace_e2_compression_v1.json')
    a=load('checkpoints/R7N-021_refine_once_v1.json')
    b=load('checkpoints/R7N-021_t_refine_once_v2.json')
    c=load('checkpoints/R7N-021_t_refine_second_cheap_v1.json')
    assert all(not x['pending_parents'] and x['complete'] for x in [a,b,c])
    safe=base['trace_terminals']+base['e2_terminals']+base['compression_terminals']+a['refined_safe']+b['refined_safe']+c['refined_safe']
    failed=c['refined_failed'];hull=tuple(tuple(map(F,p)) for p in base['root_hull'])
    trie={}
    for row in safe+failed:
        path=row['path'];node=trie;assert len(path)%2==0
        for k in range(0,len(path),2):
            assert 'leaf' not in node
            token=path[k:k+2];assert token[0] in '0123' and token[1] in 'LR'
            node=node.setdefault(token,{})
        assert not node;node['leaf']=cell(row)
    def walk(node):
        if 'leaf' in node:
            assert len(node)==1;return node['leaf']
        assert len(node)==2
        axes={int(k[0]) for k in node};assert len(axes)==1;axis=next(iter(axes))
        left=walk(node[str(axis)+'L']);right=walk(node[str(axis)+'R'])
        assert left[axis][1]==right[axis][0]
        assert all(left[k]==right[k] for k in range(4) if k!=axis)
        result=list(left);result[axis]=(left[axis][0],right[axis][1]);return tuple(result)
    assert walk(trie)==hull
    total=volume(hull);residual=sum((volume(cell(x)) for x in failed),F(0))
    assert sum((volume(cell(x)) for x in safe),F(0))+residual==total
    geometry=dict(safe_leaves=len(safe),residual_leaves=len(failed),residual_fraction=str(residual/total),
                  residual_decimal=float(residual/total),coverage='Exact path-driven binary partition; nonmidpoint split endpoints checked.')
    save('partial_geometry.json',geometry)
    return safe,geometry

GRID=10**30
def floorq(x):return F((x*GRID).__floor__(),GRID)
def ceilq(x):return F((x*GRID).__ceil__(),GRID)
@lru_cache(maxsize=50000)
def power(x,kind):
    exponent=sr.I(F(1,2)) if kind=='sqrt' else sr.I(2)+(1 if kind=='plus' else -1)*sr.iv.sqrt(3)
    a,b=sr.bounds(sr.I(x)**exponent);return floorq(a),ceilq(b)
def weights(c):
    (r0,r1),(s0,s1),(t0,t1),(y0,y1)=c
    rl=power(r0,'sqrt')[0];rh=power(r1,'sqrt')[1]
    pl=power(t0,'plus')[0];ph=power(t1,'plus')[1]
    ml=power(t0,'minus')[0];mh=power(t1,'minus')[1]
    raw=[(F(1),F(1)),(2*s0*t0**3,2*s1*t1**3),(r0*t0**4,r1*t1**4),
            (2*r0*s0*t0,2*r1*s1*t1),(2*rl*t0**2*y0,2*rh*t1**2*y1),
            (2*rl*s0*ml*y0,2*rh*s1*mh*y1),(2*rl*s0*pl*y0,2*rh*s1*ph*y1)]
    return [(floorq(a),ceilq(b)) for a,b in raw]

def pair_upper(w,i,j):
    li,ui=w[i];lj,uj=w[j];R=sum((a for k,(a,b) in enumerate(w) if k not in [i,j]),F(0))
    f=lambda x,y:x*y/(x+y+R)**2
    vals=[f(x,y) for x in [li,ui] for y in [lj,uj]]
    vals.extend(f(x,min(max(x+R,lj),uj)) for x in [li,ui])
    vals.extend(f(min(max(y+R,li),ui),y) for y in [lj,uj])
    return ceilq(max(vals))

def triple_upper(w,ids):
    R=sum((a for k,(a,b) in enumerate(w) if k not in ids),F(0));best=F(0)
    for k in ids:
        j,h=[i for i in ids if i!=k];ly,uy=w[j];lz,uz=w[h]
        for x in w[k]:
            C=R+x;f=lambda y,z:x*y*z/(C+y+z)**3
            vals=[f(y,z) for y in [ly,uy] for z in [lz,uz]]
            if ly<=C<=uy and lz<=C<=uz:vals.append(f(C,C))
            vals.extend(f(y,min(max((C+y)/2,lz),uz)) for y in [ly,uy])
            vals.extend(f(min(max((C+z)/2,ly),uy),z) for z in [lz,uz])
            best=max(best,*vals)
    return ceilq(best)

def setup():
    data=json.loads((ROOT/'fin_handoff_audit/results.json').read_text())['exact']['laplacian_intervals']
    assert data==load('inputs/FR223_20260916/inputs/fin_handoff_audit/results.json')['exact']['laplacian_intervals']
    L=[interval_pair(row) for row in data];scales=[sr.iv.sqrt(L[k]/n) for k,n in [(3,6),(4,6),(5,6),(6,12)]]
    X=[]
    for j in [0,4,6,2,3,5,1]:
        X.append([scales[a]*sr.iv.cos(2*sr.iv.pi*k*j/12) for a,k in enumerate([3,4,5])]+[scales[3]*((-1)**j)])
    dot=lambda a,b:sum((x*y for x,y in zip(a,b)),sr.I(0))
    d2={};area={}
    for i,j in combinations(range(7),2):
        v=[x-y for x,y in zip(X[i],X[j])];d2[i,j]=ceilq(sr.bounds(dot(v,v))[1])
    for i,j,k in combinations(range(7),3):
        v=[x-y for x,y in zip(X[j],X[i])];u=[x-y for x,y in zip(X[k],X[i])]
        area[i,j,k]=ceilq(sr.bounds(dot(v,v)*dot(u,u)-dot(v,u)**2)[1])
    return X,d2,area

def compression(row,w,X):
    B=[[F(int(x),int(row['basis_den'])) for x in r] for r in row['basis_num']]
    def det(A):return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
    assert any(det([B[k] for k in inds]) for inds in combinations(range(4),3))
    Z=[[sum((sr.I(B[a][b])*x[a] for a in range(4)),sr.I(0)) for b in range(3)] for x in X]
    K=[[sr.I(TAU*sum(B[a][i]*B[a][j] for a in range(4))) for j in range(3)] for i in range(3)]
    if 'center_c' in row:
        c=list(map(F,row['center_c']));D=interval_pair((sum(a for a,b in w),sum(b for a,b in w)))
        for a in range(3):
            for b in range(a,3):
                v=sum((interval_pair(w[i])*(Z[i][a]-sr.I(c[a]))*(Z[i][b]-sr.I(c[b])) for i in range(7)),sr.I(0))/D
                K[a][b]-=v;K[b][a]=K[a][b]
    else:
        denom=sum(b for a,b in w)**2
        for i,j in combinations(range(7),2):
            p=interval_pair((w[i][0]*w[j][0]/denom,pair_upper(w,i,j)))
            v=[a-b for a,b in zip(Z[i],Z[j])]
            for a in range(3):
                for b in range(a,3):K[a][b]-=p*v[a]*v[b];K[b][a]=K[a][b]
    pivots=[K[0][0],K[0][0]*K[1][1]-K[0][1]**2,det(K)]
    lower=[sr.bounds(x)[0] for x in pivots]
    gersh=[sr.bounds(K[i][i])[0]-sum(max(abs(v) for v in sr.bounds(K[i][j])) for j in range(3) if i!=j) for i in range(3)]
    return min(lower)>0 or min(gersh)>0

def run():
    if (HERE/'partial_leaf_replay.json').exists() and not (HERE/'partial_leaf_initial_progress.json').exists():
        save('partial_leaf_initial_progress.json',json.loads((HERE/'partial_leaf_replay.json').read_text()))
    safe,geometry=assembled();X,d2,areas=setup();failed=[];counts=Counter();start=time.monotonic()
    for i,row in enumerate(safe):
        w=weights(cell(row));reason=row.get('reason','')
        if reason=='SAFE_BY_TRACE':ok=sum(pair_upper(w,*ij)*d for ij,d in d2.items())<=2*TAU
        elif reason=='SAFE_BY_E2':ok=sum(triple_upper(w,ids)*a for ids,a in areas.items())<=TAU**2
        else:ok=compression(row,w,X)
        counts[reason]+=1
        if not ok:failed.append(dict(index=i,path=row['path'],reason=reason))
        if (i+1)%500==0 or i==len(safe)-1:
            out=dict(processed=i+1,total=len(safe),failed=failed,complete=i==len(safe)-1,
                     arithmetic='Outward rational 10^-30 weight/geometry bounds; 45-digit interval matrix arithmetic.',
                     geometry=geometry,method_counts=dict(counts),seconds=time.monotonic()-start,
                     scope='Partial Target-P domain only; no global ceiling.')
            save('partial_leaf_replay.json',out);print('Target P',i+1,'/',len(safe),'failures',len(failed),flush=True)

if __name__=='__main__':run()
