"""Independent R7P-054 checker for the frozen boundary-Ising cover.

This intentionally does not import boundary_cover.py.  It reads the frozen
R7P-048 polynomial specification and the saved R7P-053 certificate tree,
rebuilds interval Bernstein coefficients, checks every split/leaf, and
rechecks the R7P-052 tangent-cone Hessian certificate.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import hashlib,json,sys
import sympy as sp

class IQ:
    def __init__(self,lo,hi=None):
        self.lo=F(lo);self.hi=F(lo if hi is None else hi)
        if self.lo>self.hi:raise ValueError('reversed interval')
    @staticmethod
    def cast(x):return x if isinstance(x,IQ) else IQ(x)
    def __add__(self,o):o=IQ.cast(o);return IQ(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return IQ(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-IQ.cast(o))
    def __rsub__(self,o):return IQ.cast(o)-self
    def __mul__(self,o):
        o=IQ.cast(o);v=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)];return IQ(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=IQ.cast(o)
        if o.lo<=0<=o.hi:raise ValueError('zero denominator')
        return self*IQ(1/o.hi,1/o.lo)
    def __rtruediv__(self,o):return IQ.cast(o)/self
    def __pow__(self,n):
        if n<0:return IQ(1)/(self**(-n))
        out=IQ(1)
        for _ in range(n):out=out*self
        return out


def eval_i(expr,values):
    expr=sp.sympify(expr)
    if expr.is_Rational:return IQ(F(int(expr.p),int(expr.q)))
    if expr.is_Symbol:return values[expr]
    if expr.is_Add:
        out=IQ(0)
        for x in expr.args:out=out+eval_i(x,values)
        return out
    if expr.is_Mul:
        out=IQ(1)
        for x in expr.args:out=out*eval_i(x,values)
        return out
    if expr.is_Pow and expr.exp.is_Integer:return eval_i(expr.base,values)**int(expr.exp)
    raise ValueError(f'unsupported expression {expr}')


def power_coeffs(expr,vars_,spectral_values):
    P=sp.Poly(expr,*vars_);deg=tuple(P.degree_list());out={}
    for mon,c in P.terms():out[tuple(mon)]=eval_i(c,spectral_values)
    return deg,out


def p2b(power_coeff,degrees):
    out={};n=len(degrees)
    for idx in product(*[range(d+1) for d in degrees]):
        acc=IQ(0)
        for k in product(*[range(i+1) for i in idx]):
            a=power_coeff.get(tuple(k))
            if a is None:continue
            q=F(1)
            for j in range(n):q*=F(comb(idx[j],k[j]),comb(degrees[j],k[j]))
            acc=acc+a*q
        out[tuple(idx)]=acc
    return out


def bounds(c):return IQ(min(x.lo for x in c.values()),max(x.hi for x in c.values()))

def split_line(line):
    rows=[list(line)]
    while len(rows[-1])>1:
        p=rows[-1];rows.append([(p[i]+p[i+1])*F(1,2) for i in range(len(p)-1)])
    n=len(line)-1;L=[rows[k][0] for k in range(n+1)];R=[None]*(n+1)
    for k in range(n+1):R[n-k]=rows[k][-1]
    return L,R

def split_coeff(c,degrees,axis):
    n=len(degrees);others=[j for j in range(n) if j!=axis];L={};R={}
    for oi in product(*[range(degrees[j]+1) for j in others]):
        base=[0]*n
        for j,v in zip(others,oi):base[j]=v
        line=[]
        for a in range(degrees[axis]+1):base[axis]=a;line.append(c[tuple(base)])
        ll,rr=split_line(line)
        for a,x in enumerate(ll):base[axis]=a;L[tuple(base)]=x
        for a,x in enumerate(rr):base[axis]=a;R[tuple(base)]=x
    return L,R

def split_box(box,axis):
    box=list(box);lo,hi=box[axis];mid=(lo+hi)/2;L=list(box);R=list(box);L[axis]=(lo,mid);R[axis]=(mid,hi);return tuple(L),tuple(R)

def ij(I):return [str(I.lo),str(I.hi)]
def bj(box):return [[str(a),str(b)] for a,b in box]

def parse_path(path):
    if len(path)%2:raise ValueError('odd path length')
    out=[]
    for i in range(0,len(path),2):
        ax=int(path[i]);side=path[i+1]
        if ax not in (0,1,2) or side not in 'LR':raise ValueError('bad path token')
        out.append((ax,side))
    return out


def derivative(c,degrees,axis):
    n=degrees[axis]
    if n<=0:raise ValueError('zero degree')
    nd=list(degrees);nd[axis]-=1;nd=tuple(nd);out={}
    for idx in product(*[range(d+1) for d in nd]):
        j=list(idx);k=list(idx);k[axis]+=1
        out[idx]=(c[tuple(k)]-c[tuple(j)])*n
    return out,nd

def scale(c,q):return {k:v*q for k,v in c.items()}

def dyadic_path(lo,hi):
    lo,hi=F(lo),F(hi);w=hi-lo;den=w.denominator
    if w.numerator!=1 or den&(den-1):raise ValueError('non-dyadic cell')
    n=den.bit_length()-1;k=lo*den
    if k.denominator!=1 or hi*den!=k+1:raise ValueError('misaligned cell')
    return ''.join('R' if (int(k)>>(n-1-i))&1 else 'L' for i in range(n))

def restrict(c,degrees,box):
    out=c
    for ax,(lo,hi) in enumerate(box):
        for side in dyadic_path(lo,hi):
            L,R=split_coeff(out,degrees,ax);out=L if side=='L' else R
    return out

def hessian_bounds(A,degrees,box):
    local=restrict(A,degrees,box);w=[b-a for a,b in box]
    def D(axes):
        c=local;d=degrees;q=F(1)
        for ax in axes:c,d=derivative(c,d,ax);q/=w[ax]
        return bounds(scale(c,q))
    return {'Arr':D((0,0)),'Ass':D((1,1)),'Att':D((2,2)),'Ars':D((0,1)),'Art':D((0,2)),'Ast':D((1,2))}
def sq(I):
    if I.lo<=0<=I.hi:return IQ(0,max(I.lo*I.lo,I.hi*I.hi))
    return IQ(min(I.lo*I.lo,I.hi*I.hi),max(I.lo*I.lo,I.hi*I.hi))


def root_data(spec):
    fp=spec['frozen_polynomials'];r,s,t,l3,l4,l5=sp.symbols('r s t l3 l4 l5')
    vars_=(r,s,t);sv=(l3,l4,l5)
    S={x:IQ(*map(F,fp['spectral_intervals'][str(x)])) for x in sv}
    Aexpr=sp.sympify(fp['A'],locals={str(x):x for x in (*vars_,*sv)})
    Bexpr=sp.sympify(fp['B'],locals={str(x):x for x in (*vars_,*sv)})
    dA,pA=power_coeffs(Aexpr,vars_,S);dB,pB=power_coeffs(Bexpr,vars_,S)
    if list(dA)!=fp['A_degrees'] or list(dB)!=fp['B_degrees']:raise ValueError('degree mismatch')
    return dict(A=p2b(pA,dA),B=p2b(pB,dB),dA=dA,dB=dB,exprA=Aexpr,
                symbols=(r,s,t,l3,l4,l5),spectral=S)


def validate_local(spec,tree,root):
    c=tree['local_certificate'];box=tuple((F(a),F(b)) for a,b in c['box'])
    H=hessian_bounds(root['A'],root['dA'],box);Arr,Ass,Att,Ars,Art,Ast=[H[k] for k in ('Arr','Ass','Att','Ars','Art','Ast')]
    Nuu=Arr*Ass-sq(Ars);Nuv=Arr*Ast-(Ars*Art);Nvv=Arr*Att-sq(Art)
    if not (Arr.hi<0 and Nuu.lo>0 and Nuv.lo>0 and Nvv.lo>0):raise ValueError('local Hessian cone certificate failed')
    # Exact double-root value and gradient, recomputed from the frozen A polynomial.
    r,s,t,l3,l4,l5=root['symbols'];A=root['exprA'];tau=sp.symbols('tau',positive=True)
    tau2=(2*l3-l4)*(2*l3-l5)/(4*l3**2);r0=(1-tau)/(1+tau);mod=sp.Poly(tau**2-tau2,tau)
    for name,e in [('A',A),('Ar',sp.diff(A,r)),('As',sp.diff(A,s)),('At',sp.diff(A,t))]:
        q=sp.cancel(e.subs({r:r0,s:1,t:1}));num=sp.together(q).as_numer_denom()[0]
        rem=sp.factor(sp.rem(sp.Poly(num,tau),mod).as_expr())
        if rem!=0:raise ValueError(f'local equality identity failed: {name}')
    T=IQ(*map(F,spec['frozen_polynomials']['t_star_interval']));R=(1-T)/(1+T)
    if not (box[0][0]<=R.lo and R.hi<=box[0][1]):raise ValueError('r_star not in local box')
    # Recomputed records must match stored certificate exactly where supplied.
    if c.get('hessian_bounds')!={k:ij(v) for k,v in H.items()}:raise ValueError('stored local Hessian bounds mismatch')
    if c.get('schur_numerator_intervals')!={'Nuu':ij(Nuu),'Nuv':ij(Nuv),'Nvv':ij(Nvv)}:raise ValueError('stored local Schur bounds mismatch')
    return {'box':bj(box),'Arr':ij(Arr),'Nuu':ij(Nuu),'Nuv':ij(Nuv),'Nvv':ij(Nvv),'strict':True}


def build_trie(leaves):
    root={'children':{},'leaf':None}
    for leaf in leaves:
        node=root
        for tok in parse_path(leaf['path']):
            if node['leaf'] is not None:raise ValueError('leaf is prefix of another path')
            node=node['children'].setdefault(tok,{'children':{},'leaf':None})
        if node['leaf'] is not None or node['children']:raise ValueError('duplicate/prefix leaf path')
        node['leaf']=leaf
    return root


_PREP_CACHE={}

def _prepared(spec_bytes,spec,tree):
    # Deterministic memoization only: mutation tests often reuse the exact same
    # frozen proof specification and local certificate.  The one-shot checker
    # still rebuilds all proof data from source; repeated calls in one process
    # may reuse that identical reconstruction.
    sh=hashlib.sha256(spec_bytes).hexdigest()
    local_key=json.dumps(tree.get('local_certificate'),sort_keys=True,separators=(',',':'))
    key=(sh,hashlib.sha256(local_key.encode()).hexdigest())
    if key not in _PREP_CACHE:
        root=root_data(spec)
        local=validate_local(spec,tree,root)
        _PREP_CACHE[key]=(root,local)
    return _PREP_CACHE[key]


def structural_preflight(tree):
    """Cheap completeness/sign sanity before expensive polynomial replay.

    This never certifies a leaf: the full replay below still recomputes every
    Bernstein bound.  It only rejects malformed/mutated trees early.
    """
    trie=build_trie(tree['leaves'])
    def walk(node):
        if node['leaf'] is not None:
            leaf=node['leaf'];reason=leaf.get('reason')
            if reason=='SAFE_A_NONPOS':
                if F(leaf['A_bounds'][1])>=0:raise ValueError('stored SAFE_A sign sanity failed')
                if F(leaf['margin'])!=-F(leaf['A_bounds'][1]):raise ValueError('stored SAFE_A margin sanity failed')
            elif reason=='SAFE_B_NONNEG':
                if F(leaf['B_bounds'][0])<=0:raise ValueError('stored SAFE_B sign sanity failed')
                if F(leaf['margin'])!=F(leaf['B_bounds'][0]):raise ValueError('stored SAFE_B margin sanity failed')
            return
        ch=node['children']
        if len(ch)!=2:raise ValueError('incomplete split: coverage gap')
        axes={ax for ax,side in ch};sides={side for ax,side in ch}
        if len(axes)!=1 or sides!={'L','R'}:raise ValueError('corrupt split children')
        for child in ch.values():walk(child)
    walk(trie)
    return trie

def check(spec_path,tree_path):
    spec_path=Path(spec_path);tree_path=Path(tree_path)
    spec_bytes=spec_path.read_bytes()
    spec=json.loads(spec_bytes);tree=json.loads(tree_path.read_text())
    h=hashlib.sha256(spec_bytes).hexdigest()
    if h!=tree.get('proof_spec_sha256'):raise ValueError('proof-spec hash mismatch')
    trie=structural_preflight(tree)
    root,local=_prepared(spec_bytes,spec,tree);local_box=tuple((F(a),F(b)) for a,b in local['box'])
    counts={};minA=None;minB=None;visited=0
    def rec(node,box,A,B,path=''):
        nonlocal visited,minA,minB
        if node['leaf'] is not None:
            if node['children']:raise ValueError('leaf with children')
            leaf=node['leaf'];visited+=1
            if leaf['path']!=path or leaf['box']!=bj(box):raise ValueError('stored leaf box/path mismatch')
            a,b=bounds(A),bounds(B)
            if leaf.get('A_bounds')!=ij(a) or leaf.get('B_bounds')!=ij(b):raise ValueError('stored leaf polynomial bounds mismatch')
            reason=leaf['reason'];counts[reason]=counts.get(reason,0)+1
            if reason=='SAFE_A_NONPOS':
                if a.hi>=0:raise ValueError('SAFE_A leaf lacks strict A<0')
                if F(leaf['margin'])!=-a.hi:raise ValueError('SAFE_A margin mismatch')
                minA=-a.hi if minA is None else min(minA,-a.hi)
            elif reason=='SAFE_B_NONNEG':
                if b.lo<=0:raise ValueError('SAFE_B leaf lacks strict B>0')
                if F(leaf['margin'])!=b.lo:raise ValueError('SAFE_B margin mismatch')
                minB=b.lo if minB is None else min(minB,b.lo)
            elif reason=='LOCAL_EQUALITY_CERTIFICATE':
                if not all(box[i][0]>=local_box[i][0] and box[i][1]<=local_box[i][1] for i in range(3)):
                    raise ValueError('local-certificate leaf outside certified box')
            elif reason=='BOUNDARY_LEMMA':
                if not any(lo==hi==0 for lo,hi in box):raise ValueError('invalid boundary lemma leaf')
            else:
                raise ValueError(f'unresolved/unknown terminal reason {reason}')
            return
        ch=node['children']
        if len(ch)!=2:raise ValueError('incomplete split: coverage gap')
        axes={ax for ax,side in ch};sides={side for ax,side in ch}
        if len(axes)!=1 or sides!={'L','R'}:raise ValueError('corrupt split children')
        ax=next(iter(axes));AL,AR=split_coeff(A,root['dA'],ax);BL,BR=split_coeff(B,root['dB'],ax);boxL,boxR=split_box(box,ax)
        rec(ch[(ax,'L')],boxL,AL,BL,path+f'{ax}L');rec(ch[(ax,'R')],boxR,AR,BR,path+f'{ax}R')
    rec(trie,((F(0),F(1)),)*3,root['A'],root['B'],'')
    if visited!=tree['terminal_leaves'] or visited!=len(tree['leaves']):raise ValueError('leaf count mismatch')
    if counts!=tree['reason_counts']:raise ValueError('reason counts mismatch')
    if tree.get('unresolved_count')!=0 or not tree.get('global_cover_complete'):raise ValueError('tree is not declared complete')
    # Strict safe margins + strict local Hessian imply the only lambda2=sigma equality is the exact double root.
    if minA is None or minB is None or minA<=0 or minB<=0:raise ValueError('missing strict exterior margins')
    return {'status':'CHECK_PASS','leaves':visited,'reason_counts':counts,
            'minimum_SAFE_A_margin':str(minA),'minimum_SAFE_B_margin':str(minB),
            'local_certificate':local,'coverage':'full binary partition of [0,1]^3',
            'equality_scope':'unique target equality can occur only at the exact R7P-045 double root'}


def main(argv=None):
    argv=sys.argv[1:] if argv is None else argv
    if len(argv)!=2:
        print('usage: boundary_checker.py SPEC.json TREE.json');return 2
    try:r=check(argv[0],argv[1])
    except Exception as e:
        print('CHECK_FAIL');print(type(e).__name__+': '+str(e));return 1
    print('CHECK_PASS');print(json.dumps(r,indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
