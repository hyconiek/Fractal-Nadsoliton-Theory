"""R7P-046--051 boundary-Ising continuation.

Exact/validated layers:
- exact compactification of the physical exponential-family closure by
  r=1/X, s=1/Y, t=1/Z in [0,1]^3;
- exact zero-probability strata and rank-one boundary covariance lemma;
- interval-certified relaxed-domain negative controls;
- exact-rational interval Bernstein coefficients for the shifted characteristic
  numerators A=d^4(24 l3)^3 P(sigma), B=d^3(24 l3)^2 P'(sigma);
- bounded adaptive cover pilot.  Equality-neighborhood leaves are deliberately
  deferred to R7P-052 and are NOT counted as proof of the target inequality.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json
import sympy as sp

from intervals import QI
import boundary_ising as bi

ROOT=Path(__file__).resolve().parents[1]


def qi_json(I:QI):
    return [str(I.lo),str(I.hi)]


def eval_qi(expr, values):
    """Evaluate a rational polynomial expression with exact interval inputs."""
    expr=sp.sympify(expr)
    if expr.is_Rational:
        return QI(F(int(expr.p),int(expr.q)))
    if expr.is_Symbol:
        return values[expr]
    if expr.is_Add:
        out=QI(0)
        for x in expr.args: out=out+eval_qi(x,values)
        return out
    if expr.is_Mul:
        out=QI(1)
        for x in expr.args: out=out*eval_qi(x,values)
        return out
    if expr.is_Pow and expr.exp.is_Integer:
        return eval_qi(expr.base,values)**int(expr.exp)
    raise ValueError(f'unsupported interval expression: {expr}')


def compactified_symbolics():
    inv=bi.symbolic_invariants(); l3,l4,l5=inv['l']; pvars=inv['p']
    r,s,t=sp.symbols('r s t', nonnegative=True)
    d=1+2*s*t**3+r*t**4+2*r*s*t
    probs=[1/d,2*s*t**3/d,r*t**4/d,2*r*s*t/d]
    sub=dict(zip(pvars,probs))
    e1=sp.factor(inv['e1'].subs(sub)); e2=sp.factor(inv['e2'].subs(sub)); e3=sp.factor(inv['e3'].subs(sub))
    sigma=sp.factor((2*l3*(l4+l5)-l4*l5)/(24*l3))
    P=sp.factor(sigma**3-e1*sigma**2+e2*sigma-e3)
    Pp=sp.factor(3*sigma**2-2*e1*sigma+e2)
    # Multipliers are strictly positive on the spectral domain and d>=1.
    A=sp.factor(sp.cancel(P*d**4*(24*l3)**3))
    B=sp.factor(sp.cancel(Pp*d**3*(24*l3)**2))
    if sp.denom(A)!=1 or sp.denom(B)!=1:
        raise AssertionError('cleared characteristic numerators must be polynomials')
    return dict(l=(l3,l4,l5),vars=(r,s,t),d=d,probs=probs,sigma=sigma,A=sp.expand(A),B=sp.expand(B))


def exact_strata():
    """Complete zero-probability closure in the exact compact cube."""
    # p ∝ (1,2 s t^3,r t^4,2 r s t), r,s,t in [0,1].
    return {
        'compact_coordinates':{'r':'1/X','s':'1/Y','t':'1/Z'},
        'cube':'0<=r,s,t<=1',
        'unnormalized_weights':['1','2*s*t^3','r*t^4','2*r*s*t'],
        'normalizer':'1+2*s*t^3+r*t^4+2*r*s*t',
        'strata':[
            {'condition':'r=0, s>0, t>0','support':[1,2],
             'ratio':'p1/p2=1/(2*s*t^3)>=1/2','restriction':'p1>=1/3',
             'covariance_rank_bound':1,'lambda2':'0'},
            {'condition':'s=0, r>0, t>0','support':[1,3],
             'ratio':'p1/p3=1/(r*t^4)>=1','restriction':'p1>=1/2',
             'covariance_rank_bound':1,'lambda2':'0'},
            {'condition':'t=0 OR (r=s=0) with remaining limits','support':[1],
             'ratio':'vertex','restriction':'p1=1',
             'covariance_rank_bound':0,'lambda2':'0'},
        ],
        'excluded_supports':[[2],[3],[4],[1,4],[2,3],[2,4],[3,4],
                             [1,2,3],[1,2,4],[1,3,4],[2,3,4]],
        'boundary_curvature_conclusion':'Every zero-probability stratum has covariance rank <=1, hence lambda2=0<sigma_*.'
    }


def physical_constraints_exact(p):
    p1,p2,p3,p4=p
    return [p1*p4-p2*p3, 4*p1*p3-p2*p4, p1*p2*p2-p3*p4*p4]


def relaxed_negative_controls():
    """Interval-certify R7P-047 via the exact R7P-044 shifted-polynomial test."""
    witnesses=[
      ('drop_g1',[F(3,10),F(99,200),F(41,200),F(0)],0),
      ('drop_g2',[F(51,200),F(49,100),F(0),F(51,200)],1),
      ('drop_g3',[F(101,500),F(0),F(149,500),F(1,2)],2),
    ]
    inv=bi.symbolic_invariants(); l3,l4,l5=inv['l']; pvars=inv['p']; L=bi.strict_intervals()
    lv={l3:L[3],l4:L[4],l5:L[5]}
    sigma=(2*L[3]*(L[4]+L[5])-L[4]*L[5])/(24*L[3])
    out=[]
    for name,p,drop in witnesses:
        gs=physical_constraints_exact(p)
        vals=lv|dict(zip(pvars,map(QI,p)))
        e1=eval_qi(inv['e1'],vals); e2=eval_qi(inv['e2'],vals); e3=eval_qi(inv['e3'],vals)
        P=sigma**3-e1*sigma**2+e2*sigma-e3
        Pp=3*sigma**2-2*e1*sigma+e2
        certified=P.lo>0 and Pp.hi<0
        out.append({
            'name':name,'p':[str(x) for x in p],'dropped_constraint':drop+1,
            'g_exact':[str(x) for x in gs],
            'kept_constraints_nonnegative':all(gs[i]>=0 for i in range(3) if i!=drop),
            'dropped_constraint_negative':gs[drop]<0,
            'P_sigma_interval':qi_json(P),'Pprime_sigma_interval':qi_json(Pp),
            'lambda2_gt_sigma_certified':certified,
            'criterion':'R7P-044: P(sigma)>0 and Pprime(sigma)<0 iff lambda2>sigma, since 3sigma-trCov>0 globally.'
        })
    return out


def proof_spec():
    # Wide rational box intentionally quarantines the known double-root point;
    # R7P-052 must prove it, so these leaves are NOT safe/proved leaves.
    return {
      'id':'R7P-048-boundary-cover-spec-v1',
      'claim_id':'CLM-043',
      'domain':'exact compactified boundary-Ising closure [0,1]^3',
      'quantifiers':'Frozen proof specification for all compact coordinates and accepted strict spectral intervals.',
      'assumptions':['R7P-041--047 accepted scoped inputs','accepted strict lambda3,lambda4,lambda5 intervals'],
      'proof_type':'frozen machine-checkable proof specification; not itself a proof result',
      'inputs':['src/boundary_cover.py','inputs/fin_handoff_audit/results.json'],
      'conclusion':'Defines the exact domain, cleared characteristic criterion, equality quarantine and Bernstein cover rules used by R7P-049--055.',
      'version':1,
      'target':'lambda2(Cov)<=sigma_* on the exact four-state boundary-Ising exponential-family closure',
      'domain':{
        'coordinates':['r=1/X','s=1/Y','t=1/Z'],
        'box':[['0','1'],['0','1'],['0','1']],
        'probabilities':['1/d','2*s*t^3/d','r*t^4/d','2*r*s*t/d'],
        'd':'1+2*s*t^3+r*t^4+2*r*s*t',
        'exactness':'This closed cube maps exactly onto the physical exponential-family closure; no weak simplex relaxation is used.'
      },
      'criterion':{
        'A':'(24*lambda3)^3*d^4*P(sigma_*)',
        'B':'(24*lambda3)^2*d^3*Pprime(sigma_*)',
        'safe_disjunction':'A<=0 OR B>=0',
        'violation':'A>0 AND B<0',
        'justification':'positive clearing factors plus R7P-044 global z^2 coefficient positivity'
      },
      'boundary_lemmas':['R7P-046: zero-probability supports are {1,2},{1,3},{1}; covariance rank <=1 there.'],
      'equality_neighborhood':{
        'box':[['3/8','7/16'],['31/32','1'],['31/32','1']],
        'status':'DEFER_TO_R7P-052',
        'note':'Contains the strict-interval double-root location r_*=(1-t_*)/(1+t_*), s=t=1. Membership is not a proof of the target.'
      },
      'cover':{
        'representation':'tensor-product Bernstein boxes on [0,1]^3; exact dyadic bisection at 1/2 using de Casteljau',
        'spectral_uncertainty':'Each power coefficient is evaluated with exact rational interval arithmetic on accepted lambda3,lambda4,lambda5 intervals.',
        'reason_codes':['SAFE_A_NONPOS','SAFE_B_NONNEG','BOUNDARY_LEMMA','EQUALITY_NEIGHBORHOOD','UNRESOLVED_BUDGET','UNRESOLVED_DEPTH'],
        'pilot_leaf_budget':10000,'max_depth':18,
        'split_rule':'bisect axis maximizing physical width * max(A_degree,B_degree), ties r<s<t order by degree'
      },
      'dependency_policy':'Equality-neighborhood and unresolved leaves prevent a global PASS. Numerical midpoint/optimizer evidence is never a leaf acceptance rule.'
    }


def _power_coeff_intervals(expr, vars_, spectral_values):
    poly=sp.Poly(expr,*vars_)
    degrees=poly.degree_list()
    coeff={}
    for mon,c in poly.terms():
        coeff[tuple(mon)]=eval_qi(c,spectral_values)
    return tuple(degrees),coeff


def power_to_bernstein(power_coeff, degrees):
    """Tensor power -> Bernstein conversion on [0,1]^n, interval-safe."""
    n=len(degrees); out={}
    ranges=[range(d+1) for d in degrees]
    for idx in product(*ranges):
        acc=QI(0)
        kranges=[range(i+1) for i in idx]
        for k in product(*kranges):
            a=power_coeff.get(tuple(k))
            if a is None: continue
            rat=F(1)
            for j in range(n):
                rat*=F(comb(idx[j],k[j]),comb(degrees[j],k[j]))
            acc=acc+a*rat
        out[tuple(idx)]=acc
    return out


def bernstein_bounds(coeff):
    return QI(min(x.lo for x in coeff.values()),max(x.hi for x in coeff.values()))


def _split_line(line):
    rows=[list(line)]
    while len(rows[-1])>1:
        prev=rows[-1]
        rows.append([(prev[i]+prev[i+1])*F(1,2) for i in range(len(prev)-1)])
    n=len(line)-1
    left=[rows[k][0] for k in range(n+1)]
    right=[None]*(n+1)
    for k in range(n+1): right[n-k]=rows[k][-1]
    return left,right


def split_bernstein(coeff,degrees,axis):
    """Exact interval de Casteljau split at coordinate 1/2."""
    n=len(degrees); others=[j for j in range(n) if j!=axis]
    left={};right={}
    for oi in product(*[range(degrees[j]+1) for j in others]):
        base=[0]*n
        for j,v in zip(others,oi):base[j]=v
        line=[]
        for a in range(degrees[axis]+1):
            base[axis]=a;line.append(coeff[tuple(base)])
        L,R=_split_line(line)
        for a,x in enumerate(L):
            base[axis]=a;left[tuple(base)]=x
        for a,x in enumerate(R):
            base[axis]=a;right[tuple(base)]=x
    return left,right


def root_bernstein_data():
    sym=compactified_symbolics(); l3,l4,l5=sym['l']; L=bi.strict_intervals()
    sv={l3:L[3],l4:L[4],l5:L[5]}
    degA,powA=_power_coeff_intervals(sym['A'],sym['vars'],sv)
    degB,powB=_power_coeff_intervals(sym['B'],sym['vars'],sv)
    bernA=power_to_bernstein(powA,degA); bernB=power_to_bernstein(powB,degB)
    return {'degrees_A':degA,'degrees_B':degB,'A':bernA,'B':bernB,
            'power_terms_A':len(powA),'power_terms_B':len(powB)}


def _inside(box,outer):
    return all(box[i][0]>=outer[i][0] and box[i][1]<=outer[i][1] for i in range(3))


def classify_leaf(box,A,B,eq_box=None,local_certificate_box=None):
    """R7P-050 sound reason-code classifier. Equality is quarantine, not PASS."""
    unit=[(F(0),F(1))]*3
    if not _inside(box,unit):
        return {'reason':'OUTSIDE_COORDINATE_DOMAIN','proof_safe':False}
    a=bernstein_bounds(A); b=bernstein_bounds(B)
    if a.hi<=0:
        return {'reason':'SAFE_A_NONPOS','proof_safe':True,'margin':str(-a.hi),
                'A_bounds':qi_json(a),'B_bounds':qi_json(b)}
    if b.lo>=0:
        return {'reason':'SAFE_B_NONNEG','proof_safe':True,'margin':str(b.lo),
                'A_bounds':qi_json(a),'B_bounds':qi_json(b)}
    if local_certificate_box is not None and _inside(box,local_certificate_box):
        return {'reason':'LOCAL_EQUALITY_CERTIFICATE','proof_safe':True,
                'margin':'R7P-052 tangent-cone Hessian certificate',
                'A_bounds':qi_json(a),'B_bounds':qi_json(b)}
    # Degenerate exact boundary faces can be discharged analytically.
    if any(lo==hi==0 for lo,hi in box):
        return {'reason':'BOUNDARY_LEMMA','proof_safe':True,'margin':'sigma_*>0; rank(Cov)<=1',
                'A_bounds':qi_json(a),'B_bounds':qi_json(b)}
    if eq_box is not None and _inside(box,eq_box):
        return {'reason':'EQUALITY_NEIGHBORHOOD','proof_safe':False,
                'A_bounds':qi_json(a),'B_bounds':qi_json(b)}
    return {'reason':'UNRESOLVED','proof_safe':False,'A_bounds':qi_json(a),'B_bounds':qi_json(b)}


def split_box(box,axis):
    outL=list(box);outR=list(box);lo,hi=box[axis];mid=(lo+hi)/2
    outL[axis]=(lo,mid);outR[axis]=(mid,hi)
    return tuple(outL),tuple(outR)


def _box_json(box): return [[str(a),str(b)] for a,b in box]


def _choose_axis(box,degA,degB):
    deg=[max(degA[i],degB[i]) for i in range(3)]
    score=[(box[i][1]-box[i][0])*deg[i] for i in range(3)]
    return max(range(3),key=lambda i:(score[i],deg[i],i))


def run_cover_pilot(max_leaves=10000,max_depth=18):
    root=root_bernstein_data(); spec=proof_spec()
    eq=tuple((F(a),F(b)) for a,b in spec['equality_neighborhood']['box'])
    unit=((F(0),F(1)),)*3
    stack=[(unit,root['A'],root['B'],0,'')]
    leaves=[]; splits=0; max_stack=1
    while stack:
        box,A,B,depth,path=stack.pop()
        cls=classify_leaf(box,A,B,eq)
        if cls['reason']!='UNRESOLVED':
            leaves.append({'path':path,'depth':depth,'box':_box_json(box),**cls})
            continue
        # Splitting one current leaf into two increases eventual leaf count by one.
        if len(leaves)+len(stack)+1>=max_leaves:
            leaves.append({'path':path,'depth':depth,'box':_box_json(box),
                           **cls,'reason':'UNRESOLVED_BUDGET'})
            continue
        if depth>=max_depth:
            leaves.append({'path':path,'depth':depth,'box':_box_json(box),
                           **cls,'reason':'UNRESOLVED_DEPTH'})
            continue
        axis=_choose_axis(box,root['degrees_A'],root['degrees_B'])
        boxL,boxR=split_box(box,axis)
        AL,AR=split_bernstein(A,root['degrees_A'],axis)
        BL,BR=split_bernstein(B,root['degrees_B'],axis)
        # DFS left-first: push right then left.
        stack.append((boxR,AR,BR,depth+1,path+f'{axis}R'))
        stack.append((boxL,AL,BL,depth+1,path+f'{axis}L'))
        splits+=1; max_stack=max(max_stack,len(stack))
    counts={}
    for x in leaves: counts[x['reason']]=counts.get(x['reason'],0)+1
    unresolved=[x for x in leaves if not x.get('proof_safe',False)]
    remote=[x for x in unresolved if x['reason'] not in ('EQUALITY_NEIGHBORHOOD',)]
    return {
      'task':'R7P-051 bounded cover pilot','leaf_budget':max_leaves,'max_depth':max_depth,
      'splits':splits,'terminal_leaves':len(leaves),'max_stack':max_stack,
      'polynomials':{'A_degrees':list(root['degrees_A']),'B_degrees':list(root['degrees_B']),
                     'A_power_terms':root['power_terms_A'],'B_power_terms':root['power_terms_B'],
                     'root_A_bounds':qi_json(bernstein_bounds(root['A'])),
                     'root_B_bounds':qi_json(bernstein_bounds(root['B']))},
      'reason_counts':counts,'unresolved_count':len(unresolved),'remote_unresolved_count':len(remote),
      'global_pass':len(unresolved)==0,
      'note':'EQUALITY_NEIGHBORHOOD leaves are intentionally unresolved pending R7P-052.',
      'leaves':leaves,
    }


def build_046_050_results(include_root_coefficients=False):
    root=root_bernstein_data()
    out={
      'R7P-046':exact_strata(),
      'R7P-047':{'status':'INTERVAL_CERTIFIED_NEGATIVE_CONTROLS','witnesses':relaxed_negative_controls()},
      'R7P-048':proof_spec(),
      'R7P-049':{
        'A_degrees':list(root['degrees_A']),'B_degrees':list(root['degrees_B']),
        'A_power_terms':root['power_terms_A'],'B_power_terms':root['power_terms_B'],
        'root_A_bounds':qi_json(bernstein_bounds(root['A'])),'root_B_bounds':qi_json(bernstein_bounds(root['B'])),
        'conversion':'tensor power-to-Bernstein exact rational formula; de Casteljau dyadic subdivision',
        'spectral_coefficients':'QI exact rational intervals from accepted strict lambda3,lambda4,lambda5 intervals'
      },
      'R7P-050':{
        'safe_codes':['SAFE_A_NONPOS','SAFE_B_NONNEG','BOUNDARY_LEMMA'],
        'quarantine_codes':['EQUALITY_NEIGHBORHOOD'],
        'failure_codes':['UNRESOLVED','UNRESOLVED_BUDGET','UNRESOLVED_DEPTH','OUTSIDE_COORDINATE_DOMAIN'],
        'midpoint_or_optimizer_acceptance':False
      }
    }
    if include_root_coefficients:
        def dump(c): return {','.join(map(str,k)):qi_json(v) for k,v in c.items()}
        out['R7P-049']['root_A_coefficients']=dump(root['A']);out['R7P-049']['root_B_coefficients']=dump(root['B'])
    return out


if __name__=='__main__':
    print(json.dumps(build_046_050_results(),indent=2))

# --- R7P-052 local equality-neighborhood certificate -----------------------

def _derivative_bernstein(coeff,degrees,axis):
    """Bernstein coefficients for derivative in the LOCAL normalized coordinate."""
    n=degrees[axis]
    if n<=0: raise ValueError('zero degree derivative')
    newdeg=list(degrees);newdeg[axis]-=1;newdeg=tuple(newdeg)
    out={}
    for idx in product(*[range(d+1) for d in newdeg]):
        j=list(idx);j1=list(idx);j1[axis]+=1
        out[idx]=(coeff[tuple(j1)]-coeff[tuple(j)])*n
    return out,newdeg


def _scale_coeff(coeff,scale):
    return {k:v*scale for k,v in coeff.items()}


def _dyadic_cell_path(lo:F,hi:F):
    """Return L/R path for one exact dyadic cell [k/2^n,(k+1)/2^n]."""
    lo,hi=F(lo),F(hi);w=hi-lo
    if w<=0:return '' if w==1 else None
    den=w.denominator
    if w.numerator!=1 or den&(den-1): raise ValueError('box is not one dyadic cell')
    n=den.bit_length()-1
    k=lo*den
    if k.denominator!=1 or hi*den!=k+1: raise ValueError('misaligned dyadic cell')
    return ''.join('R' if (int(k)>>(n-1-i))&1 else 'L' for i in range(n))


def restrict_to_dyadic_box(coeff,degrees,box):
    out=coeff
    for axis,(lo,hi) in enumerate(box):
        path=_dyadic_cell_path(lo,hi)
        if path is None: raise ValueError('bad dyadic cell')
        for side in path:
            L,R=split_bernstein(out,degrees,axis); out=L if side=='L' else R
    return out


def hessian_component_bounds_on_box(Aroot,degrees,box):
    """Rigorous physical-coordinate Hessian component bounds from local Bernstein data."""
    local=restrict_to_dyadic_box(Aroot,degrees,box)
    widths=[b-a for a,b in box]
    def deriv(axes):
        c=local;d=degrees
        scale=F(1)
        for ax in axes:
            c,d=_derivative_bernstein(c,d,ax);scale/=widths[ax]
        c=_scale_coeff(c,scale)
        return bernstein_bounds(c)
    return {
      'Arr':deriv((0,0)),'Ass':deriv((1,1)),'Att':deriv((2,2)),
      'Ars':deriv((0,1)),'Art':deriv((0,2)),'Ast':deriv((1,2)),
    }


def square_interval(I:QI):
    if I.lo<=0<=I.hi:return QI(0,max(I.lo*I.lo,I.hi*I.hi))
    xs=[I.lo*I.lo,I.hi*I.hi]
    return QI(min(xs),max(xs))


def local_equality_certificate():
    """R7P-052: certify A<=0 on a dyadic box containing the double root.

    Coordinates for the tangent cone are x=r-r*, u=1-s>=0, v=1-t>=0.
    The transformed Hessian has hxx=Arr, hxu=-Ars, hxv=-Art,
    huu=Ass, huv=Ast, hvv=Att.  With hxx<0 and all three 2x2 Schur
    numerators positive, the Schur complement entries are negative.  Hence its
    quadratic form is <=0 for every x in R and u,v>=0.  Integral Taylor from
    the exact zero-gradient double root gives A<=0 throughout the box.
    """
    # One dyadic cell in every coordinate; contains r_* ~=0.40205464 and s=t=1.
    box=((F(102,256),F(103,256)),(F(255,256),F(1)),(F(255,256),F(1)))
    root=root_bernstein_data(); H=hessian_component_bounds_on_box(root['A'],root['degrees_A'],box)
    Arr,Ass,Att,Ars,Art,Ast=[H[k] for k in ('Arr','Ass','Att','Ars','Art','Ast')]
    Nuu=Arr*Ass-square_interval(Ars)
    Nvv=Arr*Att-square_interval(Art)
    Nuv=Arr*Ast-(Ars*Art)
    # Strict interval decisions.
    ok=Arr.hi<0 and Nuu.lo>0 and Nuv.lo>0 and Nvv.lo>0

    # Exact equality and zero-gradient identities modulo tau^2=tau_*^2.
    sym=compactified_symbolics();r,s,t=sym['vars'];l3,l4,l5=sym['l'];A=sym['A']
    tau=sp.symbols('tau',positive=True)
    tau2=sp.factor((2*l3-l4)*(2*l3-l5)/(4*l3**2)); rstar=(1-tau)/(1+tau)
    mod=sp.Poly(tau**2-tau2,tau)
    def rem(expr):
        q=sp.cancel(expr.subs({r:rstar,s:1,t:1}));num=sp.together(q).as_numer_denom()[0]
        return sp.factor(sp.rem(sp.Poly(num,tau),mod).as_expr())
    remainders={'A':str(rem(A)),'Ar':str(rem(sp.diff(A,r))),
                'As':str(rem(sp.diff(A,s))),'At':str(rem(sp.diff(A,t)))}
    # Verify the strict t_* interval implies r_* lies in the rational r cell.
    signs=bi.interval_signs(); T=QI(signs['t'][0],signs['t'][1]); Rstar=(1-T)/(1+T)
    contains=box[0][0]<=Rstar.lo and Rstar.hi<=box[0][1]
    return {
      'id':'R7P-052-local-equality',
      'claim_id':'CLM-042',
      'domain':'dyadic physical tangent-cone neighborhood of the R7P-045 double root',
      'quantifiers':'For all (r,s,t) in the stated dyadic box and all strict spectral tuples in the accepted rational intervals.',
      'assumptions':['R7P-045 exact double-root identities','R7P-048 cleared polynomial A and accepted spectral intervals'],
      'proof_type':'exact zero/gradient identities plus rational interval Bernstein Hessian-cone certificate',
      'inputs':['certificates/R7P-048_boundary_proof_spec.json','src/boundary_cover.py'],
      'status':'INTERVAL_CERTIFIED' if ok and contains and all(v=='0' for v in remainders.values()) else 'FAILED',
      'box':_box_json(box),'r_star_interval':qi_json(Rstar),'r_star_contained':contains,
      'equality_remainders_mod_tau2':remainders,
      'hessian_bounds':{k:qi_json(v) for k,v in H.items()},
      'schur_numerator_intervals':{'Nuu':qi_json(Nuu),'Nuv':qi_json(Nuv),'Nvv':qi_json(Nvv)},
      'conditions':{'Arr_negative':Arr.hi<0,'Nuu_positive':Nuu.lo>0,
                    'Nuv_positive':Nuv.lo>0,'Nvv_positive':Nvv.lo>0},
      'conclusion':'A<=0 on the entire local box; hence P(sigma_*)<=0 and lambda2<=sigma_* there.',
      'proof_route':'exact zero value+gradient at the double root; tangent-cone Hessian Schur negativity; integral Taylor along the segment inside the convex box.'
    }


def run_refined_cover(max_leaves=10000,max_depth=32):
    """R7P-053 refinement: use R7P-052 local theorem as a safe leaf rule."""
    cert=local_equality_certificate()
    if cert['status']!='INTERVAL_CERTIFIED': raise RuntimeError('local certificate unavailable')
    root=root_bernstein_data(); unit=((F(0),F(1)),)*3
    local_box=tuple((F(a),F(b)) for a,b in cert['box'])
    stack=[(unit,root['A'],root['B'],0,'')]; leaves=[];splits=0;max_stack=1
    while stack:
        box,A,B,depth,path=stack.pop()
        cls=classify_leaf(box,A,B,None,local_box)
        if cls['reason']!='UNRESOLVED':
            leaves.append({'path':path,'depth':depth,'box':_box_json(box),**cls});continue
        if len(leaves)+len(stack)+1>=max_leaves:
            leaves.append({'path':path,'depth':depth,'box':_box_json(box),**cls,'reason':'UNRESOLVED_BUDGET'});continue
        if depth>=max_depth:
            leaves.append({'path':path,'depth':depth,'box':_box_json(box),**cls,'reason':'UNRESOLVED_DEPTH'});continue
        axis=_choose_axis(box,root['degrees_A'],root['degrees_B'])
        boxL,boxR=split_box(box,axis);AL,AR=split_bernstein(A,root['degrees_A'],axis);BL,BR=split_bernstein(B,root['degrees_B'],axis)
        stack.append((boxR,AR,BR,depth+1,path+f'{axis}R'));stack.append((boxL,AL,BL,depth+1,path+f'{axis}L'))
        splits+=1;max_stack=max(max_stack,len(stack))
    counts={}
    for x in leaves:counts[x['reason']]=counts.get(x['reason'],0)+1
    unresolved=[x for x in leaves if not x.get('proof_safe',False)]
    return {'task':'R7P-053 refined boundary cover','leaf_budget':max_leaves,'max_depth':max_depth,
            'splits':splits,'terminal_leaves':len(leaves),'max_stack':max_stack,
            'local_certificate_box':cert['box'],'reason_counts':counts,
            'unresolved_count':len(unresolved),'global_cover_complete':len(unresolved)==0,
            'leaves':leaves}

def frozen_proof_spec():
    """Machine-checkable R7P-048 specification with exact polynomial source."""
    spec=proof_spec(); sym=compactified_symbolics();L=bi.strict_intervals(); signs=bi.interval_signs()
    spec['frozen_polynomials']={
      'variables':['r','s','t'],'spectral_variables':['l3','l4','l5'],
      'spectral_intervals':{'l3':qi_json(L[3]),'l4':qi_json(L[4]),'l5':qi_json(L[5])},
      'A':str(sym['A']),'B':str(sym['B']),
      'A_degrees':list(sp.Poly(sym['A'],*sym['vars']).degree_list()),
      'B_degrees':list(sp.Poly(sym['B'],*sym['vars']).degree_list()),
      't_star_interval':[str(signs['t'][0]),str(signs['t'][1])],
      'tau2_relation':'tau^2=((2*l3-l4)*(2*l3-l5))/(4*l3^2)',
      'r_star_relation':'r_star=(1-tau)/(1+tau)'
    }
    return spec
