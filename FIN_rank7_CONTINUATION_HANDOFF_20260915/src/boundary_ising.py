"""R7P-041--045: exact boundary-Ising algebra for the rank-seven FIN audit.

This module separates exact symbolic identities from strict-spectrum interval
sign decisions.  The four-state model is the even-parity compactified boundary
of the positive four-amplitude chart; it is not a global full-7D theorem.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp

from intervals import QI, sqrt_interval

ROOT=Path(__file__).resolve().parents[1]


def strict_intervals():
    saved=json.loads((ROOT/'inputs/fin_handoff_audit/results.json').read_text())
    vals=saved['exact']['laplacian_intervals']
    return [QI(F(x[0]),F(x[1])) for x in vals]


def four_state_table():
    """Return exact combinatorial boundary states before spectral rescaling."""
    # Six even labels j=0,2,4,6,8,10 aggregate by A=cos(pi j/2),
    # Y=(4 cos(2pi j/3)-1)/3.  C5=A(1+3Y)/4.
    rows=[]
    for j in range(0,12,2):
        A=1 if (j//2)%2==0 else -1
        # cos(2*pi*j/3) on even j is 1 or -1/2 exactly.
        B=1 if j%6==0 else F(-1,2)
        Y=F(4*B-1,3)
        C5=F(A,4)*(1+3*Y)
        rows.append((j,A,Y,B,C5))
    grouped={}
    for j,A,Y,B,C5 in rows:
        grouped.setdefault((A,int(Y)),[]).append(j)
    order=[(1,1),(1,-1),(-1,1),(-1,-1)]
    return {
        'labels':rows,
        'states':[{'A':A,'Y':Y,'multiplicity':len(grouped[(A,Y)]),'labels':grouped[(A,Y)]}
                  for A,Y in order],
        'multiplicities':[len(grouped[x]) for x in order],
    }


def exponential_weights():
    """Exact monomial weights in X=e^(2J3),Y=e^(3J4/2),Z=e^(J5/2)."""
    return ['X*Y*Z^4','2*X*Z','Y','2*Z^3']


def domain_symbolics():
    X,Y,Z=sp.symbols('X Y Z', positive=True)
    w=[X*Y*Z**4,2*X*Z,Y,2*Z**3]
    S=sum(w); p=[sp.cancel(x/S) for x in w]
    g1=sp.factor(p[0]*p[3]-p[1]*p[2])
    g2=sp.factor(4*p[0]*p[2]-p[1]*p[3])
    g3=sp.factor(p[0]*p[1]**2-p[2]*p[3]**2)
    return X,Y,Z,w,p,[g1,g2,g3]


def inverse_parameter_identities():
    p1,p2,p3,p4=sp.symbols('p1 p2 p3 p4', positive=True)
    return {
        'Z6':sp.factor(p1*p4/(p2*p3)),
        'Y2':sp.factor(4*p1*p3/(p2*p4)),
        'X3':sp.factor(p1*p2**2/(p3*p4**2)),
    }


def closure_boundary_description():
    # Tropical/exposed-support classification for nonnegative log parameters.
    # v_i are exponent vectors of the four unnormalised weights.
    return {
        'exponent_vectors':[(1,1,4),(1,0,1),(0,1,0),(0,0,3)],
        'possible_infinite_supports':[[1],[1,2],[1,3]],
        'edge_12_condition':'p3=p4=0, p1+p2=1, p1/p2>=1/2 (equiv. p1>=1/3)',
        'edge_13_condition':'p2=p4=0, p1+p3=1, p1/p3>=1 (equiv. p1>=1/2)',
        'vertex_condition':'p1=1',
        'relaxation_counterexample':'p2=1 satisfies weak polynomial inequalities but is not in the exponential-family closure',
    }


def symbolic_invariants():
    l3,l4,l5=sp.symbols('l3 l4 l5', positive=True)
    p1,p2,p3,p4=sp.symbols('p1 p2 p3 p4', nonnegative=True)
    a3=sp.sqrt(l3/6); a4=sp.sqrt(l4/6); a5=sp.sqrt(l5/6)
    V=sp.Matrix([
        [ a3, sp.Rational(3,4)*a4,  a5],
        [ a3,-sp.Rational(3,4)*a4, -sp.Rational(1,2)*a5],
        [-a3, sp.Rational(3,4)*a4, -a5],
        [-a3,-sp.Rational(3,4)*a4,  sp.Rational(1,2)*a5],
    ])
    ps=[p1,p2,p3,p4]
    distances={}
    for i in range(4):
        for j in range(i+1,4):
            d=V[i,:]-V[j,:]
            distances[f'd{i+1}{j+1}']=sp.factor((d*d.T)[0])
    e1=sp.factor(sum(ps[i-1]*ps[j-1]*val for key,val in distances.items()
                     for i,j in [(int(key[1]),int(key[2]))]))
    A123=sp.factor((l3*l4+l3*l5+l4*l5)/4)
    A124=sp.factor((4*l3*l4+4*l3*l5+l4*l5)/16)
    e2=sp.factor(A123*p1*p2*p3+A124*p1*p2*p4+A123*p1*p3*p4+A124*p2*p3*p4)
    e3=sp.factor(sp.Rational(3,8)*l3*l4*l5*p1*p2*p3*p4)
    return dict(l=(l3,l4,l5),p=(p1,p2,p3,p4),features=V,distances=distances,
                A123=A123,A124=A124,e1=e1,e2=e2,e3=e3)


def eigenvalue_count_criterion(e1,e2,e3,sigma):
    """Coefficients of P(sigma+z); positive z roots are eigenvalues > sigma."""
    P=sp.expand(sigma**3-e1*sigma**2+e2*sigma-e3)
    Pp=sp.expand(3*sigma**2-2*e1*sigma+e2)
    halfPpp=sp.expand(3*sigma-e1) # P''(sigma)/2
    return halfPpp,Pp,P


def enclosing_ball_symbolics():
    """Exact four-point trace bound proving P''(sigma)>0 for strict spectrum."""
    inv=symbolic_invariants(); l3,l4,l5=inv['l']; V=inv['features']
    # Symmetry reduces center to (0,c,0); equalize the two radius classes.
    c=sp.sqrt(6)*l5/(24*sp.sqrt(l4))
    center=sp.Matrix([0,c,0])
    radii=[]
    for i in range(4):
        d=V[i,:].T-center
        radii.append(sp.factor((d.T*d)[0]))
    assert sp.factor(radii[0]-radii[1])==0
    assert all(sp.factor(x-radii[0])==0 for x in radii)
    R2=sp.factor(radii[0])
    sigma=sp.factor((2*l3*(l4+l5)-l4*l5)/(24*l3))
    margin=sp.factor(3*sigma-R2)
    return c,R2,sigma,margin


def double_root_symbolics():
    inv=symbolic_invariants(); l3,l4,l5=inv['l']; p1,p2,p3,p4=inv['p']
    t=sp.symbols('t', positive=True)
    sigma=sp.factor((2*l3*(l4+l5)-l4*l5)/(24*l3))
    t2=sp.factor((2*l3-l4)*(2*l3-l5)/(4*l3**2))
    probs=[(1+t)/6,(1+t)/3,(1-t)/6,(1-t)/3]
    sub=dict(zip((p1,p2,p3,p4),probs))
    e1=sp.factor(inv['e1'].subs(sub));e2=sp.factor(inv['e2'].subs(sub));e3=sp.factor(inv['e3'].subs(sub))
    P=sp.factor(sigma**3-e1*sigma**2+e2*sigma-e3)
    Pp=sp.factor(3*sigma**2-2*e1*sigma+e2)
    # Reduce modulo t^2-t2 to retain exact dependence, not unrelated decimals.
    def reduce_t(expr):
        num=sp.together(expr).as_numer_denom()[0]
        den=sp.together(expr).as_numer_denom()[1]
        rem=sp.rem(sp.Poly(num,t),sp.Poly(t**2-t2,t)).as_expr()
        return sp.factor(rem),den
    P_rem,_=reduce_t(P); Pp_rem,_=reduce_t(Pp)
    rho=sp.factor((e1-2*sigma).subs(t**2,t2))
    gap=sp.factor((sigma-rho))
    g1=sp.factor((p1*p4-p2*p3).subs(sub))
    g2=sp.factor((4*p1*p3-p2*p4).subs(sub))
    g3=sp.factor((p1*p2**2-p3*p4**2).subs(sub))
    return dict(t=t,t2=t2,sigma=sigma,probs=probs,P_remainder=P_rem,Pp_remainder=Pp_rem,
                rho=rho,gap=gap,g=(g1,g2,g3))


def interval_signs():
    L=strict_intervals(); l3,l4,l5=L[3],L[4],L[5]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    t2=((2*l3-l4)*(2*l3-l5))/(4*l3*l3)
    t=sqrt_interval(t2,36)
    rho=l4*l5/(24*l3)
    gap=(l3*l4+l3*l5-l4*l5)/(12*l3)
    # Exact enclosing-ball radius without irrational center after simplification.
    R2=(16*l3*l4+9*l4*l4+10*l4*l5+l5*l5)/(96*l4)
    trace_margin=3*sigma-R2
    return {
        'sigma':(sigma.lo,sigma.hi),'t2':(t2.lo,t2.hi),'t':(t.lo,t.hi),
        'rho':(rho.lo,rho.hi),'sigma_minus_rho':(gap.lo,gap.hi),
        'trace_margin':(trace_margin.lo,trace_margin.hi),'R2':(R2.lo,R2.hi)
    }


def build_results():
    tab=four_state_table(); inv=symbolic_invariants(); clo=closure_boundary_description()
    c,R2,sigma,margin=enclosing_ball_symbolics(); dr=double_root_symbolics(); signs=interval_signs()
    def S(x):return str(sp.factor(x))
    out={
        'R7P-041':{
            'states':tab['states'],'multiplicities':tab['multiplicities'],
            'weights':exponential_weights(),
            'effective_fields':{'H_A':'J3+J5/4','H_Y':'3*J4/4-log(2)/2','K':'3*J5/4'},
            'inverse_parameter_powers':{k:S(v) for k,v in inverse_parameter_identities().items()},
        },
        'R7P-042':{
            'positive_interior_inequalities':['p1*p4>=p2*p3','4*p1*p3>=p2*p4','p1*p2^2>=p3*p4^2'],
            'exact_positive_interior_converse':True,
            'closure':clo,
            'weak_polynomial_boundary_is_outer_relaxation':True,
        },
        'R7P-043':{
            'feature_matrix':[[S(x) for x in inv['features'][i,:]] for i in range(4)],
            'pair_distance_coefficients':{k:S(v) for k,v in inv['distances'].items()},
            'A123':S(inv['A123']),'A124':S(inv['A124']),
            'e1':S(inv['e1']),'e2':S(inv['e2']),'e3':S(inv['e3']),
            'characteristic_polynomial':'t^3-e1*t^2+e2*t-e3',
        },
        'R7P-044':{
            'shifted_coefficients':{'z3':'1','z2':'3*sigma-e1','z1':'Pprime(sigma)','z0':'P(sigma)'},
            'criterion_under_z2_positive':'lambda2<=sigma iff P(sigma)<=0 OR Pprime(sigma)>=0',
            'violation_equivalence':'lambda2>sigma iff P(sigma)>0 AND Pprime(sigma)<0',
            'enclosing_center_y':S(c),'R2':S(R2),'trace_margin_symbolic':S(margin),
            'strict_trace_margin_interval':[str(signs['trace_margin'][0]),str(signs['trace_margin'][1])],
        },
        'R7P-045':{
            'sigma':S(dr['sigma']),'t_star_squared':S(dr['t2']),
            'probabilities':[S(x) for x in dr['probs']],
            'P_sigma_remainder_mod_t2':S(dr['P_remainder']),
            'Pprime_sigma_remainder_mod_t2':S(dr['Pp_remainder']),
            'remaining_eigenvalue':S(dr['rho']),'sigma_minus_remaining':S(dr['gap']),
            'domain_constraints':[S(x) for x in dr['g']],
            'strict_t_interval':[str(signs['t'][0]),str(signs['t'][1])],
            'strict_gap_interval':[str(signs['sigma_minus_rho'][0]),str(signs['sigma_minus_rho'][1])],
        }
    }
    return out


if __name__=='__main__':
    out=build_results()
    path=ROOT/'results/R7P-041_045_boundary_ising.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
