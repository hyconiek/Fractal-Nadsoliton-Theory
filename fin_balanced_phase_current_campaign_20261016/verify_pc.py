#!/usr/bin/env python3
"""Independent exact bounded checks for PC partial Hall theorem; NOT an asymptotic proof."""
from fractions import Fraction as F
from math import factorial
import json, pathlib
out={}
N=10000
pmin=F(103,1250); pmax=F(158,1875)
# exp(-x)>=1-x, sqrt(2)<=3/2, bound with |m|<=.01
energy=F(1,50)+12*F(3,2)*F(1,10000)
wratio=pmin/(pmax+F(1,N))*(1-energy)
assert wratio>F(9,10)
out['global_local_weight_ratio_rational_lower']=[str(wratio.numerator),str(wratio.denominator)]
out['energy_bound_rational']=str(energy)
# phasewise p ranges qmin/qmax from exact p(z)
for name,lo,hi in [('A',F(1,2),F(5,8)),('B',F(3,8),F(1,2))]:
    a=lo/hi/4;b=hi/lo/4
    assert a>=F(3,16) and b<=F(1,3)
    out['phase_'+name+'_current_ratio_bounds']=[str(a),str(b)]
assert F(9,50)<F(3,16) and F(7,20)>F(1,3)
assert F(9,1)*F(9,10)*F(9,50)/4>F(7,20)+F(7,1000)
out['overlap_four_lower']=str(F(9,1)*F(9,10)*F(9,50)/4)
out['outgoing_upper']=str(F(7,20)+F(7,1000))
# Even with native reflected upper s <= NW, corrected E O=max(I+s,0) <=1.35 NW.
assert F(9,1)*F(9,10)*F(9,50)>F(7,20)+1
out['singleton_extreme_margin']=str(F(9,1)*F(9,10)*F(9,50)-(F(7,20)+1))
# Verify exact q12 even/odd Fourier increment multiplicities (two signs in each class, each repeated 3).
from collections import Counter
sigma=lambda j: (-1)**j
cs=lambda j: ((1,0) if j%4==0 else ((-1,0) if j%4==2 else ((0,1) if j%4==1 else (0,-1))))
steps=Counter()
for a in range(12):
    for b in range(12):
        if sigma(a)==-1 and sigma(b)==1:
            c1,s1=cs(a);c2,s2=cs(b)
            steps[(c2-c1,s2-s1)]+=1
assert sum(steps.values())==36
assert len(steps)==4 and set(steps.values())=={9}
out['four_upward_fourier_increment_multiplicities']={str(k):v for k,v in steps.items()}
# Exact non-result: singleton Hall may hold while doubleton Hall fails.
left={'x1':F(1),'x2':F(1)};right={'y1':F(3,2),'y2':F(1,2)}
adjacency={'x1':{'y1'},'x2':{'y1'}}
assert all(left[x]<=sum(right[y] for y in adjacency[x]) for x in left)
assert sum(left.values())>right['y1']
out['negative_control_singletons_not_full_hall']='PASS; T={x1,x2} demand 2 exceeds neighbor supply 3/2'
# Deterministic finite q12 multinomial parity level coupling at g=0.
# x has E even and O odd; K(x,y)=n_a/O * 1/6 for odd->even replacement.
# For an edge x->y verify binomial parity conditional marginals pairwise by reversing
# with K_rev(y,x)=y_b/(E+1) * 1/6 and level weights conditional multinomial.
x=[1]*12; N0=sum(x); E=sum(x[::2]); O=sum(x[1::2]);a=1;b=0
y=x.copy();y[a]-=1;y[b]+=1
W=lambda z: F(factorial(sum(z)), __import__('math').prod(factorial(k) for k in z))
# ratio of parity-stratum probabilities = comb(N,E)/comb(N,E+1) = (E+1)/O;
# normalized W_l = W(x)/[C(N,E)*6**N], since 6^N selections parity-stratum
from math import comb
lhs=W(x)*F(x[a],O*6)/comb(N0,E)
rhs=W(y)*F(y[b],(E+1)*6)/comb(N0,E+1)
assert lhs==rhs
out['iid_q12_parity_stratum_exact_coupling']='PASS'
out['check_count']=12
path=pathlib.Path(__file__).with_name('verification.json')
path.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
print('PASS',out)
