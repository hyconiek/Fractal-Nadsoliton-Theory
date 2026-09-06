"""Exact rational enclosures and integer sign certificates for finite spin maps.

This checker uses no floating transcendental function for any proof decision.
The scale is 10^12; all matrix/field integer operations stay inside int64.
"""
from fractions import Fraction as F
import itertools
import math

import numpy as np


SCALE=10**12


def add(A,B):return A[0]+B[0],A[1]+B[1]
def scale(c,A):return (c*A[0],c*A[1]) if c>=0 else (c*A[1],c*A[0])
def multiply(A,B):
    values=[x*y for x in A for y in B]
    return min(values),max(values)
def divide_positive(A,B):
    if not 0<B[0]<=B[1]:raise ValueError('Nonpositive denominator')
    return multiply(A,(1/B[1],1/B[0]))


def root_interval(n,p,decimal_places=18):
    """Enclose n^(1/p) using only an integer root search."""
    q=10**decimal_places;target=n*q**p
    lo=0;hi=1<<((target.bit_length()+p-1)//p)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**p<=target:lo=mid
        else:hi=mid
    assert lo**p<=target<(lo+1)**p
    return F(lo,q),F(lo+1,q)


def cos_small_rational(x):
    # All strict arguments are in (0, 1.3); alternating terms decrease.
    x=F(x)
    if not 0<=x<=F(13,10):raise ValueError('Unsupported cosine range')
    upper=sum((-1)**k*x**(2*k)/math.factorial(2*k) for k in range(21))
    nextterm=x**42/math.factorial(42)
    return upper-nextterm,upper


def log_two():
    N=25
    lower=2*sum(F(1,3)**(2*k+1)/(2*k+1) for k in range(N))
    tail=2*F(1,3)**(2*N+1)/((2*N+1)*(1-F(1,9)))
    return lower,lower+tail


def strict_weights():
    return [divide_positive(cos_small_rational(F(18575*d+16250,100000)),
                             add((F(1),F(1)),root_interval(d**9,5)))
            for d in range(1,7)]


def legacy_weights():
    r2,r3,r6=[root_interval(n,2) for n in [2,3,6]]
    diff=scale(F(1,4),add(r6,scale(-1,r2)))
    summ=scale(F(1,4),add(r6,r2))
    angles=[diff,(F(-1,2),F(-1,2)),scale(-1,summ),scale(F(-1,2),r3),
            scale(-1,diff),(F(1,2),F(1,2)),summ,scale(F(1,2),r3)]
    amplitude=scale(4,log_two())
    return [divide_positive(multiply(amplitude,angles[(d-1)%8]),
                             (F(100+d,100),F(100+d,100))) for d in range(1,12)]


def grid_interval(I):
    lo,hi=I
    a=lo*SCALE;b=hi*SCALE
    return a.numerator//a.denominator,-((-b.numerator)//b.denominator)


def integer_matrix(intervals,cyclic=True):
    low=np.zeros((12,12),dtype=np.int64);high=low.copy()
    for i in range(12):
        for j in range(12):
            if i!=j:
                d=min(abs(i-j),12-abs(i-j)) if cyclic else abs(i-j)
                low[i,j],high[i,j]=grid_interval(intervals[d-1])
    return low,high


def sign_certificate(low,high):
    states=np.array(list(itertools.product((-1,1),repeat=12)),dtype=np.int64)
    safe_bound=12*(int(np.max(abs(low)))+int(np.max(abs(high))))+int(np.max((high-low).sum(axis=1)))
    if safe_bound>=2**63:raise OverflowError('Integer field arithmetic would not be certified')
    mid=states@(low+high).T
    radius=(high-low).sum(axis=1)
    lower,upper=mid-radius,mid+radius
    undecided=(lower<=0)&(upper>=0)
    if np.any(undecided):raise AssertionError('Some spin fields have no certified sign')
    output=np.where(lower>0,1,-1)
    successor=((output+1)//2@(2**np.arange(11,-1,-1))).astype(int)
    fixed=int(np.count_nonzero(successor==np.arange(len(states))))
    two=int(np.count_nonzero((successor[successor]==np.arange(len(states)))&
                            (successor!=np.arange(len(states))))//2)
    seen=np.zeros(len(states),bool);lengths=[]
    for start in range(len(states)):
        if seen[start]:continue
        local={};path=[];u=start
        while not seen[u] and u not in local:
            local[u]=len(path);path.append(u);u=successor[u]
        if u in local:lengths.append(len(path)-local[u])
        seen[path]=True
    assert lengths.count(1)==fixed and lengths.count(2)==two
    minimum_margin=int(np.min(np.where(lower>0,lower,-upper)))
    return dict(configurations=len(states),certified_field_signs=int(states.size),
        scale=SCALE,minimum_abs_field_lower_bound=str(F(minimum_margin,2*SCALE)),
        integer_arithmetic_absolute_bound=safe_bound,
        maximum_interval_radius=str(F(int(max(radius)),2*SCALE)),
        fixed_points=fixed,two_cycles=two,other_cycle_lengths=[k for k in lengths if k not in [1,2]],
        weight_lower_integers=low.tolist(),weight_upper_integers=high.tolist())


def run():
    intervals=legacy_weights()
    out={}
    for name,weights,cycle in [('strict_cycle',strict_weights(),True),
                               ('legacy_cycle',intervals,True),('legacy_line',intervals,False)]:
        out[name]=sign_certificate(*integer_matrix(weights,cycle))
    return out


if __name__=='__main__':
    import json
    from pathlib import Path
    result=run()
    Path(__file__).with_name('sign_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    for k,v in result.items():
        print(k,'fixed',v['fixed_points'],'2-cycles',v['two_cycles'],
              'minimum certified field magnitude',v['minimum_abs_field_lower_bound'])
