"""High precision outward source P/pi input enclosure; conditional on mpmath.iv arithmetic."""
import mpmath as mp
from fractions import Fraction
from model import *

def endpoint(z):
    s,m,e,b=z
    return Fraction(-m if s else m,1)*(2**e if e>=0 else Fraction(1,2**(-e)))
def iv_src():
    mp.iv.dps=50;I=mp.iv
    omg=I.mpf('0.18575');phi=I.mpf('0.1625');eta=I.mpf('1.8');gain=I.mpf('5.145228719489142')
    W=[I.mpf(0)]+[I.cos(omg*min(i,12-i)+phi)/(1+I.power(I.mpf(min(i,12-i)),eta)) for i in range(1,12)]
    lam={k:sum((W[d]*(1-I.cos(2*I.pi*k*d/12)) for d in range(1,12)),I.mpf(0)) for k in (3,4,5,6)}
    return gain,[sum((lam[k]/6*I.cos(2*I.pi*k*d/12) for k in (3,4,5)),I.mpf(0))+(lam[6]/12)*(-1)**d for d in range(12)]
def bound(v,target):
    dy=Fraction.from_float(float(target));lo,hi=map(endpoint,v._mpi_)
    return max(abs(lo-dy),abs(hi-dy))
def to_decimal_hi(f,scale=10**15):
    import math
    return f'{(f.numerator*scale+f.denominator-1)//f.denominator/scale:.15g}'
def run():
    m=build(A7());gain,row=iv_src(); I=mp.iv
    maxA=max(bound(row[(j-i)%12],m['A7'][i,j]) for i in range(12) for j in range(12))
    maxrate=Fraction(0)
    for a in range(12):
      for b in range(12):
        vals=[I.exp(gain/3*(row[(j-a)%12]+row[(j-b)%12])) for j in range(12)]
        su=sum(vals,I.mpf(0))
        for j in range(12):maxrate=max(maxrate,bound(vals[j]/su,m['R'][a,b,j]))
    # 1728 stationary probabilities from exact-decimal source A7/g, including normalization.
    weights=[]
    for a in range(12):
      for b in range(12):
        for c in range(12):
          ids=(a,b,c)
          weights.append(I.exp(gain/6*sum((row[(u-v)%12] for u in ids for v in ids),I.mpf(0))))
    Z=sum(weights,I.mpf(0))
    maxpi=max(bound(weights[k]/Z,m['pi'][k]) for k in range(N))
    cond=[weights[k]/sum(weights[:12],I.mpf(0)) for k in range(12)]
    p0=make_prep(m,'conditional_00')
    maxprep=max(bound(cond[k],p0[k]) for k in range(12))
    piYexact=[sum((weights[k] for k in range(y*12,(y+1)*12)),I.mpf(0))/Z for y in range(Y)]
    maxpiY=max(bound(piYexact[y],m['piY'][y]) for y in range(Y))
    result={'mathematical_input':'exact decimal omega=.18575 phi=.1625 eta=1.8 g=5.145228719489142, exact circle pi','method':'mpmath.iv at 50 decimal digits; A7 circulant Fourier, 144 categorical vectors, all 1728 normalized equilibrium weights, initial conditional prep','conditional_on':'correct outward mpmath.iv transcendental enclosures','bounds_from_exact_fraction_endpoints':{'A7_max':to_decimal_hi(maxA),'rate_max':to_decimal_hi(maxrate),'pi_max':to_decimal_hi(maxpi),'piY_max':to_decimal_hi(maxpiY),'conditional_prep_max':to_decimal_hi(maxprep)},'conservative_input_allowances':{'A7_abs':1e-12,'per_rate_abs':1e-12,'pi_abs':1e-12,'piY_abs':1e-12,'initial_conditional_abs':1e-12}, 'input_status':'INTERVAL_ENCLOSED_CONDITIONAL_ON_MPMATH_IV'}
    for k,v in (('A7',maxA),('rate',maxrate),('pi',maxpi),('piY',maxpiY),('conditional',maxprep)):
        assert v<Fraction(1,10**12),(k,float(v))
    write('RG-D02/INPUT_ENCLOSURES.json',result);print('INPUT ENCLOSURES',result['bounds_from_exact_fraction_endpoints']);return result
if __name__=='__main__':run()
