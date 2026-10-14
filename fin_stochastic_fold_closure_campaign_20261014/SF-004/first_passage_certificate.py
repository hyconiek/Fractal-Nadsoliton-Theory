#!/usr/bin/env python3
"""Exact rational interval for canonical drift=-y^2 Brownian exit-side probability.
No Monte Carlo; mathematically bounded Taylor remainder; independent of numpy/SciPy.
"""
import math,json
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path

M=25;c=F(2,3)
Iplus=sum((c**k)/F(math.factorial(k)*(3*k+1)) for k in range(M+1))
Iminus=sum(((-c)**k)/F(math.factorial(k)*(3*k+1)) for k in range(M+1))
# exp(|x|) <= exp(2/3) < exp(1) <3; integrated Taylor error <= max point error
E=3*c**(M+1)/math.factorial(M+1)
lowN=Iminus-E; hiN=Iminus+E;lowP=Iplus-E;hiP=Iplus+E
assert lowN>0 and lowP>0
lo=lowN/(lowN+hiP);hi=hiN/(hiN+lowP)
assert 0<lo<hi<1
with localcontext() as ctx:
 ctx.prec=85
 fmt=lambda x,r: str((Decimal(x.numerator)/Decimal(x.denominator)).quantize(Decimal('0.000000000000000001'),rounding=r))
 lo_print=fmt(lo,ROUND_FLOOR); hi_print=fmt(hi,ROUND_CEILING)
 mid=str(Decimal((lo+hi).numerator)/Decimal((lo+hi).denominator)/2)
cert={
 'model':'dY=(-Y**2)dt+dW; start Y0=0; bounds -1,+1; kappa*lambda=0',
 'source':'original FIN stochastic fold after exact parameter-free canonicalization conditional on stopped process theorem',
 'M':M,'c_exponent_rational':'2/3',
 'Taylor_remainder_each_integral_strict_upper':str(E.numerator)+'/'+str(E.denominator),
 'I_negative_exact_rational_truncation':str(Iminus.numerator)+'/'+str(Iminus.denominator),
 'I_positive_exact_rational_truncation':str(Iplus.numerator)+'/'+str(Iplus.denominator),
 'P_exit_upper_public_outward_interval':[lo_print,hi_print],
 'P_exit_lower_public_outward_interval':None,
 'proof':'u_plus(0)= integral_-1^0 exp((2/3)y^3)dy / integral_-1^1 exp((2/3)y^3)dy. Taylor expansion with rational integrals and exp(2/3)<3 bounds both truncated integrals by E. Exact rational monotonicity bounds quotient; public decimal endpoints outward rounded.',
 'all_operation_scope':'exact rational arithmetic and final outward Decimal quantization; fold parameters canceled by normalization; no validation of original exact folded input bytes by this computation'
}
# for lower: (1-pupper) in outward sense
with localcontext() as ctx:
 ctx.prec=85
 cert['P_exit_lower_public_outward_interval']=[fmt(1-hi,ROUND_FLOOR),fmt(1-lo,ROUND_CEILING)]
p=Path(__file__).resolve().parent
(p/'EXIT_SIDE_CERTIFICATE.json').write_text(json.dumps(cert,indent=2)+'\n')
print('P upper',lo_print,hi_print,'P lower',cert['P_exit_lower_public_outward_interval'])
