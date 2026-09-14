# R7P-057 — exact shared-field parity distributions and q bound

For the four cosine amplitudes, the parity coordinate `k=6` is constant within each parity class. Therefore its factor cancels in the conditional distributions. On odd labels the `k=3` cosine is identically zero, so `C_-` depends only on `J4,J5`; `C_+` depends on `J3,J4,J5`.

Let `E` and `O` be the even/odd partition sums with the common `J6` factors removed. Exact label enumeration gives

`E = 2 exp(J4) cosh(J3+J5) + 4 exp(-J4/2) cosh(J3-J5/2)`,

`O = 2 exp(J4) + 4 exp(-J4/2) cosh(sqrt(3) J5/2)`.

At `J6=0`, `q0=E/(E+O)`, so it suffices to prove `E>=O`. Put `a=J3`, `b=J5`, `c=J4`. Since `cosh(a+b)-1>=0` and `exp(3c/2)>=1`,

`E-O >= 2 exp(-c/2) F(a,b)`,

where

`F(a,b)=cosh(a+b)-1+2 cosh(a-b/2)-2 cosh(sqrt(3)b/2)`.

Now `F_a` is strictly increasing in `a`, and

`F_a(0,b)=sinh(b)-2sinh(b/2)=2sinh(b/2)(cosh(b/2)-1)>=0`.

Hence `F(a,b)>=F(0,b)`. Finally

`F(0,b)=sum_{n>=1} b^(2n)/(2n)! * (4^n+2-2*3^n)/4^n`.

The coefficients vanish for `n=1,2`; for every `n>=3`, `4^n>2*3^n`, so they are strictly positive. Thus `F>=0`, with equality only at `a=b=0`. Consequently

`q0>=1/2`, equality iff `J3=J5=0` (arbitrary `J4`).

For general nonnegative `J6`,

`q_even = exp(J6)E / (exp(J6)E + exp(-J6)O) >= E/(E+O)`,

so `q_even>=q0>=1/2`. This extension intentionally requires `J6>=0`; a negative-`J6` regression is included in the tests.
