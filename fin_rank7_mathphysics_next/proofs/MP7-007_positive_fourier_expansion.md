# MP7-007 — exact positive Fourier expansion of the aligned partition function

Let

`Z(phi;b)=(1/12) sum_{j=0}^{11} exp[ sum_{k=3,4,5} a_k cos(2*pi*k*j/12+phi_k) + b(-1)^j ]`

with finite `a_k>=0` and `b>=0`.

For `n>=0`, define

`I_n(a)=sum_{r>=0} (a/2)^(2r+n)/(r!(r+n)!)`, and `I_{-n}=I_n`.

Multiplying the absolutely convergent exponential series for
`exp[(a/2)e^(ix)]` and `exp[(a/2)e^(-ix)]` gives

`exp(a cos x)=sum_{m in Z} I_m(a)e^(imx)`.

Every coefficient is nonnegative; if `a>0`, every `I_m(a)` is strictly positive, while `I_m(0)=1_{m=0}`.  The coefficient sum is finite (indeed `sum_m I_m(a)=e^a`), so the product of the three series is absolutely convergent and the finite `j` average may be interchanged with the sums.

Also

`exp(b(-1)^j)=cosh(b)+sinh(b)(-1)^j`.

Put `K(m)=3m3+4m4+5m5`.  The exact character averages are

`(1/12)sum_j e^(2*pi*i*j*K/12)=1_{K=0 mod 12}`

and

`(1/12)sum_j (-1)^j e^(2*pi*i*j*K/12)=1_{K=6 mod 12}`.

Therefore

`Z(phi;b)=sum_{m in Z^3} c_m(a,b) e^(i m.phi)`

with

`c_m = prod_k I_{m_k}(a_k) [ cosh(b) 1_{K=0 mod12} + sinh(b) 1_{K=6 mod12} ]`.

Consequences:

- all `c_m>=0`;
- if all three `a_k>0` and `b>0`, the support is exactly
  `L6={m: K(m)=0 mod6}`;
- if all three `a_k>0` and `b=0`, the support is exactly
  `L12={m: K(m)=0 mod12}`;
- if some `a_k=0`, the exact support is obtained by imposing `m_k=0` on every inactive coordinate and the same congruence on the active coordinates.

The total absolute coefficient sum is bounded by
`e^(a3+a4+a5)(cosh b+sinh b)=e^(a3+a4+a5+b)`, so all regrouping above is paid.

This is an all-orders statement for `Z`.  It says nothing about coefficient signs of `log Z` or a finite cumulant truncation.
