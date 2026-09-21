# R7O3-007 transformed interval-jet audit

Status: **validated analytic contract used by the clean fixed-witness replay**.

The proof chart is `(A,u,v,y)=(sqrt(r),1-s,1-t,y)`. Thus `r=A^2`, `s=1-u`, `t=1-v`, and `z=t^sqrt(3)`. On every compact cell, `A>0`, `s,t,y>0`; therefore every aggregate weight and the common denominator are positive.

The seven weights in the frozen aggregate-state order `[0,4,6,2,3,5,1]` are

```
1,
2 s t^3,
A^2 t^4,
2 A^2 s t,
2 A t^2 y,
2 A s t^2 y / z,
2 A s t^2 z y.
```

All weights are evaluated as jets in the **same four variables** before normalization. This retains the common physical dependence and avoids the older independent-weight relaxation.

For jets `(v,g,H)`, multiplication uses
`H(fg)=H(f)g+fH(g)+grad(f)grad(g)^T+grad(g)grad(f)^T`.
For an inverse,
`grad(1/f)=-grad(f)/f^2` and
`H(1/f)=2 grad(f)grad(f)^T/f^3-H(f)/f^2`.
The denominator is the sum of the seven positive weight jets and is inverted only after its interval is proved positive.

For `z=t^sqrt(3)`, monotonicity gives the value enclosure using rational exponent brackets. The derivative enclosure uses
`|dz/dv|=sqrt(3)t^(sqrt(3)-1)` with the conservative inclusion `sqrt(3)t <= |dz/dv| <= sqrt(3)` for `0<t<=1`. The second derivative is positive and
`sqrt(3)(sqrt(3)-1)t^(sqrt(3)-2)` is enclosed conservatively by multiplying `sqrt(3)(sqrt(3)-1)` by `[1,1/t_lo]`. These are deliberately wider than the exact derivative ranges but remain enclosing.

Each normalized moment entry is `N/D`. A midpoint Taylor enclosure is formed as

`f(m)+grad f(m)*(X-m)+1/2 sum_ij sup|H_ij(X)| radius_i radius_j`.

The proof witness `(B,c)` is fixed throughout differentiation. The 3x3 matrix
`K=(67/250) B^T B - E[(B^T F-c)(B^T F-c)^T]`
is then checked by exact/outward Sylvester minors or Gershgorin. No eigensolver is called by `certify_fixed`.

Validation evidence:
- `R7O3-010_arithmetic_validation.json` compares exact rational and outward-rounded backends on a frozen corpus, including the smallest paid-margin leaf.
- the producer replay is `12425/12425 PASS`;
- the independently copied clean-directory replay is `12425/12425 PASS` and reproduces the stored rational moment and PD bounds exactly;
- hostile mutations of center, Gram, rank evidence, split and threshold are rejected.

This audit establishes the implemented enclosure contract. It does **not** promote Target S at `sigma`, full X7, a sourced physical gain, or any selector/ToE claim.
