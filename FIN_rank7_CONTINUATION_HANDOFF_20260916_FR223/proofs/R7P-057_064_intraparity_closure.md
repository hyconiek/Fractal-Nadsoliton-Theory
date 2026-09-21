# R7P-057--064 — intraparity closure from shared fields

## Scope

This note concerns only the four-amplitude parity decomposition with shared
nonnegative fields `J3,J4,J5,J6>=0`.  It proves a bound for the intraparity
covariance `W_par=q C_+ +(1-q) C_-`.  It does **not** yet bound the full
`M=W_par+b b^T` away from the extreme parity face.

## R7P-057: exact conditional laws and parity weight

Removing the constant `J6` factor inside each parity class gives

`Z_+ = 2 exp(J4) cosh(J3+J5) + 4 exp(-J4/2) cosh(J3-J5/2)`,

`Z_- = 2 exp(J4) + 4 exp(-J4/2) cosh(sqrt(3) J5/2)`.

Hence the conditional laws are independent of `J6`; the odd law is also
independent of `J3`.  With `A=J3`, `B=J5/2`, the proof in
`R7P-057_intraparity_q_bound.md` gives `Z_+>=Z_-`, equality iff
`J3=J5=0`.  Therefore at `J6=0`, `q0>=1/2`, and for `J6>=0`,
`q>=q0>=1/2`.

The same fields parameterize the two conditionals by

- even: `X=exp(2J3)`, `Y=exp(3J4/2)`, `Z=exp(J5/2)`, weights
  `(X Y Z^4, 2 X Z, Y, 2 Z^3)`;
- odd: the same `Y` and `S=exp(sqrt(3)J5/2)=Z^sqrt(3)`, weights
  `(Y,S,S^-1)`.

The relation `S=Z^sqrt(3)` is retained throughout; the two parity sectors are
not optimized independently.

## R7P-058: exact odd-sector domain and supremum

Write odd probabilities as

`p0=1-u`, `p+=(u+d)/2`, `p-=(u-d)/2`.

The exact field domain is

`0<=d<=u<=1`, `4(1-u)^2 >= u^2-d^2`,

with reconstruction

`S^2=(u+d)/(u-d)`, `Y^2=4(1-u)^2/(u^2-d^2)`.

The only nonconstant odd conditional coordinates are modes 4 and 5, giving

`C_- = (1/8) [[3 lambda4 u(1-u), -sqrt(3 lambda4 lambda5)(1-u)d],
              [same, lambda5(u-d^2)]]`.

Let `m=(3 lambda4+lambda5)/32`.  The matrix `m I-C_-` has positive
diagonal entries on the physical domain.  Its determinant is affine in
`y=d^2`; splitting at the derivative switch and at `u=2/3` reduces the proof
to three exact expressions.  The last interval has positive Bernstein
coefficients proportional to

`(lambda4+3lambda5)(9lambda4-5lambda5)`,
`(3lambda4+lambda5)(7lambda4-3lambda5)`, and positive squares.

The accepted strict spectral intervals certify all required signs.  Thus

`sup lambda1(C_-) = (3 lambda4+lambda5)/32`.

Equality is attained only in the **closure** at
`(p0,p+,p-)=(1/2,1/2,0)`, corresponding asymptotically to
`J4,J5->infinity` with `J5=sqrt(3)J4`; it is not attained at finite fields.

## R7P-059: dangerous set

The smaller eigenvalue is globally below `sigma_*`, using

`lambda_min(C_-) <= tr(C_-)/2 <= 3 lambda4/64 + lambda5/24 < sigma_*`.

Consequently `lambda1(C_-)>=sigma_*` is equivalent to
`det(sigma_* I-C_-)<=0`.  The determinant is linear in `d^2` and has the
correct strict denominator/monotonicity signs throughout the dangerous
`u` interval.  Therefore

`u in [u_-,u_+]`,

`u_±=(1 ± sqrt(1-32 sigma_*/(3lambda4+lambda5)))/2`,

and

`d^2 >= dmin^2(u)`

with

`dmin^2=((8sigma-lambda5 u)(8sigma-3lambda4 u(1-u))) /
        (lambda5(3lambda4(1-u)-8sigma))`,

intersected with the exact physical `(u,d)` domain.

## R7P-060: exact one-dimensional dominant-mass reduction

For the same fields, reconstruct

`R=(u+d)/(u-d)`, `Z=R^(1/(2sqrt(3)))`,
`Y=2(1-u)/sqrt(u^2-d^2)`.

The first even-state probability is

`p1 = 1/(1 + 2/(Y Z^3) + Z^-4 + 2/(Y Z))` when `X=1`.

For general `X>=1`, `p1` increases strictly with `X`.  For fixed `u`, both
`Y` and `Z` increase with `d`, hence `p1` increases with `d`.  A separate
strict interval cover proves `Y^2>4` throughout the dangerous interval, so
`p1` is indeed the largest even-sector mass there.  Thus the dangerous-set
minimum reduces exactly to one variable:

`X=1`, `d=dmin(u)`, `u in [u_-,u_+]`.

Numerical minimization of this exact reduced function (using spectral
midpoints only for the numerical locator) reproduces

`u≈0.5280356754`, `d≈0.4522658586`, `p_dom≈0.719071046876`.

The location/value are labelled numerical; the reduction is exact.

## R7P-061: conservative certified mass floor

Rather than certify the exact minimizer, two rational inequalities are proved
by one-dimensional Bernstein covers.  The certified dangerous interval is
strictly contained in `[2/5,3/5]`, and on that larger rational interval the
cleared denominator is positive.  The covers prove

`d/u >= 1703/2000 = 0.8515`,

`Y^2 >= 229/20 = 11.45`.

These imply `R>=(1+k)/(1-k)=3703/297`.  Since
`2sqrt(3)<693/200` and the exact rational comparison

`(3703/297)^200 > (2071/1000)^693`

holds, `Z>2071/1000`.  Also `Y>3383/1000` because
`(3383/1000)^2<229/20`.  Substitution into the monotone formula for `p1`
gives a fresh exact rational floor greater than

`p_dom > 711/1000`.

The imported decimal `0.7112098557` is not used as a proof input.

## R7P-062: covariance envelope

For the four even feature points, the second elementary covariance invariant is

`e2 = A(p1 p2 p3+p1 p3 p4)+B(p1 p2 p4+p2 p3 p4)`,

where

`A=(lambda3 lambda4+lambda3 lambda5+lambda4 lambda5)/4`,
`B=(4lambda3 lambda4+4lambda3 lambda5+lambda4 lambda5)/16`.

If one probability is at least `alpha>1/2`, AM-GM shows that the maximizing
split of the opposite pair is equal.  There are two symmetry classes of
dominant vertices; both reduce to a cubic in one variable.  Their derivatives
with respect to `alpha` are of the form
`2 C x(1-2alpha-2x)<0`, so the envelope for `p_i>=alpha` is maximized at
`p_i=alpha`.

At `alpha=711/1000`, exact-rational interval Bernstein covers for both branches
prove

`e2 < (511/2000)^2`.

Since covariance is PSD and `e2>=lambda2^2`, this yields

`lambda2(C_+) < 511/2000`.

## R7P-063: intraparity Weyl theorem

Use

`lambda2(W_par) <= q lambda2(C_+) +(1-q) lambda1(C_-)`.

Case A: if `lambda1(C_-)<=sigma_*`, R7P-055 gives
`lambda2(C_+)<=sigma_*`; therefore the right side is at most `sigma_*`.

Case B: if `lambda1(C_-)>=sigma_*`, R7P-059 identifies the dangerous set,
R7P-061 gives `p_dom>711/1000`, R7P-062 gives
`lambda2(C_+)<511/2000`, and R7P-058 gives
`lambda1(C_-)<=(3lambda4+lambda5)/32`.  Since `q>=1/2` and the latter bound
is larger than `511/2000`, the affine Weyl bound is maximized at `q=1/2`.
Strict spectral intervals certify

`2 sigma_* - (3lambda4+lambda5)/32 - 511/2000 > 0`

with lower margin about `0.00134546`.  Hence

**`lambda2(W_par)<=sigma_*` for all shared nonnegative fields
`J3,J4,J5,J6>=0`.**

## R7P-064: shortcut regressions

The theorem API explicitly checks the nonnegative-field/shared-field domain.
Two negative controls are retained:

1. If `J6<0`, the premise `q>=1/2` itself is false.  At `J3=J5=0` and
   `exp(J6)=1/2`, one has `Z_+=Z_-` and exactly `q=1/5`.
2. If the even conditional is replaced by an arbitrary simplex distribution
   rather than the exact physical boundary-Ising family, the already-certified
   R7P-047 relaxed-domain witness has `lambda2(C_+)>sigma_*`.

A bounded numerical search over independently chosen *physical* even/odd
exponential families did not find a curvature violation; that absence is not
promoted to a theorem.  The proved R7P-063 result does not require such a
relaxed theorem.
