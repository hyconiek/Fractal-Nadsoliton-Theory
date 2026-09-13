# R7P-041--R7P-045 — Exact boundary-Ising reconstruction and eigenvalue-count criterion

## Scope

This note concerns only the compactified even-parity boundary model inherited from the positive four-amplitude chart of the conditional rank-seven active-gain dual.  It neither proves the global four-field curvature ceiling nor supplies a strict FIN gain source.

## R7P-041: exact four-state model

On the even labels `j=0,2,4,6,8,10`, set

- `A = cos(pi*j/2) in {+1,-1}`,
- `B = cos(2*pi*j/3) = 1/4 + 3Y/4`, `Y in {+1,-1}`,
- `C5 = cos(5*pi*j/6) = A(1+3Y)/4`.

The four `(A,Y)` states, ordered `(++,+-,-+,--)`, have label sets

`{0}`, `{4,8}`, `{6}`, `{2,10}`

and multiplicities `(1,2,1,2)`.  Therefore a raw field

`J3*A + J4*B + J5*C5`

has the equivalent two-spin representation

`H_A A + H_Y Y + K A Y`,

with

`H_A=J3+J5/4`, `H_Y=3 J4/4-(log 2)/2`, `K=3 J5/4`.

The `-(log 2)/2` term is forced by the unequal degeneracies and cannot be dropped.  With

`X=exp(2J3)`, `Yp=exp(3J4/2)`, `Z=exp(J5/2)`,

the unnormalised probabilities are exactly

`(X Yp Z^4, 2 X Z, Yp, 2 Z^3)`.

A uniform-degeneracy two-spin model is therefore a different model.

## R7P-042: exact positive interior and the true closure

For positive probabilities `p1,...,p4`, the physical conditions `X,Yp,Z >= 1` are equivalent to

`g1=p1 p4-p2 p3 >=0`,

`g2=4 p1 p3-p2 p4 >=0`,

`g3=p1 p2^2-p3 p4^2 >=0`.

Necessity follows from the exact factorizations

`g1*S^2 = 2 X Yp Z (Z^6-1)`,

`g2*S^2 = 4 X Z^4 (Yp^2-1)`,

`g3*S^3 = 4 Yp Z^6 (X^3-1)`.

Conversely, for positive `p`, define

`Z^6=p1 p4/(p2 p3)`,

`Yp^2=4 p1 p3/(p2 p4)`,

`X^3=p1 p2^2/(p3 p4^2)`.

The inequalities make the positive roots at least one, and substitution into the three independent log-ratios reconstructs the original probability vector up to its normalization.  Hence the three inequalities exactly characterize the positive interior.

The zero-probability closure is smaller than the weak semialgebraic relaxation.  The exponent vectors of the four weights are

`v1=(1,1,4)`, `v2=(1,0,1)`, `v3=(0,1,0)`, `v4=(0,0,3)`

in the nonnegative log-parameter cone.  For every nonzero recession direction, `v1` is maximal.  The only nontrivial ties are `v1=v2` on the pure-X ray and `v1=v3` on the pure-Yp ray.  Thus the only infinite-parameter supports are subsets of `{1,2}` or `{1,3}`, plus vertex `{1}`.  The full nontrivial edge limits are

- `{1,2}` with `p1/p2 >= 1/2`, equivalently `p1>=1/3`;
- `{1,3}` with `p1/p3 >= 1`, equivalently `p1>=1/2`;
- vertex `p1=1`.

For example `p2=1` satisfies all three weak polynomial inequalities but is not in the exponential-family closure.  Therefore any Bernstein cover over only `g1,g2,g3>=0` is an outer-relaxation proof unless these extra boundary points are separately discharged.

## R7P-043: covariance invariants

After dropping the constant part of the k4 coordinate (covariance is translation invariant), the four feature vectors are

`(+sqrt(l3/6), +(3/4)sqrt(l4/6), +sqrt(l5/6))`,

`(+sqrt(l3/6), -(3/4)sqrt(l4/6), -(1/2)sqrt(l5/6))`,

`(-sqrt(l3/6), +(3/4)sqrt(l4/6), -sqrt(l5/6))`,

`(-sqrt(l3/6), -(3/4)sqrt(l4/6), +(1/2)sqrt(l5/6))`.

For the covariance matrix `M`, write

`det(t I-M)=t^3-e1 t^2+e2 t-e3`.

The trace is the pair-distance identity

`e1=sum_{i<j} p_i p_j d_ij`,

with

`d12=d34=3(l4+l5)/8`,

`d13=2(l3+l5)/3`,

`d14=d23=(16l3+9l4+l5)/24`,

`d24=(4l3+l5)/6`.

The second invariant is

`e2=A123(p1p2p3+p1p3p4)+A124(p1p2p4+p2p3p4)`,

where

`A123=(l3l4+l3l5+l4l5)/4`,

`A124=(4l3l4+4l3l5+l4l5)/16`.

The determinant is exactly

`e3=(3/8) l3 l4 l5 p1 p2 p3 p4`.

These identities were independently regenerated symbolically from the feature matrix.

## R7P-044: exact threshold-count criterion

Let the eigenvalues of the PSD covariance be `lambda_1>=lambda_2>=lambda_3>=0` and let

`P(t)=det(tI-M)`.

Set

`Q(z)=P(sigma+z)=z^3+c2 z^2+c1 z+c0`,

where

`c2=P''(sigma)/2=3 sigma-e1`, `c1=P'(sigma)`, `c0=P(sigma)`.

Positive roots of `Q` are exactly covariance eigenvalues above `sigma`.  If `c2>0`, three positive roots are impossible.  If exactly two positive roots `r1,r2` exist, the third root `r3` obeys `r3<-(r1+r2)` because the sum of roots is `-c2<0`.  Hence

`c0>0`,

`c1=r1 r2+r3(r1+r2)<r1r2-(r1+r2)^2<0`.

Conversely, if `c0>0` and `c1<0` with `c2>0`, zero positive roots are impossible because three nonpositive real roots give `c1>=0`; one positive root gives `c0<=0`; and three positive roots contradict `c2>0`.  Therefore exactly two eigenvalues exceed `sigma`.

Thus, under `c2>0`,

`lambda_2 <= sigma  iff  P(sigma)<=0 or P'(sigma)>=0`,

and equivalently

`lambda_2 > sigma  iff  P(sigma)>0 and P'(sigma)<0`.

This includes equality cases because `P(sigma)=0` covers an eigenvalue exactly at threshold.

For the strict four-state feature geometry, `c2>0` holds globally.  By symmetry the minimum enclosing ball has center

`(0, sqrt(6) l5/(24 sqrt(l4)), 0)`

and exact squared radius

`R^2=(16 l3 l4+9 l4^2+10 l4 l5+l5^2)/(96 l4)`.

For any probability distribution on the four points,

`tr Cov <= R^2`.

Strict spectral intervals certify

`3 sigma_* - R^2 > 0`

with numerical value about `0.004758950459`.  Hence the shifted-polynomial criterion is valid on the entire four-state boundary model.

This corrects an earlier symbolic expression for `R^2`; the earlier reported numerical value `~0.797570782227429` was correct.

## R7P-045: exact double root

Define

`sigma_*=[2 l3(l4+l5)-l4 l5]/(24 l3)`

and

`t_*^2=[(2l3-l4)(2l3-l5)]/(4l3^2)`.

At

`p=((1+t)/6,(1+t)/3,(1-t)/6,(1-t)/3)`

with `t^2=t_*^2`, exact polynomial reduction gives

`P(sigma_*)=0`, `P'(sigma_*)=0`.

Therefore two covariance eigenvalues equal `sigma_*`.  The remaining eigenvalue is

`rho=l4 l5/(24 l3)`

and

`sigma_*-rho=(l3l4+l3l5-l4l5)/(12l3)>0`

by the accepted strict spectral intervals.  The same intervals certify `0<t_*<1`.

The physical-domain constraints are

`g1=0`, `g2=0`, `g3=t(t^2+3)/27>0`.

Equivalently `Yp=Z=1` and `X=(1+t)/(1-t)>1`, so this is a genuine positive-probability point on the `J4=J5=0` field boundary, not a zero-probability simplex point.

## Nonclaims

Nothing here proves `lambda_2<=sigma_*` over the full finite four-field interior, proves global rank-seven minimality, derives active gain from strict FIN, or identifies a physical clock, selector, Standard Model, gravity sector, or ToE closure.
