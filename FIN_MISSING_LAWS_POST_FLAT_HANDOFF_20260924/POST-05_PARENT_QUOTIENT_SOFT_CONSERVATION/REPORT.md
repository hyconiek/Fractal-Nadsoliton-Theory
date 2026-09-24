# POST-05 — parent quotient and soft conservation

Status: **EFFECTIVE_EQUIVALENCE_ONLY**

The strict parent has

`lambda1=0.754121154207`, `lambda2=1.577049514428`,

while A7 retains only lambda3..lambda6.

A sparse nonnegative parent

`[0.53385079, 0.19757896, 0.05169472, 0.00175987, 0.0, 0.0]`

reproduces lambda3..lambda6 with error `1.471e-08`
while changing lambda1,lambda2 to
`0.449293, 1.338646`.

Its embedded-jump shell distribution differs from strict by TV
`0.134417`.  Thus equality after hard moment projection is only an
**effective-lane equivalence**, not a gauge equivalence of the whole theory.

For soft moment penalties, a constrained eigenmode returns as

`lambda_eff=lambda/(1+tau lambda)`.

Strict values:
- tau=1: [0.4299139500133143, 0.6119593378390671];
- tau=10: [0.08829205909402459, 0.09403714683796138];
- tau=100: [0.009869130700741077, 0.009936989993638562].

Controlled moment leakage is therefore a natural parent-sensitive held-out
observable.
