# POST-02 — exact composition law

Status: **EXACT_THEOREM**

For a quadratic edge cost

`E_l(Delta)=1/2 Delta^T K(l) Delta`

series composition with exact elimination of the intermediate state gives

`K(l1+l2)^(-1)=K(l1)^(-1)+K(l2)^(-1)`.

With continuity, `R(l)=K(l)^(-1)` solves the additive Cauchy equation, hence

`K(l)=G/l`.

Parallel independent channels add precisions, so an additive cross-sectional
channel count gives the finite-volume form `Area/l`.

An 8-segment random 7D SPD replay closed with relative Frobenius error
`1.366e-15`.

For Gaussian self-similar increments `Var Delta_l ~ l^(2H)`, exact independent
concatenation requires `(l1+l2)^(2H)=l1^(2H)+l2^(2H)` for all positive lengths,
so H=1/2.  Persistent/telegraph dynamics can approach this law asymptotically,
but do not satisfy it exactly at all scales.
