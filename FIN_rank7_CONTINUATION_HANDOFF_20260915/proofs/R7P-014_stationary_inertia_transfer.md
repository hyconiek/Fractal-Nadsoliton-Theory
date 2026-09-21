# R7P-014 — exact stationary Hessian-index transfer

Let `p>0`, `g>0`, and let `Q` be any orthonormal basis of the probability
tangent space `1^perp`. At a common stationary point of the joint functional
from R7P-011, its Hessian in `(delta p, delta theta)` is

```
J = [ Q^T diag(1/p) Q    -Q^T X ]
    [ -X^T Q              I/g   ].
```

The first diagonal block is strictly positive. Its Schur complement is

`H7 = I/g - X^T [diag(p)-p p^T] X`,

the full seven-coordinate dual Hessian. Hence Sylvester inertia additivity gives

`n_-(J)=n_-(H7)`, `n_0(J)=n_0(H7)`.

The second diagonal block `I/g` is also strictly positive. Its Schur complement
is the primal tangent Hessian

`Hp = Q^T [diag(1/p)-g X X^T] Q`.

Therefore

**`n_-(Hp)=n_-(H7)` and `n_0(Hp)=n_0(H7)`**, while `Hp` has exactly four more
positive directions than `H7` because their dimensions are 11 and 7.

The equality is a theorem at a common interior stationary state. It does not
license transferring Hessian signatures through arbitrary nonlinear coordinate
restrictions at noncritical points, and it does not identify the four-amplitude
H4 with H7.
