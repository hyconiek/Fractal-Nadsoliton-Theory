# R7P-011 — primal/dual correspondence

For `g>0`, let `u=1/12`, `X=X7` with `X^T 1=0`, and

`F(p,theta)=D(p||u)+||theta||^2/(2g)-theta^T X^T p`.

For fixed `p`, the unique minimizer is `theta=g X^T p`, giving

`D(p||u)-g||X^T p||^2/2 = V_g(p)`.

For fixed finite `theta`, strict convexity of relative entropy gives the unique
interior minimizer

`p_j=exp((X theta)_j)/sum_i exp((X theta)_i)`

and the reduced value

`Phi_g(theta)=||theta||^2/(2g)-log[(1/12)sum_j exp((X theta)_j)]`.

The joint problem is coercive in `theta` and compact in `p`; hence its minimum
exists. Partial minimization in either order therefore gives the same global
infimum. In particular, a global primal minimizer for `g>0` cannot be a boundary
simplex point: paired with its finite theta-minimizer it would be a joint global
minimum, while fixed-theta minimization has a unique strictly positive softmax.
At a common stationary point,

`theta=g X^T p`, `p=softmax(X theta)`.

This is joint minimization, not a minimax theorem, and it does not assert
`p-u in Range(A7)`. At `g=0` the theta representation with `1/g` is undefined;
the primal problem reduces to `D(p||u)` with unique minimizer `p=u`.
