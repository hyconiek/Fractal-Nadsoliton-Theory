# R7P-015 — compact domains for stationary roots and negative sublevels

Let `R=max_j ||X7_j||`. Stationarity gives `theta=g X7^T p`, and `X7^T p` is a
convex combination of the row feature vectors. Hence every stationary point
satisfies

`||theta|| <= g R`.

For the dual objective, `log mean exp(X theta) <= R ||theta||`, so

`Phi_g(theta) >= ||theta||^2/(2g)-R||theta||`.

Therefore every point with `Phi_g(theta)<=0` satisfies

`||theta|| <= 2 g R`.

More generally `Phi_g<=c` lies inside the positive root of
`r^2/(2g)-R r-c<=0`, namely

`r <= gR + sqrt(g^2 R^2 + 2 g c)`

whenever the stated sublevel is nonempty. The stationary ball and the sublevel
ball are different proof domains and must not be interchanged.
