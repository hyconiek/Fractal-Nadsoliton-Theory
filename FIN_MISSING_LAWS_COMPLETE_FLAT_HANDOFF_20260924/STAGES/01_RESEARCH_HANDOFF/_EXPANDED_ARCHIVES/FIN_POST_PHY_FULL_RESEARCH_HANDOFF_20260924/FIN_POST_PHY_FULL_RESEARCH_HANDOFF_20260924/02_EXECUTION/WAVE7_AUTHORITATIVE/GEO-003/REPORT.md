# GEO-003 — reconcile operator geometry with path geometry

Status: **CONDITIONAL_RESPONSE_GEOMETRY_PASS_WITH_ORDER_OF_LIMITS_GATE**.

The structured-grid control already shows the central mismatch: on `(Z_m)^4` RNG is the 8-axis graph and its scaled operator converges to the flat Laplacian, while graph shortest-path distance converges to the L1 norm, not Euclidean distance (e.g. `(h,h,h,h)` has graph distance `4h` and Euclidean distance `2h`).

The correct operator-associated continuum metric is instead the intrinsic metric of the limiting Dirichlet form. For `E(f)=int grad f^T K grad f`, the energy-density constraint `grad f^T K grad f <=1` gives the dual metric

`d_K(x,y)=sqrt((x-y)^T K^-1 (x-y))`

(up to torus geodesic wrapping). For K=I this is Euclidean. For the GEO-002 rectangular component `K_a=diag(1/a,a)`, it becomes `sqrt(a dx^2 + dy^2/a)`, matching the operator anisotropy rather than the raw graph L1 metric.

A response-distance estimator also has a strict order-of-limits gate. At a fixed finite graph, an off-diagonal heat kernel has leading behavior `p_t(i,j)=Theta(t^k)` where k is graph distance in jumps, so `sqrt(-4t log p_t(i,j)) -> 0` as `t->0`. It therefore cannot yield the continuum distance if the time limit is taken first. After mesh refinement at fixed macroscopic t, the heat response converges on the declared band-limited family; only then can the continuum small-time asymptotic recover `d_K`.

Thus operator geometry and path geometry are reconcilable, but only by changing the observable to a response/energy metric and paying the mesh-before-time limit plus the REF-006 conditional refinement law. This is not a derivation of spacetime.
