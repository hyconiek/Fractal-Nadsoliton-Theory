# MP7-026 — controlled fold remainder and local scaling theorem

Scientific state: **PROVED_INTERVAL_ASSISTED_LOCAL_FOLD_SCALING**.

## Inputs and scope

This result imports the immutable validated R7P-031 simple-fold certificate and the same accepted outward spectral intervals.  It does not replace the certified spectral tuple by floating midpoints.  The source certificate SHA-256 is

`5537174281996415a20ea3dffce23c33258c51396b3850aedf4168d9280e352e`.

Let `(s_f,g_f,v)` be the simple fold supplied by that certificate, with `||v||=1`, `F=grad Phi`, `F(s_f,g_f)=0`, `H_f v=0`, and

- `a=v^T F_g in [-0.212571640084014,-0.212571626020459]`,
- `b=D^3 Phi[v,v,v] in [0.118982527450499,0.118982667103224]`.

Use the local chart

` s = s_f + xi v + E y,   epsilon = g-g_f,   E=(e1,e2,e3). `

The leading 3-by-3 principal block of `H_f` is positive definite, so this is a valid transverse chart because the fourth component of `v` is bounded away from zero.

The theorem below is **local**.  It does not identify a global transition or a global escape barrier.

## 1. Uniform transverse elimination

On the interval box containing every point needed below,

- each Cartesian `s_i` stays within radius `0.002` of the numerical fold center,
- `0 <= epsilon <= 1e-7`,
- the accepted spectral intervals are retained.

Direct interval evaluation of the centered feature moments gives conservative tensor bounds

- `||T3|| <= 0.222833872594371`,
- `||T4|| <= 0.190171798000148`.

An interval LDL factorization of the transverse block `A=E^T H E`, followed by

`lambda_min(A) >= min(D)/||L^{-1}||_F^2`,

gives the uniform lower bound

`lambda_min(A) >= 0.00930411916943287`.

Therefore the three transverse equations `E^T F=0` have at most one solution `y` for fixed `(xi,epsilon)`.  Strong monotonicity plus the outward boundary estimate below gives existence as well.

Taylor expansion about the *exact* fold uses the identities `F_f=0` and `H_f v=0`.  For `x=|xi|`, the transverse residual at `y=0` is bounded by

`||T(xi,0,epsilon)|| <= ||E^T F_g|| epsilon + 1/2 ||E^T T3[v,v]|| x^2`

plus the explicit mixed `x epsilon/g_f^2`, `epsilon^2`, and third-order Taylor remainder terms.  At the final outer branch radius and `epsilon=1e-7` this implies

`||y|| <= 4.3070e-6`.

The total amplitude displacement is at most

`||delta s|| <= 6.1104e-4`,

well inside the `0.002` box used to obtain the derivative bounds.  Thus the local estimate is self-consistent.

## 2. Reduced stationarity and paid remainder

After transverse elimination define

`r(xi,epsilon) = v^T F(s_f+xi v+E y(xi,epsilon), g_f+epsilon)`.

The second-order fold part is

`r_0 = a epsilon + (b/2) xi^2`.

Using the tensor bounds above, the exact second-order mixed `g` derivatives, and the third-order Taylor integral remainder gives, throughout the final branch bracket,

`|r-r_0| <= 0.005965415221 epsilon`.

This is the missing uniform remainder payment.  It is not obtained by fitting a numerical exponent.

The formal coefficient

`sqrt(-2a/b)`

lies in

`[1.89027850228933, 1.89027967415187]`.

Choose a 1.5% outward bracket.  For every `0<epsilon<=1e-7`, the inner and outer test points satisfy strict opposite signs because

- inner leading sign margin: `0.006329320165 epsilon`,
- outer leading sign margin: `0.006424977822 epsilon`,
- full paid remainder: at most `0.005965415221 epsilon`.

Hence one root exists on each side with

`1.86192432475499 <= |xi|/sqrt(epsilon) <= 1.91863386926415`.

## 3. Exactly two local branches

Differentiate the transverse equations.  With

`q = v + E y_xi`,

one has

`E^T H q=0`

and, because the eliminated coordinates have zero gradient,

`r_xixi = T3[q,q,q]`.

The local variation of `T3`, together with the inverse transverse bound, gives

`|r_xixi-b| <= 0.010035734194`.

Since

`b_lower - 0.010035734194 >= 0.108946793257 > 0`,

`r` is strictly convex throughout the certified local interval.  The two sign-bracketed roots are therefore the **only two** reduced roots in that interval.  The negative-`xi` branch is the index-one side and the positive-`xi` branch is the local-minimum side.

## 4. Controlled local energy difference

On the transverse solution manifold,

`d Phi_reduced / d xi = r`.

Therefore the saddle-minus-minimum local energy difference is obtained by integrating `r` between the two roots.  Integrating the leading cubic normal form and the paid uniform remainder gives

`0.512686799589 <= [Phi_saddle-Phi_min]/epsilon^(3/2) <= 0.558650696110`.

The formal coefficient from MP7-025 was

`[0.535759433151, 0.535759800736]`,

which lies inside this controlled interval.

This is a **local branch energy difference**.  It is not promoted to a global escape barrier without a separating-path/global variational theorem.

## 5. Actual soft Hessian eigenvalue

The original fold has one zero H4 eigenvalue.  The positive leading 3-by-3 principal block and Cauchy interlacing give a certified lower bound on the second H4 eigenvalue at the fold,

`lambda_2(H_f) >= 0.0172259553883460`.

Across the branch tube,

`||H-H_f|| <= 1.36166781195e-4`.

In the orthogonal decomposition generated by the fold null vector and its complement, the small eigenvalue obeys the Schur formula.  Bounding the off-diagonal correction by

`||H-H_f||^2 / (gap_f - 2||H-H_f||)`

gives the controlled result

`0.214905205154 <= |lambda_soft|/sqrt(epsilon) <= 0.234915431887`.

Its sign is negative on the `xi<0` branch and positive on the `xi>0` branch.  The formal coefficient

`[0.224910245779,0.224910385210]`

is again enclosed.

## 6. Theorem statement

For every supplied strict spectral tuple covered by the accepted outward intervals and for the R7P-031 simple fold associated with that tuple, every

`0 < epsilon = g-g_f <= 1e-7`

has exactly two stationary branches in the declared local chart.  They satisfy the controlled bounds

- `|xi| = C_xi sqrt(epsilon)`, `C_xi in [1.861924324755,1.918633869264]`;
- `Phi_saddle-Phi_min = C_E epsilon^(3/2)`, `C_E in [0.512686799589,0.558650696110]`;
- `|lambda_soft| = C_H sqrt(epsilon)`, `C_H in [0.214905205154,0.234915431887]`.

Thus the three MP7-026 fold exponents are now controlled on an explicit punctured interval, with certified finite-error coefficients.

## 7. Nonconclusions

This theorem does not prove that this fold is the first fold globally, that the local saddle is the global escape saddle, that the local equal-energy crossing is the first global transition, or that any particular mobility/noise law supplies a physical time scale.  Those remain separate obligations.
