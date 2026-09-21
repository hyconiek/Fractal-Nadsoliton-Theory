# MP7-020 — normalized spectral margin

## Result

Let `tau0=67/250`. On the complete declared shared-field domain of the four-amplitude model,

`lambda2(M4) <= tau0 - m_global`

with the rigorous bound

`m_global >= 6500000000000000000000000000/121203210760018863485060548125326941`

so numerically `m_global >= 5.36289423295059e-08` and
`lambda2(M4) <= 0.267999946371058`.

This is a normalized covariance spectral margin. It is **not** the raw R7O3 PD minor margin `13/500000000`.

## Compact cover

The refined compact cover has 25,656 terminal cells: 13,231 inherited R7N-safe terminals plus 12,425 R7O3 repair leaves. Every terminal receives a quantitative margin. The certificate classes are converted as follows.

1. `TRACE`: for PSD `M`, `trace(M) >= lambda1+lambda2 >= 2 lambda2`, hence `lambda2 <= trace(M)/2` and `m >= tau0-trace_upper/2`.
2. `E2`: `e2(M) >= lambda1 lambda2 >= lambda2^2`, hence `lambda2 <= sqrt(e2_upper)` and `m >= tau0-sqrt(e2_upper)`.
3. Three-dimensional witness compression: `K=tau0 G-B^T M B`. If `K>0`, then `K>=m G` implies by min-max `lambda2(M)<=tau0-m`. For a 3x3 SPD `K`,
   `lambda_min(K) >= 4 det(K)/trace(K)^2`. Since `B^T M B>=0`, `trace(K)<=tau0 trace(G)`, and `lambda_max(G)<=trace(G)`. Therefore
   `m >= 4 det_lower(K)/(tau0^2 trace(G)^3)`.
   A positive Gershgorin lower bound is used when stronger.

For old R7N records whose exact interval endpoints were serialized as Python floats, the checker first replaces each saved lower endpoint `f` by `nextafter(f,-infinity)`, converted exactly to a rational. Because `f` was the correctly rounded float of the original exact rational endpoint, this predecessor float is a rigorous lower bound. R7O3 stores exact determinant and Gram endpoints directly.

The compact minimum occurs at R7O3 leaf 10878 / original parent 4884. There

- `det(K) >= 13/500000000`,
- `trace(G) = 30000017389/10000000000`,
- the determinant/trace inequality gives exactly the displayed `m_global`.

No compact terminal failed the margin extraction.

## Tails

The accepted FR1 J3/J4/J5 tail theorems bound `lambda2` by quantities below the sharper threshold `sigma_*`; FR42 gives `lambda2(M4)<=sigma_*` on the global large-J6 tail. The supplied enclosure has

`sigma_* <= 131141253604923881958598857/490351715494107500000000000`,

so on those tails

`tau0-lambda2 >= tau0-sigma_* >= 273006147496928041401143/490351715494107500000000000 ~= 0.000556755771155`.

This is much larger than the compact minimum. The compact hull retains the tail threshold boundaries, so the compact/tail union covers the full declared domain without a gap. Therefore the global margin is the compact leaf-10878 margin.

## Nonconclusions

This margin is tied to the supplied C4 feature normalization and the declared shared-field model. It is invariant only after carrying the corresponding quadratic metric as in MP7-019. It is not Target S, not a full-X7 covariance gap, and not a physical constant.
