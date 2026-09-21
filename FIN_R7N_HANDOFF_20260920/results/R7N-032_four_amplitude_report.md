# R7N-032 — Four-amplitude lane closeout

## Target P — tau0 = 67/250

Target P remains **open globally**. A rigorous bounded campaign produced a certified partial cover of the compact residual hull. After the 10,000-leaf pilot, tighter trace/e2/compression primitives and two directed `t` refinements, the final unresolved set occupies `0.3631049472406795` of the compact hull volume. Thus `0.6368950527593205` of that hull is certified by the current proof tree, in addition to the accepted unbounded tails.

The final residual contains 5,432 cells. Numerical values at cell centers are navigation only: the largest observed center value was `lambda2(M4) ~= 0.25392497687445353 < 0.268`, and no center exceeded `tau0`. This is not a global bound and does not eliminate the residual.

The two directed `t` passes reduced residual volume by about 23.23% and then 13.16% relative to their incoming residuals. A further identical split is disallowed by the anti-loop policy without a new analytic/enclosure primitive.

## Target S — sharp sigma

Target S remains **open**. It was not re-entered after Target P retained a 36.31% compact-hull residual. No tau0-only leaf is promoted to sigma, and the campaign makes no sharp-ceiling or equality claim beyond the accepted baseline sigma-safe regions.

## Gain consequence

The exact implication `sigma < 67/250` and `g <= 250/67` was checked. If Target P were later closed, the four-amplitude Cartesian Hessian would have at most one strictly negative eigenvalue for supplied `0 < g <= 250/67`, with the endpoint caveat from the campaign plan. Since Target P is not globally closed, this remains a conditional global consequence; it is not promoted here.

## Scope and nonconclusions

Nothing in this partial result transfers to the full X7 model. It is not a stationary-point census, a global minimizer theorem, a physical gain/source law, or a proof of the sharp sigma ceiling.

## Next mathematical atom

The highest-value new atom is a physical-coupled matrix enclosure preserving the common `(r,s,t,y)` dependence—especially the coupled `t`, `t^3`, `t^4`, and `t^(2+/-sqrt(3))` terms—rather than a third global `t` split.
