# R7P-091--095 phase continuation status

## R7P-091 quartic complement exhaustion
A bounded interval cover was executed on the full periodic three-torus, with
certified quartic root neighborhoods of radius 0.05 removed only when a cell was
fully contained in such a neighborhood.  Under a hard budget of 2000 processed
cells, 396 cells were excluded by a gradient component whose interval omitted
zero.  The pass retained 1272 explicit unresolved leaves. Therefore quartic
exhaustion is **not proved**. The leaf list is saved; no cell was discarded as
"numerical noise".

## R7P-092 direct full-root isolation
Every one of the 60 quartic roots was used only as a starting locator for a
direct full log-mgf solve. The resulting 60 full roots are distinct, have maximum
torus displacement `0.0318902462465` and zero Morse-index changes. Each full
root is then independently enclosed by an interval Krawczyk box (canonical
radius `1e-5`) with certified Hessian inertia. The smallest full Hessian absolute
eigenvalue is `5.2923037774e-5`. This proves a local 60-to-60 correspondence;
it does **not** exclude additional full roots elsewhere.

## R7P-093 extra full roots
The task terminates unresolved. R7P-086 did not supply a sharp uniform gradient
remainder and R7P-091 did not exhaust the quartic complement. Consequently the
local 60-to-60 correspondence cannot be promoted to a full-function census.
The smallest missing atom is a certified complement exclusion, either directly
for the full function or by a quartic gradient gap strictly larger than a
uniform gradient remainder.

## R7P-094 amplitude robustness
No upstream interval certificate for the angular coexistence amplitudes was
supplied or reconstructed. All phase results in R7P-081--093 therefore apply to
the explicit decimal fixture only. This is the terminal scoped result required
by the task; no transfer is made to radial coexistence amplitudes or to an
unspecified amplitude neighborhood.

## R7P-095 first pure-k6 angular instability
For the pure alternating direction at fixed Cartesian dual radius, the
constrained spherical Hessian is the ambient log-mgf Hessian minus the Lagrange
multiplier. The tangent eigenvalues are analytic. The k3-cosine sector crosses
when

`lambda3 (1+tanh x) = lambda6 tanh(x)/x`, `x=sqrt(lambda6/12) r`.

The left-minus-right function is strictly increasing for x>0 because
`x sech^2(x)-tanh(x)<0`. Strict spectral intervals certify the root in
`r in [0.41421132290, 0.41421132293]`. At the upper endpoint all other tangent
sectors remain strictly negative, with the nearest competitor the k5 pair at
about `-1.48e-3`. Thus k3 is the first angular instability. This is an angular
sphere theorem, not the radial localization transition.
