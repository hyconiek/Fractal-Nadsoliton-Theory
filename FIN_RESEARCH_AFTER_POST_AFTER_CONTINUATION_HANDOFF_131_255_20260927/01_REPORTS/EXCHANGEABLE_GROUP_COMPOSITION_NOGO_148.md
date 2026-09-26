# EXCHANGEABLE-GROUP-COMPOSITION-NOGO-148
## Arbitrarily tagging the finite-copy ensemble into groups does not create physical subsystems

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, commit `ad15a909...`.

Status:
exact large-N algebra for the already declared exchangeable-copy ensemble.

## 1. Split one population into two tagged groups

Let the fixed group fractions be

    alpha,
    1-alpha,

with empirical compositions

    p1,
    p2.

The total composition is

    pbar
      =
      alpha p1+(1-alpha)p2.

Because the existing finite-copy Hamiltonian depends only on the TOTAL
composition, the tagged large-deviation potential is

    V_alpha(p1,p2)
      =
      alpha D(p1||u0)
      +(1-alpha)D(p2||u0)
      -(g/2) pbar^T A7 pbar.

The interaction term contains no additional information about which copy
belongs to which tagged group.

## 2. Equal groups: relative mode has no A7 interaction

For alpha=1/2 write

    p1=pbar+delta,
    p2=pbar-delta.

Then

    (g/2) pbar^T A7 pbar

is exactly independent of delta.

So at fixed pbar the relative coordinate delta is controlled only by the
entropy terms.

To quadratic order:

    V(pbar+delta,pbar-delta)
      =
      V(pbar,pbar)
      +(1/2) delta^T diag(1/pbar) delta
      +O(delta^4).

There is no A7 stiffness/softening in the relative mode.

## 3. Mean-field dynamics gives exact bare decay of group differences

The same heat-bath target acts on every tagged group:

    q = softmax(g A7 pbar).

Hence

    p1_dot=q-p1,
    p2_dot=q-p2.

Subtracting,

    boxed:
    d/dt (p1-p2)
      =
      -(p1-p2).

For any number of tagged groups with fixed weights,

    r_a=p_a-pbar

satisfies

    boxed:
    r_a_dot=-r_a.

Thus arbitrary subgroup differences decay at the microscopic refresh rate.

They do not inherit:
- the slow three-sector rate;
- a new metastable scale;
- a spatial coupling;
- a new FIN phase structure.

## 4. Composition verdict

The existing exchangeable-copy ensemble supplies a thermodynamic population,
not a decomposition into interacting spatial/local cells.

Tagging copies creates bookkeeping subgroups only.

Therefore a nontrivial multi-unit FIN theory requires a new JOINT interaction
law that distinguishes subsystem incidence.

This result prevents a false shortcut from

    many exchangeable copies

to

    many physical FIN cells.
