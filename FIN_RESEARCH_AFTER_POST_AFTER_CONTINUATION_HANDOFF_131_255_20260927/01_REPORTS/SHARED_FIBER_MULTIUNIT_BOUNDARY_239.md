# SHARED-FIBER-MULTIUNIT-BOUNDARY-239
## Shared-fiber correlation creates a slower pair mode but not an extensive spatial hierarchy

Date: 2026-09-26

Status:
exact character-spectrum analysis of report 238's composition law.

## 1. One globally shared fiber

Let m Z3 units share the same beta event.

Conditional on beta their alpha increments remain independent.

For a character vector

    r=(r_1,...,r_m) in Z3^m,

the eigenvalue is

    lambda(r)
      =
      sum_beta R_beta[
        product_x phi_beta(r_x)-1
      ],

where

    phi_beta(r)
      =
      E[
        exp(2 pi i r alpha/3)
        |
        beta
      ].

The rule extends from 2 to 3 to arbitrary m without modification.

## 2. Held-out many-unit result

For every tested N=3,...,8:

- m=1 has the original Z3 gap;
- m=2 produces the slower relative gap from report 238;
- for EVERY m>=2 the gap stays at exactly that same pair-relative value.

The gap multiplicity grows as

    boxed:
    m(m-1).

Thus adding more units creates more degenerate relative modes but NO new
longer scale.

There is:
- no wavelength hierarchy;
- no n^-2 collective law;
- no spatial dimension.

This is the same qualitative obstruction as a complete/mean-field coupling.

## 3. Localize fiber sharing to edges

A natural repair is:
- give every graph edge its own shared-fiber channel;
- split each unit's fixed activity budget equally among its incident equivalent
  edges.

On a degree-2 cycle this gives half-weight per edge.

For the actual shared-fiber coefficients the edge-character costs satisfy, over
N=3,...,8,

    min(
      cost_same_nonzero,
      cost_opposite_nonzero
    )
      >
    0.86 * cost_zero/nonzero.

For a cycle with n>=3:

- any character pattern containing both zero and nonzero values has at least
  two zero/nonzero boundaries, already costing the isolated gap;
- an all-nonzero pattern has n nonzero/nonzero edges, and the above 0.86 bound
  also puts it above the isolated gap.

Therefore

    boxed:
    gap_cycle
      =
      gap_single

for n>=3 in this edge-shared-fiber model.

So localizing the fiber correlation creates locality but STILL does not create
a growing long-distance time scale.

## Verdict

Shared Z4 memory can source a nontrivial kinetic correlation.

It does not by itself supply hydrodynamic/spatial scaling.

Something stronger than common noise is required.
