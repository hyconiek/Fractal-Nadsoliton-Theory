# FR8/FR9 — diagonal local staircase around the compactified equality point

These post-handoff certificates reuse the exact second-order interval-AD and
Schur-cone checker of R7P-068.  They do not enlarge the theorem by convexity;
they add two explicitly certified boxes to the union of safe local regions.

Coordinates are

- `x=r-r_*`, `r=exp(-2 J3)`,
- `u=1-exp(-3 J4/2)`,
- `v=1-exp(-J5/2)`,
- `e=1-q_even`.

## FR8: diagonal r-v arm

The checker certifies the full box

`|x| <= 1/6500`, `0<=u<=1/8192`, `0<=v<=1/4600`, `0<=e<=1/100000`.

Both the boundary and endpoint Schur-cone sufficient conditions pass with
strict interval margins.  Hence `lambda2(M4)<=sigma_*` throughout this box.

## FR9: diagonal r-u arm

The checker certifies the full box

`|x| <= 1/6500`, `0<=u<=1/2432`, `0<=v<=1/8192`, `0<=e<=1/100000`.

The exploratory checker continued to pass through `u=1/2406` and failed at
`u=1/2404` with the other declared exploratory radii fixed.  The theorem uses
the deliberately safer rational threshold `1/2432`; the pass/fail transition
is only a method diagnostic, not a physical boundary.

## Scope

FR8/FR9 are four-amplitude local certificates.  Their union with R7P-068,
FR6, FR7, FR3 and the global tails is legitimate; their convex hull is not
certified.  No statement transfers to the full seven-coordinate Hessian.
