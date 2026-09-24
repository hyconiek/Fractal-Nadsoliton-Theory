# NL-03 — complete transverse-mode reduction

## Result

`PASS_COMPLETE_PRIMITIVE_FOUR_INVARIANT_PHASE_QUOTIENT`.

For phases `(phi1,...,phi5)` the common carrier translation is

`phi_k -> phi_k + k alpha`, k=1,...,5.

The integer invariant lattice therefore is

`K = { n in Z^5 : n dot (1,2,3,4,5)=0 }`,

which has rank four.  A primitive basis is

```
r2 = 2 phi1 - phi2
r3 = phi1 + phi2 - phi3
r4 = phi1 + phi3 - phi4
r5 = phi1 + phi4 - phi5
```

The 4x5 row matrix R has rank four, `R k=0`, and the gcd of its 4x4 minors is
one, so these rows generate the full integer invariant lattice, not merely a
finite-index sublattice.

The previous proposed variables are projections of this complete sector:

`beta = phi3 - 2 phi4 + phi5 = r4-r5`,

`gamma = 4 phi3 - 3 phi4 = -r2-r3+3r4`.

Thus beta alone misses three invariant directions and `(beta,gamma)` still
misses two.  On the all-orders locked manifold `phi_k=k alpha`, all four r's
vanish.  Combined with the accepted PHA-001 transverse-gap result, the r-sector
is a gapped/internal-deformation sector on the declared positive-amplitude
window, not four additional Goldstone coordinates.

The frame `R^T R` has spectrum
`0, 1, 2.1225962543, 3.5235479603, 7.3538557854`; the only zero direction is
carrier translation.
