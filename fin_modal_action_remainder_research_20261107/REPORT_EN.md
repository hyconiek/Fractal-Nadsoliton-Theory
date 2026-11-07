# Actual modal-action identities with an O(epsilon^3) remainder

Using CM's exact canonical chart and TP2's critical local-continuum bound,
the corrected modal actions now satisfy their actual three-channel work
identities with a remainder bounded uniformly by C_S epsilon^3 on
0<=t<=S epsilon^-2. This remainder is smaller than the required epsilon^2
scale as epsilon tends to zero.

The leading work retains the two odd-field profiles K f^2 and K g^2,
and the mixed even-field profile K f g, with their original factors of
12 in the action equations. The canonical momenta are kept throughout
the derivation. Explicit cubic corrections reconstruct the action
observables in physical variables, including CM's initial shape slip.

| Remainder contribution | Integrated upper bound |
|---|---|
| Cubic terms with at least two continuum factors | O(epsilon^3) |
| Mixed quartic terms, including the canonical transformation terms | O(epsilon^3) |
| Hamiltonian Taylor remainder of degree at least five | O(epsilon^3) |

The pure quartic normal form commutes with both actions. The proof also
gives L2-in-time bounds O(epsilon^2) on the canonical interaction-amplitude
derivatives and O(epsilon^3) on derivatives of their quadratic products.
These provide quantitative inputs for the next causal-response transfer.

The leading actual radiation work has size O(epsilon^2). Its replacement
by CD's dissipative action polynomial still needs an o(epsilon^2) error
estimate. MA15 states that obligation for both actions, retaining all
three channels. This is the next task before an actual damping law can
be admitted.

## Verification and scope

All 53 recorded checks passed using SymPy 1.14.0. They verify the canonical
flow sign, both explicit cubic action corrections, the three-channel work
coefficients, the full cubic density expansion and the vanishing brackets
with the quartic normal form. The mixed quartic audit finds, in each
action bracket, eight monomials with one continuum factor, four with two,
and two with three. Every coefficient retains a localized bound profile.
Four exact rational norm-history controls confirm the integrated powers,
and predecessor hashes are checked.

These are conditional analytical results under the recorded global and
local-continuum bounds. Independent proof review and formal verification
remain outstanding. The finite norm-history controls test exponent
bookkeeping rather than simulate the prepared nonlinear trajectory.
