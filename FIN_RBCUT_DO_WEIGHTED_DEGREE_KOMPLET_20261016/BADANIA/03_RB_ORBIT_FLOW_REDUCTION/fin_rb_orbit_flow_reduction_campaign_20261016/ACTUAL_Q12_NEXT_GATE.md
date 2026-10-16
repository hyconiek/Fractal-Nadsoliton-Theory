# Actual q12 all-N gate after the orbit theorem

Scientific status: `EXACT_ORBIT_REDUCTION; ACTUAL_MT_QUOTIENT_INEQUALITY_OPEN`.

For the frozen actual q12 MT count graph, construct the 12-element parity-preserving label action `j -> j+2k` and `j -> 2k-j`. Each orbit is a set of genuine admissible count vectors, not a tagged independent copy or a physical spatial cell. For each orbit O calculate **symbolically with the actual Gibbs weights and fully reflected native generator**:

- `s_O=sum_(n in O)(N/2)W_n L0_D m(n)`;
- `W_O=sum_(n in O) W_n`;
- `u_OP=sum_(n->n' legal on H_N, n in O,n' in P) W_n q0(n,n')`.

The strict native-relative route is equivalent, not just sufficient at fixed N, to

`sup_(J subset orbits) [sum_(O in J)s_O]_+/sum_(O in J,P notin J)u_OP < infinity`.

Proving a bound uniform over large N yields a generator via explicit edge-capacity-proportional lifting; it would then allow the previous conditional rare-leak bracket argument. But **divergence of this ratio need not refute** the weaker original class.

To target the original vertex cap directly, solve or analytically bound the quotient problem

`Abar f=sbar, f>=0, sum_incident f<= C N Wbar`.

A negative verdict requires an exact orbit-level dual pair `(phi,beta)` meeting `phi_O-phi_P<=beta_O+beta_P`, `beta>=0`, and `sum phi s >CN sum beta W` at the relevant N/cap quantifiers; or with beta0, an actual negative lower ideal. This is a **complete** finite feasibility certificate, not a numerical sample or floating threshold test.

**Suggested next intrinsic step, bounded:** attempt a native monotone orbit-level transport estimate for the irregular-shell cuts that carry potential obstructions, retaining exact membership in D_N and E_N. If this remains unavailable, report the one missing weighted cut/LP inequality and stop; no new threshold fitting, detector or independent physics source.

The parity-preserving orbit count is generically reduced by only a factor near 12 in the number of vertices. The quotient is exact and structurally useful, but remains far too large for unrestricted enumeration at the true asymptotic MT N. There is no derived cross-N analytical inequality in this package.
