# Coherent superposition preserves the cubic comparison budget

PS extends EF from one Gaussian packet to absolutely summable collections
of even exact-potential Gaussian packets with bounded relative widths,
arbitrary phases and a common mass. All quadratic and cubic cross terms
are retained. Constants are independent of the number of packets, provided
the sum M of absolute coefficients remains bounded.

The integrated squared cubic trace is bounded by C M^6 epsilon^2 on the
critical window. Changing coefficients by total absolute amount delta,
with both coefficient masses at most M, changes the response by a squared
trace bounded by C M^4 delta^2 epsilon^2. These are analytical comparison
bounds using the predecessor estimates, not formal proof certificates.

This does not cover arbitrary energy-small errors or the actual forced
nonlinear continuum. It does not bound the leading linear packet trace
by O(epsilon^2). The actual nonlinear cubic and modal-action gates remain
open. Next: build a causal packet interface for source histories and check
whether its coefficient norm is controlled by existing actual estimates.

## Mixed-term numerical diagnostic

| dt | Full mixture trace norm | Omitted cross-term norm | Diagonal-only norm |
|---|---|---|---|
| 0.02 | 1.368225659 | 0.810693808 | 0.586998934 |
| 0.01 | 1.368098283 | 0.810606672 | 0.586960614 |

These are L2-in-time trace norms, not squared integrals. The omitted cross
response has about 59.25% of the full response norm; summing diagonal
single-packet responses is therefore inaccurate in this example. The
norms do not add because cross and diagonal traces can interfere.
Halving dt changes the full norm by about 0.0093%. The maximum relative
energy/source-work defect decreases from 0.0005212 to 0.0001303.
Doubling both input amplitudes multiplies the computed cubic response
by eight (zero difference at recorded floating-point precision).
No asymptotic scaling claim is inferred from this one R and mass.
