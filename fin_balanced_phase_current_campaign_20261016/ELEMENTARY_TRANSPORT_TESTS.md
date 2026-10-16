# PC-001 — elementary actual-fibre checks

The frozen row/column identities are kept **exact**: every phase P at level l has incoming total `J_{P,l-1}`, outgoing total `J_{P,l}`, and their difference is the exact sum of microscopic sources `s_x`. Corrections `a_x` merely redistribute incoming marginals within an individual phase fibre, preserving their sum. Edge support is full H, including E_N-incident phase crossings. Zeros and `R=T(z)` phase-A ties are retained. This accounting is algebraic and does not by itself supply c_l.

**Actual geometric check.** Every upward odd→even count replacement changes the integer Fourier pair by exactly one of `(±1,±1)`; 9 distinct odd/even label pairs realize each. At x∉E_N, the SC shell argument gives at least 9 retained same-phase H successors. Exact Gibbs ratios on these edges exceed .9 for all sufficiently large N. See `PARTIAL_HALL_PROOF.md`.

**Marginal window check (dependent on PC03/PC04).** The supplied compact phase-prefix asymptotic and an analytically bounded quartic integral imply the incoming bulk interval [.18,.35]×NW and all-column lower .18×NW. E corrections cannot spoil the column lower, and their compact-fibre subtraction is negligible by the actual MT departure exponential estimate. No phase is discarded for small mass.

**Negative check.** Local degree, singleton Hall or the positivity of the marginals does not establish full weighted Hall. `verify_pc.py` supplies a small exact abstract example with individually feasible singletons and a failing doubleton; it is *not* a FIN generator.

The user-requested original PC01–PC05 file bytes were not mounted in this execution, so their hashes, original numerical margins and independently checked PC04 source proof are not asserted. The new result is conditional on those precise inherited local-limit inputs.
