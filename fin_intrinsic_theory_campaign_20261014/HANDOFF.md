# FIN intrinsic-theory campaign handoff — TH-001–TH-005

Date: 2026-10-14.

## Disposition

The first intrinsic-theory package has produced one theorem-level advance and one sharply isolated critical-limit obligation:

1. **TH-002 path LDP:** Kraaij Theorem 2.3 applies to the actual leave-one-out q12 empirical chain for initial sequences with a good initial LDP. Deterministic starts and the declared Gibbs equilibrium starts are explicitly paid examples. The resulting action is the non-quadratic edge-flux relative-entropy action.
2. **Local quasipotential:** the exact reversal symmetry gives a basin-local `F(p)-F(min)` quasipotential along reversible relaxation orbits. This does not pay the global Gamma/capacity gate.
3. **TH-003 information result:** seven retained means are not path-sufficient; an exact N=120 count-state witness separates their fluctuation generator. Seven means plus four existing second moments form an exact 11D occupation lift.
4. **TH-004 microscopic fold bridge:** the actual jump generator gives soft noise coefficient exactly `2g_f`, not an assumed auxiliary Langevin covariance. The supplied interval gives variance coefficient `[7.0312894136791835,7.031289453679183]`.
5. **TH-005 stochastic fold:** the critical scaling and candidate SDE are derived at generator-characteristic level. Full stopped process convergence is **not** claimed. One precise lemma remains: fast-mode/hidden-mode averaging and containment on the accelerated time scale.

## Strongest justified scientific statement

The supplied FIN heat-bath model now has a rigorous interacting path large-deviation theory in its declared mathematical scope, plus a source-derived critical diffusion candidate whose drift/noise coefficients are fixed by the microscopic generator and certified fold. The next proof-grade move, if the architect continues this lane, is **only** the fast-mode averaging lemma needed to turn the stochastic-fold generator balance into a stopped convergence theorem.

## Architecture boundary

TH-006–TH-009 were not executed. TH-010 is reserved for architect review; this handoff does not choose a new physical programme. Empirical promotion remains stopped without independent input.

## Verification and replay

- `independent_verify.py` rebuilds the kernel, rank-seven factor, fold, fluctuation witness and Hamiltonian identities without importing `theory_core.py`: **41/41 checks PASS**.
- A clean-directory replay regenerated `NUMERICAL_REPLAY.json`, `TH-005/FOLD_GENERATOR_DIAGNOSTICS.json` and `INDEPENDENT_VERIFICATION.json` byte-identically.
- `PACKAGE_VERIFICATION.json` checks the required task contract and JSON parseability.
- The independent verifier is an algebraic/numerical replay. It does not substitute for Kraaij's published proof and does not prove the remaining TH-005 fast-mode lemma.

## Source/provenance note

The user supplied the campaign, derived-bridge synthesis and theoretical literature database in this conversation. Selected repository sources at baseline `1b7cd447` and the full statement of Kraaij Theorem 2.3 were inspected through their public raw/HTML sources. Web-inspected bytes were not materialized locally, so this handoff records their exact URLs rather than inventing SHA-256 values. Local supplied inputs are hash-bound in `SOURCE_INPUTS.sha256`.
