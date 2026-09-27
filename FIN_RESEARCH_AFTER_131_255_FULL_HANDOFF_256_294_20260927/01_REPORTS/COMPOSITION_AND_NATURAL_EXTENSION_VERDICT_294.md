# COMPOSITION-AND-NATURAL-EXTENSION-VERDICT-294
## Finite reversible convergence is controlled; minimal closed-content selects SWAP conditionally; FIN retains measurable internal fingerprints beyond the coarse clock

Date: 2026-09-27

This report synthesizes 290-293 and the current repository guardrail.

## 1. Finite-to-natural-extension convergence is solved at process level

For any primitive finite Markov kernel P, the finite cyclic bridge

    mu_L
      proportional to
    product P(x_i,x_(i+1))

with deterministic cyclic shift converges on every fixed cylinder to the stationary two-sided Markov natural extension.

For the Q3 full-reset event kernel U, every non-wrapping finite window is already EXACT.

Thus the finite cyclic carrier has a rigorous process-convergence interpretation.

## 2. This does not source geometry by itself

Base cylinder statistics do not expose all internal recodings/factorizations.

Operational subsystem access is still required.

So convergence to natural extension is not equivalent to derivation of physical space.

## 3. Closed-system minimality gives a conditional alpha selector

In the mixed reset/SWAP family, conditioned on the same event schedule, M hidden reset events require at least

    M log_2 q

fresh hidden-content bits.

Pure visible SWAP requires no new target-symbol entropy.

Therefore if:
- the visible records are declared complete;
- extra hidden content is minimized;

then

    alpha=0

is uniquely selected.

This is the strongest current source argument for conservative closure.

It remains conditional on the definition of the complete system.

## 4. Internal symmetry still does not give higher spatial rank

Any canonically selected internal symmetry must lie in the center of the internal automorphism group.

For Q3:

    Z(S3)={e}.

For Q12:

    Z(D12)={e,r^6},

and r^6 is still an internal antipodal state transformation.

Therefore no second subsystem-translation direction is sourced.

## 5. FIN-specific information survives above the coarse level

Although pure long-wavelength SWAP diffusion is universal after rho calibration, the resolved 12-state effective chain has dimensionless mode ratios

    -lambda_k/rho

that depend on the shell structure q_d.

At N=8, representative hidden ratios are:

    R1≈1.4197
    R2≈1.6265
    R3≈0.9308
    R5≈1.2205
    R6≈1.2067.

A comparator with the same:
- rho;
- Z3 generator;
- total exit rate

still differs by about 0.12 in hidden-mode response at t=1/rho.

So FIN-specific tests should target the pre-hydrodynamic resolved dynamics, not universal diffusion alone.

## 6. Current strongest chain

    exact finite-N FIN process
      ->
    memory-aware 12-state effective chain
      ->
    exact Z3 quotient
      ->
    canonical natural extension
      ->
    finite cyclic Markov-bridge approximation
      ->
    [complete-system + minimal hidden-content premise]
      ->
    pure visible content transport
      ->
    conservative SWAP hydrodynamics.

Every remaining conditional arrow is now explicitly named.

## Next proof-grade programme

### P0 — COMPLETE-SYSTEM-CRITERION-295

Can FIN define operationally when all relevant degrees of freedom have been included, rather than declaring closure by hand?

Use:
- preparation/intervention/readout completeness;
- absence of residual memory after enlarging the state;
- information-balance closure.

### P0 — MICROSCOPIC-PERIODIC-BRIDGE-296

Apply the finite cyclic-bridge construction directly to a sampled skeleton of the exact finite-N leave-one-out Gibbs process and quantify cylinder convergence before the metastable reduction.

This tests information continuity at the actual microscopic lane.

### P1 — FINGERPRINT-OPTIMAL-PROBE-297

Find the lowest-complexity observable maximizing discrimination between:
- actual FIN q_d;
- same-rho/same-exit countermodels.

Include measurement noise.

### Parallel — B4 proof lane

Continue capacity asymptotics only where it contributes to a rigorous large-N communication exponent.

## Main update

The locality programme is no longer blocked by a vague notion of reversible history.

The remaining fundamental ambiguity is now:

    boxed:
    what criterion tells FIN that a chosen set of variables is the COMPLETE closed system?

Once that is answered, the minimal-hidden-content theorem has a concrete route to selecting conservative transport.
