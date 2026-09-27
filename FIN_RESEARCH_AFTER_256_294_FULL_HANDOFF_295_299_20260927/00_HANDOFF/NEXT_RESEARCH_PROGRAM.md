# NEXT RESEARCH PROGRAM AFTER REPORT 299

Date: 2026-09-27

## Immediate numbered task

# 300 — COUPLED-MEMORY-GAP-LARGEN

### Scientific question

Does the symmetry sector that actually carries Mori–Zwanzig memory retain an O(1) hidden decay rate while the resolved metastable Z3 clock rho_N collapses with N?

The correct asymptotic object is **not** the global microscopic gap. Reports 133, 234 and 299 show that extremely slow hidden modes can be symmetry-decoupled from the memory source.

### Required quantities

For each feasible N, ideally N=9,10,11,12:

1. exact/symmetry-reduced finite-N leave-one-out Gibbs generator;
2. localized/Z3 resolved projection consistent with reports 210–234;
3. memory source vector/operator C;
4. hidden generator H on the unresolved subspace;
5. all hidden eigenmodes or certified intervals needed near the coupled edge;
6. coupling weight of each candidate slow hidden mode;
7. coupled gap gamma_c(N): smallest decay rate with rigorously nonzero/admitted coupling;
8. rho_N from the same effective lane;
9. M0(N), M1(N), epsilon_mem(N)=rho_N*M1/M0;
10. initial-slip residue Z_N and MZ relative eigenvalue error.

### Primary acceptance routes

A. Strong route:

```text
gamma_c(N) >= gamma_* > 0
```

with a proved N-uniform lower bound over the relevant coupled symmetry sector, together with an independent proof that `rho_N -> 0`.

Then for fixed theta>0:

```text
memory tail after theta/rho_N
  <= exp[-gamma_* theta/rho_N]
  -> 0.
```

B. Moment route:

Directly prove

```text
rho_N M1(N)/M0(N) -> 0.
```

By report 299's exact Stieltjes tail theorem this gives, for every fixed theta>0,

```text
[int_{theta/rho_N}^infinity K_N(t)dt]/M0(N)
  <= epsilon_mem(N)/theta
  -> 0.
```

### Failure criteria

Do not claim asymptotic Markov closure if:

- gamma_c collapses on the same scale as rho;
- a previously tiny coupling weight becomes asymptotically relevant;
- M1/M0 grows like 1/rho or faster;
- results are available only as an extrapolation from N<=8;
- numerical eigensolvers fail to separate near-zero symmetry sectors reliably.

If any of these occurs, keep an explicit memory kernel/auxiliary-variable theory.

### Computational discipline

- exploit D12 / quotient symmetry before attempting large sparse diagonalization;
- distinguish exact zero coupling from numerical near-zero by symmetry/projector arguments where possible;
- record solver residuals and tolerances;
- do not identify the smallest global hidden eigenvalue with gamma_c unless its coupling weight is certified nonzero;
- preserve finite-N source data and independent replays.

---

## Parallel foundational P0

# AMBIENT-SIMULTANEOUS-CARRIER-SOURCE

Report 295 proves only relative reversible closure. Absolute physical completeness remains unsourced.

Required target:

```text
FIN internal structure
  -> typed contemporaneous carrier C_now
  -> exclusion/non-equivalence of natural-extension history coordinates as contemporaneous resources
  -> theorem 295-A on C_now
  + record-content continuity 269/278
  -> re-test compulsory SWAP.
```

Do not solve this by naming an arbitrary time slice, by assuming factor-specific controls, or by calling the natural-extension coordinates spatial sites.

---

## Secondary fingerprint lane

Reports 297–298 give a compact operational test suite:

- seven-bin scalar sensor `Y=cos(2*pi J/12)`;
- one-time shell fingerprint at `rho*t ~= 0.542706`;
- multi-time log-curvature/rate-drift test for memory;
- late-time shell-generator reconstruction.

Next useful refinements:

1. joint nuisance treatment for rho and symmetric/asymmetric readout error;
2. finite-sample uncertainty on reconstructed q_d;
3. preregistered train/calibration/test split;
4. comparator families beyond the equalized same-rho/same-exit model;
5. N-scaling of optimal sampling times once 300 is available.

These are tests of the effective FIN lane, not fundamental-ontology proofs.
