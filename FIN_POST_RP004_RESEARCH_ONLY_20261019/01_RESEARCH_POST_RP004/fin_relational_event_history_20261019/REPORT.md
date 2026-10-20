# FIN — Relational event history and anonymous return obstruction

Date: 2026-10-19. Status: **SCOPED_EVENT_INDEXED_RECONSTRUCTION_AND_EXACT_RETURN_NONIDENTIFIABILITY; F01_OPEN_SOURCE**.

## What was actually done

- Derived the event-indexed one-step Fourier correlation rule `m_k=1-lambda_k/s` from the **assumed** constant-exit-rate radial Z12 CTMC and gave the complete positive-kernel + normalisation inverse test `a=B^-1(1-m)>0`, `2sum(a1..a5)+a6=1`.
- Kept a sharp division between a relational count `n` of actual jumps and a model or physical clock. Poissonization is a consequence of the chosen CTMC, **not an independent source of stochastic dynamics**.
- Gave an explicit exact nonisomorphic pair of positive Z12 generators whose two-time equality/return correlations agree for **all continuous lags and all jump counts**; proved nonisomorphism by `11/24000` difference in the cubic edge-weight invariant.
- Finite diagnostic on decimal-frozen *strict* weights found only identity and symmetry-equivalent `k1↔k5` among 120 positive reconstructed spectral permutations. The exact counterpair is for the **wider radial family**, not the naturally selected FIN coefficient vector.

## Why it matters

A scalar, unlabelled, all-lag return correlation is not rich enough to infer a unique law, even after clock calibration. To recover normalized rates from purely event-indexed data one must provide a mathematically legitimate observation of labelled modal directions or shell-resolved relative displacements. The source/availability of that observation is itself unproved. This is a concrete F01/F06 interface constraint, not fundamental genesis.

## Boundaries

- Assumed finite reversible D12-symmetric radial single-label Markov family, **not** the native many-copy N2 heat bath; no model promotion.
- No generator/physical time derived from intrinsic nadsoliton ontology, no unique state law, no empirical data, no physical claims.
- No full-repository checkout: AGENTS is a mounted supplied historical file; the original user-owned files are unchanged.

Reproduce using `python verify.py` then independently inspect PAPER.md and the explicit exact algebra. See HANDOFF.md for the next source gate.
