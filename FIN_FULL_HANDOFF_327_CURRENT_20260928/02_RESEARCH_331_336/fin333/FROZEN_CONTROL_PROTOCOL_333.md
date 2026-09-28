# FIN 333 — FROZEN CONTROL PROTOCOL
## Freeze before a genuinely new g/protocol holdout

Date: 2026-09-28

Status: **FROZEN BEFORE NEW HOLDOUT**

Based only on already-opened N<=10 calculations from tasks 329–332, freeze the following preparation protocol:

1. initial microscopic state: `n0=N`, all other counts zero;
2. preparation dynamics: unrestricted leave-one-out heat-bath at the declared production gain `g`;
3. control field: `theta=2` on label 0, equivalently `kappa_N=2N`;
4. preparation duration: `Tprep=4` generator-time units;
5. switch field off exactly after Tprep;
6. continue with the unmodified production generator;
7. no hard J=0 wall;
8. no postselection;
9. retain all localized outcomes plus an explicit unlocalized/other-phase outcome.

No future holdout may retune `theta` or `Tprep` and still be called a test of this frozen protocol.

Frozen JSON SHA-256:

`b7ae70e9b8e482e7a180a960a97f099d01b9f50403596f17817c5f194f10d3ff`

## Why theta=2

Task 329 found that fixed `theta=kappa/N` is materially more stable across the few-defect regime than fixed kappa. Task 332 then showed that theta=2 gives small escape and sub-0.1% conservative preparation error by t=4 for N>=7, without a basin oracle.

This is an **operational model choice**, not a fundamental law of FIN.

## Next valid holdout

Choose a gain or control condition not used to select this protocol, freeze the full observed-process predictor before opening it, then evaluate all outcomes without postselection.
