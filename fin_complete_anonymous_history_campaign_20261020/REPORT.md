# FIN complete anonymous history — terminal research report

**Date:** 2026-10-20. **Tasks:** AH-001, AH-002 (activated after positive AH-001), AH-003. **Scientific status:** mathematical theorems conditional on the frozen **one-label 12-state positive reversible radial jump Markov chain**. Not an independent journal review, experimental result, intrinsic FIN law, native multi-copy leave-one-out kinetics or proof of ToE.

## New results

**AH-001: positive full-history isomorphism theorem.** For two positive twelve-state radial jump matrices `J` and `K`, equality of the full collection of first-appearance equality-pattern prefix laws forces `K=PJP^T` for some permutation matrix `P`. All states are covered almost surely; empirical transition frequencies of the anonymously renamed states converge to a random permutation conjugate of J. Equality of finite-cylinder distributions implies equality of the law of that estimator, whose support lies on a finite permutation orbit. The orbits must intersect and hence coincide.

**AH-002: explicit finite observation length.** Let `b=10a^2`, `0<a<=1/11`, `d_iso(J,K)>=delta>0`, `0<eta<1`. A sufficient integer is

```
m >= max(2/b, (8/b^2) log(24/eta),
                (384/(delta^2*b^2)) log(528/eta))
L=2m+2.
```

There is a deterministic rule using only a single equality prefix of length L with conditional error <=eta under either alternative. This is a conservative theorem, not a practical resource estimate. It uses the two-step lower bound `J²(i,j)>=10a²`, one submartingale lower-tail estimate for visits, 132 edge martingale Azuma bounds, and a triangle inequality on the quotient space by permutations. No independence of path increments was presumed.

**AH-003: exact scope boundary.** Complete anonymous path law identifies J only up to arbitrary relabelling, and cannot infer the overall continuous-time activity `s` in `Q=s(J-I)`, the physical meaning of `Z12` coordinates, or an independently sourced evolution mechanism. Full mean bracket with known pi is a distinct, labelled one-hot data type; it was not treated as anonymous.

## Relation to earlier studies

The existing E2 all-lag return nonidentifiability and E3 four-event witness are preserved: a particular pair separates at length four, but general uniqueness requires *all* finite prefixes. Alternating moment recovery remains a weaker observation hierarchy. The six-mode labelled inverse and mean-bracket correction are not rerun. New theorems do not claim literature priority; identifying a Markov chain from sufficiently rich path statistics is a standard broad phenomenon, whereas this packet pays the precise frozen anonymous contract and one explicit conservative finite-prefix bound.

## Independent local diagnostics

The included `verify.py` uses exact rational arithmetic to check positive radial normalization, two-step minorization and the alternating-moment example on `w=(1,2,3,4,5,6)`; it performs a **short, seeded 3200-transition diagnostic** to verify first-appearance code reconstruction up to the induced random permutation. It also evaluates the sufficient analytic sample-size formula without attempting a simulation of that astronomical length. Checks pass after correction of an initial diagnostic implementation error (failure to update the sampled current state). A numerical script's PASS is not a proof of ergodicity, Azuma-Hoeffding, or the path-measure theorem.

## Provenance and stops

The original live path `/media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison` is not mounted here. Available and hashed are the user's current exact AH assignment, archived `AGENTS(1).md`, GEN archive and predecessor local history proofs. The newer upstream intake, audit and correction are available only as verbatim excerpts in the current supplied assignment, not as independently authenticated live files. Duplicate check was restricted to these available sources; no priority claim.

No source files, historical strategies, AGENTS, physical premises or Git checkout were modified. No hardware, physical experiment, spending, externally sent message, commit or push.

**Final status:** `AH001_PROVED_CONDITIONALLY; AH002_EXPLICIT_FINITE_BOUND; AH003_SOURCE_GATE_UNCHANGED; F01_F03_F06_OPEN_SOURCE`. Stop at architecture.
