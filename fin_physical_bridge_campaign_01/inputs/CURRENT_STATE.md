# FIN — CURRENT STATE AFTER 327–336 + 337-WIP

Date: 2026-09-28

## Barrier lane

- 327 repairs the proof-status weakness identified in the external review of 326.
- `Gamma=B4` status inside this packet: **PASS_DIRECTED_REPLAY**.
- External/repository acceptance is not claimed.
- Full Arb/MPFR replay remains an optional formal-audit lane.

## Preparation lane

- 328: exact discrete defect weights; D<=6 preparation tail <=0.1603% over studied N=3..12 for every kappa>=0.
- 329: exact control variable `theta=kappa/N`; fixed kappa, fixed theta and fixed KL are different protocols.
- 330: unrestricted biased heat-bath has exact biased-Gibbs detailed balance; deep-seed preparation mixes rapidly even where worst-state gap is poor.
- 331: exact defect preparation removes most of the former 1–2% empirical preparation-transfer error.
- 332: unrestricted theta=2 preparation acts as soft confinement, making a hard basin oracle unnecessary operationally for larger N if escaped/other outcomes remain explicit.
- 333: frozen future protocol = deep seed, theta=2, Tprep=4, unrestricted biased heat-bath, field off after Tprep, no wall, no postselection.
- 334: not executed; optional full MPFR/Arb barrier replay.
- 335: not executed; new-g holdout deliberately postponed until the full dynamics/response predictor is frozen.
- 336: D<=6 defect-to-post-burn response libraries constructed for N=7..10.

## Current WIP

For frozen theta=2, only 78 D<=2 states carry:

- N=7: 99.98806636%;
- N=8: 99.99094605%;
- N=9: 99.99233793%;
- N=10: 99.99307156%.

Thus the direct D>2 TV tail is only 0.01193%, 0.00905%, 0.00766%, 0.00693% respectively.

This strongly motivates a D=0/1/2 cluster-response law, but that law is **not yet completed**.

## Highest-value next task

Complete the low-defect cross-N/cross-g response law before opening a new gain. After it is frozen, execute the task335 new-g holdout under unchanged task333 control.
