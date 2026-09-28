# FIN — CURRENT STATE AFTER RESEARCH 327–330

Date: 2026-09-28

## 1. Barrier / metastability lane

The external review of 300–326 correctly challenged the proof status of original report 326 because not every exclusion decision in the global branch-and-bound had a surfaced propagated rounding certificate.

Research 327 independently replayed the separator calculation with directed interval decisions. The robust six-dimensional half-wall cover processed 3,402,039 boxes with zero unresolved cells; the triple-junction five-dimensional cover processed 59,271 boxes with zero unresolved cells. Critical boxes were independently replayed at 70-digit interval precision on a high-precision mathematical kernel. The d4 root, local Hessian, and explicit path were also replayed at high precision.

Current status:

`Gamma = B4` is **PASS_DIRECTED_REPLAY** inside this research packet, but external/repository acceptance is not claimed. A full Arb/MPFR replay remains a formal-audit option.

Consequently the report-324 implication to the capacity exponent can be used conditionally on acceptance of the 327 proof standard. The Eyring–Kramers prefactor remains open.

## 2. Preparation law lane

Research 328 replaces an empirical preparation correction by an exact finite-defect representation. For defect count `D=N-n0`, the pinning family obeys

`mu_{N,kappa}/mu_{N,0} proportional to exp[-(kappa/N)D]`.

The exact relative configuration weights include the finite multinomial factor, exact one-defect costs and exact defect-defect interaction term. Truncating to `D<=6` retains only 12,376 defect configurations and has a worst upper preparation-TV tail below 0.1603% over all studied N=3..12 and every `kappa>=0`.

This is currently the strongest positive, probability-preserving preparation representation in the finite-N lane.

## 3. Control parameter lane

Research 329 establishes that the direct defect field is

`theta=kappa/N`.

Fixed kappa, fixed theta and fixed KL cost are different experiments. Over the current few-defect regime, fixed theta produces much more stable preparation statistics than fixed kappa. This does not select a physical protocol; it only identifies the correct conjugate parameter of the declared statistical family.

## 4. Preparation mechanism lane

Research 330 proves detailed balance of the unrestricted biased leave-one-out heat-bath chain with the full biased Gibbs distribution.

Hard reflection at the declared `J=0` basin preserves conditional weights on each connected component, but the basin graph is generally disconnected. Therefore hard reflection + arbitrary initialization does not realize the full conditional distribution exactly.

For deep-seed initialization, however, the seed component carries nearly all target mass. At fixed `theta=2`, its missing mass is <=1.47e-6 at N=4 and falls to ~4.9e-14 by N=10.

Worst-state spectral mixing is not uniform: N=9 has a slow pair with gap ~0.03725. But those modes have negligible overlap with the deep seed. Direct seed-to-target semigroup propagation gives TV<1% by t≈1.5–1.75 and TV<0.1% by t≈4–4.25 for every N=3..10 tested.

Thus the operational preparation statement is now:

`deep seed + fixed control theta + reflecting biased heat-bath + finite preparation time`

approximates the declared basin-conditioned target with two explicit errors:

1. missing connected-component mass;
2. finite-time mixing error.

The remaining physical gap is sourcing the reflecting/feedback controller without an oracle-like basin label.

## 5. Current highest-value next task

P0 is **331 — DEFECT-BASED-OBSERVED-PROCESS-PREDICTOR**: insert the exact D<=6 preparation law into the observed-process pipeline and test whether it explains/reduces `eps_prep` without PCA or per-N refitting.

A new blind test should preferably change `g` or the declared preparation protocol after the model is frozen, rather than only increment N again.
