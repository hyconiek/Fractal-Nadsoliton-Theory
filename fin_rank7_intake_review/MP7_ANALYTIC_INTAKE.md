# MP7 intake: phase alignment and selected interval certificates

Date: 2026-09-22. Source: `../fin_rank7_mathphysics_next/`.

## Accepted statements in the supplied finite model

1. **MP7-007–011 (analytic).** For fixed nonnegative amplitudes and alternating amplitude `b>=0`, the twelve-label partition function obeys `Z(phi;b)<=Z(0;b)`. Equality fields, after discarding phases of absent modes, are label translations of the aligned field. An odd translation handles `b<0`. Every global minimizer of the specified full X7 primal/dual model is therefore translation-equivalent to an aligned nonnegative C4 representative. This says nothing about arbitrary stationary points.
2. **MP7-016 (interval assisted, at `g=37/10`).** The aligned nonnegative stationary system has exactly three roots: uniform, localized and saddle. Both nonzero roots have strictly positive certified energy. With MP7-011, the uniform state is the unique full-X7 global minimizer at this gain.
3. **MP7-017 (interval assisted, supplied spectral tuples).** The certified localized/uniform equal-energy event in the R7P-026 gain box is the first global transition: below it the uniform state is uniquely global; at it the uniform state and the twelve localized D12 images are the complete minimizing set. The event lies in `[3.7183448971203875, 3.7183448991203876]`. No claim about every larger gain follows.
4. **MP7-020 (interval assisted, declared shared-field C4 domain).** With decreasing eigenvalue order, `lambda2(M4)<=0.267999946371058`; the rigorous normalized separation from `67/250` is at least `6500000000000000000000000000/121203210760018863485060548125326941` (about `5.36289423295e-8`). This strengthens the separately accepted Target P result in its stated normalization. The old smallest PD minor is not itself a covariance spectral gap.

## Checks and provenance

The analytic proof was checked directly: an absolutely convergent positive Fourier expansion of `exp(a cos x)` gives nonnegative coefficients after the twelve-label character average. Pairing opposite frequencies proves the inequality. Equality makes the phase character trivial on the support kernel; the kernel/image calculation handles zero-amplitude strata. The mediator quadratic cost is phase independent, and the finite Gibbs variational identity connects dual and primal minima.

The absent ZIPs are packaging artifacts, **not missing scientific inputs**. All 54 predecessor files named in `replay/MP7-045_INPUT_MANIFEST.sha256` match SHA-256 against the unpacked `FIN_R7N_HANDOFF_20260920`, `FIN_R7O3_TARGETP_HANDOFF_20260920`, `fin_rank7_followup`, and `fin_rank7_intake_review` directories. All 165 present entries of the MP7 package manifest match; the four missing entries are the three ZIPs and their archive checksum list.

MP7-016, MP7-017 and MP7-020 were rerun with those unpacked inputs and output directed to `/tmp/mp7_intake_replay_20260922`. Their scientific JSON fields matched the package exactly; MP7-020 differed only in `elapsed_seconds`. This is a fresh local replay of the same checker implementations, not independent reimplementation or proof-assistant formalization. One `work/scripts/mp7_022_023_response_quant.py` hash in the older clean-replay input manifest differs from the current package script; it does not enter these three reruns and remains a separate provenance item for MP7-022/023.

Other MP7 interval-assisted, finite-N and dynamical claims remain package results awaiting targeted intake. Existing R7O3 Target P acceptance is not reopened. No selector, physical gain/time/temperature, laboratory realization, kernel bridge, role transfer, `L_total`, SM/GR or ToE closure follows.
