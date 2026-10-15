# FIN SRN-001–SRN-005 ONLY — Native angular observation with common-scale invariance

Date: 2026-10-15. Baseline named in the source plan `380e9dc0` with preserved LI/KO/DKO/NOC intakes; no claim of local git checkout of that commit. Only SRN research appears in this package; older packages are input hashes, not copied results.

## Scientific verdict

**CONDITIONAL_NATIVE_RAW_COMMON_SCALE_INVARIANT_CLASSIFIER_CONSISTENCY**, under the inherited exact q12 A7/control spectral and MT/KO/DKO/NOC local-Gibbs, moments, positive-cylinder and labelled parent-coupling premises. Four distinct pieces have been paid in this submission: exact common-gain event invariance and positive angular mean gap; a native two-type hypergeometric conditional anti-concentration argument; a four-dimensional smooth-generator plus mollifier/Poisson/Duhamel estimate for the *raw nonsmooth* event; and asymptotic misclassification o(1) for the specified original iid-start parents with k_N=floor(N^(63/64)), M_N=floor(N^(1/4)). The proof is analytic and conditional; finite checks are bounded diagnostics, not a mechanically verified theorem.

## New SRN results only

- SRN-001: exact pointwise H(a_iW1,a_iW3)=H(W1,W3) for any finite a_i!=0, all ties preserved. p_gamma(alpha)=(2/pi)atan sqrt(1+alpha(1/gamma-1)), p_control-p_A7=d_ang alpha+O(alpha²), d_ang=(1/gamma_control-1/gamma_A7)/(2pi)>0. Separate gain static-Gaussian overlap explicitly blocks broader invariance.
- SRN-002: native conditional P(|W3±W1-z|<=eta|n)<=C_c(eta+k^-1/2) on the accepted positive MT cylinder using two hypergeometric label pairs (0,6)/(2,8), factorial maximum-atom bounds and all discrete ties.
- SRN-003: the genuine 4D chart (F1,R1,F3,R3) has conditional stationary L1 generator error C epsilon_N,k, with all off-diagonal brackets paid. Analysis-only smoothing and discrete OU Poisson corrector yield the raw observation tail bound C[eta^-2/(M alpha²)+eta^-3 epsilon_N,k/alpha+(eta+k^-1/2)/alpha]. This is **not** a growing-time Gaussian process approximation.
- SRN-004: the frozen raw angular midpoint classifier, with k_N=floor(N^(63/64)), M_N=floor(N^(1/4)), eta_N=N^-1/32 used only in proof, has both misclassification probabilities o(1) under the inherited iid-parent/local coupling. Its slowest explicit polynomial term is N^-1/64; no calibrated finite-N N0 or optimality follows.
- SRN-005: limits and one physical-interface open question, no unrequested next study.

## Scope and limitation

Observed data consist ONLY of simultaneous pairwise block Fourier aggregates at M times (or their common nonzero samplewise gain multiples). There is no robustness to separate channel gains, additive offsets/noise/erasures, nor source selection or physical apparatus validation. Local-Gibbs/MT analytical prerequisites are inherited as *conditional* source premises and not independently audited from unavailable review files. This package does not reproduce previous KO/DKO/NOC theorems as discoveries.

## Integrity and replay

Run `python verify_srn.py` to replay bounded exact projector, generator, hypergeometric, algebraic invariance, ideal-angle and exponent checks; see `INDEPENDENT_VERIFICATION.json`. `MANIFEST.sha256` binds all SRN deliverables excluding its own line and `PACKAGE_VERIFICATION.json`. The output is a submitted scoped analytical proof, ready for a separate mathematical intake — not peer review or machine proof.

**Stop.** Preserve AGENTS/history and do not automatically activate hardware, physics, lower-bound or other campaigns.

**Asymptotic feasibility caveat:** `k_N<=N/2` is guaranteed only for extremely large `N` (e.g. `N>=2^64`), not at laboratory sizes. The slow `N^-1/64` convergence and unknown constants rule out practical N0 or resource advice. Integer floors are verified with exact integer roots, not floating-point powers.
