# SRN-only scientific report

Question: can the already supplied FIN q12 mean-field parent distinguish A7 versus the frozen +0.01(P3−P4) comparator using a raw **common-scale-invariant** bounded event on a shrinking block observed at a growing number of dependent times?

Result: conditional yes, at the fixed analytical sequence k_N=floor(N^(63/64)), M_N=floor(N^(1/4)), eta_N=N^(-1/32) *only for proof*, with unknown per-time nonzero common gain allowed. The method uses a native hypergeometric anti-concentration lemma and a four-coordinate Poisson corrector, not previous cosine or sign-only theorems in disguise. Under the existing MT domain/moment/coupling premises, the misclassification error tends to zero for each of the unchanged iid-start parents. The slowest certified asymptotic polynomial from the submitted estimates is N^(-1/64), with undisplayed constants and an exponentially small boundary error under MT.

An amplitude-independent event is not detector-independent: independent gains, offsets/noise and channel mapping are unpaid. More importantly, the result is not a source for A7 or a physical prediction. All earlier campaign results are dependencies, not copied into this deliverable.

**Practical limitation:** the declared asymptotic sequence may have `k_N>N/2` at every ordinary finite `N`; the mathematical guarantee `k_N<=N/2` only becomes automatic around `N>=2^64`. These rates are proof devices, not hardware instructions.
