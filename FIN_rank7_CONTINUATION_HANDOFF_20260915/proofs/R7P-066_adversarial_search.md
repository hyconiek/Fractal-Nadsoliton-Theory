# R7P-066 — bounded adversarial off-face search ledger

The search uses the exact R7P-067 compactification, so its numerical domain is
the whole closed probability-law closure `[0,1]^4`, not a finite `J` cutoff.
The objective is `lambda2(M4)-sigma_*`.

Three fixed DE seeds (`66066,66166,66266`) were run, together with the exact
extreme-face slice, finite `J6=2,4,8,12,20` slices, and one-sided `J4/J5`
perturbations around the known equality locator. All DE runs returned the same
boundary point to floating precision. The largest positive display residual is
~`2.8e-16`, consistent with roundoff at the exact double root.

Every finite candidate in the ledger is recomputed through two implementations:
(1) compact weights and (2) direct exponentials in `J`. Their probability and
covariance discrepancies are recorded. No search result is promoted to a proof.
If a future run produces a stable positive gap, R7P-070 must certify it before
any ceiling revision.
