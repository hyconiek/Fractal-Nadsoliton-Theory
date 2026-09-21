# R7O3-035 — Target P theorem candidate

**Status:** interval-certified theorem candidate for supervisory review; not automatically merged into the authoritative accepted-result register.

For the declared rank-seven C4 shared-field model with nonnegative fields, let `M4` be the four-feature covariance matrix. Then

`lambda_2(M4) <= 67/250`.

Here eigenvalues are ordered increasingly or, equivalently, the assertion says that **at most one eigenvalue of M4 is strictly larger than `67/250`**.

## Certified chain

1. Accepted R7N tail lemmas reduce the unbounded shared nonnegative field domain to the compact hull
   `r∈[1/900,1]`, `s∈[1/128,1]`, `t∈[1/9,1]`, `y∈[10^-6,1]`, plus the certified tails.
2. The final R7N compact partition has 18,663 exact leaves: 13,231 already accepted safe leaves and the 5,432 residual parents addressed here.
3. The 5,432 R7O3 parents have exact nonoverlapping tree covers with zero unresolved leaves and 12,425 active SAFE terminals.
4. Every one of those 12,425 terminals passes the fixed-witness centered-moment certificate; a clean-directory replay independently recomputes all 12,425 and reproduces the stored rational bounds exactly.
5. The exact compact-tree join and accepted tail hashes pass independently; hostile mutations are rejected.

The centered-moment implication used on each terminal is
`Cov(B^T F) <= E[(B^T F-c)(B^T F-c)^T]`.
If the saved full-rank rational `B` satisfies
`(67/250) B^T B - E[(B^T F-c)(B^T F-c)^T] > 0`,
then the min-max principle excludes two covariance directions above the threshold.

## Gain consequence

For the four-amplitude Cartesian Hessian `H4 = I/g - M4`, the theorem candidate implies
`index_negative(H4) <= 1` for every supplied `0 < g <= 250/67`.
At the endpoint this is a non-strict spectral ceiling: it does **not** by itself exclude a zero eigenvalue.

## Nonconclusions

This does not prove the sharper `sigma` ceiling (Target S), full-X7 globality, gain provenance, physical-role/source claims, selector closure, laboratory realization, Standard Model/GR closure, `L_total`, or a theory of everything.
