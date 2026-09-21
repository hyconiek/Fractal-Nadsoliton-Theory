# MP7-010 — complete zero-amplitude equality strata

Let `S` be the active subset of `{3,4,5}`.  Inactive phases are redundant coordinates because `I_m(0)=0` unless `m=0`.

For `b>0`, set `n=6`; for `b=0`, set `n=12`.  Define

`f_S: Z^S -> Z_n`,  `f_S(m)=sum_{k in S} k m_k mod n`.

The exact positive Fourier support is `ker f_S`.  Equality therefore means the active phase character is trivial on `ker f_S`, hence factors through `im f_S`.  Write

`d=gcd(n,{k:k in S})`, with `d=n` when `S` is empty.

Then `im f_S` has size `n/d`, and every distinct equality **field** is a label translate of the aligned field:

`phi_k = 2*pi*q*k/n` for active `k`,

with `q` taken modulo the stabilizer.  The number of distinct equality fields is `n/d`; phases of inactive modes are discarded rather than counted.

| active S | `b>0`, n=6 | `b=0`, n=12 |
|---|---:|---:|
| empty | 1 | 1 |
| {3} | 2 | 4 |
| {4} | 3 | 3 |
| {5} | 6 | 12 |
| {3,4} | 6 | 12 |
| {3,5} | 6 | 12 |
| {4,5} | 6 | 12 |
| {3,4,5} | 6 | 12 |

For the empty set, the fixed-`b>0` field is the pure alternating field and the `b=0` field is identically zero.  There is no continuum of distinct states coming from absent-mode phases.

The machine-readable table and determinant checks are in `results/MP7-007_010_phase_alignment_diagnostics.json`.
