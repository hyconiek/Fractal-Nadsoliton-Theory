# MP7-024 — reusable response theorem interface

Scientific state: **PROVED_INTERFACE_SCOPE**.

For a nondegenerate stationary solution of `s=g mu(s)` with `H=I/g-M`, the following
interfaces are distinct.

| Experiment | Perturbation | Response | Required gate |
|---|---|---|---|
| gain continuation | change `g` | `ds/dg = H^{-1}s/g^2` | `H` invertible |
| dual source | add `-f^T s` to `Phi` | `ds/df = H^{-1}` | nondegenerate root; stability if used as susceptibility |
| microscopic feature field | `s=g mu(s+h)` | `d mu/dh=(I-gM)^{-1}M` | `I-gM` invertible |

The formulas are coordinate statements relative to the supplied Euclidean dual metric.
After a linear coordinate change, transform both covariance and quadratic metric as in
MP7-019. No formula is to be continued through a fold by ordinary inversion.

For numerical exports, a response record must include: branch/root certificate, gain
interval, coordinate metric, source convention, denominator/inverse margin and the
observable contraction actually reported. Missing global acceptance must remain visible
on every global corollary.

## Certified fixtures now available

The interface has two quantitative fixtures:

1. `R7P-026` localized coexistence root: direct interval inversion certifies all four
   gain-response components, recorded in `results/MP7-022_023_quantitative_response.json`.
2. `R7P-031` simple fold: MP7-020 plus `H4 v=0` yields a certified leading-covariance
   spectral separation greater than `0.0164429020`.

These fixtures provide regression targets for future response calculations without
turning the supplied feature normalization into a physical unit convention.
