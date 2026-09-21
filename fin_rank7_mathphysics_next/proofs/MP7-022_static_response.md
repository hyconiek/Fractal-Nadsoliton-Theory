# MP7-022 — static response with explicit source conventions

For an aligned stationary branch let

`F(s,g)=s/g-mu(s)=0`, `M=d mu/ds`, `H=I/g-M`.

## Gain derivative

Differentiate `F=0` along a nondegenerate branch:

`H ds/dg - s/g^2 = 0`,

so

`ds/dg = H^{-1}s/g^2`.

The formula is valid only while `H` is invertible; it cannot be continued through a fold by ordinary inversion.

## Source conjugate directly to mediator amplitude

Adding `-f^T s` to the dual potential changes stationarity to `grad Phi=f`. Therefore at a stable/nondegenerate root

`ds/df = H^{-1}`.

This is the susceptibility to the **dual mediator source** `f`.

## Microscopic feature field inside the exponential family

Instead insert `h` in the Gibbs field and solve

`s = g mu(s+h)`.

Linearizing gives

`(I-gM) ds = g M dh`.

Because `s=g mu` on the branch,

`d mu/dh = (1/g) ds/dh = (I-gM)^{-1}M`.

The inverse exists iff `I-gM` is nonsingular. As `g->0`, the response tends to the covariance `M` of the supplied exponential family.

These two source experiments are different. `H^{-1}` is not interchangeable with `(I-gM)^{-1}M`. The covariance is the Fisher information matrix for the finite exponential family in the declared feature coordinates; no spacetime interpretation follows.

## Quantitative certified branch response at the coexistence root

Using the full interval root box from `R7P-026_equal_energy_event.json`, the accepted
spectral intervals, and direct interval inversion of the C4 Hessian box, the gain
continuation response is enclosed by

```
ds3/dg in [1.2015128668728788, 1.2015129485905275]
ds4/dg in [1.2790600695187322, 1.2790601580950933]
ds5/dg in [1.3500919102461282, 1.3500920049771038]
ds6/dg in [0.9633378506440589, 0.9633379203717149].
```

Hence all four aligned amplitudes increase with gain on this certified localized
branch at coexistence.  This is a local branch statement, not a global monotonicity
theorem for arbitrary stationary solutions.

The calculation is replayed by `scripts/mp7_022_023_response_quant.py` and recorded
in `results/MP7-022_023_quantitative_response.json`.
