# FIN ST3522–ST3671 — Ten-Round Dimension-Bridge Falsification

## Executive conclusion

Ten adaptive rounds tested whether strict heat/wave dynamics derive
\(d_w=2\) and whether the strict damping exponent \(\eta=1.8\) derives a
Hausdorff dimension. Neither identification is unavoidable.

The same finite \(C_{12}\) kernel admits inequivalent continuum classes:

- a long-range continuation with fractional symbol order
  \(\alpha=\eta-1=0.8\);
- a finite-range continuation with quadratic order \(\alpha=2\).

The finite repository data select neither.

## Rounds 1–3: heat and dual dynamics

Existence, positivity and compactness of a heat semigroup do not imply normal
diffusion. In the Cantor hierarchy,

\[
d_w=\frac{\log(2\rho)}{\log r},
\]

and every \(\rho>1/2\) yields a valid compact heat operator. The condition
\(d_w=2\) is the additional relation \(2\rho=r^2\), still leaving \(r\) free.

For a one-dimensional tail \(J(n)\sim n^{-\eta}\), \(1<\eta<3\),

\[
\lambda(k)=2\sum_{n\ge1}J(n)(1-\cos kn)
\sim C|k|^{\eta-1}.

\]

For \(\eta=1.8\):

\[
\boxed{\alpha=0.8,qquad d_w=0.8,qquad d_s=\frac2{0.8}=2.5.}
\]

A 200,000-term numerical replay gives slope \(0.8002084570\).

The dual dispersions are

\[
\omega_U(k)\sim|k|^\alpha,
\qquad
\omega_W(k)\sim|k|^{\alpha/2}.

\]

For \(\alpha=0.8\), both relevant group velocities diverge at small \(k\).
A finite nonzero wave speed requires \(\alpha=2\). Thus dual dynamics do not
derive a Lorentz cone from the literal long-range extension.

## Round 4: type obstruction

A pair-kernel tail exponent and a Hausdorff dimension are different objects.
For a standard fractional kernel in dimension \(D\), the tail behaves as

\[
|x-y|^{-(D+\alpha)}.

\]

The tail exponent combines geometry and operator order; it is not generally
equal to either. Therefore \(D_H=\eta\) remains an independent semantic axiom.

## Round 5: status of the earlier closure

The economical branch

\[
D_H=\eta=1.8,qquad d_w=2

\]

still gives the coherent values

\[
r=1.4697344923,qquad
\rho=1.0800597389,qquad
d_s=1.8.

\]

It survives as an axiom-augmented model, not a strict consequence. It is
compatible with a finite-range continuation and conflicts with the literal
long-range walk exponent.

## Rounds 6–8: product dimensions and the four-dimensional target

For a Kronecker-sum product,

\[
\operatorname{Tr}e^{-t(L_b\otimes I+I\otimes L_f)}
=\operatorname{Tr}e^{-tL_b}\operatorname{Tr}e^{-tL_f},

\]

so asymptotic spectral dimensions add:

\[
d_{s,\mathrm{total}}=d_{s,b}+d_{s,f}.

\]

The long-range base value \(2.5\) plus fibre value \(1.8\) gives \(4.3\), not
four. Imposing total dimension four would instead require fibre dimension
\(1.5\), corresponding to

\[
\rho=2^{1/3}=1.259921\ldots.

\]

A finite-range one-dimensional base gives another solution:

\[
1+3=4,qquad \rho=2^{-1/3}=0.793700\ldots.

\]

More generally, every positive split \(d_b+d_f=4\) is possible. The target
four and the decomposition are additional assumptions.

A relative product scale \(\gamma\) moves finite-window crossovers without
changing the ideal exponent. Finite cutoffs, observation windows and the order
of \(N\to\infty\), \(t\to0\) can engineer apparent plateaus.

## Round 9: falsification requirements

Distinguishing the continuum branches requires multiple sizes, multiple
calibrated times, heat traces or mode counts, explicit small-wave-number tests,
fibre-depth scaling and cutoff variation. A single finite \(C_{12}\) spectrum
cannot establish an asymptotic dimension. No external multiscale record is
present.

## Round 10: gate

| Requirement | Result |
|---|---|
| fractional-order theorem | PASS |
| dual-dispersion theorem | PASS |
| product-dimension theorem | PASS |
| strict \(d_w=2\) | FAIL |
| strict \(D_H=\eta\) | FAIL |
| unique base continuation | FAIL |
| unique four-dimensional split | FAIL |
| relative product scale | FAIL |
| physical units | FAIL |
| OA evidence | FAIL |

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

## Deepest surviving interpretation

FIN supports several mathematically coherent continuum universality classes,
including fractional long-range and normal finite-range diffusion. The strict
finite operator does not select one, and its damping exponent cannot be
identified directly with geometric dimension.

The decisive next theorem must select the continuum class and assign \(\eta\)
to one typed role. Without it, dimensional physics remains conditional.

No physical \(3+1\) spacetime, Lorentz symmetry, clock, units, laboratory
evidence, Standard Model, gravity, \(L_{\rm total}\), or Theory of Everything
closure follows.
