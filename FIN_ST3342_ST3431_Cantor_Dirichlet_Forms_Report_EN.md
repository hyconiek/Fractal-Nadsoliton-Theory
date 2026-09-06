# FIN ST3342–ST3431 — Cantor Dirichlet Forms and Metric-Energy No-Go

## Dirichlet family

On the depth-\(N\) binary tree, assign conductance

\[
c_\ell=\rho^\ell,
\qquad \rho>0,
\]

and energy

\[
\mathcal E_\rho(f)=\frac12\sum_e c_e|df(e)|^2.
\]

Every \(\rho\) gives a positive Dirichlet form with one constant zero mode.
The Cantor boundary topology and uniform Bernoulli measure are unchanged.

At depth six, the same 127-node tree gives:

| \(\rho\) | spectral gap | largest eigenvalue |
|---:|---:|---:|
| 0.5 | 0.0016936039 | 1.75936535 |
| 1 | 0.0085151004 | 5.54832478 |
| 2 | 0.0241407463 | 225.19876449 |

Thus metric topology and measure do not determine diffusion dynamics.

## Competing normalization axioms

The total unit-jump energy at level \(\ell\) scales as

\[
(2\rho)^\ell.
\]

Requiring constant total layer energy selects \(\rho=1/2\). Requiring constant
conductance per edge selects \(\rho=1\). Fine-scale stiffening selects
\(\rho>1\). These natural conditions are mutually inequivalent.

## Metric-measure-energy exponents

With \(\delta_\ell=r^{-\ell}\),

\[
\mu_\ell=2^{-\ell}=\delta_\ell^D,
\qquad
D=\frac{\log2}{\log r},
\]

and

\[
c_\ell=\rho^\ell=\delta_\ell^{-\sigma},
\qquad
\sigma=\frac{\log\rho}{\log r}.
\]

Equivalently \(\rho=r^\sigma\). This relation changes coordinates but does not
select \(r\) and \(\sigma\). An overall energy/time multiplier also remains.

## Eta ambiguity

Strict \(\eta=1.8\) could be assigned to Hausdorff dimension, energy exponent,
or walk dimension. These produce different metrics and spectra. No current
theorem assigns its semantic role on the refinement tree.

## Verdict

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

FIN supplies a natural fractal state space and a large family of compatible
Dirichlet forms, but not a unique resistance form, diffusion, spectral
dimension, clock or physical geometry.

The next decisive theorem must be a variational or renormalization fixed-point
principle jointly selecting \(r\) and \(\rho\). Without it, the metric-energy
family must remain explicit.

No selector, physical scale, apparatus evidence, Standard Model, gravity,
\(L_{\rm total}\), or Theory of Everything closure follows.
