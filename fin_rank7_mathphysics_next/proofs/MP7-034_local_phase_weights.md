# MP7-034 — local phase weights and symmetry multiplicity

Scientific state: **PROVED_INTERVAL_ASSISTED_LOCAL_CAP_WEIGHTS**.

## Scope and inputs

This task uses the accepted local equal-energy root `R7P-026`, its full-X7
positive-stability certificate `R7P-027`, and the explicit finite-N auxiliary
field model of MP7-031/032.  It is a local-cap theorem.  No global complement
mass bound is assumed or concluded.

The source hashes used by the checker are serialized in
`results/MP7-034_local_phase_weights.json`.

## 1. Full seven-dimensional Laplace determinant

For a stable mediator minimum define

`G7 = I7 - g M7 = g H7`.

At the reflection-even localized root, the exact parity splitting gives

`H7 = H4 direct_sum Hsin`.

The accepted interval LDL pivots for these two blocks therefore give a direct
interval enclosure of `det H7` by multiplying the seven pivots. Multiplication
by the certified event interval `g^7` gives

`det G7_localized in [0.1725430228889216, 0.1725430740419904]`.

Hence the single-minimum determinant prefactor is

`det(G7_localized)^(-1/2) in [2.407416508999, 2.407416865857]`.

At the uniform root, orthogonality of the supplied Fourier columns gives

`M7_uniform = diag(lambda3/12,lambda3/12,lambda4/12,lambda4/12,
                   lambda5/12,lambda5/12,lambda6/12)`.

Thus

`det G7_uniform = prod_{k=3,4,5}(1-g lambda_k/12)^2 (1-g lambda6/12)`

and direct outward evaluation gives

`det G7_uniform in [0.000354248229942486, 0.000354248232796172]`,

so

`det(G7_uniform)^(-1/2) in [53.130775359865, 53.130775573865]`.

The uniform root is therefore much broader in the Gaussian sense than one
localized minimum.

## 2. Exact localized orbit size

The equal-energy localized root is aligned and all four active amplitudes
`k=3,4,5,6` are strictly nonzero.

A translation by `a` labels can stabilize the field only if every active
Fourier phase is fixed.  In particular the nonzero `k=5` component requires

`exp(2 pi i 5 a/12)=1`.

Because `gcd(5,12)=1`, this forces `a=0 mod 12`.  Thus there is no nontrivial
translation stabilizer.  Reflection `j -> -j` fixes an aligned cosine field and
the alternating component, so the stabilizer contains exactly that reflection
and the identity.  A second reflection would compose with the first to produce
a forbidden nontrivial translation.

Therefore the stabilizer has order 2 and the D12 orbit has exactly

`24/2 = 12`

distinct localized minima.  The numerical group-action check finds twelve
translation images and a minimum center-to-center distance about `4.62949`;
this numerical value is only a separation check, not the orbit-count proof.

## 3. Explicit disjoint local caps

On the certified root box, interval evaluation of the full Hessian followed by
a midpoint-plus-radius spectral bound gives

`lambda_min(H7) >= 0.105825424562247`.

The Krawczyk image lies at least about `9.99865e-10` from every corresponding
root-box face.  Hence a mediator-space ball of radius

`r = 4.99932761855e-10`

about the exact root lies inside the validated box.  Strong convexity gives the
explicit boundary gap

`Phi(theta)-Phi(theta_*) >= (lambda_min/2) r^2
                         >= 1.32246e-20`.

Translated copies have the same gap and are disjoint by the orbit separation.
The radius is deliberately tiny because it reuses the existing local
certificate; it is sufficient for a rigorous fixed-cap Laplace asymptotic, but
it does not provide a practical moderate-N error bound.

## 4. Total localized-family versus uniform prefactor

For the twelve translated localized minima, the leading local-cap ratio at the
exact local equal-energy event is

`R_pref = 12 sqrt(det G7_uniform / det G7_localized)`

with the certified enclosure

`R_pref in [0.543733792626906, 0.543733875416134]`.

Equivalently,

`A := log 12 + 0.5 log(det G7_uniform/det G7_localized)`

satisfies

`A in [-0.609295503686632, -0.609295351426053]`.

Thus the twelvefold multiplicity does **not** make the localized family heavier
at equal energy: at Gaussian order its total local mass is only about 54.37% of
the uniform local mass.

## Nonconclusions

This is a statement about explicitly defined mediator-space local caps.  It
does not show that these caps contain essentially all finite-N equilibrium
mass, does not prove the local equal-energy event is the first global
transition, and does not supply a physical value of N or temperature.
