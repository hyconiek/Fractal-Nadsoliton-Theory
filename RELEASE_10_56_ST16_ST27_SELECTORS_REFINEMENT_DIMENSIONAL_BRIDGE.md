# Release 10.56 — ST16–ST27 Selectors, Refinement, and the Dimensional Bridge

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and accompanying research software
- **Version:** 1.0.0
- **Publication date:** 2026-08-10
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Abstract

Release 10.56 executes Programs ST16–ST27, the second local
shadow-to-physics research batch for the finite strict FIN operator. The work
investigates the Hamiltonian-memory transition, selector-free observable
algebras, exact positive refinements, cross-channel predictions, interval
branch control, propagation bounds, operational minimality, gauge-connection
sources, nonlinear saturation, formal verification, dimensional scaling, and
adversarial sector refits.

The strongest constructive result is an analytic reduction of the sampled
Hamiltonian-memory transition to a one-dimensional zero-mode Jordan-chain
condition. The strongest obstructions show that stabilizer-invariant functions
of the strict operator cannot generate either the missing full observable
algebra or a nonzero oriented continuous U(1) connection, and that exact
positive translation-covariant refinement still has an infinite unresolved
fiber.

The release does not establish a physical interpretation or Theory of
Everything. It identifies more precisely which missing structures cannot be
obtained from the symmetric finite operator alone.

## Principal results

1. **Gauge-Jordan collision reduction.** For the declared Hamiltonian-memory
   realization, the zero-mode chain extends exactly when
   `(J0 v0)^T v1 = 0`, with `J v1 = v0`. Hidden-block elimination reduces this
   condition to a quadratic polynomial in inverse mediator speed. The single
   positive numerical root is `1.278014751820449`; the reduced and full-block
   invariants agree within `2.9e-15`. The root is not interval certified.

2. **Selector-free algebra no-go.** The real symmetric dihedral-invariant
   algebra on `C12` has dimension seven and already equals `C*(A)`. Every
   deterministic stabilizer-natural self-adjoint generator remains in this
   algebra and cannot generate `M_12(C)`. State-dependent or explicitly
   symmetry-broken generators remain outside the theorem.

3. **Infinite exact positive refinement fiber.** A family `A24(q)`, `q >= 0`,
   consists of positive translation-invariant graph Laplacians satisfying
   `A24(q) J = J A12` exactly. The free fine antipodal weight changes unresolved
   odd modes. Positivity, translation symmetry, and exact coarse intertwining
   therefore do not select a continuum or refinement functor.

4. **Derived cross-channel preregistration.** The mixed identity
   `U_t P_tau = exp[-(tau+i t)A]` replaces the definitional strict-kernel
   fingerprint. A hashed 400+400 synthetic holdout has no false decisions under
   its frozen threshold. It remains local synthetic evidence without
   independent custody.

5. **Interval branch control.** Fifty-five-decimal interval evaluation gives a
   minimum distinct spectral gap above `0.0435757`, a unitary principal-branch
   margin above `2.0407`, and a wave arccos branch margin above `2.2998`. This
   certifies scalar projector matching conditional on the `mpmath.iv`
   implementation.

6. **Locality bounds.** Finite-range interactions satisfy an exact path-order
   exponential-series bound. Hermitian truncations satisfy a Duhamel bound. The
   literal long-range continuation has a fitted absolute-tail exponent about
   `-0.843`, but is signed from distance eight and does not define a Markov heat
   channel.

7. **Operational minimality after packaging.** The nine ST01 maps compress to
   five obligation groups: sector, dimensional calibration, refinement,
   operational process structure, and external record/custody. Finite
   countermodels show why none of these five can be deleted within the declared
   operational semantics.

8. **Continuous connection-source no-go.** The five-dimensional space of
   translation-invariant antisymmetric one-forms has zero reflection-invariant
   subspace. A deterministic stabilizer-equivariant function of strict `A`
   cannot produce a nonzero continuous oriented U(1) connection.

9. **Saturation universality classes.** For
   `q_a(rho)=rho/(1+rho/R)^(1-a)`, all members share the same attractive local
   fourth-order field jet. With a positive linear density term, the effective
   potential is coercive for `a<1/2`, has a sharp coefficient threshold at
   `a=1/2`, and is unbounded below for `a>1/2`.

10. **Formal replay remains blocked.** A Lean/Mathlib source is exported and
    hashed, but the installed launcher has no configured toolchain. The source
    is not described as machine checked; no network installation was attempted.

11. **Multiweight dimensional equivariance.** Unitary/heat, wave, Green, and
    Gibbs channels transform with different positive-scaling weights. A single
    dimension assignment to `A` cannot make it both a first-order time
    generator and a second-order wave stiffness.

12. **Adversarial no-refit test.** Isospectral rotations pass every
    eigenvalue-only comparison and every sector-specific refit. A frozen
    matrix-level common-generator protocol rejects all sampled rotations of
    size `0.01` or larger at threshold `1e-3`.

## Evidence inventory

- 12 executed research programs;
- one analytic zero-mode/Jordan-chain reduction;
- four explicit exact positive `C24` lifts and an infinite parameter theorem;
- 800 preregistered mixed-channel synthetic holdouts;
- 55-decimal scalar interval evaluation;
- 540 isospectral adversarial channel trials;
- 6 generated figures;
- 1 exported but unreplayed Lean source;
- 15 live acceptance tests, all passing.

## Files included

- `FIN_ST16_ST27_Selectors_Refinement_and_Dimensional_Bridge_Report_EN.pdf` —
  complete English research report;
- `FIN_ST16_ST27_Selectors_Refinement_and_Dimensional_Bridge_Report_EN.tex` —
  report source;
- `fin_st16_st27_research.py` — deterministic research implementation;
- `test_fin_st16_st27.py` — live acceptance tests;
- `FIN_ST16_ST27_Results.json` — machine-readable results;
- `FIN_ST16_ST27_Summary.csv` — concise program ledger;
- `FIN_ST19_Cross_Sector_Preregistration.json` — frozen mixed-channel record;
- `FIN_ST20_Scalar_Interval_Certificate.json` — interval branch certificate;
- `FIN_ST25_Formal_Core.lean` — unreplayed proof-assistant source;
- `FIN_ST16_ST27_Figures/` — six figures;
- `FIN_ST16_ST27_INPUTS.sha256` — principal input hashes;
- `FIN_ST16_ST27_RELEASE_MANIFEST.sha256` — release integrity manifest.

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-st16-mpl python3 fin_st16_st27_research.py
MPLCONFIGDIR=/tmp/fin-st16-test-mpl python3 -m unittest -v test_fin_st16_st27.py
pdflatex -interaction=nonstopmode -halt-on-error FIN_ST16_ST27_Selectors_Refinement_and_Dimensional_Bridge_Report_EN.tex
sha256sum -c FIN_ST16_ST27_RELEASE_MANIFEST.sha256
```

The collision location and adversarial curves retain binary64 numerical trust
boundaries. The interval certificate trusts `mpmath.iv`. The Lean source is not
machine checked in the recorded environment.

## Recommended next programs

ST28–ST39 are ranked in the report. The leading studies are:

1. interval certification of the ST16 stationary branch and collision root;
2. classification of state-dependent and orbit-valued symmetry breaking;
3. classification of associative dyadic refinements inside the ST18 fiber;
4. finite-count likelihood testing for the derived mixed-channel observable.

## Scope and non-claims

Release 10.56 supplies no canonical selector or orientation, dimensional
constant, physical clock, unique continuum, derived gauge field, FIN-internal
nonlinear law, legacy-to-strict completion, legacy-role transfer theorem,
apparatus, laboratory data, independent custody, Standard Model sector,
gravitational sector, total physical Lagrangian, or Theory-of-Everything
closure.
