# KINETIC-NATURALITY-43
## Refinement fixes kinetic scaling shape but not relative kinetic ratios

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- **EXACT refinement/representation no-go** in the declared local quadratic
  kinetic classes;
- exact parameter count after quotienting the common clock scale;
- consistent with the repository's existing `m(ell)=mu0 ell` conditional wave
  bridge and `same statics/different kinetics` theorem;
- no physical kinetic coefficients are derived.

This continues `PHASE-HIDDEN-RECIPROCAL-DYNAMICS-42`.

---

## 1. Question

Can the same arbitrary-split refinement principle that fixed

    c(ell)=kappa0/ell

also collapse the independent kinetic coefficients of

    theta, k=1 hidden plane, k=2 hidden plane

to one common scale?

The answer is **no** under the currently admitted assumptions.

---

## 2. Additive inertial storage

For any kinetic species r, let its segment storage/mass be

    m_r(ell).

Exact subdivision additivity requires

    m_r(a+b)=m_r(a)+m_r(b)

for all positive a,b.

Under the same regularity assumption already used in the repository, the
Cauchy equation gives

    boxed:
    m_r(ell)=mu_r ell.

So refinement fixes the **length dependence** of every additive inertial
storage law.

This reproduces the existing conditional result `m(ell)=mu0 ell`.

---

## 3. Three inequivalent kinetic sectors

For the autonomous phase-hidden system there are three symmetry-inequivalent
quadratic velocity sectors:

1. phase tangent theta;
2. hidden real Fourier plane k=1;
3. hidden real Fourier plane k=2.

D12 symmetry enforces isotropy *inside* each real two-dimensional hidden plane,
but it does not identify inequivalent irreducible representations.

Therefore refinement gives

    m_theta(ell)=mu_theta ell,
    m_1(ell)=mu_1 ell,
    m_2(ell)=mu_2 ell,

with three independent positive constants.

No arbitrary-split equation relates

    mu_theta, mu_1, mu_2.

---

## 4. Why symmetry cannot equate them

The k=1 and k=2 hidden planes are inequivalent real irreducible
representations of D12.

A D12-invariant constant quadratic form on

    R_theta + H_(k=1) + H_(k=2)

is block diagonal and scalar on each irreducible block.

Hence its general form is

    diag(
      mu_theta,
      mu_1 I_2,
      mu_2 I_2
    ).

Schur-type representation rigidity removes off-block mixing but leaves one
positive coefficient per inequivalent sector.

Thus symmetry reduces tensorial freedom but does not produce universal
equipartition of kinetic storage.

---

## 5. Clock quotient still leaves two observables

A global time rescaling multiplies all inertial coefficients by one common
factor relative to the potential normalization.

Use that freedom to set, for example,

    mu_theta=1.

Then the two ratios

    rho_1=mu_1/mu_theta,
    rho_2=mu_2/mu_theta

remain.

These are dimensionless and cannot be removed by changing the unit of time.

So the obstruction is stronger than the absolute-clock gauge.

---

## 6. Distinct observable consequences

For small uncoupled quadratic curvatures

    K_theta, K_1, K_2,

the inertial frequencies scale as

    omega_theta^2 = K_theta/mu_theta,
    omega_1^2     = K_1/mu_1,
    omega_2^2     = K_2/mu_2.

Changing rho_1 or rho_2 changes frequency ratios while leaving the static
potential unchanged.

Therefore the surviving parameters are operationally meaningful in the
conditional inertial model.

---

## 7. Relaxational/Onsager analogue

Let a local reversible/gradient kinetic coefficient per segment be

    gamma_r(ell).

If one imposes the analogous additive local dissipation/storage measure,

    gamma_r(a+b)=gamma_r(a)+gamma_r(b),

regularity again gives

    gamma_r(ell)=Gamma_r ell.

D12 invariance gives the same block structure

    diag(
      Gamma_theta,
      Gamma_1 I_2,
      Gamma_2 I_2
    ).

After quotienting one overall rate, the two ratios

    Gamma_1/Gamma_theta,
    Gamma_2/Gamma_theta

survive.

Thus the no-go is not specific to inertia.

---

## 8. Refinement does not select a universal speed

The existing intercell wave bridge combines

    c(ell)=kappa0/ell
and
    m(ell)=mu0 ell,

giving a conditional speed

    v=sqrt(kappa0/mu0).

For multiple inequivalent kinetic sectors, the same logic yields distinct
conditional speeds/rates unless an additional law relates their mu_r or
Gamma_r.

Refinement by itself provides no such relation.

Therefore neither a single universal propagation speed nor a single universal
relaxation rate follows from the current refinement principle.

---

## 9. Stronger source premise that WOULD close the ratios

Any of the following would be sufficient but is currently extra input:

1. one universal kinetic density shared by all irreducible sectors;
2. a microscopic update law whose projected mobility is proven proportional
   to one common metric on theta + H4;
3. an information/Fisher metric theorem fixing the relative block weights;
4. an enlarged symmetry mixing theta, k=1 and k=2 sectors.

None is presently derived from the accepted FIN static/refinement core.

---

## 10. Exact disposition

`KINETIC-NATURALITY-43` closes with:

    REFINEMENT_FIXES_SCALING_SHAPE
    BUT_NOT_RELATIVE_KINETIC_METRIC.

The current source hierarchy is now:

    static edge law:
        shape fixed, one kappa0

    hidden->geometry interaction:
        functional form fixed conditionally, one beta

    kinetic law:
        temporal category not selected;
        even inside one category, two nontrivial relative ratios survive.

This identifies a very specific missing physical principle.

---

## 11. Next research atom

### FISHER-KINETIC-METRIC-44

The most economical remaining route is to test whether the already present
multinomial/Fisher geometry fixes the relative metric on

    theta + H4

rather than introducing a new kinetic axiom.

Required steps:

1. compute the pullback Fisher metric onto the localized phase tangent and the
   four hidden Fourier directions at the same declared state;
2. transform it into the phase-locked basis;
3. compare its three irreducible block coefficients;
4. test refinement/state dependence.

Acceptance:
- a canonical relative metric (up to one scalar) stable on the required state
  class; or
- a no-go showing state dependence / preparation dependence prevents it from
  being a universal kinetic metric.

Important boundary:
even a successful Fisher result would provide a natural metric, not by itself
prove that physical dynamics is gradient, inertial, or quantum.
