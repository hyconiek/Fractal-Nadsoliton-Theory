# PHASE-LOCKED-HIDDEN-SCALAR-37
## Internal symmetry forbids a bare linear scalar, but the localized FIN phase supplies state-relative scalar channels

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- **EXACT representation-theoretic no-go** for a bare D12-invariant linear scalar from H4;
- **EXACT construction/classification** of state-relative linear scalars once the
  localized retained phase is supplied;
- **CONDITIONAL one-scalar candidate** under an additional point-evaluation
  principle;
- not a sourced physical backreaction law.

This continues `REFINEMENT-FORCES-BRIDGE-36`.

## 1. Hidden representation

The discarded hidden sector is

    H4 = k=1 real plane + k=2 real plane.

Write its complex amplitudes as

    z1, z2.

Under a cyclic translation by angle alpha=2 pi a/12,

    z1 -> exp(i alpha) z1,
    z2 -> exp(i 2 alpha) z2.

Reflection complex-conjugates the Fourier amplitudes.

## 2. No bare linear scalar

A D12-invariant linear scalar would be a nonzero invariant covector on H4.

But neither the k=1 nor k=2 plane contains the trivial representation.
Equivalently, one-step rotation acts by nonzero angles pi/6 and pi/3, so it has
no fixed nonzero vector in either plane.

Therefore

    Hom_D12(H4,R)=0.

So, without another state-dependent object, the four hidden coordinates cannot
linearly produce an internal-symmetry-preserving scalar stiffness/source.

This is an exact no-go.

## 3. First symmetry-preserving bare scalars are quadratic

At quadratic order the D12-invariant scalar space has dimension two:

    I1 = |z1|^2,
    I2 = |z2|^2.

Thus symmetry alone permits

    alpha I1 + beta I2.

If the physical hidden fluctuation is eps*z with eps=N^-1/2, such bare
quadratic scalars arise naturally at order eps^2=1/N rather than eps.

This is structurally consistent with the previously established fact that
stationary hidden feedback first appears at O(1/N), although it does not by
itself derive that Edgeworth result.

## 4. Localized retained phase supplies a relational selector

The accepted state-sourced geometry construction uses the intrinsic phase

    chi = m4 conj(m3) / |m4 m3| in S1,

whenever m3 and m4 are nonzero.

Because m_k transforms with Fourier character k,

    chi -> exp(i alpha) chi

under cyclic translation, while reflection sends chi to conj(chi).

Therefore chi itself transforms in the same k=1 representation carried by z1.

## 5. Two exact state-relative linear scalar channels

The combinations

    s1 = Re[z1 conj(chi)],
    s2 = Re[z2 conj(chi)^2]

are invariant under all D12 translations and reflections.

Proof under translation:

    z1 conj(chi)
      -> exp(i alpha) z1 * exp(-i alpha) conj(chi)
      = z1 conj(chi),

and similarly for k=2.

Under reflection each product is complex-conjugated, so its real part is
unchanged.

Thus the localized order parameter converts the hidden k=1 and k=2 covariants
into two internal-symmetry-preserving scalars without an externally supplied
orientation.

## 6. Classification

Among real scalars that are

- linear in the hidden amplitudes z1,z2,
- allowed to depend on the phase chi only through its Fourier characters,
- D12-invariant,
- reflection-even,

the general form is

    u_rel = alpha s1 + beta s2.

So the state-relative linear scalar space has dimension exactly two.

The phase order parameter removes the zero-dimensional no-go, but does not by
symmetry alone choose the relative weights alpha:beta.

## 7. Point-evaluation candidate

The hidden sector is already a truncated trigonometric field on the internal
circle:

    h(theta)
      = A1 cos(theta)+B1 sin(theta)
        +A2 cos(2 theta)+B2 sin(2 theta),

with canonical Fourier normalization inherited from the H4 basis.

Evaluating this actual hidden field at the cell's own intrinsic phase gives,
up to the fixed Fourier normalization convention,

    h(theta_chi)
      proportional to
      s1+s2.

Thus the additional principle

    "the emergent cell samples its hidden field locally at its own phase"

collapses the two-parameter state-relative family to one overall coupling
constant.

This is a natural, target-blind candidate because it uses only:
- the existing hidden Fourier field,
- the existing intrinsic retained phase,
- point evaluation.

But point evaluation itself is an additional constitutive/locality premise.
D12 covariance alone allows independent filtering of k=1 and k=2.

## 8. Coupling to the refinement-rigid spatial edge law

Combining reports 36 and 37 yields the conditional chain

    retained localized phase chi_i
        ->
    state-relative hidden scalar u_i
        ->
    metric edge conductance modulation

with

    delta c_ij / c_ij
      = beta0 (u_i+u_j)/2.

If point evaluation is adopted,

    u_i = h_i(theta_chi_i)

up to normalization.

Then arbitrary-split refinement fixes the functional form of the edge
backreaction; only one global coupling beta0 remains.

This is the first internally coherent route found in this continuation from

    hidden memory
      -> retained order parameter
      -> emergent geometry response.

It remains conditional.

## 9. Why this is conceptually important

Before localization:
- H4 cannot linearly source a scalar without breaking D12.

After localization:
- the retained phase itself provides the reference needed to turn hidden
  covariants into relational scalars.

So the symmetry-breaking/order parameter is not merely a label for geometry.
It can also act as the missing relational reference that makes hidden memory
capable of backreacting on that geometry.

This matches a more precise emergence picture:

    relation
      -> ordered relational phase
      -> state-relative hidden scalar
      -> change of relational geometry
      -> new relations.

No external coordinate direction is inserted at the scalarization step.

## 10. Remaining ambiguity

Two distinct source questions remain:

1. Why point evaluation rather than a D12-equivariant spectral filter
   alpha s1 + beta s2?
2. What fixes the single overall coupling between the resulting scalar and the
   intercell Dirichlet stiffness?

Until those are answered, no GR-like source law is derived.

## 11. Next research atom

### PHASE-LOCALITY-SELECTOR-38

Test whether existing FIN locality/refinement principles select point evaluation
over the two-parameter filtered family.

A useful formal target is:

Given the hidden trigonometric polynomial h(theta) in modes 1 and 2 and a cell
phase theta_chi, classify all linear functionals F[h,theta_chi] satisfying:

1. D12 covariance/invariance;
2. locality under phase refinement;
3. compatibility with restriction/refinement of the trigonometric carrier;
4. no new spectral scale or fitted mode weight.

Acceptance:
- uniqueness of point evaluation up to one scalar; or
- an explicit surviving nontrivial filter family.

This is now the shortest path toward a sourced hidden-memory backreaction law.
