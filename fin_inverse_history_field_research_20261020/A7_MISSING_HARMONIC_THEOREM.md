# IF-002: one mean interaction is insufficient on the actual N2 reference

Date: 2026-10-20. Conditional on the versioned A7 Fourier/spectral inputs.
All g>0 statements below are finite algebra/analysis, not a critical MT claim.

## 1. Actual two-slot inverse field, not a bare twelve-state prior

Reuse the exact HYP quotient. With q=12, S=sum_j psi_j,

    p_j=psi_j(2S-psi_j), psi_j>0,
    s_j=psi_j^2/p_j=psi_j/(2S-psi_j), f_j=log s_j.    (1)

There is one self continuation and TWO distinct labelled continuations
to every different separation. These multiplicities are essential.
Each psi_j<S, sum p=1 gives 1/2<S^2<1, and p_j<S^2.

Set x=p_j/S^2 and t=psi_j/(2S)=(1-sqrt(1-x))/2. Direct algebra gives

    f_j=log p_j-log(4S^2)+H(p_j/S^2),
    H(x)=2log[2/(1+sqrt(1-x))]
        =sum_(n>=1) [binom(2n,n)/(n 4^n)] x^n.       (2)

For 0<=x<1 the series is absolutely convergent with positive coefficients.
To verify it, differentiate: H'(x)=((1-x)^(-1/2)-1)/x, then integrate
the binomial series from zero. Finite q and max p/S^2<1 pay uniform
convergence over all j; Fourier sums can therefore be interchanged.

As a separate check, for fixed S, f as a function of w=log p has

    df/dw=1/(1-2t),
    d^2f/dw^2=2t(1-t)/(1-2t)^3>0.                   (3)

Three distinct p levels already preclude an affine relation f=theta log p+c.

## 2. Positivity in the exact Fourier convention

Use hat u_r=(1/12)sum_j u_j exp(-2pi i rj/12). For real even u, define
C_r(u)=2 hat u_r for r=1,...,5. The actual pair interaction E=A7 row/2
has only the following nonzero Fourier coefficients:

    hat E_(+/-3)=Lambda3/24,
    hat E_(+/-4)=Lambda4/24,
    hat E_(+/-5)=Lambda5/24,
    hat E_6=Lambda6/24,                              (4)

all positive. In particular C_1(E)=C_2(E)=0. Products multiply by cyclic
convolution, so every coefficient of exp(a E) is nonnegative for a>=0.
For a>0 the r=1,2 coefficients are strictly positive already in E^2:

    C_1(E^2)=(Lambda3 Lambda4+Lambda4 Lambda5
                         +Lambda5 Lambda6)/144>0,
    C_2(E^2)=(2Lambda3 Lambda5+2Lambda4 Lambda6
                         +Lambda5^2)/288>0.         (5)

The exponential series is uniformly absolutely convergent on this finite
set, justifying the convolution calculation even where E_j is negative.
Since p_j=exp(g E_j)/Z, (2) now implies

    C_1(f)>0, C_2(f)>0 for EVERY g>0.                 (6)

Indeed log p=gE-log Z contributes zero to these modes, and each H term
has nonnegative coefficients, with its n=1 term already strictly positive.
There can be no cancellation from unobserved higher powers.

Thus f is not theta E+c for ANY real theta. By IF-001 the actual Gibbs
occupation cannot be the history-entropy maximizer under only its mean
interaction. This also excludes a fitted comparison-temperature multiplier;
it does not exclude independently justified additional history constraints.
The new harmonic is a DUAL STATE WEIGHT, not a new physical mediator,
kernel term, mass, field or proof that A7 itself has changed rank.

## 3. A conservative public bound at one frozen finite-model gain

The first H coefficient is 1/2 and exp(gE)'s r=1 contribution includes
g^2 E^2/2. Consequently

    C_1(f)>=g^2 C_1(E^2)/(4 Z S^2).                  (7)

At g=12/Lambda6, the admitted exact spectral boxes give

    g>5, C_1(E^2)>1/10, g E_0<10/3.

The maximum absolute E_j is at most E_0 because (4) is nonnegative.
An exact rational exponential Taylor sum and geometric tail give
exp(10/3)<29. Hence Z<12*29, S^2<1, and

    C_1(f)>5/2784>1/600.                             (8)

The [certificate](CERTIFICATE.json) recomputes these strict rational
payments. No numerical solve for the psi vector is needed for (8).

## 4. Reference and parameter controls

At g=0 the target is uniform, f is constant and the obstruction disappears.
If one instead discards the labelled multiplicities and uses the complete
q-state reference K_jk=1, the selected process is iid p and its inverse
field is log p=gE-log Z. Then one mean does suffice. That is an explicit
different-reference countercontrol, NOT a repair of the frozen model.
The same complete-reference control applies to a single labelled slot.

The theorem does not derive Gibbs, the entropy principle or the natural
reference. It disproves one specific shortcut for removing a supplied
premise; natural F01/F03/F06 remain OPEN_SOURCE.
