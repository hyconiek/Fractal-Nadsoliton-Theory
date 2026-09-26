# ACCEPTANCE-NATURALITY-51
## Multiplicative field-difference composition uniquely selects the half-density rule

Date: 2026-09-26

Status:
- exact functional-equation theorem;
- conditional source result only.

## 1. Starting class

Detailed balance gives

    f(x)=exp(x/2) psi(x),

with positive even psi.

## 2. Cocycle premise

Assume independent additive field increments compose multiplicatively:

    f(x+y)=f(x)f(y),
    f(0)=1,

with continuity (or the standard weaker regularity hypotheses).

Then the multiplicative Cauchy equation gives

    f(x)=exp(c x).

Detailed balance requires

    exp(2 c x)=exp(x),

so

    c=1/2.

Therefore uniquely

    f(x)=exp(x/2),
    psi(x)=1.

## 3. Equivalent psi proof

The cocycle gives

    psi(x+y)=psi(x)psi(y).

A continuous positive solution is exp(a x).  Evenness forces a=0.

## 4. Common reversible rules fail this stronger premise

Barker and Metropolis satisfy detailed balance but fail the multiplicative
cocycle.  Their exclusion would therefore come from the new composition
principle, not from reversibility itself.

## 5. Boundary

The cocycle is a plausible kinetic naturality condition, but it is not
currently an admitted FIN law.

If it is supplied, the localized generator becomes

    q_ij =
      (gamma/12+eta W_ij) exp[(h_j-h_i)/2].

After quotienting the overall clock scale, the remaining kinetic parameter is

    rho=eta/gamma.
