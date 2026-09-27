# RESET-SWAP-COMPOSITION-BUDGET-NONUNIQUENESS-288
## The one-unit Q3 law does not determine what fraction of relaxation becomes transport after closing the environment

Date: 2026-09-27

Status:
exact finite-generator family and full n=4 replay.

This report strengthens the open/closed guardrail.

The isolated effective unit has

    Q3
      =
    rho(U-I),

so its nontrivial local relaxation rate is rho.

Suppose a multiunit completion allocates this same rate budget between:

1. local fresh-bath reset;
2. nearest-neighbor record exchange.

Let

    alpha in [0,1]

be the local-reset fraction.

Then:
- local reset rate = alpha rho;
- total swap participation rate per site = (1-alpha)rho;
- on a degree-two cycle each edge swaps at rate (1-alpha)rho/2.

## 1. Same initial single-site calibration

If the neighboring visible sites are initially independent and uniform, a swap replaces
the site label by a uniform neighbor label.

Therefore:
- the reset part contributes alpha rho(U-I);
- the swap part contributes (1-alpha)rho(U-I)

to the initial one-site marginal drift.

Their sum is exactly:

    boxed:
    rho(U-I)

for EVERY alpha.

So the accepted one-unit generator and its clock do not determine alpha.

## 2. Closed collective dynamics differs

For a color Fourier mode with spatial wave number q:

    boxed:
    lambda_alpha(q)
      =
      rho[
        alpha
        +
        (1-alpha)(
          1-cos q
        )
      ].

This was replayed on the complete 3^4 state generator with relative residual below
1.4e-15.

At q=0:

    boxed:
    lambda_alpha(0)
      =
      alpha rho.

The spatial diffusion coefficient is:

    boxed:
    D_alpha
      =
      (1-alpha)rho/2.

## 3. Universality-class endpoints

### alpha=0

Pure visible-site SWAP:

    lambda(0)=0.

Color counts are conserved.

There is a gapless diffusive hydrodynamic branch.

### alpha=1

Pure local fresh bath:

    lambda(q)=rho

for every q.

No spatial transport scale exists.

### 0<alpha<1

Reaction-diffusion:

    finite color gap
      alpha rho

plus diffusion

    D=(1-alpha)rho/2.

## 4. Composition no-go

All alpha values:
- agree with the same one-unit Q3 under a fresh initially uniform environment;
- use the same calibrated rho.

Yet they predict different:
- conservation laws;
- long-time memory;
- hydrodynamic gap;
- spatial response.

Therefore:

    boxed:
    one-unit FIN dynamics does NOT determine the closed multiunit completion.

This remains true even after the SWAP gate itself has been selected.

## 5. What would select alpha=0?

Pure SWAP requires an additional global closure principle, for example:

    every event exchanges content only among the declared visible closed-system records,
    with no continuing fresh hidden reset bath.

That principle is stronger than the isolated Q3 law.

It may be motivated by:
- minimal additional environment;
- globally closed information accounting;
- content-preserving transport.

But it has not yet been derived from A7.

## Verdict

The earlier "same local clock -> unique conservative transport" narrative must be read
conditionally.

The local rate rho is sourced.

The allocation of that rate between:
- visible transport;
- hidden/open reset

is another composition datum unless a closed-system principle fixes it.
