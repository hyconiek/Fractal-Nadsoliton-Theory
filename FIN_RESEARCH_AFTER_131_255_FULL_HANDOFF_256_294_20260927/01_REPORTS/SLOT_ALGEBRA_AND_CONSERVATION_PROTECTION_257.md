# SLOT-ALGEBRA-AND-CONSERVATION-PROTECTION-257
## Information preservation and detailed balance do not protect transport; exact additive conservation selects SWAP

Date: 2026-09-27

Status:
exact finite-state theorem and exact generator replay.

Consider the three exchange-symmetric Z3 reset dilations

    F_c(x,y)
      =
      (y+c,x+c),

    c in Z3.

The gates c=+1 and c=-1 are inverses of one another.
F_0 is SWAP.

All three are bijections.
With fresh uniform environment trits all three realize the same one-unit full reset.

Thus one-unit heat-bath behavior does not distinguish them.

## 1. Additive-charge theorem

Take any real one-site observable f:Z3->R and define the global additive charge

    Q_f
      =
      sum_i f(X_i).

For a local gate F_c to conserve Q_f for EVERY pair input x,y we require

    f(y+c)+f(x+c)
      =
    f(x)+f(y)

for all x,y.

Set x=y. Then

    f(x+c)=f(x)

for all x.

For c=+1 or c=-1, the shift acts transitively on Z3.

Therefore f must be constant.

Hence:

    boxed:
    F_+ and F_- preserve no nonconstant additive one-site charge.

By contrast SWAP satisfies

    {x,y} after = {x,y} before,

so it preserves EVERY additive composition observable.

Thus among this exact reset-dilation family:

    boxed:
    existence of any nontrivial additive conserved density selects F_0 = SWAP.

## 2. Reversible perturbation that destroys diffusion

Define the symmetric gate mixture

    F_0 with weight 1-epsilon,
    F_+ with weight epsilon/2,
    F_- with weight epsilon/2.

Because F_+^(-1)=F_-, the continuous-time process remains reversible under the uniform measure.

It is still built entirely from bijections.

For the color Fourier observable

    chi(x)=exp(2 pi i x/3),

and spatial wave number q, the exact decay rate is

    boxed:
    lambda_epsilon(q)
      =
      rho[
        1-
        (1-3 epsilon/2) cos q
      ].

This identity was replayed on the full 3^4 state generator to numerical precision better than 1.3e-14 for nonzero target modes.

At q=0:

    boxed:
    lambda_epsilon(0)
      =
      3 rho epsilon / 2.

So ANY fixed epsilon>0 opens a nonzero density gap.

The pure n^(-2) hydrodynamic mode disappears.

## 3. Small-q form

For epsilon<2/3,

    lambda_epsilon(q)
      =
      3 rho epsilon/2
      +
      [rho/2 (1-3epsilon/2)] q^2
      +
      O(q^4).

Thus the perturbation converts diffusion into reaction-diffusion.

To keep the reaction gap parametrically below the cycle diffusion gap

    ~ 2 pi^2 rho/n^2,

one would need approximately

    epsilon << 4 pi^2/(3 n^2).

A size-dependent tuning of this kind cannot be accepted as a universal transport law unless FIN derives it.

## Verdict

The protected slow mode comes from a SPECIFIC conservation law.

It does not follow from:
- bijectivity;
- global information preservation;
- detailed balance;
- exchange symmetry;
- correct one-unit reset behavior.

This identifies the missing law precisely:

    why are the physically allowed elementary transformations restricted to those that preserve the relevant composition charges?

SWAP is protected once that premise is supplied.
The premise itself is not yet derived.
