# HYDRODYNAMIC-MEMORY-CLOSURE-259
## Closing the SWAP environment converts short cell memory into a hydrodynamic power-law tail

Date: 2026-09-27

Status:
exact stochastic calculation for the SWAP cycle / infinite line.

This report follows the external audit requirement that memory decay be re-evaluated after composing many units.

Let edge swaps occur on a cycle at rate

    rho/2

per edge.

Track one color-Fourier label initially at one site.

The tagged information performs the continuous-time nearest-neighbor random walk with total jump rate rho.

## 1. Exact return-memory coefficient

On an n-cycle:

    boxed:
    R_n(t)
      =
      (1/n)
      sum_(m=0)^(n-1)
      exp{
        -rho t[
          1-cos(2 pi m/n)
        ]
      }.

On the infinite line:

    boxed:
    R_inf(t)
      =
      exp(-rho t) I_0(rho t).

Therefore

    boxed:
    R_inf(t)
      ~
      1/sqrt(2 pi rho t).

The retained information has a t^(-1/2) tail, not exponential decay.

## 2. Exact scalar Mori-Zwanzig memory kernel

Write the one-site reduced equation as

    u_dot(t)
      =
      -rho u(t)
      +
      integral_0^t K(t-s)u(s) ds.

Since

    Rhat_inf(s)
      =
      1/sqrt{s(s+2rho)},

the memory transform is

    boxed:
    Khat(s)
      =
      s+rho-sqrt{s(s+2rho)}.

Inverse Laplace transform:

    boxed:
    K(t)
      =
      rho exp(-rho t) I_1(rho t)/t.

At large t:

    boxed:
    K(t)
      ~
      sqrt(rho/(2 pi)) t^(-3/2).

The numerical Laplace replay agrees with the closed form to better than 6e-15 relative error.

## 3. Memory moments

For the infinite carrier:

    boxed:
    M0
      =
      integral K(t)dt
      =
      rho,

but

    boxed:
    M1
      =
      integral t K(t)dt
      =
      infinity.

Thus the first-memory-moment local closure used successfully inside one finite FIN cell does NOT survive the infinite hydrodynamic carrier.

## 4. New exact finite-size theorem

For a finite n-cycle, retain the conserved global mode.

The return resolvent is

    Rhat_n(s)
      =
      (1/n)
      sum_m
      1/[s+rho(1-cos q_m)].

Near s=0:

    Rhat_n(s)
      =
      1/(n s)
      +
      O(1).

Hence

    1/Rhat_n(s)
      =
      n s
      +
      O(s^2).

Using

    Khat_n(s)
      =
      s+rho-1/Rhat_n(s),

we obtain EXACTLY:

    boxed:
    M0(n)=rho,

    boxed:
    M1(n)=n-1.

This was replayed numerically for

    n=3,4,6,9,12,24,64,128

with relative finite-difference errors from about 3e-8 to 4e-5.

## Scientific consequence

There is no n-uniform "short memory" theorem for the local one-site projection.

As the network grows:

    M1(n)
      ->
    infinity.

The correct effective-theory move is therefore NOT to integrate the diffusive conserved modes into a faster memory kernel.

They must remain explicit slow fields.

This is the precise mathematical bridge from the earlier finite-cell memory programme to hydrodynamics.
