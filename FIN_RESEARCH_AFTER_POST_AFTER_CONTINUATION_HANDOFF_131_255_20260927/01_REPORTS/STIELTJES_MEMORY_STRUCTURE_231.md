# STIELTJES-MEMORY-STRUCTURE-231
## Reversibility makes each sector memory kernel a positive mixture of exponentials; the one-pole model is the canonical two-moment Stieltjes approximation

Date: 2026-09-26

Status:
proof-grade linear-algebra statement for the reversible microscopic generator.

Let

    H=-QSQ >= 0

on the hidden subspace and let

    c=QSP b

for one resolved symmetry mode.

Then the exact Mori-Zwanzig memory kernel is

    boxed:
    K(t)
      =
      <c,exp(-tH)c>.

By the spectral theorem there is a positive measure mu such that

    boxed:
    K(t)
      =
      integral exp(-gamma t) dmu(gamma).

Therefore K(t) is completely monotone.

Its Laplace transform

    Khat(z)
      =
      <c,(z+H)^(-1)c>

is a Stieltjes function.

## 1. Moment meaning

The first moments are

    M0
      =
      integral gamma^(-1) dmu,

    M1
      =
      integral gamma^(-2) dmu.

The one-pole approximation

    a/(z+gamma_eff)

that matches M0,M1 has uniquely

    boxed:
    gamma_eff=M0/M1,
    a=M0^2/M1.

Thus report 230's model is the natural one-pole Stieltjes/Padé approximation,
not an arbitrary exponential fit.

## 2. Effective memory rate is a true spectral average

Because

    gamma_eff
      =
      [integral gamma * gamma^(-2)dmu]
      /
      [integral gamma^(-2)dmu],

it is a positive weighted average of hidden decay rates.

Hence it lies between the smallest and largest hidden rates carrying that
memory sector.

## 3. Cauchy bound

Let

    K(0)=integral dmu=||c||^2.

Cauchy-Schwarz gives

    M0^2
      <=
      K(0) M1.

Therefore

    boxed:
    a=M0^2/M1
      <=
    K(0).

Equality holds only for a single exact hidden decay rate.

At N=6 the ratio a/K(0) is:


    k=1:
      0.677440

    k=2:
      0.676885

    k=3:
      0.706156

    k=4:
      0.705208

    k=5:
      0.698333

    k=6:
      0.708609

So the hidden memory spectrum has nonzero width, but roughly 68-71% of the
instantaneous kernel strength is captured by the moment-matched single pole.

Despite that spectral width, report 230 shows that the resulting resolved
dynamics is already sub-percent accurate.

## 4. Consequence

The observed rapid memory decay is not merely numerical curve shape.

It follows from a positive spectral representation enforced by reversibility.

This provides a rigorous mathematical basis for replacing the memory tail by
a small number of auxiliary relaxing variables.
