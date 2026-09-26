# NONSTATIONARY-EDGEWORTH-INITIAL-SLIP-32

Date: 2026-09-26

Status:
**PROOF-GRADE CONDITIONAL AT THE LEADING GAUSSIAN-PROJECTED NONSTATIONARY LEVEL;
PARTIAL, NOT A FULL NONSTATIONARY EDGEWORTH CLOSURE.**

Scope:
the declared finite-N empirical heat-bath + ME7 model already used by
EDGEWORTH-STATIONARY-15. No physical clock or laboratory interpretation is
promoted.

## 1. Starting point

Use the accepted stationary coordinates

    eps = N^(-1/2),
    z = y - H F^(-1) x,

for which the leading generator factorizes

    L0 = Lx + Lz,

    Lz = -z.grad_z + Sigma_z : Hess_z.

The first hidden coupling in the full system-size expansion is

    V = (1/2) sum_a z_a O_a,

    O_a = T_a : Hess_x,
    T_a = X^T diag(Y_a) X.

At second order the additional hidden local term is

    W = sum_a z_a W_a,

    W_a = -(1/6) sum_i Y_ia xi_i^(tensor 3) : third_x.

In stationary equilibrium, E[z]=0 and Cov(z)=Sigma_z, so PVP=PWP=0 and the
accepted stationary memory starts at eps^2=1/N.

## 2. Nonstationary leading hidden law

Take a leading-order initial hidden Gaussian law independent of the visible
initial variable x, with

    E[z(0)] = m0,
    Cov[z(0)] = C0.

Under Lz the hidden residual is an Ornstein-Uhlenbeck process. Therefore

    m(t) = exp(-t) m0,                                           (1)

    C(t) = Sigma_z + exp(-2t)(C0-Sigma_z).                      (2)

For t>=s,

    Cov[z(t),z(s)] = exp(-(t-s)) C(s).                          (3)

These identities are exact for the leading OU residual.

## 3. First nonstationary reduced correction

Averaging the O(eps) hidden coupling V over the evolving hidden law gives

    A1_hidden(t)
      = (1/2) sum_a m_a(t) O_a
      = (exp(-t)/2) sum_a m0_a O_a.                             (4)

Hence the visible reduced generator contains

    eps A1_hidden(t)
      = [exp(-t)/(2 sqrt(N))] sum_a m0_a O_a.                   (5)

This term is absent only if the initial hidden mean vanishes (or lies in an
operator-null direction).

Therefore the statement "genuine hidden feedback starts at 1/N" is strictly a
stationary/mean-zero statement. Away from that preparation class an
N^(-1/2) initial-slip term survives.

Important structural point:
A1_hidden is a diffusion modulation, not a visible drift. It vanishes on
linear visible observables but is generically visible on quadratic and higher
observables.

## 4. Explicit hidden-preparation no-go

Take two ensembles with the same visible initial marginal and hidden means

    m0^(+) = +delta e_a,
    m0^(-) = -delta e_a,

with the same hidden covariance.

Their O(eps) visible generator difference is

    Delta L_vis(t)
      = eps delta exp(-t) O_a.                                 (6)

Choose the quadratic visible observable

    f_a(x) = x^T T_a x.

Since Hess(f_a)=2 T_a,

    O_a f_a = 2 ||T_a||_F^2.

Thus at t=0

    Delta [d/dt E f_a]
      = 2 delta ||T_a||_F^2 / sqrt(N) + O(1/N).                (7)

Whenever T_a != 0, two states with identical visible preparation but different
hidden preparation give observably different short-time quadratic evolution.

Therefore the visible initial marginal alone does NOT determine the
nonstationary reduced dynamics to order N^(-1/2). A hidden reset/stationarity
law or hidden-preparation information is an additional required datum.

This is a conditional mathematical no-go for autonomous visible closure away
from stationarity.

## 5. Strict FIN hidden channels are non-null

At the uniform strict state the four tensors have Frobenius norms

    ||T_1||_F = 1.5681170511980722
    ||T_2||_F = 1.5681170511980718
    ||T_3||_F = 1.4319036378909369
    ||T_4||_F = 1.4319036378909370.

So none of the four hidden Fourier channels is an operator-null direction.
The preparation no-go is therefore active in the frozen strict model.

For delta=1, the coefficient multiplying 1/sqrt(N) in (7) is approximately

    channel 1: 4.9180
    channel 2: 4.9180
    channel 3: 4.1007
    channel 4: 4.1007.

These numbers are dimensionless in the declared heat-bath clock convention.

## 6. Nonstationary covariance produces an aging memory kernel

For a zero-mean but nonstationary hidden Gaussian, m0=0 and C0!=Sigma_z.
Then the O(eps) term (5) is absent, but the connected two-V contraction gives

    K_NS(t,s)
      = exp(-(t-s))/4
        sum_ab C_ab(s)
        O_b exp((t-s)Lx) O_a,          t>=s.                   (8)

Using (2),

    K_NS(t,s) = K_stat(t-s) + K_age(t,s),                      (9)

where

    K_age(t,s)
      = exp(-(t+s))/4
        sum_ab (C0-Sigma_z)_ab
        O_b exp((t-s)Lx) O_a.                                 (10)

Thus a covariance mismatch produces explicit dependence on t+s, not merely on
the lag t-s. Time-translation invariance of the stationary memory kernel is
lost.

This is a precise mathematical form of "aging" / preparation memory in the
declared reduced model.

## 7. Second-order mean term

Because W is linear in z,

    E[W(t)] = sum_a exp(-t) m0_a W_a.                          (11)

It supplies an additional local O(1/N) preparation-dependent correction.

At O(1/N), time-ordered products of the O(N^-1/2) mean term must also be kept.
Equivalently one may formulate a cumulant/TCL generator using connected hidden
covariances. Equation (8) is the connected memory kernel; deterministic
mean-products are generated by iterating (5).

## 8. What is closed and what remains open

Closed at leading Gaussian-projected nonstationary order, under a factorized
initial visible/hidden law:
- the exact N^-1/2 initial-slip operator;
- the exact hidden-mean decay;
- the covariance-aging correction to the O(1/N) memory kernel;
- the explicit no-go showing the visible marginal alone is insufficient.

Still open for the full P0 NONSTATIONARY EDGEWORTH task:
1. initial x-z correlations, where E[z0|x0] is x-dependent;
2. the finite-N correction to the nonstationary conditional projector;
3. non-Gaussian initial hidden cumulants and their first surviving orders;
4. a uniform remainder theorem in a joint (t,N) window;
5. matching to an actual operational probability/readout law.

Therefore this report should not be labelled a complete nonstationary
Edgeworth closure.

## 9. Consequence for the relational-clock programme

The result sharpens the previous valuation theorem.

The hidden mean produces an O(N^-1/2) history term, but its direct operator is

    O_a = T_a:Hess.

So:
- linear visible observables are protected at this order;
- quadratic and higher visible observables are not;
- any escape-profile observable whose leading coefficients depend on the
  diffusion/quadratic sector can inherit preparation memory.

The next decisive question is therefore not simply whether hidden memory
exists. It is whether the chosen operational escape/readout map couples to
the four T_a directions at an order lower than its intrinsic transverse
profile correction.

That is the exact bridge required before interpreting gamma as a robust
emergent-time diagnostic.

## 10. Next research atom

NONSTATIONARY-EDGEWORTH-CORRELATED-33

Allow a leading conditional hidden mean

    E[z0 | x0=x] = m0 + M x + higher terms.

Derive the resulting O(N^-1/2) visible operator and decide whether the Mx term
creates effective drift after projection or only state-dependent diffusion.
Then classify which low-degree visible observables can distinguish it.

Acceptance:
an explicit operator formula plus a separation/null-space theorem for the
seven-dimensional visible polynomial space through degree four.
