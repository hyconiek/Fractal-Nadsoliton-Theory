# RI-003 — Causal nonlinear filter of the limiting classical record

**Scientific disposition:** self-contained filtering derivation for the *limiting* diffusion, at stationary quartic initial law. No finite-N filter limit, Kalman closure or physical detection assertion.

Let `b(z)=chi*lambda*z-z^3/3`, `Gf=b f'+f''`, `sigma_obs²=2/k`, `dZ=b(Z)ds+sqrt(2)dB`, `dY=Z ds+sigma_obs dW`, B,W independent, and `Z0~rho_lambda`. Set `F_s^Y=sigma(Y_r:0<=r<=s)` (usual completion), `pi_s(f)=E[f(Z_s)|F_s^Y]`. State martingale and observation noise are independent, but Y and Z are not.

**Stationary moment payment, including unbounded observation h(z)=z.** Quartic rho_lambda has `E exp(a Z0²)<infinity` for every finite a and polynomial moments of all orders. Stationarity of the scalar diffusion and Jensen give for each finite S and a>0

`E exp(a int_0^S Z_s² ds) <= (1/S)int_0^S E exp(a S Z_s²)ds = E exp(a S Z0²)<infinity`.

In particular the Novikov quadratic expectation with coefficient `k/4` is finite. Polynomial/cubic generator terms are integrable after localization, and tail events are controlled by the same moments. The signal diffusion is nonexplosive since `G z²=2chi*lambda*z²-(2/3)z⁴+2` has a Lyapunov upper bound.

**Reference measure and likelihood.** Under Q, Z has its supplied diffusion law and Y is an independent Brownian with variance `2/k` per time. For finite S define

`L_S=exp[(k/2)int_0^S Z_s dY_s -(k/4)int_0^S Z_s² ds]`.

Conditional on Z under Q this is a genuine Brownian exponential martingale in Y; stationary moment bound also pays an inverse/change-of-measure localization if required. Define physical P on F_S via L_S dQ. Under P, Y has the required drift Z; no new law is postulated. The unnormalized conditional law is

`rho_s(f)=E_Q[L_s f(Z_s)|F_s^Y]`, `pi_s(f)=rho_s(f)/rho_s(1)`, with `rho_s(1)>0` Q-a.s.

By Ito's product rule, conditional projection and independence of signal and reference Y, the unnormalized Zakai equation is

`d rho_s(f)=rho_s(Gf)ds+(k/2)rho_s(zf)dY_s`.

Quotient Ito then yields for smooth bounded f with bounded derivatives (extended by localization to polynomial test functions satisfying moment payments):

`d pi_s(f)=pi_s(Gf)ds + (k/2)[pi_s(zf)-pi_s(z)pi_s(f)] [dY_s-pi_s(z)ds]`.

With `I_s=sqrt(k/2)(Y_s-int_0^s pi_r(z)dr)`, tower conditioning shows I is a continuous local F^Y-martingale with bracket s; the L² moment payment makes it true, so Levy's characterization makes I Brownian **in the observation filtration**. It need not be the hidden W nor independent of Z.

Existence of the normalized filter is by the actual conditional distribution. Uniqueness is asserted in the **Kallianpur–Striebel admissible class** of nonnegative normalized conditional laws obtained from this frozen signal semigroup, prior, observation path and likelihood, not uniqueness of arbitrary formal measure-valued stochastic PDE solutions without constraints.

Non-Gaussianity: for f=z, the first conditional equation couples to pi(z²), and for f=z² the generator yields z⁴; there is no exact finite Gaussian/Kalman moment closure. The frozen RI-002 conditional probabilities are obtained from this *same joint (Z,Y) law* by `E[1_{C=c}1_{Q(V1)=r1,Q(V2)=r2}]` (or an appropriate backward smoothing functional), not by treating terminal C as observable at t<=1/2.

Literature benchmark: Kallianpur–Striebel / Kushner–Stratonovich filtering, checked against the primary model equations in Campillo–Corso, *A Guided Tour of the Equations of Nonlinear Filtering for Diffusion Processes*, arXiv:2606.09272 (2026). The finite-horizon stationary quartic-moment argument above pays this model's unbounded signal, rather than assuming a bounded h theorem.
