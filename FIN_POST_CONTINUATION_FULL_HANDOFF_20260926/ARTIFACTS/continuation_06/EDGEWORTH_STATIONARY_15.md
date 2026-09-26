# EDGEWORTH-STATIONARY-15 — complete stationary reduced generator through O(1/N)

Status: **PROOF_GRADE_GAUSSIAN-PROJECTED_WITH_INDEPENDENT_GENERATOR_REPLAY**.

Scope: the explicitly declared finite-N heat-bath refresh model and its system-size
expansion about an interior stationary state `p*=q(mu*)`.  This is not a FIN-sourced
physical clock or a universal stochastic law.

## 1. Coordinates and assumptions

Let `Q=12`, `X` be the strict seven-column retained feature matrix, and `Y` an
Euclidean-orthonormal basis for the discarded real Fourier modes `k=1,2`.
Write

    Gamma = X^T X,       R = X Gamma^{-1}.

For an interior stationary state `p*`, put

    mu* = X^T p*,
    S   = diag(p*) - p* p*^T,
    F   = X^T S X,
    H   = Y^T S X,
    G   = Y^T S Y.

Assume `F>0`.  Define

    K = H F^{-1},
    Sigma_z = G - H F^{-1} H^T.

With `eps=N^{-1/2}`, visible and hidden fluctuations are

    x = eps^{-1} X^T(p-p*),
    y = eps^{-1} Y^T(p-p*),
    z = y-Kx.

Hence the exact tangent reconstruction is

    p = p* + eps [ d_x + Y z ],
    d_x = (R+Y K)x.

The accepted stationary decoupling theorem gives, at Gaussian order,

    L_0 = L_x + L_z,
    L_x = ((gF-I)x).grad + F:Hess,
    L_z = -z.grad_z + Sigma_z:Hess_z.

Thus `z` is an independent stationary OU residual with covariance `Sigma_z`.

## 2. Expansion of the heat-bath target

Let

    a_i = g (X x)_i,
    abar = sum_i p*_i a_i,
    atilde_i = a_i-abar,
    kappa_r = sum_i p*_i atilde_i^r.

Then

    q(mu*+eps x) = p* + eps q1 + eps^2 q2 + eps^3 q3 + O(eps^4),

with

    q1_i = p*_i atilde_i,
    q2_i = (p*_i/2)(atilde_i^2-kappa_2),
    q3_i = (p*_i/6)(atilde_i^3-3 kappa_2 atilde_i-kappa_3).

Put `m=X^T q1=gF x`.  Let `xi_i=X_i-mu*` be the centered retained row vector.

## 3. Full finite-N backward generator through O(eps^2)

For smooth visible observables, the exact jump generator has

    L_N = L_0 + eps (L_1^vis + V) + eps^2 (L_2^loc + W) + O(eps^3).

The visible local first correction is

    L_1^vis f = b1.grad f + (1/2) D1:Hess f,

where

    b1 = X^T q2,
    D1 = sum_i (d_x+q1)_i xi_i xi_i^T.

The only hidden term at this order is the diffusion modulation

    V = (1/2) sum_{a=1}^4 z_a O_a,
    O_a := T_a:Hess,
    T_a := sum_i Y_{ia} xi_i xi_i^T.

Because `sum_i Y_{ia}=0` and `X^T Y=0`, in fact

    T_a = X^T diag(Y_a) X,

so the four coupling tensors are state-independent.

The local second correction is

    L_2^loc f = b2.grad f
              + (1/2) D2:Hess f
              + (1/6) C31:third(f)
              + (1/24) C40:fourth(f),

with

    b2 = X^T q3,

    D2 = sum_i q2_i xi_i xi_i^T - x m^T - m x^T,

    C31 = sum_i (q1-d_x)_i xi_i^{tensor 3}
          + Sym_3[(m-x) tensor F],

    C40 = 2 sum_i p*_i xi_i^{tensor 4}
          + 2 [F_ab F_cd + F_ac F_bd + F_ad F_bc].

The remaining hidden second-order term is linear in `z`,

    W = -(1/6) sum_i (Yz)_i xi_i^{tensor 3}:third,

and therefore `P W P=0` under stationary hidden averaging.

## 4. Projected reduced equation and the compact memory kernel

Let `P` average over the stationary Gaussian law of `z`.  Since

    P V P = 0,

and visible `L_1^vis` preserves the `P` subspace, the Mori--Zwanzig/Duhamel
expansion through `eps^2` gives

    dU/dt = [L_x + eps L_1^vis + eps^2 L_2^loc] U
            + eps^2 integral_0^t K(t-s) U(s) ds
            + O(eps^3),

for projected visible observables (with the stationary hidden projection and no
initial-slip term).

The complete hidden-memory kernel at this order is

    K(tau) = (e^{-tau}/4) sum_{a,b}
             (Sigma_z)_{ab} O_b exp(tau L_x) O_a.              (1)

Equivalently in Laplace form,

    Khat(s) = (1/4) sum_{a,b} (Sigma_z)_{ab}
              O_b (s+1-L_x)^{-1} O_a.                         (2)

Thus the first genuine stationary hidden-memory feedback is `O(1/N)` and has
exactly four OU channels.  It is not a Gaussian friction term.

## 5. Quartic observables and the previously isolated k=5 kernel

For `f_c(x)=(c.x)^4`, define

    h_a(c)=c^T T_a c,
    c_tau = exp(A tau)c,       A=gF-I.

Then (1) gives exactly

    K(tau) f_c = 6 e^{-tau} h(c)^T Sigma_z h(c_tau).            (3)

At uniform equilibrium, `Sigma_z=I_4/12`, so

    K(tau) f_c = (1/2)e^{-tau} sum_a h_a(c)h_a(c_tau).

In label space `phi=Xc`, `h_a=Y_a^T phi^2`, hence

    K(tau) f_c
      = 6 e^{-tau} <P_H(phi^2),P_H(phi_tau^2)>_u.              (4)

For the pure retained `k=5` feature,

    K(tau) f_5
      = (lambda5^2/48) exp[-(1+2 gamma5)tau],
    gamma5 = 1-g lambda5/12.

At `g_eq=3.7183448981203875` this is

    0.11007480820919507 * exp(-1.5754981825711944 tau).

After restoring the prefactor `1/N`, this exactly reproduces the older isolated
memory-kernel component `0.1100748/N * exp(-1.575498 t)`.  It is therefore a
corollary of the single four-channel operator (1), not an unrelated effect.

## 6. Uniform channel spectrum

At `g_eq`, the retained OU decay rates are

    gamma3 = 0.3922344001359638,
    gamma4 = 0.3184370325847773,
    gamma5 = 0.2877490912855972,
    gamma6 = 0.2742466130695381.

The allowed hidden pair channels from the exact Z12 selection rule decay at
`1+gamma_r+gamma_s`:

hidden k=1:
- (3,4): 1.710671432720741
- (4,5): 1.6061861238703745
- (5,6): 1.5619957043551351

hidden k=2:
- (3,5): 1.6799834914215608
- (4,6): 1.5926836456543154
- (5,5): 1.5754981825711944

All rates are in the dimensionless declared heat-bath clock.

## 7. Independent generator replay

`edgeworth_uniform_local_check.py` compares the exact finite-N jump generator
against the displayed `L0,L1,L2` coefficients for random polynomial probes of
degrees 1--4.  Symmetric +/-eps evaluation followed by Richardson cancellation
gives maximum errors about

    2.4e-11 in L1,
    2.0e-9  in L2.

`edgeworth_uniform_hidden_check.py` repeats the check with nonzero hidden `z` and
verifies the `V` and `W` terms; maximum errors are about

    1.1e-11 in L1,
    5.5e-10 in L2.

`edgeworth_stationary_general_check.py` independently repeats the full expansion
around uniform, saddle and localized stationary states.  Stationarity residuals
are <=1.1e-16; across all three states the largest coefficient errors are about

    5.4e-11 in L1,
    1.8e-9  in L2.

The hidden covariance eigenvalues reproduced are:

uniform:
    0.08333333, 0.08333333, 0.08333333, 0.08333333
saddle:
    0.03124955, 0.03869532, 0.05683540, 0.06468681
localized:
    0.00395660, 0.00460135, 0.01352956, 0.01992300.

## 8. Relation to GENERAL-QUARTIC-THEOREM-13

The four tensors `T_a` are the operator representation of the same quadratic
hidden map

    B: Sym^2(V7) -> H4,     B(phi,psi)=P_H(phi psi),

whose exact rank is four.  Therefore the quartic closure fingerprint and the
first stationary reduced memory are controlled by the same four hidden Fourier
channels.  This makes precise the statement that linearly discarded modes are
not nonlinearly invisible.

## 9. Methodological boundary

Promoted within the declared model: complete local system-size backward-generator expansion and the memory kernel obtained by projection onto the leading stationary Gaussian hidden residual through `O(1/N)`, conditional on an interior stationary heat-bath state with invertible visible Fisher block `F`.  The finite-N correction to the true stationary conditional projector is treated separately in `FIBER_EDGEWORTH_BIAS_17.md`.

Not promoted: physical time, SI memory scale, nonstationary Edgeworth closure,
arbitrary activity rules, laboratory evidence, a spacetime dimension, a
particle interpretation, QW-2191, legacy-to-strict role transfer, `L_total`,
SM/GR, quantum-gravity or ToE closure.
