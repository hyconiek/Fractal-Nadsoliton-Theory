# IF-001: what field would actually supply a given occupation law?

Date: 2026-10-20. Finite symmetric-reference mathematics; the general
Perron/entropy variational method is known, not a new FIN discovery.

## 1. Exact variational identity, without presupposing Markov histories

Let K be finite, connected, symmetric and 0/1, K_xx=1. For any real state
field f define T_f(x,y)=K_xy exp[(f_x+f_y)/2]. Finite Perron--Frobenius
gives its simple largest eigenvalue rho>0 and positive real eigenvector
psi normalized by sum psi_x^2=1. Positivity on a connected support with
self loops ensures primitivity. Put

    P_f(x,y)=T_f(x,y) psi_y/(rho psi_x),
    pi_f(x)=psi_x^2.                                 (1)

Rows sum to one; symmetry of T_f gives detailed balance and stationarity.
For any stationary legal history with marginal nu, flow d and entropy
rate h, its expected log P_f is nu(f)-log rho: the log psi endpoint terms
cancel by the two equal margins. Therefore

    log rho-[h+nu(f)]
       =KL(d || nu_x P_f(x,y))
          +[H(X_1|X_0)-h]>=0.                        (2)

The formula also allows zero nu entries, with the usual zero-mass convention.
Equality forces d=nu P_f, then stationarity/irreducibility force nu=pi_f.
Equality of finite conditional entropies with the one-step entropy forces
the actual full history to be Markov. Thus (1) is the unique stationary
maximizer of h+nu(f), not merely an assumed Markov ansatz.

## 2. Inverse field and its exact gauge

For supplied pi>0 reuse the unique z>0 with pi_x=z_x(Kz)_x. Write

    P*_xy=K_xy z_x z_y/pi_x,
    s_x=P*_xx=z_x^2/pi_x,
    f*_x=log s_x.                                    (3)

Then T_(f*) sqrt(pi)=sqrt(pi), so rho=1 and (1) equals P*. Conversely,
if a field f has pi_f=pi, use psi=sqrt(pi) and set
z'_x=sqrt(pi_x) exp(f_x/2)/sqrt(rho). The stationary flow in (1) is
K_xy z'_x z'_y. Uniqueness of the positive scaling forces z'=z. Hence

    f_x=log s_x+log rho.                             (4)

The required state field is unique up to an additive constant, which
changes rho but not P_f. It is computed FROM the target pi; (4) is not
an independent natural provider. It merely makes the missing input precise.

## 3. A single expectation really would require a scalar field

For a fixed nonconstant E let H(u) be the largest stationary history
entropy subject to nu(E)=u and legal support. Feasible stationary flow
matrices form a compact convex set; conditional entropy is continuous and
concave in the flow. Every such flow has a stationary Markov realization,
so the maximum is attained and H is concave. Self loops realize every
u in [min E,max E], including the point-mass endpoints.

For any interior u, a finite concave function has a finite supergradient
ell. The inequality H(v)<=H(u)+ell(v-u) shows that every constrained
maximizer also maximizes h+theta nu(E), theta=-ell. By (2) its occupation
law must be pi_(theta E). Thus the target Gibbs law can solve this
one-mean problem only if its inverse field (3) is affine in E.

This concerns a stationary EXPECTATION constraint. An exact finite-cycle
scalar energy equality can have additional arithmetic content for irrational
energies; no finite microcanonical energy-conditioning theorem is asserted.

## 4. Smooth entropy gradient, used in the explicit alternative

Let h*(pi) be the optimum with the full positive pi supplied. The coercive
dual already proved in the HYP review makes z(pi) smooth: in log coordinates
its derivative is diag(pi)(I+P*), with positive definite symmetric Hessian.
For any mass-zero tangent v, differentiation of sum pi=sum_xy K_xy z_x z_y
gives sum_x pi_x (log z_x)'=0. Differentiating

    h*(pi)=sum pi log pi-2sum pi log z

therefore yields

    d h*(pi)[v]=-sum_x v_x f*_x,
    (I+P*) (log z)'=v/pi,
    (f*)'=2(I+P*)^(-1)(v/pi)-v/pi.                   (5)

If every self probability is at least 1/B, write
P*=B^(-1)I+(1-B^(-1))R with R stochastic. A convergent alternating
Neumann series gives ||(I+P*)^(-1)||_infinity<=B/2. Consequently along
a straight line of positive occupations,

    |h*''|<=||v||_1 (B+1)||v/pi||_infinity.           (6)

These are finite sensitivity bounds, not postulated mobility or noise laws.
