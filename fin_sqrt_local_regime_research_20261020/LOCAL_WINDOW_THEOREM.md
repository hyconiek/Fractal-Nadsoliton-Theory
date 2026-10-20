# SR-002: the actual SQRT parent has a conditional local window

Date: 2026-10-20. This theorem concerns Q, not native FIN and not a
reflected physical process. All MT dependencies are specified in
[CONTRACT.md](CONTRACT.md) and the actual input ledger.

## 1. Exit probability from the stipulated iid diagnostic

Let tau_D=0 on initial outside states and otherwise be the actual first
exit. PSD A and eventually g_N>=0 imply mu0<=W pointwise. Detailed balance
from SR-001 makes W invariant for exp(tQ). Positivity of the finite Markov
kernel gives, as finite measures,

    mu0 exp(tQ)<=W exp(tQ)=W.                         (1)

There is no division by the exponentially large full Gibbs partition
function. Put delta_N=mu0(D^c)=o(1). MT's true departure layer B_N obeys
F_g(n/N)>=b there. Method-of-types bounds give W(n)<=exp(-NF_g(n/N))
and #types<=(N+1)^q. The outward-jump compensator and the actual rate cap
therefore yield, for every finite T>=0,

    Pr_iid(tau_D<=T)
      <=delta_N+RN T (N+1)^q exp(-Nb).                (2)

The first term includes all initial exceptions. No exit is removed or
postselected. For each fixed 0<zeta<b, residence through exp(zeta N) has
probability tending to one. This is one-sided trapping, not a mean escape
law, matching barrier, sign lifetime or physical duration.

## 2. Initial density and polynomial local relaxation

Let nu=mu0(.|D), m_D=mu0(D). The actual density is

    f=dnu/dpi_D=[W(D)/m_D]exp(-V_N).

Since V_N>=0, f<=B_N^0=W(D)/m_D. MT coercivity F>=-C_Lambda/N and type
bounds imply W(D)<=exp(C_Lambda)(N+1)^q. Also integral f dpi_D=1, so

    ||f-1||_2^2=integral f^2 dpi_D-1<=B_N^0-1,
    ||f-1||_2<=C(N+1)^(q/2)                         (3)

eventually. This avoids inferring an L2 bound from TV closeness. Reversibility
and (2) of SR-001 give for the reflected comparison

    TV(nu exp(tQ_D),pi_D)
      <=(1/2)sqrt(B_N^0-1)
           exp[-a_* gamma_* t/(N+1)^K].              (4)

## 3. Couple before actual exit and retain the exception only once

Couple a reflected chain started from nu with the parent started from
mu0. On the event of an initial inside state their initial values agree;
they use identical internal clocks until the parent exits. On the initial
outside event, draw the reflected initial state from nu separately. Its
marginal is still nu. Thus the coupling failure is at most

    delta_N+Pr(initial inside and exit by t)=Pr(tau_D<=t).

Combining (2) and (4), for all t>=0,

    TV(Law_iid(n(t)),pi_D)
      <=delta_N+RN t (N+1)^q exp(-Nb)
          +(1/2)sqrt(B_N^0-1)
             exp[-a_* gamma_* t/(N+1)^K].             (5)

The measure pi_D is embedded in the full count space with zero mass outside
D. The older MT bound allowed a redundant second delta_N; this sharpened
coupling does not change its validity. In particular the unpostselected
parent obeys

    sup_(t_loc(N)<=t<=exp(zeta N))
      TV(Law_iid(n(t)),pi_D,N)->0,                    (6)
    t_loc(N)=1+(N+1)^(K+2), 0<zeta<b.

The constants are existential and may be enormous. No useful finite N0,
numerical burn-in or probability at an illustrative small N is supplied.
The same pi_D gives the inherited joint static quartic/Gaussian limit and
its stationary polynomial moments; no unbounded iid-start moment is
transferred through an unspecified TV error. No labelled-history mixing,
critical trajectory theorem, global mixing upper bound or source selection
is inferred from (6).
