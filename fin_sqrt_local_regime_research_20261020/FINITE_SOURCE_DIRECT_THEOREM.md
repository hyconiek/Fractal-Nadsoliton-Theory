# SR-003: direct local window for a finite entropy source

Date: 2026-10-20. This is a direct discrete-time theorem in the already
authorized hypothesis family. It does not use full-path TV convergence to
SQRT, which SL correctly disproves for ideal continuous timestamps.

## 1. Exact source, not a Poissonized replacement

For each visible N let pi be its full Gibbs count law. SC/RT give the actual
visible count source from N visible and L uniform hidden labelled slots:
reference K+cI, c=(q-1)L, K_nn=1 and
K_(n,n-e_a+e_b)=sqrt(n_a(n_b+1)). Write its positive scaling as z and

    pi_x=z_x[(Kz)_x+c z_x],
    P_c(x,y)=(K_xy+c 1_(x=y))z_x z_y/pi_x.

The visible source is exactly strongly lumpable, hence Markov for every
fine preparation, not merely stationary coarse averages. Its supplied
fine stationary law is visible labelled Gibbs times uniform hidden law.
Uniform iid visible/hidden initialization is an additional diagnostic of
this transition kernel, not the stationary entropy-maximizing history.

Set s_x=sqrt(pi_x), r_x=sqrt(c)z_x/s_x and

    A_N=max_x (Ks)_x/s_x=1+max_x sum_(y!=x)Q(x,y)
         <=1+RN=:Abar_N.

SL's exact positive scaling estimate gives 1-A_N/c<=r_x<=1. Therefore,
when c>=2 Abar_N,

    1/2<=r_x<=1,
    c P_c(x,y)=Q(x,y)r_x r_y, x!=y,                 (1)
    (a_*/4) B(x,y)<=c P_c(x,y)<=a^* B(x,y).

These statements use the true scaling, not replacing it by sqrt(pi/c).
No numerical solution or large lattice enumeration is required.

## 2. Reflected discrete mixing and its negative-eigenvalue gate

Define P_D by retaining every internal P_c edge and adding each forbidden
outgoing probability to its self entry. It is reversible at pi_D. From
(1), its associated form satisfies

    gap(c(P_D-I))>=a_* gamma_*(N+1)^(-K)/4=:g_N^0.

For the unreflected kernel,
1-P_c(x,x)<=A_N/c<=1/2; reflection only increases its diagonal.
Hence P_D=(I+S_D)/2 for a stochastic reversible S_D. All eigenvalues of
P_D are nonnegative. This explicit laziness payment is needed: a continuous
Dirichlet gap alone would not exclude near-minus-one discrete eigenvalues.
Its nonconstant spectral contraction is consequently at most

    (1-g_N^0/c)^k<=exp(-g_N^0 k/c).

Using exactly the density B_N^0 from SR-002,

    TV(nu P_D^k,pi_D)
      <=(1/2)sqrt(B_N^0-1) exp(-g_N^0 k/c).          (2)

This proves mixing of the count projection only, not of every fine hidden
coordinate or an arbitrary tagged record.

## 3. Actual finite-source exit and local window

Detailed balance gives W P_c=W. Thus mu0 P_c^j<=W for every integer j.
The probability of an outward step is at most RN/c from an inside state.
Summing the actual discrete outward-event probabilities, including the
initial exception, gives

    Pr_iid(tau_D^step<=k)
      <=delta_N+(RN k/c)(N+1)^q exp(-Nb).            (3)

There is no CTMC clock or path-equivalence assumption in this inequality.
The same pre-exit coupling as SR-002 and (2) yield

    TV(Law_iid(Y_k),pi_D)
      <=delta_N+(RN k/c)(N+1)^q exp(-Nb)
          +(1/2)sqrt(B_N^0-1)exp(-g_N^0 k/c).        (4)

For X_N(t)=Y_floor(c_N t), any sequence c_N=(q-1)L_N>=2(1+RN) satisfies
the local-equilibrium window (6) of SR-002 with X_N in place of n(t).
Indeed floor(c_N t)/c_N>=t-1/c_N; the chosen t_loc has ample fixed slack.
Exit through floor(c_N exp(zeta N)) also tends to zero for zeta<b.

A concrete sufficient choice is

    L_N=ceil(2[1+RN]/(q-1)), c_N=(q-1)L_N.           (5)

It uses O(N) uniform hidden slots and polynomially many intrinsic steps for
burn-in, since c_N t_loc(N) is polynomial. No practical constant, numeric
N0 or physical interpretation of those slots is asserted. The cap and
mixing constants are uniform over every larger c_N as well.

## 4. What this resolves, and what it does not

The former double-limit warning is answered for this **local-equilibrium/
residence property** by a direct finite-source proof. No exponentially many
spectators are needed for that particular conclusion. This is not agreement
of finite-source and SQRT histories through exponential times.

For that different target the SL sampled-record bound still costs
3 A_N^2 T_N/c_N (and its J1 coupling has its own endpoint terms). With the
uniform Abar_N=1+RN, sufficient worst-case vanishing error on
T_N=exp(zeta N) would require c_N/(N^2 exp(zeta N))->infinity. O(N) does
not pay that sufficient budget; failure of an upper bound is not a lower
bound proving such huge resources necessary. Full continuous-path TV stays
inappropriate because finite-grid jumps have atomic timestamps.

Here only the internal comparison generator c(P_c-I) is used to express
forms and spectral gaps. The actual finite source remains discrete, and
the proof never silently substitutes a Poissonized process for its paths.
Its total size is N+L_N, with a supplied product static law, not an enlarged
fully interacting native FIN system. No physical autonomy/bath/clock or
strict F01 conclusion follows.
