# Exact chiral feedback sector: accepted results and remaining gate

ST8651, 9 September 2026. Research checkpoint, not the final discovery report.
This is mathematics of the **supplied mixed-density projected learning law**;
it is not a derivation of that law, its initial data, or physical observables.

## 1. Model and exact embedding

Let Π delete the diagonal. Use

    rho_dot = i[K,rho],   K_dot = eta (Π Re rho - gamma K),

with positive eta,gamma, real symmetric zero-diagonal K and a density rho.
Complex conjugation maps this convention to the previously audited minus-i
law, leaving the real kernel path unchanged. No force has been added.

W is the supplied twelve-site real symmetric circulant strict kernel.
Write its doublet eigenvalues as lambda_k and c_k=1/12+gamma lambda_k.
For i,j=0,...,11, set E_ij=cos(pi(i+j)/6)/3 and
F_ij=sin(pi(i+j)/6)/3 when i+j is odd, and zero otherwise.
Set G=[E,F]/(2i). Exact arithmetic over Q(sqrt(3),i) verifies:

    [E,F]=2iG, [F,G]=2iE, [G,E]=2iF,
    E²=F²=G²=P,  P²=P, Tr P=4,
    {E,F}={F,G}={G,E}=0.

All three have zero diagonal and row sums. E,F are real symmetric;
G is purely imaginary Hermitian. They commute with **every** real symmetric
C12 circulant W, not only its rounded strict instance. P=P_1+P_5.
Equivalently these are two Pauli blocks with opposite handed embeddings.

Let R>0, a=R/gamma, kappa=eta gamma²/R and tau=a t. The ansatz

    K=W+a(x E+y F),
    rho=I/12+gamma W+R(u E+v F+z G)

is exactly invariant in the full twelve-site equations, and reduces them to

    x'=kappa(u-x),  y'=kappa(v-y),
    u'=-2yz,        v'=2xz,        z'=-2(xv-yu).                 (1)

Primes denote tau. The algebra above, Pi E=E, Pi F=F and Re G=0
prove the assertion by direct substitution. Nothing is projected onto this
sector during evolution; only the initial condition is restricted to it.

## 2. Admissibility is paid separately

Assume u²+v²+z²=1, |w(0)|<=1 where w=x+iy, and put b=u+iv.
Then |b|<=1 is invariant, and variation of constants in
w'=kappa(b-w) gives |w|<=1 for all positive time.

The two affected density blocks have eigenvalues c_1 +/- R and c_5 +/- R.
All other eigenvalues are unchanged. Thus c_j>=0 in the other sectors and
R<=min(c_1,c_5) suffice for rho>=0. This condition is independent of seed
size and phase. The vertex populations remain exactly 1/12.

Only odd-distance edges of K change, by at most a/3. Thus

    a/3 < min(W_1,W_3,W_5)

preserves strict positive edges. All rows retain their original sum s.
This is a forward-invariant family for the original law, not a general
positivity theorem for that law. The prior counterexample to unrestricted
edge-positivity preservation remains valid.

For gamma=1/20, R=1/2000, a=1/100 and eta=1/5, kappa=1.
The inherited exact strict weight/eigenvalue enclosures can pay these strict
inequalities; numerical margins alone are not used as a universal proof.

## 3. Global convergence from a nonzero real-coherence seed

**Theorem.** For every kappa>0, 0<epsilon<=1 and phi in R, the solution
of (1) with w(0)=0, b(0)=epsilon exp(i phi),
z(0)=+sqrt(1-epsilon²) converges to exactly one point

    w=b=exp(i Theta), z=0.                                    (2)

The same statement holds with the negative initial z. Here convergence
to one point does not mean a common or canonical point for all inputs.

**Proof.** The spin sphere and |w|<=1 give a compact invariant set, so the
polynomial vector field exists for all positive time. Define

    V=|w|²-2 Re(conj(w)b),
    V+1=|w-b|²+z²,   V'=-2 kappa |w-b|².                     (3)

The largest invariant subset of V'=0 consists of the two poles w=b=0,
z=+/-1, and the equilibrium circle (2): w=b implies w'=0, while
b'=2izw forces zw=0.

Let q=|w|². Direct differentiation gives

    q''+kappa q'=2|w'|²,   q'(0)=0.                          (4)

Since w'(0)=kappa epsilon exp(i phi) is nonzero, the integrating-factor
formula proves q'(tau)>0 for every tau>0. In particular its limiting
value is positive. LaSalle's invariance argument applied to (3) excludes
the poles and gives |w|->1, w-b->0 and z->0.

It remains to exclude wandering around the circle. At any tau_0>0 choose
a continuous lift theta of arg w. From (1),

    theta'=kappa Im(conj(w)b)/q=-kappa z'/(2q).

Integration by parts yields

    theta(T)-theta(tau_0)
      = -kappa/2 [z/q]_(tau_0)^T
        -kappa/2 integral_(tau_0)^T z q'/q² d tau.            (5)

The boundary term converges. The integral is absolutely convergent because
|z|<=1, q'>=0 and integral q'/q²=1/q(tau_0)-1/q(T).
Thus theta has a finite limit and the entire state tends to (2). QED.

This proof rules out a periodic orbit or chaotic limiting attractor for
the stated seed family in this exact sector. It does not classify arbitrary
initial conditions in the full twelve-site system.

## 4. Exact instability and failure of continuous selection

The zero-seed chiral state w=b=0,z=1 is itself stationary. On the spin leaf
its linearized complex eigenvalues obey

    mu²+kappa mu-2i kappa=0.

The unstable root mu_u=alpha+i beta=(-kappa+sqrt(kappa²+8i kappa))/2
has alpha>0: the real part of the square root is strictly larger than
kappa. There is one unstable complex pair and one stable pair
mu_s=-kappa-mu_u. The opposite pole reverses the rotation direction.

More strongly, the theorem above gives an exact nonlinear instability:
arbitrarily small nonzero seeds approach |w|=1, while the zero seed stays
at |w|=0. The final kernel displacement has Frobenius norm 2R/gamma,
independent of epsilon. The time to see it diverges as epsilon->0; there
is no discontinuity of finite-time ODE solutions.

Equation (1) is equivariant under (w,b)->exp(i phi)(w,b). Therefore its
endpoint satisfies Theta(epsilon,phi)=Theta(epsilon,0)+phi modulo 2pi.
For every desired endpoint on (2) and every epsilon>0, a seed phase gives
that endpoint. Taking epsilon->0 proves:

**No-continuous-selector corollary.** The cluster set of the endpoint map
at the chiral pole is the **whole** equilibrium circle. No continuous
extension, and no unique geometry independent of vanishing seed phase,
exists there. This uses global convergence and exact symmetry, not a fit
to a logarithmic spiral or an unproved generic instability principle.

All these initial conditions have the same density spectrum, same kernel
W, same uniform populations and same initial learning energy. The selected
angle is coherence data, not identifiable from those quantities.

The circle is normally attracting **inside the four-dimensional spin leaf**:
at w=b=1,z=0 its normal eigenvalues are -kappa and the roots of
lambda²+kappa lambda+4=0; the fourth eigenvalue is the tangent zero.
This is not a full-system stability certificate.

## 5. A proved local logarithmic law, and a separate unproved extension

On the two-real-dimensional unstable manifold of the positive pole, choose
an equivariant smooth complex coordinate zeta tangent to the unstable
eigenspace. The S1 action is scalar rotation, so its vector field is

    zeta'=zeta (alpha+i beta+f(|zeta|²)),  f(0)=0.

Smoothness and symmetry give f=O(|zeta|²). For sufficiently small radius,
outward radial velocity is positive. The unstable manifold has w=C zeta
+O(|zeta|³), C!=0. Formula (4), integrated from negative infinity, gives
q'>0 on every nonstationary outgoing trajectory; the global convergence
proof above therefore also applies to those trajectories.

Let h(r) be the endpoint phase for the initial unstable coordinate r>0.
The endpoint map is smooth away from the pole by the normally attracting
circle and finite-time smooth dependence. It is constant along an orbit
and equivariant under phase rotations. Hence

    h'(r)=-(beta+Im f(r²))/(r(alpha+Re f(r²)))
          =-beta/(alpha r)+O(r).

Integration proves h(r)=-(beta/alpha) log r+C_0+O(r²).
This is ordinary saddle-focus logarithmic winding, not evidence for a
new universal scaling law. The asymptotic scale ratio exp(2pi alpha/|beta|)
depends on the supplied kappa.

For the physically simpler curve **K(0)=W exactly**, b(0)=epsilon>0,
the same formula is strongly supported by the saved ODE runs, but the
incoming stable-component passage estimate has not yet been proved here.
Do not transfer the unstable-manifold theorem silently to that curve.
The accepted no-continuous-selector corollary does not need this extension.

## 6. What is converted, and what is not derived

On the stated seed family the density eigenvalues and von Neumann entropy
are constant. Vertex Shannon entropy is always log12. Nevertheless K
changes, and in the limit the two affected doublets split by +/-R/gamma.
The state purity excess over I/12+gamma W is 4R², while the final squared
kernel Frobenius norm increases by 4R²/gamma². The learning functional
F=gamma||K||²/2-Tr(K rho) decreases by 2R²/gamma.

Thus this is an exact conversion of supplied chiral/coherence structure
into real geometry within a candidate feedback law. It refutes a claim
that such feedback must continuously erase seed choices or uniquely select
strict. It does not derive the chiral resource R, its sign, the seed,
the background strict spectrum, an initial density, a clock or an apparatus.
It does not discharge QW-2191 or any legacy-to-strict/physical-role gate.

Scientific priority and significance beyond this FIN-specific construction
remain under audit. No theorem of fundamental physics or global optimizer
selection is claimed.
