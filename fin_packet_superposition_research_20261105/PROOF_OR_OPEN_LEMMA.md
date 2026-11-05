# PS: coherent packet superposition and cubic response stability

Date: 2026-11-05. Base: 6fb252e1c84cc2f1b1251c91a0ab2e9dc61ea53b.
Conditional mathematical comparison theorem using TP, PC and EF. The
argument is self-audited, not formally verified or independently reviewed.

## 1. A specified class, with uniform constants

Fix 0<a_-<=a_+<infinity, S>0 and one common m^2=4+w in the
original band. Let R>=max(1,1/a_-), T=S R^2. An atom is

    p_(a,theta)(t)=Re[e^(i theta) Wcal exp(-it Omega_m) h_(aR)],
    a in [a_-,a_+], theta in [0,2pi),
    h_b(x)=b^-1/2 pi^-1/4 exp[-x^2/(2b^2)].             (PS1)

All atoms solve the same exact-potential linear equation, are even and
continuous, and have finite energy. A nonzero theta generally gives
nonzero initial velocity; these are comparison data only. No incoming
field is inserted into the actual RP preparation.

Define p=sum_j alpha_j p_j with real coefficients and
M=sum_j |alpha_j|<infinity. This is either a finite sum or an absolutely
convergent series. Constants below depend on a_-,a_+,S and the fixed mass
band but not on R, the number of atoms, phases or coefficients. In
particular a_- cannot silently tend to zero as epsilon tends to zero.
The estimate is stated for any such representation; no claim that energy
alone bounds M is made. A signed measure variant is not needed here.

## 2. Uniform atom estimates and their bilinear polarization

Put tau=t/R^2. TP's local expansion becomes

    p_j=R^-1/2 r0 Re[e^-imt B_j(tau)]+e_j,
    B_j(tau)=a_j^-1/2 e^(i theta_j) A(tau/a_j^2),
    ||e_j||H1_-6<=C R^-5/2.                           (PS2)

B_j and its first tau derivative are uniformly bounded. All PC3 bounds
hold uniformly with R in place of the atom width: free half-wave H1 norm
O(1), A_1 norm O(R^-1/2), and positive-weight potential forcing
O(R^-1/2). To see this, repeat PC's rescaled Gaussian formula with width
aR and rescaled time t/(aR)^2<=S/a_-^2. The compact interval of a gives
uniform Gaussian moments and Schwartz bounds. Multiplication by a unit
complex phase and taking a real part preserve these estimates. The
velocity estimates follow from the same Fourier formula, not an imposed
zero-velocity condition.

For any two atoms define the symmetric bilinear zero-data responses

    Y[j,k]''+H Y[j,k]=-6 P_o(K p_j p_k),
    c[j,k]''+3 c[j,k]=-6<g,K p_j p_k>.                (PS3)

The PC normal form is bilinear in its two inputs. Replace its diagonal
correction by the symmetric polarization
Q[j,k]=(Q[p_j+p_k]-Q[p_j]-Q[p_k])/2, and similarly for F.
Equivalently keep both input orderings in PC4--PC6. Each bilinear bound
PC8--PC10 applies with one norm from each atom; nothing requires equal
widths or equal phases. Thus

    ||Q[j,k]||X0<=C R^-1/2, ||F[j,k]||X6<=C R^-1.

The exact projected identity retains -U_H,c(t)Q[j,k](0) and the full
causal integral. EF's dual odd local energy argument gives

    sup_[0,T] ||(Y[j,k],Y[j,k]')||X0<=C.              (PS4)

This is a bilinear proof, not an inference from checking diagonal atoms.
The continuous projection is the physical phase projection from CN.

For the scalar response, (PS2) yields a source with carriers 0,+/-2m,
slow envelopes B_j B_k and B_j conjugate(B_k), amplitude O(R^-1),
and scalar error O(R^-3). Pairing the local error with g is permitted
because g is exponentially localized. At the smallest mass,
4m^2-3>=18. Thus all these carriers are uniformly separated from the
shape oscillator frequency sqrt(3). EF's integration by parts for the
C1 envelopes, including switch-on and error accumulation, proves

    sup_[0,T] (|c[j,k]|+|c[j,k]'|)<=C R^-1.          (PS5)

Common mass is essential here. Packets from different operators are not
covered by this proof. A spatial translation outside this even family is
also not covered. Products have no discarded beat-frequency term.

## 3. Keep every cubic cross term

Let E_L(F) denote the zero-data retarded even continuous response to F.
Define ordered trilinear maps

    C[i;j,k]=E_L(-12 P_e(K p_i Y[j,k])),
    B[i;j,k]=E_L(-12 P_e(K p_i c[j,k] g)),
    D[i,j,k]=E_L(-P_e(p_i p_j p_k)).                 (PS6)

C and B are symmetric only in the last two slots; no false complete
symmetry is required. Write ||Z||Q=||q(Wcal*Z)||L2(0,T).
The first slot has sup norm O(R^-1/2). By PS4 the C source is O(R^-1/2)
in global L2. CT gives ||C[i;j,k]||Q<=C R^2. For B, PS5 and localization
of g give source L2_6 size O(R^-3/2), so LE gives ||B||Q<=C R^1/2.
For D use two atom sup norms and one atom L2 norm to obtain source
L2 size O(R^-1); CT gives ||D||Q<=C R^3/2. Therefore

    ||C[i;j,k]||Q<=C R^2,
    ||B[i;j,k]||Q<=C R^1/2,
    ||D[i,j,k]||Q<=C R^3/2.                         (PS7)

For p=sum alpha_j p_j, its full cubic continuous comparison response is

    Z[p]=sum_(i,j,k) alpha_i alpha_j alpha_k
                   (C[i;j,k]+B[i;j,k]+g4 D[i,j,k]). (PS8)

The quadratic response is sum_(j,k) alpha_j alpha_k Y[j,k], and likewise
for c. Ordered sums automatically include factors of two in distinct
quadratic indices and all cubic cross terms. Summing only alpha_j^3
Z[p_j] does not compute (PS8). By absolute convergence and the triangle
inequality,

    ||Z[p]||Q<=C M^3(R^2+R^1/2+|g4|R^3/2).          (PS9)

For countable sums, first take finite truncations. Inputs converge in
continuous-in-time energy on the fixed interval, quadratic responses in
energy, and cubic responses in Q; the estimates make them Cauchy.
Standard bounded multiplication H1 x H1 -> L2 in one dimension and
finite-time energy Duhamel continuity identify the limit with the
retarded response of the limiting source. This also proves the result
independent of a chosen representation, though the upper bound uses M.

At input amplitude epsilon, R=epsilon^-1 and fixed M,g4, (PS9) gives

    integral_0^(S epsilon^-2)|q(Wcal*epsilon^3 Z[p])|^2dt
                                                    <=C_(S,g4) M^6 epsilon^2.

This concerns the cubic response only, not the order-epsilon input trace,
which TP already showed can have squared time integral of order epsilon.

## 4. Quantified stability in the same atomic class

For p and p_tilde use a common union of their atoms. Let their coefficient
masses be at most M and let delta=sum_j |alpha_j-beta_j|. Trilinearity and

    T(p,p,p)-T(v,v,v)
      =T(p-v,p,p)+T(v,p-v,p)+T(v,v,p-v)

give, separately for the three maps, a bound by 3 M^2 delta times the
corresponding power in PS7. Hence

    ||Z[p]-Z[p_tilde]||Q<=C_(g4) M^2 delta R^2,      (PS10)

and the physically scaled squared trace difference is bounded by
C_(g4) M^4 delta^2 epsilon^2 at R=epsilon^-1. The statement holds for
any common representation, including appending the atoms of a remainder.
It is useful for truncations whose omitted coefficient mass tends to zero.
It does not control arbitrary small-energy remainders, new mass bands,
uneven translations, or nonlinear forced histories. Energy closeness
alone does not supply the weighted Fourier and slow-envelope conditions
used above; no implication from energy closeness to small delta is claimed.

## 5. Remaining research obligation

The previous single-packet estimate now covers an absolutely summable,
coherent, common-scale near-threshold family and its cubic cross terms.
It remains necessary to construct such a representation, with uniform M,
for an appropriate component of the actual generated continuum, or to
control the component outside this class in a propagator-admissible norm.
The actual continuum is forced and starts at zero; it is not a free sum
of the atoms PS1. The next step is a causal packet representation/interface
for forced histories, explicitly paying the source-time coefficient norm
and identifying which available source estimates do or do not bound it.
No matching theorem, full nonlinear approximation, modal action-error
closure, empirical FIN validation or new physical law follows here.
