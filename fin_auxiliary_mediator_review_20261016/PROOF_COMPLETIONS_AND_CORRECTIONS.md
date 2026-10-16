# HM proof audit and stationary path completion

Date:2026-10-16. Analytical model theorem under inherited exact spectral and
local-Gibbs premises, not mechanized or external peer review. Primary FIN
is unchanged; AUX-HB-OU-kappa remains a diagnostic comparison.

## 1. Reuse and source binding

The actual MP7-032 Gaussian identity, NT intake and MT proof completions are
available locally and bound in verification.json. This pays the current
dependency, not the producer's historical access to missing uploaded texts,
AGENTS snapshot or archive wrappers.

The cavity identity, reset-lift construction and formal r_kappa already
appeared in the preceding MEDIATOR_BRIDGE_SEED.md. They are reused and
specialized here, not newly discovered general laws. The main new process
result is HM-004's stationary mixed jump/diffusion critical limit, reviewed
below. Generic positive-memory and hidden-realization principles are old.

## 2. Cavity and persistent sampler

The cavity integral is correct for g>0 and constant A_jj, with the actual
parent denominator N. Full proposal followed by departure rejection gives
the native reflected rule. A conditional sampler renormalized only over
allowed labels would change it.

The persistent joint law is also correct. Globally, label refresh and the
conditional OU generator are separately reversible for Pi. Under restriction
n in D_N, the conditional label law is no longer an unrestricted product;
however pairwise detailed balance survives **rejection** of outgoing jumps.
The field conditional rho_n remains the same Gaussian.

Finite labelled/count states, bounded total label rate N and a linear OU
drift with bounded state-dependent means give a nonexplosive construction.
The OU part is a direct sum of closed Gaussian forms. The symmetric label
jump operator has L2 norm at most2N, by bounded exit rates and reversibility.
It is a bounded symmetric perturbation. This supplies a conservative
self-adjoint mixed generator rather than only a formal invariant density.

For a count observable m, Lm is a finite linear combination of bounded
smooth softmax functions. Their first/second derivatives are bounded; the
OU action grows at most linearly in theta and is Gaussian-square-integrable.
Constants in theta belong to each OU domain, and Gaussian cutoffs justify
the smooth tests. Consequently m is in Dom(L^2); the compressed first and
second derivatives used by HM-002 are paid. The same holds with rejection,
where the finite count-dependent indicators introduce no theta singularity.

## 3. Correct reflected memory witness

Globally, with P=E[.|n],

    ||QLPm||^2 = sum_n pi(n) Var_rho_n(phi(theta)),
    phi(theta)=sigma.softmax(X theta).

Gaussian full support and phi(t e6)=tanh(sqrt(Lambda6/12)t) make this strictly
positive. For reflection write b_D(n,theta) for the omitted parity drift:

    L_D m=phi(theta)-m-b_D(n,theta), |b_D|<=2 1_BD.

The exact full-domain reflected identity is therefore

    ||Q L_D Pm||^2
       =sum_n pi_D(n) Var_rho_n(phi(theta)-b_D(n,theta)),

**not** the uncorrected average of Var(phi) on the entire D_N. At an interior
type n0 all neighbours are retained, so it is at least

    pi_D(n0) Var_rho_n0(phi)>0.

The balanced-uniform type and all its one-step neighbours belong to the
fixed MT local domain for all sufficiently large N in12N: their free-energy
and count displacements vanish while the centre has a fixed sublevel buffer.
Thus the source's local nonsemigroup conclusion survives. Its short witness
formula needs the boundary correction above when interpreted globally on D_N.

The positive identity P L^2 P-(PLP)^2=(QLP)^*(QLP) then proves that the
stationary count compression is not a semigroup. This is stronger than a
strong-lumpability counterexample for one nonstationary preparation. It
does not exclude every mediator or an asymptotic approximation.

The producer's two-state compression fixture has PL=0, so its displayed
defect is0=0. Its True flags for the full-posterior caveat/reset witness are
descriptive, not computational verification. The intake replaces those
diagnostic gaps with a nonzero exact four-state defect and finite reset
intertwining/covariance fixtures. The analytical Gaussian proof remains
the evidence for the actual continuous-state claim.

## 4. MT moments and the actual departure layer pay HM-004's inputs

MT's admitted polynomial moment domination gives, under pi_D,N,

    E|m|^q=O(N^(-q/4)), E||w||^q=O(N^(-q/2))

for every fixed positive integer q. The conditional Gaussian law supplies
the field moments used by HM, including

    E||xi||^q=O(N^(-q/2)), E|h-m|^q=O(N^(-q/2)).

Here xi is X theta minus its parity-line part. The actual one-step departure
layer has pi_D,N(B_D)<=poly(N)exp(-bN) from the inherited positive boundary
free energy and polynomial local partition normalization. This is the
**upper** bound paid in the MT completion, not the original MT-003 sentence
with its reversed inequality. In particular N^(3/4)pi_D,N(B_D)->0.

These bounds concern the unchanged equilibrium measure. They are reusable
for the **stationary diagnostic** because its count marginal and rho_n are
exactly these measures. The old primary-model iid mixing/trapping theorem
is not transferred to its changed generator.

## 5. Complete tightness and the joint martingale identification

The source's drift, bracket and initial-law calculations for
s=(kappa m+h)/(1+kappa), Z_N(t)=N^(1/4)s(sqrt(N)t) are correct. Fix kappa1
or4 and bounded fixed lambda. In the reflected stationary diagnostic,

    Z_N(t)=Z_N(0)+integral_0^t b_kappa(Z_N(u))du+M_N(t)+R_N(t),
    b_kappa(z)=r_kappa(chi lambda z-z^3/3),
    E sup_(t<=S)|R_N(t)| ->0.

The last assertion follows by bounding R's total variation with the source's
stationary L1 ledger. Its four polynomial powers are negative, and the actual
departure layer is exponentially small. Taylor bounds are global: softmax
derivatives are bounded, and tanh's polynomial remainder is globally dominated.
No uniform-in-theta finite tube is silently assumed.

The critically scaled predictable bracket density has a deterministic upper
bound C: |m|,|phi|<=1 and the h-noise coefficient is constant at fixed N.
Reflection only deletes jump-bracket terms. Jumps are O(N^(-3/4)); the
remaining noise is continuous. Thus E sup|M_N|^2<=C_S. Ito's fourth-moment
formula/BDG, including the positive fourth-jump compensator, gives
E sup|M_N|^4<=C_S. The martingale and its square-minus-bracket have the
uniform integrability needed below.

Stationary moments at fixed times alone would not establish maximum tightness.
Here use the exact equation. For w=M_N+R_N, the function Y=Z_N-w is absolutely
continuous and obeys Y'=b_kappa(Y+w). On an event sup|w|<=K, the confining
drift points inwards outside a fixed radius depending on K and lambda.
A deterministic comparison bounds sup|Y| by its initial magnitude and that
radius. Tight Z_N(0), BDG and sup|R_N|->0 therefore give compact containment.
After localization the main drift is bounded. Its increments, the vanishing
residual and bounded-bracket martingale imply Aldous tightness; vanishing
jumps make every limit continuous.

Pass the martingale and square identities jointly with the initial variable
in the limiting natural filtration. The bracket converges to2r_kappa t by
the paid stationary L1 calculation. The continuous martingale is a Brownian
driver with that bracket, with future increments independent of Z0. The
quartic SDE is nonexplosive and uniquely solvable. This identifies the whole
sequence on each fixed scaled horizon, not merely generator balance:

    dZ=r_kappa(chi lambda Z-Z^3/3)dt+sqrt(2r_kappa)dB,
    Z0~rho_lambda, r_kappa=kappa/(1+kappa).

Admit this stationary reflected diagnostic theorem under the inherited
exact model/spectral/MT premises. No numerical N0, kappa_N limit, arbitrary
preparation, physical seconds or global capacity theorem follows.

## 6. Additional corollary: the actual count parity has the same path limit

The source proves its theorem for s; fixed-time closeness to m alone would
not transfer an entire path. The actual exponential boundary bound permits
this completion. Put d=h-m and c_N=g_N/g_c. Off the departure layer,

    Ld=-(1+kappa)d+F_N,
    F_N=kappa(c_N-1)m+h-phi(theta).

On the layer add a forcing of magnitude<=2. From the stationary high moments,
E|d|^4=O(N^-2), E|F_N|^4=O(N^-3), with added C pi_D(B_D) for reflection.
The unscaled d-martingale has bracket density<=C/N and jumps<=2/N. Variation
of constants and fourth-moment BDG on **one model-time unit** yield

    E sup_(j<=u<=j+1)|d(u)|^4 <=C(N^-2+pi_D(B_D)).

Stationarity makes this uniform in integer j. A union bound over O(sqrt(N))
intervals and Markov's inequality give, for fixed epsilon>0,

    P(sup_(u<=sqrt(N)S) N^(1/4)|h(u)-m(u)|>epsilon)
       <= C_S epsilon^-4 [N^-1/2+N^(3/2)pi_D(B_D)] ->0.

Hence N^(1/4)m(sqrt(N)t) and Z_N have the same continuous path limit.
This extra corollary uses the actual stronger MT exponential boundary
payment, not merely the source's weaker N^(3/4) boundary condition.

## 7. Clock language, exact lifts and physical limits

Source COEFFICIENT_NEGATIVE_CONTROLS.json says same_native_clock=false.
The actual **label attempt convention is1 in every model**. What differs
is collective relaxation/mobility, and the extra hidden rate kappa is supplied.
Do not silently change clock units by fitting r_kappa after the comparison.

The exact reset-lift proof is correct. Sum only over n'!=n. Its no-count-jump
residual covariance formula is valid; contributions with a count jump vanish
after the fresh centred independent residual draw, so it also equals the
full covariance in stationary initialization. Gaussian moments justify it;
a bounded truncation can be used if only bounded readouts are licensed.
The hidden reset eta is a separate analysis rate, not the strict kernel's
compression exponent eta=9/5.

This distinguishes invisible eta in exact-native lifts from kappa in feedback
models whose count laws change. It does not imply that kappa is unidentifiable
**within a frozen feedback family**. Neither a mathematically identifiable
parameter nor a chosen sampler supplies nature's hidden kinetic law.
The next plan addresses visible feedback versus observation-induced memory,
without adopting a new physical field or claiming empirical evidence.
