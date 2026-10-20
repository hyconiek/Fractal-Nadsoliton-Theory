# CF-001: derive the predictive state from a specified finite source

Date: 2026-10-20. Conditional mathematical continuation of the ES source
probe, not an independently sourced FIN dynamics or a new observer law.
Filtering and lumpability are established methods; the contribution here
is their explicit, source-first use with the ES preparation and null marks.

## 1. Inputs before the target generator

Let Z be a **supplied finite CTMC** with conservative generator L and
supplied initial law mu. Let `R:S->E` be a deterministic readout and
`D_t=R(Z_t)`. These are actual input objects, not things inferred from a
desired B. The observed filtration is `F^D_t=sigma(D_s:0<=s<=t)`.
For `F_i=R^-1(i)`, use the blocks L_ij of the *original* generator and set

    a_ij(z)=sum_(z' in F_j) L_zz', z in F_i, j!=i.

No physical measurement disturbance or output noise is added. The process
D may be non-Markov. The analysis state

    eta_t(z)=Pr(Z_t=z | F^D_t)

is a conditional probability on the source state, not an extra ontic
substrate or a claim of dimension-minimal predictive storage.

## 2. Exact update, including unobserved internal transitions

Given D_0=i with positive probability,
`eta_0(z)=mu(z)1_{F_i}(z)/mu(F_i)`. If there is no **observed** change for
duration u, the unnormalized posterior row is

    beta(u)=eta_0 exp(u L_ii),
    eta(u)=beta(u)/(beta(u)1).                         (1)

Here L_ii is the killed subgenerator: retain all internal transitions and
the original diagonal, so the loss of row mass is precisely observable
departure. It is not the conservative reflected generator, and it is not
obtained by deleting all hidden jumps. Finite total source rates imply
positive no-observed-jump survival for finite u (no source jump is already
a positive-probability subevent), so the normalizer in (1) is positive.

At an observed change i->j,

    alpha_ij(t)=eta_(t-) L_ij 1,
    eta_t=eta_(t-) L_ij / alpha_ij(t).                 (2)

If alpha_ij=0 the event has zero instantaneous density under the source;
there is no valid conditional update by division. Keep such impossible
events explicit rather than fabricating a posterior.

Proof: finite-state Markov propagation killed at first exit from F_i solves
the linear forward equation with L_ii, by conditioning on one short source
step. Its mass is the survival probability. An observed jump i->j has
first-order transition matrix L_ij du. Multiplication followed by Bayes
normalization gives (1),(2). Iterating gives the density for any finite
observed history, with no observed attempt/internal marks supplied.

Between observed jumps the normalized equation is

    eta_dot=eta L_ii + alpha_out eta,
    alpha_out=sum_(j!=i)eta L_ij 1.                   (3)

The pair `(D,eta)` is therefore an exact sufficient predictive state. Given
eta, future source laws are the mixture of the original Markov source laws
from eta. This proves predictive sufficiency, not physical implementation
or minimality. Equations (2),(3) define a derived piecewise-deterministic
filter process; its bounded jump rates and the killed-semigroup construction
pay existence on every finite horizon.

## 3. The correct projected intensity and exact closure

For the observed i->j counting process, testing conditional expectations
against bounded F^D-predictable functions gives the compensator (2), namely

    alpha_ij(t)=E[a_ij(Z_(t-)) | F^D_(t-)]
                on {D_(t-)=i}.                       (4)

This is a **history** posterior, not conditioning only on D_t=i. For
example, the stationary conditional average over a fibre may be correct
while (4) depends on the observed past.

For a target generator Q, exact observable Markov closure for a particular
source preparation is paid if `alpha_ij(t)=q_ij` almost surely dt on active
rows. The coordinate-martingale argument in ES section 8 then determines
the complete law. It is a preparation-specific condition.

The stronger pointwise property `a_ij(z)=q_ij` on every source fibre is
sufficient for every initial source law. Conversely, if the same Q is
claimed from every deterministic source start, a short-time derivative
forces this property. This is the existing strong-lumpability criterion,
reused rather than presented as a new source theorem. It is **not
necessary** for one specially prepared projected law.

## 4. Why this is progress, but not F01 closure

This derives a causal state map and actual predictable coefficients from
an already specified `(L,R,mu)`, instead of demanding that observed labels
be Markov by fiat. The later comparison theorem can evaluate these objects
without using Q to define the filter. But the independently grounded FIN
source L, its carrier, readout and preparation are still separate inputs.
Constructing L as a desired-HB dilation would not source HB. Continuous or
deterministic field providers need their own domain/readout/filter proof;
the finite source theorem is not silently applied to them.
