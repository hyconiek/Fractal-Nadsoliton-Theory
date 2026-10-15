# A bounded new proof tool: removing one tag's feedback in path TV

Date: 2026-10-15. Status: analytical model lemma. This is an auxiliary
reference construction, **not a replacement microscopic FIN law**.
It pays one obstruction in the next likelihood proof, not full sufficiency.

## 1. Two processes with identical initial laws

Distinguish one tagged slot from its N-1 background slots, and let n_bg be
their count vector. In the actual parent P:

- the tag refresh distribution is `softmax[(g/N) A n_bg]`;
- background slot i refreshes using `softmax[(g/N) A(n_bg-e_i+e_tag)]`.

Define an auxiliary comparison R using the **same initial joint law**, gain,
per-slot attempt clock and denominator N. Leave the tag rate unchanged, but
let each background slot use `softmax[(g/N) A(n_bg-e_i)]`. Its background is
autonomous after time zero; its initial state is not yet assumed independent
of the tag. No denominator N-1 renormalization is introduced.

## 2. A global categorical entropy bound

For P_h=softmax(h), F(h)=log sum exp(h_j) and any delta,

    D(P_(h+delta) || P_h)
      =integral_0^1 t Var_(P_(h+t delta))(delta_j) dt
      <= (max delta-min delta)^2/8.

The Hessian of F is the categorical covariance; a variable in an interval
of length d has variance <=d^2/4. Integrating t gives the displayed factor.
For delta=(g/N)A[:,tag], put M_A=max|A_ij|. Thus each changed refresh pays
at most `g^2 M_A^2/(2N^2)`, uniformly over all configurations.

Construct both CTMCs using the same standard rate-one attempt representation,
including null self refreshes as *unobserved analysis marks*. The chain rule
for conditional mark likelihoods and compensation of the Poisson attempts
give, on a model-time interval [0,T],

    D(P_marked || R_marked)
      <= (N-1) T g^2 M_A^2/(2N^2).

This can also be obtained by expanding the product of conditional mark
probabilities between the common Poisson times. Both initial laws are equal,
so there is no unpaid initial entropy term. Forgetting marks contracts
relative entropy. Pinsker therefore gives for the **complete unmarked paths**

    TV(P_path,R_path) <= |g| M_A sqrt((N-1)T)/(2N).

At bounded g and T=sqrt(N)S, this is O(N^(-1/4)). It is global in the finite
state space; it needs no practical MT N0, numerical spectral gap or rare-path
truncation. Null marks are not licensed additional observational data.

## 3. Immediate useful consequences and unpaid obligations

Every bounded target/record decision experiment is changed by at most this
path-TV bound. In particular the optimal 0-1 Bayes risks differ by at most
it, since every candidate decision rule has bounded loss. Take T through the
**future target time**, not only through the observation time.

Under R the pair (background, tag parity) is Markov: the background is
autonomous and the tag flips at `(1-S_tag q(n_bg))/2`, where
`q(n_bg)=sigma^T softmax[(g/N)A n_bg]`. The tag's discarded full q12 label no
longer affects background rates. This does not justify replacing a parity
observer by a full-label or attempt-mark observer.

Still unpaid: replacing correlated local-Gibbs initial data by a suitable
product reference; the fair-tag change of measure; the joint stochastic
integral/log-likelihood limit; and approximation of conditional likelihoods
by finite observed cylinders. Background autonomy alone is not initial
independence, and the likelihood comparison alone is not statistical
sufficiency. Those are precisely the next proof gates.

The lemma is a proof tool within the existing mean-field model, not physical
locality, an independently sourced detector or a change in FIN ontology.
