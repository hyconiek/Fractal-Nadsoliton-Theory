# COMPLETE-SYSTEM-CRITERION-295
## Exact relative reversible-closure criterion, and a no-go for absolute completeness from projected path data alone

Date: 2026-09-27

Status:
- exact finite-state relative closure theorem;
- exact entropy lower bounds for hidden forward randomness and hidden reverse memory;
- exact finite reversible-tape counterexample for any prescribed finite observation horizon;
- exact all-time counterexample supplied by the canonical two-sided/natural-extension construction;
- **no unconditional FIN theorem selecting the simultaneous visible records as the absolute complete carrier.**

This report attacks the central open premise left by reports 288, 291 and 294:

> When can FIN say that the selected degrees of freedom are the complete closed system, rather than a subsystem coupled to an unmodelled reservoir?

The answer has two parts.

1. Once an ambient reversible FIN carrier is independently specified, there is an exact internal criterion for whether a candidate state is self-contained.
2. From the projected path law of the candidate state alone, absolute completeness is not identifiable: a hidden reversible tape can reproduce a fresh reset bath for an arbitrarily long finite horizon, and the natural extension reproduces it for all time.

So task 295 yields a **positive relative theorem plus an absolute no-go**.

---

## 1. Setup

Let

    Omega

be an ambient state space and let every primitive event/controller label

    e in E

act by a bijection

    T_e : Omega -> Omega.

The controller/event record is conditioned on explicitly. This preserves the guardrail from reports 257, 291 and the post-255 audit: stochastic averaging must not be silently identified with one invertible microscopic event.

Let the proposed complete state be a projection

    pi : Omega -> X.

The question is whether the omitted fiber coordinates in

    pi^(-1)(x)

contain information that can influence future or past visible evolution.

---

## 2. Exact bidirectional fiber-closure criterion

For a primitive event e define two conditions.

### Forward fiber closure

For all omega,omega' in Omega,

    pi(omega)=pi(omega')

implies

    pi(T_e omega)=pi(T_e omega').

Equivalently, the next visible state depends only on the present visible state, not on which hidden point in the fiber was present.

### Backward fiber closure

For all omega,omega' in Omega,

    pi(omega)=pi(omega')

implies

    pi(T_e^(-1) omega)=pi(T_e^(-1) omega').

Equivalently, the previous visible state is also determined by the present visible state.

### Theorem 295-A — reversible factor theorem

If forward and backward fiber closure hold for e, then there exists a unique bijection

    f_e : X -> X

such that

    boxed:
    pi o T_e = f_e o pi.

Conversely, if such a bijection f_e exists, both fiber-closure conditions hold.

### Proof

Forward closure makes

    f_e(x)=pi(T_e omega)

well-defined for any omega with pi(omega)=x.

Backward closure similarly defines

    g_e(x)=pi(T_e^(-1) omega).

Then

    g_e(f_e(pi(omega)))
      =
    pi(T_e^(-1) T_e omega)
      =
    pi(omega),

and similarly

    f_e(g_e(x))=x.

Therefore

    g_e=f_e^(-1),

so f_e is bijective.

QED.

### Interpretation

A candidate state X is event-closed relative to Omega exactly when every primitive reversible ambient event descends to a reversible event on X.

If forward closure fails, omitted coordinates contain information that changes the future.

If backward closure fails, visible evolution erases information that must be stored somewhere outside X in any reversible realization.

This is distribution-free and stronger than a fitted Markov test.

---

## 3. Entropic operational version

When only the induced event kernel is available, define on the reachable support, using a full-support test distribution,

    F_e = H(X' | X,e)

and

    B_e = H(X | X',e).

Call:
- F_e the **fresh-output defect**;
- B_e the **reverse-memory defect**.

For a finite state set:

    boxed:
    F_e=B_e=0

iff the event kernel is a permutation on its support.

Thus

    D_e=F_e+B_e

is a statistical reversibility/self-containment defect.

It is not a replacement for the exact fiber test: a test distribution can miss measure-zero or unreachable fibers. But with full support on a finite reachable class it is equivalent to the permutation test.

---

## 4. Hidden-information lower bounds

Suppose the visible stochastic event is realized by a deterministic reversible dilation on

    (X,Z).

Because X' is a deterministic function of (X,Z,e),

    H(X'|X,e)
      <=
    H(Z|X,e).

Therefore

    boxed:
    H(Z|X,e) >= F_e.

So F_e is a lower bound on hidden information that must be available to create the visible branching.

By applying the same argument to the inverse dilation,

    boxed:
    H(Z'|X',e) >= B_e.

Thus B_e lower-bounds hidden output memory required to preserve information apparently erased from X.

These are exact information inequalities.

---

## 5. Reset versus SWAP

Let |S|=q.

### Full one-record reset

For

    K(x,x')=1/q,

under full-support/uniform testing,

    boxed:
    F_reset=B_reset=log_2 q.

For q=3,

    F_reset=B_reset
      =
    log_2 3
      approximately
    1.5849625007 bits/event.

So the visible record alone is not a reversible closed event state.

### Full two-record SWAP

For

    (x,y)->(y,x),

both directions are deterministic and bijective:

    boxed:
    F_SWAP=B_SWAP=0.

### One-site projection of SWAP

If the neighbor is hidden and initially uniform, the one-site projected kernel is again the full reset kernel.

Hence

    F_one-site=B_one-site=log_2 q.

This is the desired behavior of a completeness test:

- the **pair/global SWAP state** is closed;
- one site by itself is not.

Therefore the criterion diagnoses missing simultaneous degrees of freedom rather than merely labeling SWAP as intrinsically closed.

---

## 6. Recovery of the report-291 entropy rate

In the alpha reset/SWAP family, condition on primitive event identity.

A hidden reset event contributes

    F=log_2 q

fresh bits.

A visible SWAP event contributes

    F=0.

With n sites, local reset rate alpha rho and time T,

    E[M_reset]
      =
    alpha n rho T.

Therefore

    boxed:
    E[hidden fresh information]
      >=
    alpha n rho T log_2 q.

This is exactly the lower bound of report 291.

The two-sided reversibility defect additionally gives

    boxed:
    E[D]
      =
    2 alpha n rho T log_2 q

for the uniform reset channel.

Thus report 291 is recovered as the forward half of a more general reversible-closure theorem.

---

## 7. Finite reversible-tape no-go

The previous section detects that X alone is not closed.

It does **not** prove that no larger closed carrier exists.

For any alphabet S of size q and any desired finite horizon L, define the enlarged state

    (x,r_0,...,r_(L-1),p)

with pointer

    p in Z_L.

One event:

1. swaps x with r_p;
2. increments p mod L.

Explicitly,

    x' = r_p,
    r'_p = x,
    p' = p+1 mod L,

with all other tape cells unchanged.

This map is a bijection. Its inverse decrements the pointer and swaps back.

Initialize the tape symbols independently and uniformly.

Then, for a fixed initial visible x, the first L visible outputs are exactly iid uniform reset symbols.

Therefore:

    boxed:
    for every finite observation horizon H,
    a finite closed reversible system can reproduce the fresh-bath reset process exactly through H

by taking L>=H.

### q=3, L=5 exact replay

The supplied machine replay enumerates all 3^5 tapes and all enlarged states.

- enlarged state count: 3645;
- distinct images: 3645;
- inverse replay: exact;
- total-variation distance to iid reset trajectories:

      h=1: 0
      h=2: 0
      h=3: 0
      h=4: 0
      h=5: 0.

At the first wraparound,

      h=6: TV=2/3,

because the finite tape begins to return stored information.

This proves that **no finite-horizon projected experiment can certify absence of a hidden reservoir in general.**

---

## 8. Infinite/all-time no-go from the natural extension

The obstruction is stronger than a finite-horizon trick.

For an iid q-symbol reset process, take the two-sided sequence space

    Omega=S^Z

with Bernoulli product measure and shift

    sigma(...,x_-1,x_0,x_1,...)
      =
    (...,x_0,x_1,x_2,...).

The shift is invertible.

Read out only coordinate 0.

Then the observed process is exactly iid reset for all time.

This is precisely the canonical natural-extension pattern already accepted in reports 278, 281, 286 and 290.

Hence:

    boxed:
    an apparently open fresh-reset process has an exact globally invertible closed realization.

Even forward interventions on the current visible coordinate do not reveal the future tape: after shifting, the modified present record moves into the past while untouched future coordinates continue to supply fresh symbols.

Therefore complete forward preparation/intervention/readout access to the projected coordinate does not by itself prove absolute openness or absolute completeness.

---

## 9. Audit of the candidate tests proposed for 295

### A. Preparation/intervention/readout completeness

Useful for defining the operational algebra, but not sufficient for absolute completeness.

A hidden future tape can remain outside the admitted visible intervention algebra while reproducing all visible forward statistics.

### B. Vanishing residual Mori-Zwanzig memory

Not sufficient.

The full-reset projected process is already exactly Markov and memoryless, yet requires hidden fresh symbols in a reversible realization.

Therefore:

    boxed:
    zero residual memory != complete closed system.

### C. Stability under state enlargement

Not sufficient by itself.

A larger exact dilation can leave every projected path statistic unchanged.

### D. Information-balance closure

This is the decisive test for whether the **candidate variables themselves** are self-contained.

Positive F_e or B_e proves that hidden information is required.

But zeroing the defect by enlarging the state does not say that the enlargement is the final ontology.

### E. Target-symbol entropy accounting

This is exactly captured by F_e.

It supplies a rigorous lower bound and recovers report 291.

---

## 10. FIN-relative complete-system criterion

The strongest noncircular criterion available is therefore **relative to an independently sourced ambient FIN carrier**.

Let:
- Omega_FIN be the accepted ambient reversible carrier;
- E_FIN be the accepted primitive event/controller algebra;
- pi be the proposed state description.

Call pi **FIN-relatively complete** if:

1. the operational factorization/access conditions of report 256 identify the candidate factors without using external coordinates;
2. for every primitive e in E_FIN, forward fiber closure holds;
3. for every primitive e in E_FIN, backward fiber closure holds;
4. the controller/event record needed to specify e is explicitly accounted for rather than hidden in a stochastic average.

Then every primitive event descends to a permutation of the candidate state.

No omitted ambient degree of freedom can affect its event-conditioned future or past.

This is an exact complete-system certificate **relative to Omega_FIN and E_FIN**.

---

## 11. Why the criterion does not yet make alpha=0 unconditional

Suppose the n simultaneous visible records are independently established as the ambient complete carrier.

Then a local fresh reset has

    F=log_2 q>0

and fails the criterion.

Pure visible SWAP has

    F=B=0

and passes.

Within the alpha reset/SWAP family this would force

    boxed:
    alpha=0.

However current FIN does not yet prove that those simultaneous visible records are the exhaustive ambient carrier.

The canonical natural extension instead supplies a larger reversible history carrier capable of storing/supplying reset information.

Thus 295 cannot promote

    "visible simultaneous records are the complete system"

to an unconditional theorem.

The missing source has become sharper:

    boxed:
    FIN needs an independently derived ambient-carrier / simultaneous-slice boundary
    that determines which record factors are admissible as contemporaneous system degrees of freedom
    and which history/future records may not be silently used as a fresh spatial reservoir.

This is a typed composition problem, not a memory-decay problem.

---

## 12. Relation to content preservation

Reversible closure alone still does not force additive-content conservation.

Reports 257 and 269 already show closed bijective exact-reset gates that can rewrite record content.

Therefore the correct chain is:

    operational subsystem carrier
      +
    relative complete-system / reversible-factor closure (295)
      +
    record-content continuity / zero-rewrite principle (269,278)
      ->
    SWAP within the declared two-record reset-dilation class.

The complete-system criterion removes hidden fresh-bath events.

The content-continuity criterion removes closed content-rewriting bijections.

They solve different obstructions and both are required.

---

## 13. Verdict

P0-295 returns a mixed but proof-grade result.

### Positive theorem

There is an exact FIN-constructible criterion for **relative** completeness once an ambient reversible carrier and primitive event algebra are specified:

    boxed:
    bidirectional fiber closure
    <=>
    every primitive ambient event descends to a visible permutation.

Its entropic form gives explicit hidden-information lower bounds.

### No-go

There is no criterion based only on projected path statistics, finite residual memory, finite-horizon interventions or invariance under enlargement that can prove **absolute** completeness.

A finite reversible tape defeats every finite horizon, and the natural extension gives an exact all-time reversible completion.

### Consequence for report 291

The alpha=0 selector remains conditional.

295 replaces the vague word "complete" by an exact mathematical test, but the test still needs FIN to source the ambient simultaneous carrier on which it is applied.

The new foundational blocker is therefore not:

    "how do we test closure?"

That part is solved.

It is:

    boxed:
    what FIN-internal law sources the ambient simultaneous carrier / slice boundary
    rather than allowing the canonical history reservoir to count as additional system content?

---

## 14. Next proof-grade move

Keep the previously scheduled

    MICROSCOPIC-PERIODIC-BRIDGE-296

as the next numbered task, because it tests whether the reversible information-continuity construction survives directly at the exact microscopic leave-one-out Gibbs lane.

In parallel, add a new P0 source task after 296:

    AMBIENT-SIMULTANEOUS-CARRIER-SOURCE

with the acceptance test:

- derive the admissible contemporaneous carrier from existing FIN structure;
- distinguish it from future/past natural-extension records without semantic relabeling;
- then reapply theorem 295-A plus reports 269/278 to test whether conservative SWAP becomes compulsory.

Do not call this source solved merely because a convenient visible slice has been chosen.
