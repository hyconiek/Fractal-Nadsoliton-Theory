# FIN: Event-Indexed Relational Histories and an All-Lag Return-Correlation Obstruction on Z12

**19 October 2026 — bounded independent mathematical research, not a physical theory or independent peer review.**

## Abstract

The previous FIN six-mode certificate reconstructs a positive radial continuous-time generator on a twelve-label state space from six **labelled** stationary Fourier autocorrelations at a supplied common time lag. Here we remove the calibrated model-time coordinate and consider the embedded chain of **actual state changes**, and separately weaken observations to a label-independent equality/return record. Two contrasting statements result. (E1) Six Fourier correlations after **one actual jump** recover all six *relative* shell rates, with an exact polyhedral compatibility condition, without a time unit. The absolute jump frequency remains absent. (E2) Even a complete continuum of two-time return probabilities at all lags, including the correct common model clock, fails to identify a radial generator **up to arbitrary state permutation** in the admissible positive family. We construct an explicit pair with rates in Q(sqrt(3)), identical return records for all continuous times and all jump counts, but an unequal permutation-invariant cubic edge-weight moment. Hence a scalar, anonymous correlation history cannot in general serve as a complete source object for FIN. This is an exact class obstruction, not a theorem ruling out intrinsic correlated-process ontologies. All forward probabilities below still depend on an assumed Markov generator.

## 1. Source guardrails and precise mathematical object

- Inputs: archived `AGENTS(1).md`, `GEN001-gen956(1).md`, the previously generated six-mode `PAPER.md` and `HANDOFF.md`, FIN strategy/blocker snapshots. The files are **uploaded historical snapshots**, not the live local Git checkout. The six-mode inverse matrix and all numerical constants are rechecked independently here.
- The 12 labels are *alternative states of one finite random process*, not derived locations in physical 3D space. This class is the reversible radial single-label comparison model `Q=-A`, **not** the native many-copy leave-one-out heat-bath process on the N2 label-pair state space. No native HB, PDE, source of coefficients, physical clock, apparatus, or 3D geometry is claimed.
- For 12 labels `i∈Z12`, let `d(i,j)=min(|i-j|,12-|i-j|)`. Take `w_d>0`, d=1..6, and set `Q_ij=w_d` for i≠j, `Q_ii=-s`, where `s=2Σ_{d=1}^5 w_d+w_6`. Uniform `π=1/12` is reversible. Let `A=-Q` and `B` be the exact Fourier matrix in the previous H1 paper: `λ=Bw`, with `det B=-3456 sqrt(3)` and the explicit inverse therein. Angular modes k=1..5 each have multiplicity two; k=6 has multiplicity one.
- All rates and times are dimensionless model objects. Continuity and the Markov property are part of this **comparison class**; they are not deduced from the primordial informational nadsoliton.

## 2. E1 — exact actual-jump history and source/clock separation

**Theorem E1 (relational event-index reconstruction).** For the given finite CTMC, record the ordered states `Y_n` immediately following the n-th **actual state change**. The chain `Y_n` has transition matrix

`J_ij=w_{d(i,j)}/s` for `i≠j`, `J_ii=0`.

Define real stationary normalised Fourier functions `f_k(j)=sqrt(2) cos(2πkj/12)` for k=1..5 and `f_6(j)=(-1)^j`. Then

`m_k = Eπ[f_k(Y_1) f_k(Y_0)] = 1 - λ_k/s`.

The six event-indexed one-step correlations uniquely determine the **normalised** edge rates `a=w/s`:

`a = B^{-1}(1-m)` (six-component vectors).

Conversely, a real six-component vector `m` is the event-indexed modal signature of a strictly positive radial jump chain with no self-jumps **if and only if**

`a=B^{-1}(1-m)>0` componentwise, and `2Σ_{d=1}^5 a_d+a_6=1`.

Equivalently, the second condition is `1+2Σ_{k=1}^5 m_k+m_6=0` (the zero trace of J). This statement is on a fixed *recognisable labelled modal basis*: without it, see Theorem E2.

**Proof.** Because every state has the same true exit intensity s, holding times are iid exponential(s), independent of jump destinations, and the first next-state probabilities are qij/s, with no fictitious self-transition. Thus `Q=s(J-I)`. Fourier diagonalisation gives `J f_k=(1-λ_k/s)f_k`; stationary inner products give m. Apply the previously audited inverse B. Conversely strict positivity and normalisation define a zero-diagonal symmetric stochastic J; it has modal data m. The trace identity follows from `tr J=0`, eigenvalues 1 (once), m_k twice for k=1..5 and m6 once. □

**Corollary E1a (Poissonization).** The true number `N(t)` of state changes is a Poisson process of rate s and independent of the discrete skeleton. Consequently

`e^{tQ} = Σ_{n≥0} e^{-st}(st)^n/n! J^n`

and `Eπ[f_k(X_t)f_k(X_0)] = E[m_k^{N(t)}]=exp[-s(1-m_k)t]`.

This is **known finite-state Markov-chain/uniformization mathematics**, specialised exactly to the radial Z12 family. It is not a theorem of endogenous physical seconds: the continuum time t and exponentiation were assumptions of the original comparison CTMC. If all holding-time records are erased, `w→c w`, `s→c s` leaves the entire jump-indexed path law unchanged. The ratio `w_d/s` is identifiable, the overall activity is not. No special stochastic or finite-H mechanism is thereby sourced.

## 3. E2 — explicit all-lag anonymous return counterexample

For stationary `X_t`, consider only the label-independent event `X_t=X_0`, with return-correlation function

`R(t)=Prπ(X_t=X_0)=tr(e^{tQ})/12 = [1+2Σ_{k=1}^5 e^{-λ_k t}+e^{-λ6 t}]/12.`

This quantity is legitimate and invariant under arbitrary permutations of the 12 state labels. It is only a **two-time equality statistic**. It is not the complete law of all multi-time partitions/equality-pattern records.

**Theorem E2 (two nonisomorphic laws, identical returns at every time).** Fix the two positive shell-rate vectors:

`w = (16/15+sqrt(3)/30, 1, 61/60, 1, 16/15-sqrt(3)/30, 1)`,

`v = (43/40+sqrt(3)/40, 59/60, 1, 1, 43/40-sqrt(3)/40, 31/30)`.

Their common exit rate is `s=113/10`. Their 6 modal eigenvalue vectors are respectively

`λ=(121/10,122/10,123/10,124/10,125/10,126/10)` and

`λ'=(122/10,121/10,123/10,124/10,125/10,126/10)`.

Thus for **every** t≥0,

`R_w(t)=R_v(t)=[1+2Σ_{j=1}^5 exp(-(12+j/10)t)+exp(-(63/5)t)]/12`.

The embedded jump kernels Jw and Jv also have identical probabilities of returning to the original state after **every integer number of actual jumps**:

`R_w^jump(n)=R_v^jump(n) = [1+2Σ_{j=1}^5 (1-(12+j/10)/(113/10))^n + (1-(63/5)/(113/10))^n]/12`.

Nevertheless **there exists no permutation of the 12 vertices conjugating Qw to Qv**. The permutation-invariant directed-edge entrywise cubic statistic satisfies

`(1/12)Σ_{i≠j} (Qw_ij)^3 - (1/12)Σ_{i≠j}(Qv_ij)^3
 = (2Σ_{d=1}^5 w_d^3+w6^3)-(2Σ_{d=1}^5 v_d^3+v6^3)
 =11/24000 ≠ 0`.

**Proof.** Multiply each vector by the exact B from H1 (or the independently checked exact inverse). The modal spectra differ only by a swap of k=1 and k=2, both double eigenvalues; the simple k=6 eigenvalue is unaltered. Their full eigenvalue **multisets** are therefore identical, so the heat-trace/return functions coincide. Equal s gives identical Poissonized and event-indexed trace moments. For every vertex relabelling, the multiset of off-diagonal transition-rate entries, hence their third power sum, remains unchanged. But direct expansion in Q(sqrt(3)) gives the exact nonzero difference 11/24000. Therefore no such permutation exists. Positivity is immediate from the explicit rational/radical expressions. □

This is more informative than merely swapping k=1↔k=5, which corresponds to the internal relabelling `j→5j mod 12` and is a gauge/isomorphism. Here the swap between k=1 and k=2 is **not** an affine `Z12` automorphism, and the cubic invariant rules out **every** arbitrary state permutation.

## 4. Model-specific control at the displayed strict decimal weights

Using the six weights printed in the FIN strict kernel compendium as **decimal-frozen inputs**, a bounded numerical enumeration of all 5! permutations of the five double Fourier eigenvalues, while fixing the simple k=6 eigenvalue, finds exactly two strictly positive inverse images at a threshold `1e-11`: identity and the expected `k1↔k5` automorphism. This is a diagnostic, not an outward-rational proof about the exact transcendental strict kernel and not a theorem that nonisomorphic return counterexamples exist in any arbitrarily small neighbourhood of that actual tuple. The exact E2 counterpair concerns the **positive radial class**, not specifically the strict numerical FIN weights.

## 5. What this does and does not pay toward the missing intrinsic history source

The positive result E1 clarifies a real **relational clock interface**: ordered actual state changes can parameterize histories without any external unit, and six properly typed relational-mode records determine the normalised rates within an already supplied Markov family.

The negative result E2 shows that scalar all-lag equality correlations are not a universal complete history provider, even if their time scale is independently known and even if all rates are positive and D12 symmetric. Thus a provider which exports only the trace/return function cannot eliminate all same-static alternatives or recover the six labelled modal rates. A stronger provider could supply **distance-resolved one-step jump probabilities** `Pr(d(Y1,Y0)=d)=2 w_d/s` for d=1..5 or `w6/s` for d=6. But the interpretation of d as an intrinsically sourced observable is a separate F06 carrier/observation problem: internal Z12 label distance is **not** a proved spatial distance.

Even full event-indexed histories do not infer a physical duration per event. If a candidate intrinsic state law exports the labelled jump destination distribution but no waiting-time/activity law, an entire positive ray `cQ` remains. If it exports a Poisson event law, the rate is an input unless derived from another independently grounded object. Repackaging e^{tQ}, return probabilities, or the relative jump probabilities computed from assumed Q as a fundamental provider would remain circular.

**Stops:** F01 source law/open; F03 activity/open; F06 carrier, composition and geometry/open; strict/legacy split, selector, model units and physical/EV bridges unchanged. This work **does not** adopt a new law, clock, space, medium, or physical instrument.

## 6. Scope, priority and reproducibility

The jump-chain decomposition is standard; the 6×6 inverse was already paid in the prior local FIN H1. The exact radical all-lag, nonisomorphic Z12 counterpair, the trace/moment obstruction and the jump-history identifiability interface are the scoped new outputs of this research turn. No priority over the literature is claimed. The symbolic checker verifies determinant, inverse, exact counterpair, and exact cubic invariant. Finite numerical controls check heat kernels, event-index returns, and Poissonization but are **not replacements** for analytical arguments. Neither analysis licenses QM, GR, SM, light, particle masses, thermodynamics, or ToE claims.
