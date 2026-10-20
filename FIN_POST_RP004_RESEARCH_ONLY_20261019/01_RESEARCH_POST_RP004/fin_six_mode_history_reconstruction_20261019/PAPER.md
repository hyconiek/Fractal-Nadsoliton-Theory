# FIN after the Fisher–Rao obstruction: a minimal six-mode temporal-history certificate for strict Z12

**Date:** 2026-10-19. **Research status:** scoped finite-dimensional analytical theorem, exact symbolic replay, numerical controls. **Not:** a law derived from primordial information, a new physical-clock source, a proof of emergence, independent peer review, or empirical verification.

## Abstract

FIN's central source-law problem F01 remains unresolved: the strict static graph kernel determines mathematical semigroups only **once their operational meaning and time law are supplied**. Historical least-change and KL-history proposals, including GEN-532–554 and GEN-654, do not select a nontrivial kinetic law. Subsequent studies ST8397–ST8547 already established the exact logarithmic-mean gradient representation, refinement limits and a source-level nonuniqueness obstruction; we do not re-prove them as discoveries. Instead, we ask a sharply typed question: **within the finite continuous-time, stationary, reversible, D12-symmetric one-label Markov family on strict Z12, what *temporal correlation object* is sufficient and minimal to identify all graph-jump rates?**

We prove that the six labelled Fourier-mode stationary autocorrelations at any **one common positive model-time lag** recover all six radial rates via an explicitly invertible real matrix, with determinant `-3456 sqrt(3)` and an explicit inverse. If one modal record is absent, an open, positive one-parameter class of different generators shares every retained modal autocorrelation **at all lags**. For missing mode 6 an exact vector is `(-1,+1,-1,+1,-1,+1)`. We derive a finite-error reconstruction bound, exact necessary-and-sufficient log-correlation admissibility inequalities, and identify the still-unfixed comparison-clock gauge. This is a *minimal conditional history-provider interface*, not its independent production from the nadsoliton.

## 1. Provenance and scope

Inputs inspected: current-conversation `AGENTS(1).md`, `GEN001-gen956(1).md`, prior `PAPER_INTRINSIC_DYNAMICS_FIN.md`, the local strategy/blocker documents. AGENTS §§ST6897–ST8547 already paid generic reversible activity freedom and entropy/mobility nonuniqueness. Older GEN-003 already tested relative operator algebra under marginalization, but this neither fixes temporal rates nor proves a unique source. The existing strict weights are quoted below. External background spectral reconstruction and reversible-chain theory is **known mathematics** (see literature links in SOURCES.md). Novelty is restricted to the explicit Z12 six-mode certificate/counterfamily/robustness result in this particular contract; source-independent literature priority is not claimed.

**Carrier:** `Z12 = {0,…,11}` as a set of alternative *states of one label*, not a derived physical 12-site spacetime. Let radial strictly positive weights `w_d` be indexed by cyclic separation `d=1,…,6`:

```
w1=0.46998567264502017
w2=0.19204355169010282
w3=0.09142861427792497
w4=0.04702916874565040
w5=0.02413122336363006
w6=0.011070817321442113
```

Let `Q_ij=w_distance(i,j)` for `i != j`, and `Q_ii=-sum_{j!=i}Q_ij`. Thus `Q=-A`, with the **supplied** symmetric graph Laplacian `A`; this Markov interpretation and dimensionless model time are premises. The invariant probability `π_j=1/12` is uniform. The six unknown positive `w` parametrize the entire specified reversible radial class. No assumption that nature selects this class is exported.

**History input:** one stationary path law under `Q`, with six *independently available*, labelled, time-resolved modal autocorrelations at a common model-time lag `τ>0`. Their labels and normalization are part of the observation contract. The spectral rates are not manufactured by plugging target `Q` into a proposed independent provider.

## 2. Theorem H1: a complete temporal correlation fingerprint

Define normalized real observables, for `k=1,…,5`,

\[
f_k(j)=\sqrt2\cos(2\pi kj/12),\qquad f_6(j)=(-1)^j.
\]

Each has stationary mean zero and variance one. They are orthogonal under π. Let

\[
C_k(\tau)=\mathbb E_\pi[f_k(X_\tau)f_k(X_0)] .
\]

**Theorem H1.** For the frozen Markov family and `τ>0`,

\[
C_k(\tau)=\exp(-\lambda_k\tau),\qquad
\lambda_k=\sum_{d=1}^5 2w_d(1-\cos(2\pi kd/12))+w_6(1-(-1)^k).
\]

Consequently `λ_k=-log(C_k(τ))/τ`, and the six radial weights are uniquely determined by the six correlations.

**Proof.** A real symmetric circulant `A` commutes with translations, so each `f_k` is its eigenfunction. For distances `d<6`, each site has two neighbours; the antipodal shell `d=6` has only one. Hence `A f_k=λ_k f_k` with the stated expression. `P_τ=e^{τQ}=e^{-τA}` is selfadjoint in `L²(π)`; therefore `C_k(τ)=<f_k,P_τ f_k>_π=e^{-τλ_k}`. The invertibility of the exact 6×6 coefficient matrix below closes the reconstruction. □

The matrix `λ=B w` is **exactly**

\[
B=\begin{pmatrix}
2-\sqrt3&1&2&3&2+\sqrt3&2\\
1&3&4&3&1&0\\
2&4&2&0&2&2\\
3&3&0&3&3&0\\
2+\sqrt3&1&2&3&2-\sqrt3&2\\
4&0&4&0&4&0
\end{pmatrix},
\quad \det B=-3456\sqrt3\neq0.
\]

With `r=sqrt(3)`, an **explicit left and right inverse** is

\[
B^{-1}=\begin{pmatrix}
-r/12&-1/12&0&1/12&r/12&1/12\\
-1/12&1/12&1/6&1/12&-1/12&-1/12\\
0&1/6&0&-1/6&0&1/12\\
1/12&1/12&-1/6&1/12&1/12&-1/12\\
r/12&-1/12&0&1/12&-r/12&1/12\\
1/6&-1/6&1/6&-1/6&1/6&-1/12
\end{pmatrix}.
\]

`verify.py` checks `B B^-1=I` and the determinant **symbolically**. This proves injectivity independent of numerical conditioning. In particular, the input spectral rates are the *output* of a chosen path law; their availability is not a theorem that a natural intrinsic observer can access them.

## 3. Theorem H2: five autocorrelations are insufficient, even for all times

Let `v=(-1,1,-1,1,-1,1)^T` and `w(η)=w+η v`. Then

\[
Bv=(0,0,0,0,0,-12)^T.
\]

**Theorem H2.** For each `|η|<min_d w_d`, all entries `w(η)` are strictly positive. The generator `Q(η)` has the same uniform equilibrium, reversibility, dihedral symmetry and **exactly the same autocorrelation functions `C_k(τ)` for every τ≥0 and `k=1,…,5`**, but `λ_6(η)=λ_6(0)-12η`. Hence the five-mode data cannot identify the full radial generator, even with infinitely many perfect timestamps. Any omitted individual `k` leaves a one-dimensional positive local alternative by deleting one row of an invertible `B`.

**Proof.** Positivity is open in six-dimensional weight space. The displayed exact multiplication leaves five eigenvalues unchanged and moves one. Since each `C_k(τ)` is a single exponential, all five entire functions agree. Invertibility of `B` implies any five rows have full row rank, leaving a one-dimensional kernel. Perturbing any interior positive baseline along that kernel yields distinct admissible rates. □

Numerical control `η=0.004` yields `max_{k≤5}|λ'_k-λ_k|<2e-15` and `λ'_6-λ_6=-0.048`. This counterfamily does **not** preserve *all* possible mixed observables; an observation with a nonzero projection onto mode 6 may detect the difference. The theorem concerns the declared five labelled modal histories.

## 4. Theorem H3: quantitative stability and the clock gauge

Let `Ĉ_k=C_k(τ)+δ_k`, with `|δ_k|≤η`, where `C_min=min_k C_k(τ)>η`. Define `ŵ=B^-1[-log(Ĉ_k)/τ]`.

**Theorem H3.**

\[
\boxed{\|\widehat w-w\|_\infty \leq
\frac{11}{12}\frac{\eta}{\tau(C_{\min}-\eta)}.}
\]

**Proof.** For `c,c+δ≥C_min-η>0`, the mean-value theorem gives `|log(c+δ)-log c|≤η/(C_min-η)`. Apply `B^-1` and use the exact maximum absolute row-sum `||B^-1||_∞=11/12`. □

This establishes finite-record-error sensitivity **conditional on** an absolute deterministic upper bound on the six correlation errors. It does not provide a statistical confidence interval, a sample count or a calibrated apparatus. In a controlled synthetic perturbation at `τ=0.6`, `η=10^-5`, the numerical weight error is `3.16×10^-5` against the rigorous deterministic bound `6.23×10^-5`.

**Clock ambiguity:** If the lag `τ` is not calibrated independently, the six correlations identify `τ w` and hence all **ratios** `w_d/w_e`, but not an absolute rate scale: `w→c w` with `τ→τ/c` leaves the data unchanged. Even `τ` provided as a model clock is not a physical second. An arbitrary nonlinear reparametrization of histories lies outside the time-homogeneous CTMC contract and is not addressed by the scalar gauge theorem.

## 5. Theorem H4: a coefficient-free falsification criterion for six history records

Suppose six normalized positive numbers `c_k in (0,1)` are asserted to equal the modal stationary autocorrelations at the **same declared lag** `τ>0` of **some** strictly positive radial reversible Markov law on `Z12`. Define `b_k=-log c_k` and `a=B^-1 b` using the **fixed exact** inverse in H1.

**Theorem H4.** Such a generator exists **if and only if** all six components `a_d` are strictly positive. When this holds, it is unique at fixed `τ`, with `w_d=a_d/τ`. Equivalently, the admissible six-component logarithmic correlation vector lies in the *open simplicial cone* `B R_{>0}^6`; its closure corresponds to nonnegative weights and may include reducible generators. **Six linear inequalities in `-log c_k`**, not six arbitrary fitted constraints, characterize model compatibility at one lag.

**Proof.** Necessity follows by applying the inverse in H1 to `b=τ B w` with all positive `w`. For sufficiency set `w=a/τ`. All off-diagonal CTMC rates are strictly positive, and the resulting symmetric circulant Markov generator has correlation eigenvalues `exp(-τ B w)=c` by H1. Uniqueness again follows from `det B≠0`. □

This test is not vacuous: for `τ=1` and the individually decaying positive modal targets `c=(e^-3,e^-1,e^-1,e^-1,e^-1,e^-1)`, inversion produces `w_1=1/12-sqrt(3)/6<0` and `w_2=-1/12<0`. Thus even though all six `c_k` are in `(0,1)`, **no** strictly positive radial reversible Z12 Markov generator realizes them at the given lag. The counterexample does not assert that these target correlations are generated by some other physically admissible common process.

For an independently supplied **single-lag** dataset, H4 tests compatibility with the declared radial family. To conclude the whole two-time process is a stationary homogeneous CTMC, one must *also* test the semigroup pattern `C_k(t+s)=C_k(t) C_k(s)` and the other operational hypotheses. The test does not select the frozen original strict weights unless the actual six values match that frozen prediction.

## 6. A4 source gate: what these theorems do NOT solve

H1–H4 solve an **inverse problem after a trajectory law is given**. They do not create that path law, select the uniform stationary ensemble from the nadsoliton, justify the observed state labels, prove Markovianity, source D12 symmetry or radial positivity, or select the required intrinsic clock. Constructing the modal correlations `C_k(τ)` from a postulated `Q` and then reconstructing `Q` is formally correct but **circular as a theory of physical origin**.

Exactly one typed missing primitive can be stated:

> **HISTORY-PROVIDER-Z12:** an independently justified intrinsic history object or set of correlation functionals, with a state/readout map and time-ordering, which produces six labelled `C_k(τ)` without using the target generator; or a source theorem deriving its equivalent nonlinear microscopic law. It must discriminate from the explicit H2 counterfamily (especially the k6 signature), and, separately, justify any transformation from strict Z12 to the interacting A7 and conditional 1+1 field models.

This is a **definition of a missing interface**, not a discovery of that source. A static graph does not supply it. If it cannot be derived from the presently admitted equations, F01 remains `OPEN_SOURCE`. Never promote reconstruction of a prescribed process to fundamental generation of the universe.

## 7. Relation to the historical front, originality, and stop rule

- Earlier GEN-001 marginalization and ST8547 Schur/Onsager refinement concern changing descriptions and mobility; they do not provide the independent correlations defined here.
- Past FIN RG-E, KO and AC campaigns already established Fourier identifiability and generator recovery **for other conditional observation types and interacting models**; we reuse their conceptual distinction and do not assert that inverse spectral reconstruction is new in the literature.
- The strict operator `A` here is *not* the rank-seven `A7`, and the `Z12` single-label Markov evolution is *not* the labelled interacting `N2` heat-bath, the 1+1 legacy field action or physical spacetime.
- Even six perfect time correlations do not choose which class is fundamental: another non-Markov realization or hidden-state model can have the same selected correlations.
- The target is closed once H1–H4 and the positive H2 counterfamily are checked. Do **not** run generic source inventories, coarse-reduction sweeps, B4, RP004, QM/GR/SM or claim `L_total`.

**Disposition:** `SIX_MODE_HISTORY_INVERSE_CERTIFICATE_PAID_CONDITIONALLY; FUNDAMENTAL_HISTORY_SOURCE_UNPAID`. This is a model theorem, not an EV claim.
