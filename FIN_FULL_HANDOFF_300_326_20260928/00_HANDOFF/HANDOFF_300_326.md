# FIN — FULL RESEARCH HANDOFF 300–326

Date: **2026-09-28**

Scope: every research task **300 through 326 inclusive**, continuing the predecessor handoff through 299.

Predecessor bundle:

`FIN_RESEARCH_AFTER_256_294_FULL_HANDOFF_295_299_20260927.zip`

The predecessor is not duplicated here. Task 300 starts from its accepted boundary.

---

# 1. Executive state after 326

The 300–326 campaign transformed the strongest one-unit FIN lane from a finite-N effective-generator calculation into a controlled **observed-process prediction program** and then closed the large-N **metastable exponential barrier**.

The strongest current chain is

```text
exact finite-N leave-one-out Gibbs chain
    -> exact D12-equivariant metastable basin map
    -> controlled operational preparation family
    -> short microscopic memory / initial-slip layer
    -> memory-renormalized 12-state D12-circulant generator
    -> exact Z3 quotient
    -> q_d-free clock + response-shape law across N
    -> finite-window observed-process prediction with explicit error budget
    -> blind/direct holdouts through N=12
    -> analytic preparation boundary-layer law
    -> global mod-3 separator certificate
    -> Gamma = B4
    -> capacity exponent = B4.
```

The two most important mathematical conclusions at the current boundary are:

### A. Observed-process control

For the declared preparation/readout protocol, frozen cross-N predictions survive direct microscopic holdouts through `N=12`. The current process-level safety envelope from report 313 is

\[
\boxed{B_{313}=0.05582267457195117\;TV}.
\]

At large opened N the microscopic reduction/history error becomes tiny; the dominant finite-N uncertainty is preparation transfer.

### B. Metastable exponent

Report 324 proves, for the declared leave-one-out chain, that the capacity exponent equals the global continuous communication height `Gamma`. Report 326 now certifies

\[
\boxed{\Gamma=B_4}.
\]

Therefore

\[
\boxed{
\lim_{N\to\infty}
-\frac1N\log\operatorname{cap}_N(A,B)
=B_4
}
\]

for the declared mod-3 metastable source/target sets under the conditions of 324.

The high-precision barrier value used by the final certificate is

\[
B_4=0.6622191371274619267\ldots
\]

(the earlier stored `0.6622191371274597` is the same numerical barrier at ordinary precision).

**The prefactor and controlled finite-N remainder remain open.**

---

# 2. Epistemic boundary — do not overclaim

The campaign does **not** derive a fundamental physical theory.

Still not derived:

- a physical source of the dimensionless gain `g`;
- SI time, length, temperature or energy calibration;
- spacetime, Lorentz symmetry, quantum mechanics, QFT, Standard Model or GR;
- a fundamental source of simultaneous physical units / spatial incidence (the blocker already isolated in 295 and earlier composition work);
- an Eyring–Kramers prefactor with rigorous finite-N remainder;
- a laboratory preparation/apparatus/raw-data validation.

The current results are rigorous / controlled mathematics for the **declared finite FIN statistical model** and its exact leave-one-out Gibbs heat-bath dynamics.

The user's philosophical working hypothesis “one FIN from which everything arises” remains a hypothesis, not a theorem exported by these tasks.

---

# 3. Frozen microscopic convention

Unless a report explicitly says otherwise, the 300–326 lane uses:

- rank-7 FIN potential;
- working gain
  \[
  g=5.145228719489142;
  \]
- exact occupation-count representation of `N` labelled copies;
- exact leave-one-out Gibbs heat-bath generator;
- localized basin definition obtained from
  \[
  \theta_0=(g/N)X_7^T n
  \]
  followed by full `V_g` descent and exact `D12` orbit propagation;
- 12 localized labels `J=0,...,11` plus residual/nonlocalized states;
- seven-bin reflection compression only after checking the relevant symmetry;
- symmetric categorical readout error when an explicit readout model is used.

Do not silently substitute the earlier empirical-refresh finite-N law.

---

# 4. Research chronology 300–326

## 300 — CONTROLLED-MICRO-TO-HISTOGRAM

First direct microscopic observable pipeline at `N=6`.

Key results:

- `12,376` count states;
- localized mass `0.992496493817772`;
- mass per localized basin `0.082708041151481`;
- exact microscopic histograms regenerated from the leave-one-out generator;
- preparation kill-test: the macro statement “prepare J=0” is insufficient because different microscopic distributions inside that basin produce different later histograms.

This task establishes the end-to-end microscopic -> observable route.

## 301 — WEIGHTED-MEMORY-PREPARATION-CLOSURE

Strengthens report 299. Short memory alone is insufficient; one must also control the unresolved spectral weight / initial slip.

For `N=6` the exact slow residue and fast gap give a bound of the form

\[
|C(t)-Ze^{-rt}|\le(1-Z)e^{-\gamma_{fast}t}.
\]

The exact initial-slip deficit decreases across `N=3..8`; at `N=8` it is below one percent. This supports a short boundary layer followed by a controlled slow generator, but is not an N->infinity proof by itself.

## 302 — JOINT-FINGERPRINT-PROTOCOL

Combines the earlier late-time shell fingerprint and early-time memory fingerprint into one frozen `N=6` protocol.

Important split:

- early times test multi-rate memory using slope curvature;
- late times test D12 response geometry after memory has mostly died.

The protocol includes explicit preparation and readout assumptions rather than treating `J=0` as a complete experiment.

## 303 — HELD-OUT-N7-MICROSCOPIC-PREDICTION

First direct heldout `N=7` microscopic histogram test.

Results:

- pre-existing `Q12(N=7)` predicts the new microscopic late histogram to about `0.8565% TV`;
- a simple `N=3..6 -> N=7` cross-N extrapolation predicts the observable histogram to about `1.29% TV`;
- individual rare shell rate `q1` can be wrong by ~21% while the observable remains accurate.

Conclusion: observable combinations generalize substantially better than every individual shell rate.

## 304 — OBSERVABLE-LEVEL-CROSS-N-LAW

Removes the six `q_d` rates from the predictive interface.

For any D12-circulant generator define

\[
\lambda_k=-\rho R_k,
\qquad R_4=1,
\qquad \tau=\rho t.
\]

Then the full localized response is determined by the dimensionless `R_k` and `tau`. Clock and shape are predicted separately.

Frozen `N=3..6` fits predict `N=7` and effective `N=8` histograms well. The preparation uncertainty becomes larger than the remaining clock/shape error and is identified as the next bottleneck.

## 305 — PREPARATION-CONTRACT-PREDICTION-SET

Propagates **all** microscopic distributions supported inside `J=0` at `N=7`.

Exact convex-hull optimization over all `2616` microstates gives a minimum distance from the declared comparator of

\[
6.328\%\;TV
\]

at the frozen late time.

Introduces the operational pinning family

\[
\mu_{N,\kappa}(n)\propto\pi_N(n)e^{\kappa n_0/N},
\]

and a preparation-information cost `D_KL(mu||pi)/ln2`.

## 306 — PREPARATION-UNIFORM-CROSS-N-BOUND

Shows that **macro-only J=0 fails**, while a declared operational preparation family can be controlled.

For arbitrary microdistributions in `J=0`, the later prediction set has diameter

\[
53.0271\%\;TV,
\]

so a single effective initialization must incur at least

\[
26.5136\%\;TV
\]

worst-case error.

For the pinning family `0<=kappa<=12`, a low-dimensional initial-slip map

\[
\alpha_k=\alpha_k(N,m),
\qquad m=E[n_0/N],
\]

trained on `N=3..6` gives heldout `N=7` max error ~`1.58% TV`.

Continuous-in-kappa certificate:

\[
TV(P_{pred},P_{micro})\le2.08145\%,
\]

with predicted FIN-comparator separation at least `6.89124%`, leaving certified margin

\[
\boxed{4.80979\%\;TV}.
\]

## 307 — FINITE-WINDOW-PROCESS-PREDICTION

Moves from one-time histograms to the two-time law

\[
P(Y_{t_1},Y_{t_2}).
\]

Heldout `N=7`:

- for `Delta tau=0.5`, max joint error ~`1.994% TV`;
- conditional-transition error ~`1.581%`;
- training envelope `5.129%`;
- predicted FIN-comparator separation ~`6.339%`.

This is the first process-level discriminator rather than a marginal-only test.

## 308 — THREE-TIME-PROCESS-CONSISTENCY

Tests

\[
P(Y_{t_1},Y_{t_2},Y_{t_3}).
\]

Heldout `N=7` is evaluated with a reversible spectral representation. Rank-120 and rank-140 agree to ~`1e-10 TV`.

Certified heldout error:

\[
\le2.5937\%\;TV
\]

versus training envelope `6.0955%`.

Predicted three-time FIN-comparator separation is ~`10.243%`; after all certified errors a positive margin of ~`3.99 p.p.` remains.

## 309 — FINITE-WINDOW-PROCESS-ERROR-THEOREM

Proves an exact path-law composition inequality. For process laws `P,Q`,

\[
TV(P_{1:m},Q_{1:m})
\le TV(P_1,Q_1)
+\sum_r E_{H_r\sim P}TV(P(Y_{r+1}|H_r),Q(Y_{r+1}|H_r)).
\]

Also derives an explicit rare-history / high-probability version.

Four-time `N=7` heldout passes, but strict supremum over every history can be much larger than the mass-weighted error. This motivates separating true history defect from cross-N preparation error.

## 310 — HISTORY-UNIFORM-LATE-BURN-CLOSURE

Separates same-N history closure from cross-N extrapolation.

Training-only selection gives earliest burn-in

\[
\boxed{t_{burn}=24}
\]

for the criterion `sup_h delta(h)<=5%` across `N=3..6`.

Heldout `N=7` same-N history defect at `t=24`:

\[
\boxed{2.533\%\;strict},
\qquad 0.765\%\;weighted.
\]

Therefore the earlier ~10% full-contract strict error was not evidence for a new persistent memory variable. A one-dimensional transition residual exists and can reduce worst-history error, but is not Pareto-superior in average/joint error and is not made canonical.

## 311 — POSTERIOR-ROBUST-CROSS-N-CLOSURE

Shows why strict posterior sets are fragile.

A prior uncertainty of only ~`1.61% TV` can be amplified by Bayes conditioning on rare readouts into a very broad posterior prediction set. The worst case is large but low probability.

Heldout `N=7` decomposition:

- transition cross-N error ~`0.573%`;
- same-N history defect ~`2.517% strict / 0.759% weighted`;
- preparation/posterior is the dominant strict error source.

Verdict:

- **FAIL** for a uniform set-valued posterior theorem;
- observed two-time joint law remains good.

## 312 — STRUCTURED-PREPARATION-PREDICTION-SET + DIRECT N8

First direct microscopic `N=8` holdout after a frozen six-dimensional preparation hull.

The structured latent preparation set **fails** to cover all `N=8` preparations / posteriors exactly.

But the observed joint law passes strongly:

\[
\boxed{TV_{joint,max}=1.74418\%}.
\]

Direct `N=8` state count:

\[
75,582.
\]

Clock prediction error is only ~`0.437%`.

Error decomposition shows preparation dominates:

- reduction/history ~`0.383%`;
- transition cross-N ~`0.253%`;
- preparation ~`1.430%`.

This establishes the policy: **joint-law is the primary certification object, not strict rare-history posterior coverage.**

## 313 — JOINT-LAW-STABILITY-BOUND

Proves the additivity theorem

\[
TV(P_{micro},P_{pred})
\le
\epsilon_{red}+\epsilon_{trans}+\epsilon_{prep}.
\]

Grouped `N=3..6` maximum of the correlated component sum gives the frozen process envelope

\[
\boxed{B_{313}=5.582267457195117\%\;TV}.
\]

This becomes the primary pre-holdout safety envelope for later direct N tests.

## 314 — PREPARATION-RESIDUAL-DIRECTION-LAW

Preparation residuals are strongly low-dimensional:

- first principal component initially dominates;
- first two components explain about `98.7%` in the N<=8 analysis.

However the residual amplitude is not a smooth enough function of N for an aggressive cross-N correction. A low-rank correction is frozen only as optional diagnostic, not as a replacement for envelope 313.

## 315 — DIRECT N9 HOLDOUT

Direct `N=9` state count:

\[
167,960.
\]

Frozen prediction passes the report-313 envelope. A conservative classification ambiguity budget initially gives a certified error below ~`4.04% TV`.

Subsequent exactification in 316 shows the actual joint error after enforcing a valid probability simplex is even smaller (~`1.80% TV`).

Clock:

\[
\rho_9^{exact}\approx0.00381767405,
\]

with pre-holdout prediction error ~`2.28%`.

## 316 — SIMPLEX-NATIVE-PREPARATION-CORRECTION

Important formal repair: raw cross-N preparation predictions can leave the probability simplex. Projection onto the simplex is therefore made **mandatory** before interpretation as a probability law.

The raw `N=9` prediction had a negative component near `-0.00666`; after projection:

- preparation error ~`1.534% TV`;
- joint error ~`1.805% TV`.

A probability-preserving low-rank residual correction is tested under a frozen Pareto rule. Best acceptable model:

```text
rank = 1
quadratic z(m)
shrinkage s = 0.20
```

It yields a modest improvement without materially degrading any grouped holdout. It remains optional; envelope 313 remains primary.

## 317 — DIRECT N10 HOLDOUT

Frozen-before-open direct `N=10` test.

Microscopic size:

\[
352,716\text{ states},
\qquad22,523,436\text{ nonzero generator entries}.
\]

Exact slow clock:

\[
\rho_{10}=0.00206014358214.
\]

Original clock-law prediction error grows to ~`4.31%`, identifying the next bottleneck.

Certified process error:

\[
\boxed{3.2563\%\;TV}<B_{313}.
\]

Reduction/history is already tiny (~`0.094%`).

## 318 — BARRIER-AWARE-CLOCK-LAW

Replaces the simple log-linear clock with

\[
\beta_N=-\log(\rho_N/\rho_{N-1})
=B_4-\frac{c}{N}-\frac{d}{N^2}.
\]

`B4` is frozen from the barrier lane, not fitted freely.

Rolling holdouts through `N=10` reduce maximum clock error from ~`4.31%` to ~`2.03%` and maximum transition-row error from ~`1.01%` to ~`0.48%`.

At this stage `B4` was still a barrier candidate, so the law was explicitly conditional on the barrier lane.

## 319 — FROZEN N11 PREDICTOR

Full N=11 prediction frozen before microscopic construction using the 318 clock architecture and the existing shape/preparation rules.

Predicted clock:

\[
\rho_{11}^{pred}=0.001125516559.
\]

Predicted FIN-comparator separation ~`8.493% TV`.

## 320 — DIRECT N11 HOLDOUT

Direct state count:

\[
705,432.
\]

Introduces an exact `C12` representation reduction of the reversible generator. The method was first validated on `N=10`, reproducing known eigenvalues to ~`1e-15`.

Exact `N=11` clock:

\[
\rho_{11}=0.00109960859450,
\]

clock error ~`2.356%`.

Certified process error:

\[
\boxed{2.2405\%\;TV}.
\]

Reduction/history ~`0.050%`; preparation remains dominant.

## 321 — FROZEN N12 PREDICTOR

Predictor frozen after N11 but before N12.

\[
\rho_{12}^{pred}=0.00059654031813.
\]

Predicted FIN-comparator separation ~`9.040% TV`.

## 322 — DIRECT N12 HOLDOUT

Direct state count:

\[
1,352,078.
\]

Exact C12 orbit reduction:

\[
112,720\text{ cyclic orbits},
\]

including shorter stabilizer-compatible periodic orbits rather than discarding them.

Exact clock:

\[
\rho_{12}=0.000581379607710,
\]

clock error ~`2.608%`.

Certified process error:

\[
\boxed{2.2884\%\;TV}.
\]

Error decomposition:

- reduction/history ~`0.0327%`;
- transition ~`0.543%`;
- preparation ~`1.834%`.

Observed local exponent

\[
\beta_{12}=0.6373056592
\]

moves closer to `B4`.

## 323 — PREPARATION-BOUNDARY-LAYER-LAW

Turns preparation transfer from a mostly empirical error into an analytic finite-N object.

Let

\[
D=N-n_0.
\]

The pinning family satisfies exactly

\[
\frac{d\mu_{N,\kappa}}{d\mu_{N,0}}
\propto e^{-(\kappa/N)D}.
\]

Thus the true small parameter is `kappa/N`, and

\[
\frac d{d\kappa}E[n_0/N]
=\frac{\operatorname{Var}(D)}{N^2}.
\]

For any later observable the preparation response is `O(1/N)`.

Local Laplace geometry gives

\[
p_0^*=0.9806923477,
\]

\[
N\operatorname{Var}(p_0)\to0.02123797358,
\]

and

\[
\boxed{
m_N(\kappa)
=0.9806923477+
\frac{-0.1358252077+0.02123797358\,\kappa}{N}
+O(N^{-2}).
}
\]

The analogous six-component initial-slip coefficients are structurally identified but not yet fully derived.

## 324 — B4-CAPACITY/PREFACTOR

Proves the key conditional large-N theorem.

For the declared exact leave-one-out Gibbs chain:

- number of count states is polynomial in N;
- type weights satisfy the same large-deviation potential `Phi` up to `O(log N)`;
- allowed single-step rates are subexponential;
- Thomson path and Dirichlet cut bounds squeeze the capacity exponent to the continuous communication height.

Therefore, if

\[
\Gamma=B_4,
\]

then

\[
\boxed{
\lim_{N\to\infty}-\frac1N\log\operatorname{cap}_N=B_4.
}
\]

At 324 the missing blocker was **only** the global proof that `Gamma=B4`; the prefactor was also still open.

## 325 — GLOBAL-B4-COMMUNICATION-HEIGHT

Exact topological reduction of the remaining barrier problem.

Define mod-3 masses `P_a`. Any path changing sector must cross

\[
P_0=P_1\ge P_2
\]

or a symmetry-equivalent wall.

On the representative wall, unrestricted full-simplex minimizations and a constrained KKT census strongly identify the `d=4` orbit as the wall minimum.

Most important constructive result: the equality-wall stationary problem reduces exactly from the 11-simplex to a **six-dimensional** fixed-point system with the normal multiplier eliminated analytically:

\[
z(y)=\frac12\log\frac{S_1(y)}{S_0(y)},
\qquad
y=gB^Tp(y).
\]

Status remained PARTIAL because complete interval exhaustion was still missing.

## 326 — CONSTRAINED-SEPARATOR-INTERVAL-EXHAUSTION

Closes the 325 blocker.

### Interior half-wall certificate

Exact representative half-wall:

\[
P_0=P_1\ge P_2.
\]

Internal `C4` symmetry reduces to wedge

\[
y_0\ge|y_1|.
\]

A boundary-safe six-dimensional exhaustive cover processed

\[
\boxed{1,548,513\text{ boxes}}
\]

with zero unresolved boxes. Boxes that can touch `P0=P1=P2=1/3` are not improperly eliminated by interior stationarity.

Safety policy:

- every Taylor lower bound is inflated downward by `1e-8`;
- an energy exclusion is accepted only with margin `>1e-6`.

### Triple-junction certificate

The boundary

\[
P_0=P_1=P_2=1/3
\]

has an exact five-dimensional dual reduction. A separate cover processed

\[
\boxed{54,203\text{ boxes}}
\]

with zero unresolved boxes and excludes any triple-junction point at or below `V_d4`.

### Local d4 certificate

High-precision interval Krawczyk isolates the representative root. On the full local cube `|y-y*|_inf<=0.02`,

\[
\boxed{
\lambda_{min}(\nabla^2\Phi_{sep})
\ge0.0551088593996>0.
}
\]

### Upper route

High-precision interval subdivision proves monotonicity of the explicit straight path

```text
localized minimum -> d4 saddle -> neighboring localized minimum.
```

Therefore the path maximum is exactly the d4 saddle.

### Final theorem

\[
\boxed{\Gamma=B_4.}
\]

Combining with 324:

\[
\boxed{
\lim_{N\to\infty}
-\frac1N\log\operatorname{cap}_N(A,B)
=B_4.
}
\]

This closes the exponential metastable scale. **The Eyring–Kramers prefactor remains open.**

Arithmetic qualification: global covers use conservative long-double enclosure formulas with large explicit safety guards; local/path certificates use high-precision `mpmath.iv`. An independent directed-rounding MPFI/Arb replay is recommended as a formal audit upgrade.

---

# 5. Direct holdout ledger

The central process-level holdout sequence after the prediction architecture stabilized is:

| N | direct status | exact/observed slow clock | principal process result |
|---:|---|---:|---|
| 7 | PASS | `rho≈0.0126404` | multi-time process and preparation-uniform tests pass |
| 8 | PASS | `0.00698967187` | two-time joint max error `1.744% TV` |
| 9 | PASS | `0.00381767405` | valid-simplex joint error later tightened to `~1.805% TV` |
| 10 | PASS | `0.00206014358` | certified joint error `3.256% TV` under conservative ambiguity budget |
| 11 | PASS | `0.00109960859` | certified joint error `2.2405% TV` |
| 12 | PASS | `0.000581379608` | certified joint error `2.2884% TV` |

Do not reinterpret this table as an asymptotic proof by itself. The asymptotic barrier proof comes from 324+326, not from fitting the holdout sequence.

---

# 6. Current error hierarchy

At the largest opened N, the hierarchy is consistently:

```text
microscopic history/reduction error  << preparation-transfer error,
```

with cross-N transition/clock error in between.

Representative N=12 decomposition:

```text
reduction/history  ~0.0327% TV
transition         ~0.5433% TV
preparation        ~1.8339% TV
```

Therefore the next predictive improvement should target the analytic preparation map, not add arbitrary hidden-memory variables.

---

# 7. Closed / downgraded ideas — do not restart blindly

1. **Macro-only `J=0` as a complete preparation** — false; report 306 gives a 53% prediction-set diameter.
2. **One new six-rate fit at every N** — rejected as predictive strategy; observable-level clock/shape law is preferred.
3. **Strict posterior-set coverage as primary validation** — too unstable on rare outcomes; use joint-law certification.
4. **Calling preparation-transfer error “persistent microscopic memory”** — wrong; report 310 separates them.
5. **State-only scalar memory correction** — failed grouped cross-validation.
6. **Aggressive low-rank preparation correction** — violates Pareto or probability simplex; only the small optional 316 correction survives.
7. **Raw cross-N prior without simplex projection** — invalid; projection is mandatory since 316.
8. **Simple log-linear clock** — superseded by barrier-aware law 318.
9. **Generic 12D branch-and-bound for B4** — too loose; replaced by the exact 6D separator reduction and 5D triple-boundary reduction.
10. **Pinsker / coordinatewise KL quadratic shortcut for B4** — too weak; documented in 325.
11. **B4 as only a mapped-saddle candidate** — superseded by the global separator certificate 326 in the declared model.

---

# 8. Current theorem / claim register

## Exact / theorem-level within declared model

- exact leave-one-out finite-N microscopic generator;
- exact D12 basin/orbit symmetry operations used by the declared basin construction;
- exact algebraic D12 response representation and exact effective Z3 lumping once Q12 is declared;
- exact process-TV composition theorem of 309;
- exact joint-law stability theorem of 313;
- exact preparation defect-tilt identity and derivative identities of 323;
- conditional capacity-exponent theorem of 324;
- exact topological mod-3 separator reduction of 325;
- computer-assisted global separator result `Gamma=B4` of 326;
- therefore capacity exponential rate `B4` when 324 and 326 are combined.

## Strong finite-N validated results

- operational preparation prediction set through N=7;
- multi-time observed-process predictions;
- direct blind/heldout process tests through N=12;
- barrier-aware finite-N clock law;
- exact C12 representation reduction for large-N slow spectra.

## Still conditional / open

- Eyring–Kramers prefactor;
- controlled moderate-N asymptotic remainder;
- six analytic post-burn preparation coefficients;
- directed-rounding independent replay of the global 326 cover;
- physical interpretation beyond the declared dimensionless statistical model;
- simultaneous-carrier / spatial-incidence source problem from the earlier composition program.

---

# 9. Next research program

## P0 — 327 EYRING-KRAMERS-PREFACTOR-AND-REMAINDER

Now that `Gamma=B4` is closed, derive the actual capacity / exit-time prefactor.

Acceptance tests:

1. derive the local saddle and valley Hessian factors in the **correct discrete leave-one-out chain geometry**;
2. account for the multiplicity of d4 gates and the C4 intra-sector quasi-equilibrium;
3. derive a no-refit asymptotic formula, not `A,alpha` fitted to N;
4. compare it to exact D4-quotient capacities already available through N=12;
5. give a controlled remainder bound or a clearly stated conditional theorem.

Do not promote the currently drifting fitted power `N^alpha` to a theorem.

## P1 — 328 PREPARATION-OBSERVABLE-LAPLACE-COEFFICIENTS

Complete task 323 by computing the six analytic post-burn Fourier coefficients

\[
\alpha_k(N,\kappa)
=
\alpha_k^{(0)}
+\frac{a_k+\kappa b_k}{N}
+O(N^{-2}).
\]

This attacks the dominant finite-N observed-process error.

## P2 — DIRECTED-ROUNDING-REPLAY-326

Independent arithmetic audit:

- replay the already-defined 6D and 5D covers with MPFI/Arb or explicit directed rounding;
- do **not** alter boxes, thresholds or mathematical reductions;
- compare exclusion-ledger hashes.

This is an audit/verification upgrade, not a new scientific hypothesis.

## Parallel fundamental lane

The one-unit/metastability successes here do not solve the older composition blocker: FIN still does not derive the exhaustive simultaneous physical carrier / incidence structure needed for an emergent-space claim. Continue that lane only from the explicit 295-era source problem, not by relabelling the present count-space or barrier graph as physical space.

---

# 10. Handoff instructions for the next agent

1. Read this file first.
2. Treat the predecessor through 299 as accepted baseline unless auditing it explicitly.
3. Read reports 323–326 before attempting new asymptotics.
4. Do **not** open another N merely to extend the numerical sequence before addressing 327/328 unless a specific falsification test requires it.
5. Preserve the distinction:
   - theorem/exact algebra;
   - computer-assisted certificate;
   - finite-N validation;
   - conditional theorem;
   - hypothesis/physical interpretation.
6. Keep the raw prior -> simplex projection rule from 316.
7. Keep `B313=0.05582267457195117` as the established process-level empirical safety envelope unless a new, predeclared recalibration program replaces it.
8. Do not weaken the 326 conclusion back to “mapped saddle compatibility”: within the declared model, the global separator has now been exhausted and `Gamma=B4` is the current accepted result, subject to the arithmetic-audit note.
9. Do not strengthen 326 into an Eyring–Kramers prefactor claim; that is 327.
10. No automatic GitHub commit/push is authorized by this handoff.

---

# 11. Package organization

`01_REPORTS/` contains reports 300–326.

`02_PACKAGES/` in the FULL bundle contains the original continuation ZIPs, including large exact finite-state caches.

`03_KEY_ARTIFACTS/` contains frozen prediction/certificate JSONs needed to orient the next program quickly.

`04_MANIFESTS/` contains SHA-256 inventories and predecessor reference hashes.
