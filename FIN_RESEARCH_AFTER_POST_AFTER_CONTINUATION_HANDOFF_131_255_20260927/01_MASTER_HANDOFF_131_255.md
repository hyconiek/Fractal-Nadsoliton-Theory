# FIN — FULL CONTINUATION HANDOFF AFTER `FIN_POST_AFTER_CONTINUATION_FULL_HANDOFF_20260926`

## Scope

This handoff contains the research produced **after** the committed baseline handoff `FIN_POST_AFTER_CONTINUATION_FULL_HANDOFF_20260926`, i.e. reports **131–255** and their replay/data artifacts.

Baseline repository commit: `ad15a9098ecc5e1282f964ea8b159a8ec608d7c5`.

The dominant finite-N microscopic convention in this continuation is the **exact leave-one-out Gibbs heat-bath**. Earlier empirical-refresh results are not silently transferred into this lane.

## Executive verdict

The post-baseline campaign achieved one genuinely controlled effective level:

```text
exact finite-N leave-one-out Gibbs chain
    -> non-exact projection with short positive Mori-Zwanzig memory
    -> memory-renormalized 12-localized-state effective dynamics
    -> exact semigroup lumping to a Z3 metastable process.
```

The main result is therefore no longer merely that FIN has a rich stationary landscape. It now has a concrete example where **one microscopic stochastic law produces a nontrivial effective stochastic unit without fitting a separate coarse transition law or a separate coarse clock**.

At the same time, the campaign produced several important negative results:

- the representative localized minimum does **not** recursively split by the simple binary pitchfork mechanism;
- the k6 parity variable is not exactly Markov and is a poor economical first emergent unit;
- one-unit dynamics alone does not determine two-unit interaction;
- MaxEnt at fixed marginal laws selects independence rather than interaction;
- exchangeability, fixed valence, conserved edge number, holonomy, pair identity, and shared common-noise/fiber constructions do not by themselves derive spatial locality;
- microscopic copy labels are mean-field all-to-all variables, not spatial sites;
- state-space communication geometry is not physical space;
- information-record factorization is not yet spatial factorization.

The strongest new constructive multi-unit result is conditional but notable:

```text
Z3 MaxEnt reset
 + reversible deterministic dilation
 + exchange symmetry
 + Z3 equivariance
 + involution
 -> unique local gate = SWAP.
```

Closing the environment by placing these swaps on a cycle gives dimensionless diffusion with

`lambda_1 ~ 2*pi^2*rho/n^2`

without introducing a new continuous coupling coefficient. The remaining central blocker is not the value of a coupling `kappa`, but the **typed source of simultaneous spatial units and their elementary incidence**.

---

# 1. Campaign 131–137 — strengthened metastable-unit test

## 131 — recursive child route

The strengthened test froze the exact leave-one-out microscopic process, corrected the full-X7 k6 stability window, included competing escape channels and directly tested the representative main localized minimum.

Key result: the main localized branch was continued from just above the localization fold to `g=50`; no reflection-odd soft mode was found. The simple proposal

`stable localized parent -> two stable children -> repeat`

fails over the tested range. The unrestricted recursive-pitchfork route was closed as a P0 strategy.

The k6 pair survives only as a useful endogenous binary-incidence prototype, not as an economical first physical unit.

## 132 — k6 autonomy and projection nonlumpability

Two parallel reports carry number 132.

`K6_METASTABLE_ISOLATION_132` narrowed the candidate k6 autonomy interval to approximately

`g6 < g < 5.1503747508...`

based on known direct-pair and localized-escape saddles. Global autonomy was not proved.

`PARITY_PROJECTION_NONLUMPABILITY_132` gave an exact finite-N counterexample: the leave-one-out chain does not close exactly on parity count.

## 133–134 — correct process gate and memory

The weak matrix condition `P L J = L_coarse` was shown insufficient for full process equivalence. The proper alternatives are:

- exact intertwining/lumpability;
- approximate semigroup control on a declared time window;
- or an explicit non-Markov memory equation.

The parity projection has rapidly decaying Mori-Zwanzig memory at small N, but fast decay does **not** mean zero integrated effect.

## 135–137 — switch to non-binary metastable sectors

The existing 12 localized minima admit an exact three-sector quotient at the effective metastable level. Small-N exact committor tests showed the k6 pair is not yet a clean autonomous two-state unit. The campaign therefore switched from imposed binary hierarchy to the dynamically selected three-sector architecture.

---

# 2. Campaign 138–147 — one microscopic contract -> one controlled effective level

## Microscopic contract

Report 138 froze the finite-N convention. New finite-N statements in this lane are computed directly from the leave-one-out Gibbs generator.

## Three-sector variable

The 12 localized minima are grouped as

`C0={0,3,6,9}`, `C1={1,4,7,10}`, `C2={2,5,8,11}`.

Compared with parity, the mod-3 projection better resolves the populated metastable structure. Exact capacities and hitting times show a genuine slow coarse process.

## Short memory, large self-energy

For the three-sector projection, hidden memory decays rapidly, but its integrated self-energy strongly renormalizes the naive projected drift. This establishes the important pattern

`microscopic Markov -> short non-Markov memory -> renormalized local kinetic law`.

The memory-controlling hidden subspace is selected by symmetry coupling to `QSP`; the globally slowest hidden eigenmode need not control the observable memory.

## Parameter-free Z3 effective generator

After the memory boundary layer, the effective three-state generator has the symmetric form

```text
Q3 = k * [ -2  1  1
            1 -2  1
            1  1 -2 ]
```

with `k` derived from the same microscopic generator rather than fitted independently.

At N=6 the first-moment memory closure reproduces the exact slow rate to roughly 0.03–0.04%.

## Large-N trend through N=8

The coarse transition rate becomes rapidly slower with N while the coupled memory time remains O(1). The separation of memory and metastable time improves.

Whole deterministic-basin capacity was found to overcount recrossing flux and to have the wrong finite-N exponent. This motivated metastable cores rather than whole basins.

## Second level no-go

The isolated three-state cell has no further nontrivial autonomous quotient; the next exact reduction is equilibrium. Therefore a genuine second scale needs **multiple simultaneous units**, not more internal branch searching.

---

# 3. Campaign 148–177 — composition and incidence source audit

## Exchangeable grouping no-go

Arbitrarily tagging one exchangeable microscopic population into groups does not create physical subsystems. The grouping is bookkeeping unless a new relational law distinguishes units.

## Two-cell interaction family

Within the minimal symmetric quadratic A7 grammar, two-cell diagonal consistency and exchange symmetry still leave an unsourced dimensionless interaction ratio. This made the composition problem explicit rather than semantic.

## Incidence no-gos

For identical candidate units, permutation-equivariant incidence without pair-specific relational data yields only empty or complete graphs. Intrinsic state types yield at best block-complete graphs.

The dynamically generated conductance between alternative states of **one** unit cannot be silently reused as incidence among **simultaneous** units.

## Dynamic relations and conservation

A minimal binary relation model is dense unless system-size-dependent tuning is introduced. Fixing total edge number can enforce sparsity, but leaves edge identity exchangeable. Fixing local valence yields a random regular ensemble, not unique locality. Even conserved pair number/valence does not preserve neighbor identity under natural reversible swaps.

## Holonomy / edge-action tests

The existing FIN phase connection is flat and therefore cannot select incidence cycles. A cycle observable presupposes edges and is not itself an edge source. Minimal node-edge-holonomy actions did not jointly derive sparse persistent geometry.

## Pair identity

Quenched pair identity can preserve geometry only if it is declared. This moves, rather than solves, the source problem.

Main conclusion of this campaign: FIN needs a **pair-specific relational degree of freedom or transformation structure whose identity is dynamically sourced**, not merely a scalar global edge charge.

---

# 4. Campaign 178–192 — from two-port topology to information-preserving transformations

## Degree-two candidate

If every unit has degree two and the graph is connected, topology is uniquely a cycle `C_n`. This gives exact 1D graph spectrum/resistance conditionally, but fixed degree and connectedness are not yet sourced.

## Permutation transformation route

Exact deterministic information preservation on a finite set forces a bijection/permutation. A single invariant orbit of that permutation gives one cycle. This is a more principled route to two-port topology than imposing degree=2 directly.

The strict C12 operator can be written as a polynomial in one cyclic permutation; the strongest strict shell also recovers the internal C12 skeleton in the chosen basis.

However, a cycle of **states** or **history positions** is not automatically a cycle of simultaneous physical sites.

## Reversible dilation

Every stationary leave-one-out Markov chain admits an exact invertible natural extension on two-sided history space. Thus subsystem stochasticity is compatible with global information preservation.

The price is an explicit history-bearing environment whose information capacity grows with path history. This keeps the working idea “information is transferred rather than destroyed” mathematically viable without claiming a finite fundamental reservoir.

---

# 5. Campaign 193–203 — dimensionless time lane and copy-locality no-go

## Dimensionless history time

Three internal clocks were constructed and compared:

1. reversible history-shift depth;
2. Poisson refresh counts;
3. accumulated transition surprisal.

They agree asymptotically on the same dimensionless history depth. The memory-renormalized coarse law reproduces the microscopic dimensionless clock without fitting a separate coarse rate; the N=8 clock distortion is about `0.00358%`.

Continuous time is the exact Poissonization

`L_N = N*rho*(K_N-I)`.

The global positive `rho` remains a clock gauge. No SI second is derived.

## Arrow of time

At equilibrium the reversible path law has no thermodynamic arrow. Relative entropy to equilibrium decreases after a special nonequilibrium preparation. Therefore temporal ordering and the thermodynamic arrow are distinct.

## Copy incidence no-go

Two operational definitions — intervention influence and noncommutation of local refresh operators — both give a complete graph on microscopic copies. Pair influence is O(1/N), but there are O(N) partners. This is mean-field scaling, not finite-degree locality.

Thus internal copy number N must not be reinterpreted as the number of spatial sites.

---

# 6. Campaign 204–226 — CRT factorization, 12-state dynamics, quotient selection

## Exact internal factorization

Chinese remainder coordinates give

`Z12 ~= Z3 x Z4`

with exact shift factorization. The strict operator preserves pure Z3 and Z4 quotient subspaces.

The important transition shells have a clean typing:

- d=3: pure Z4/fiber move;
- d=4: pure Z3/base move;
- d=5: mixed move.

The d3+d4 graph is exactly `C3 square C4` in state space.

## Memory-sector correction

The three-sector Mori-Zwanzig memory source has only k=0,4,8 symmetry content. Earlier H4 hidden modes k=1,2 do **not** generate this particular coarse memory. H4 preparation/frame memory and metastable coarse-graining memory are distinct mechanisms.

At N=6 in the slow k=4,8 sector:

`A_k ~= -0.1138873027`

`M0_k ~= +0.0906339246`

`M1_k ~= 0.02805023885`

so M0 cancels about 79.6% of the naive instantaneous projected rate. The M0+M1 closure gives

`lambda_MZ ~= -0.02261891215`

versus exact

`lambda_exact ~= -0.02261127519`.

## 12-state Markov effective chain

A memory-aware D12-circulant Markov generator was reconstructed on all 12 localized basins. For N=3..8 all six shell rates are positive. At N=6 the Fourier eigenvalues agree with full microscopic slow modes below ~0.083% error.

## Exact Z3 effective intertwining

Once the 12-state effective generator is obtained, projection to `j mod 3` obeys the strong identity

`Q12 R = R Q3`,

hence

`exp(t Q12) R = R exp(t Q3)`

for all t. This is an exact process-level lumping at the effective level.

## Z3 vs Z4 correction

Z4 is also an accurate quotient, and for N=3..8 it has a slightly smaller global spectral gap. Therefore “Z3 is simply the slowest quotient” is false.

Z3 is preferred specifically as a **metastable identity carrier** because:

- its escape/conductance is lower;
- raw microscopic capacity also ranks it as more isolated;
- the barrier filtration selects exactly the three mod-3 components between B3 and B4.

## Barrier filtration

At the working gain:

`B3 = 0.6448515873278782`

`B4 = 0.6622191371274597`

`B5 = 0.7826547097738674`.

Below B3 the 12 minima are separate. At B3 they merge into exactly three C4 components `j mod 3`. At B4 those three components merge into one connected graph.

A critical correction: direct d5 saddle barrier B5 is not the d5 **communication height**; a lower multi-step route through d3/d4 reaches height only B4. Effective q_d shell rates must not be naively paired one-to-one with same-distance direct saddles.

---

# 7. Campaign 227–235 — metastable cores, large-N trend, controlled memory tail

## Deep cores

Broad committor cores such as q>1/2 still overcount transition-layer flux. A threshold such as q>0.95 may work at one N but is not universal and must not be promoted to a law.

A threshold-free construction was used: one deepest/highest-stationary-probability discrete state per localized minimum, four seeds per Z3 sector.

The ratio of deep-core capacity rate to independently derived effective exit rate improves monotonically for N=3..8:

`0.790, 0.833, 0.847, 0.897, 0.928, 0.956`.

This strongly supports the recrossing interpretation.

## Symmetry-reduced capacity to N=11

An exact D4 symmetry quotient extended the deep-core capacity calculation to N=11. New rates `3*cap`:

- N=9: `0.00247627654014`;
- N=10: `0.00135199815227`;
- N=11: `0.000726643165759`.

Local exponents rise toward B4 from below; latest steps are about `0.5858`, `0.6052`, `0.6209`. This materially strengthens B4 as the candidate large-N communication exponent, but does not prove it.

## Semigroup / initial slip

Pure rate-only Markov closure misses a small initial-slip amplitude. The same M1 predicts the pole residue

`Z_k = 1/(1+M1_k)`

without time-domain fitting. After the short memory layer, sampled errors fall below ~0.2%, then ~0.1% at later times.

## One-pole auxiliary memory

Matching M0,M1 with

`gamma = M0/M1`, `a = M0^2/M1`

and adding one auxiliary relaxation variable reproduces the N=6 projected transient from t=0 at sub-percent sampled error. Typical memory time is ~0.30–0.34 microscopic clock units.

## Stieltjes structure and tail bound

Reversibility gives

`K(t)=<c, exp(-Ht)c> = integral exp(-gamma t) dmu(gamma)`

with positive spectral measure. Thus K is completely monotone and its Laplace transform is Stieltjes.

Numerically resolved coupled hidden gaps at N=6 are roughly 0.677–0.825. This gives the rigorous form

`int_T^infty K(t) dt <= [K(0)/gamma] exp(-gamma T)`.

At T=8 the remaining tail is bounded below ~3% of M0 in every tested sector.

## P0 verdict

One controlled effective level is now functioning. Remaining work in this lane is proof-grade large-N asymptotics and uniform error bounds, not discovery of the basic architecture.

---

# 8. Campaign 236–244 — multiunit re-entry, coupling simplex, unique swap

## Two-unit coupling nonuniqueness

For two symmetric Z3 units with exact marginal Q3, the most symmetric reversible joint continuous-time law still has three increment classes A,B,C with

`A+B+C=k`.

Collective modes depend on A,B,C even though every one-unit marginal is identical. One-unit data therefore do not determine interaction.

## MaxEnt composition

Maximizing joint transition entropy at fixed marginals selects the independent product law / Kronecker sum. MaxEnt alone does not create interaction.

## Shared Z4 fiber

If two units are **declared** to share the same Z4 fiber event, existing q_d rates generate a parameter-free correlated Z3 law. At N=8:

`A/k ~= 0.2528`, `B/k ~= 0.3901`, `C/k ~= 0.3572`.

A relative collective mode is ~13–14% slower than the isolated Z3 gap. However global/shared-fiber and edge-local shared-fiber constructions do not create an extensive n^-2 hierarchy. Common noise is not enough for spatial scaling.

## Unique reversible reset dilation

For one system trit plus one uniform environment trit, requiring:

- bijection;
- exact reset;
- system-environment exchange symmetry;
- diagonal Z3 equivariance;
- involution;

leaves exactly one deterministic gate:

`F(x,e)=(e,x)` — SWAP.

Exhaustive 9! enumeration confirms one survivor.

## Closed swap cycle

Place retained Z3 units on a cycle and swap neighbors at rate `rho/2`. One-point density obeys the exact cycle heat equation. Fourier mode m has

`lambda_m = rho[1-cos(2*pi*m/n)]`,

so

`lambda_1 ~ 2*pi^2*rho/n^2`.

This is the first current multi-unit construction with:

- globally information-preserving elementary transformations;
- no new continuous interaction coefficient;
- a collective time scale growing with system size.

But cycle incidence remains conditional.

---

# 9. Campaign 245–255 — transformation geometry, simultaneous carrier, record factorization

## Commuting transformations -> torus geometry

If d commuting bijections act freely and transitively with finite orders n_mu, their elementary Cayley graph is

`C_n1 square ... square C_nd`.

With equal local activity budget the low-k spectrum is the d-dimensional discrete diffusion spectrum. Heat-kernel replay on L=64 gives effective spectral dimensions close to 1,2,3 for d=1,2,3.

## Generator-set nonuniqueness

The same group can admit different elementary generator sets and hence different graph dimension. Z12 with generator +1 is one cycle; with pure CRT generators +3,+4 it is `C3 square C4`. Group factorization alone does not source dimension.

## Barrier-selected generator basis

Conditioned on mapped d3,d4,d5 saddle families, minimum-bottleneck connectivity uniquely selects `{3,4}` because B4 < B5 and neither 3 nor 4 alone connects all 12 states.

Thus transition costs can select elementary directions in state space. This is a source strategy for dimension, not yet physical space.

## Simultaneous finite carrier

A finite cyclic register of L simultaneous Z3 slots with global cyclic permutation is bijective. If all environment slots are initially independent uniform trits, the observed slot receives exact fresh resets for L-1 update events and recurs at event L.

After Poissonization with tau=rho t:

`p0^(L)(tau)=exp(-tau) sum_q tau^(qL)/(qL)!`

and point-state total-variation error from the ideal infinite fresh bath is

`TV=(2/3)[p0^(L)(tau)-exp(-tau)]`.

The error is O(tau^L/L!) at small tau. This is the first explicit carrier whose slots coexist in one global state rather than merely representing history positions.

## One-orbit source candidate

For any permutation of L finite slots, the observed fresh-reset horizon equals the length of its permutation cycle minus one. Therefore one L-cycle uniquely maximizes exact fresh-reset duration for fixed environment capacity.

This gives an operational motivation for one orbit: maximal reversible use of finite environment memory.

## Minimal environment factorization

To produce m exact independent uniform trit resets from a closed deterministic environment E:

`H(E) >= m log 3`, hence `|E| >= 3^m`.

At minimum capacity `|E|=3^m`, the deterministic map from E to the m output trits must be bijective. Therefore

`E ~= Z3^m`.

So a multi-factor information carrier is forced at minimum capacity.

## Current typed boundary

The factors are naturally labelled by future reset records `(Y1,...,Ym)`. They coexist mathematically, but their canonical meaning is still temporal/record-like. Entropy, bijectivity and cyclic scheduling do not yet prove that these factors are simultaneous **physical spatial sites**.

This is the current sharp blocker.

---

# 10. Critical corrections that must survive into the next campaign

1. **Do not restore the unrestricted recursive binary pitchfork programme.** It failed on the representative localized branch.
2. **Do not call the parity bit exactly Markov.** It is not.
3. **Do not use `PLJ=Lc` as a process-level closure criterion.** Use intertwining, semigroup error, or explicit memory.
4. **Do not mix empirical-refresh and leave-one-out finite-N formulas without rederivation.**
5. **Do not discard short memory solely because it decays fast.** Its integrated self-energy is O(1) relative to PLP.
6. **Do not identify H4 preparation memory with the Z3 metastable Mori-Zwanzig memory.** The latter lives in k=4,8 (plus invariant k=0 sector).
7. **Do not say Z3 is uniquely the slowest quotient.** Z4 has a slightly smaller global gap over N=3..8. Z3 is preferred as a metastable identity carrier by escape/capacity/barrier filtration.
8. **Do not equate q5 with direct barrier B5.** d5 communication can use a lower B4 route.
9. **Do not claim B4 is already a proven large-N Eyring-Kramers exponent.** N<=11 supports it increasingly strongly, but proof is open.
10. **Do not reinterpret microscopic copy labels as spatial sites.** Their influence graph is complete mean-field.
11. **Do not reuse one-unit state-space conductance as inter-unit incidence.** These are different typed objects.
12. **Do not claim MaxEnt derives interactions.** At fixed marginals it selects independence.
13. **Do not claim shared Z4 fiber already yields space.** It gives correlation, not n^-2 hydrodynamic scaling.
14. **Do not claim the SWAP cycle is physical space.** The local gate is strongly sourced; the simultaneous incidence interpretation remains conditional.
15. **Do not claim dimension is a group invariant.** It depends on elementary generator selection and dynamical cost scale.
16. **Do not identify minimum-capacity record factors with physical spatial factors.** Record-to-spatial role transfer remains open.
17. **Do not claim an SI clock.** The global rate rho is still a gauge/calibration scale.
18. **Do not claim a fundamental arrow from detailed balance.** The thermodynamic arrow appears after nonequilibrium preparation.
19. **Do not claim QM, GR, particles or physical spacetime have been derived.** The strongest controlled physics remains statistical/effective dynamics plus conditional transport geometry.

---

# 11. Key numerical ledger

Working gain in the main metastability lane:

`g = 5.145228719489142`.

Mapped localized barriers:

- `B3 = 0.6448515873278782`
- `B4 = 0.6622191371274597`
- `B5 = 0.7826547097738674`

N=6 slow Z3 memory sector:

- `A_k ~= -0.1138873027`
- `M0_k ~= +0.0906339246`
- `M1_k ~= 0.02805023885`
- `lambda_MZ ~= -0.02261891215`
- `lambda_exact ~= -0.02261127519`

Deep-core capacity continuation:

- N=9: `3 cap ~= 0.00247627654014`
- N=10: `3 cap ~= 0.00135199815227`
- N=11: `3 cap ~= 0.000726643165759`

Latest local exponent trend toward B4:

`~0.5858, 0.6052, 0.6209`.

Dimensionless clock:

- exact Poissonization `L_N=N rho(K_N-I)`;
- coarse clock N=8 distortion ~`0.00358%`;
- absolute rho remains free.

Swap-cycle diffusion:

`lambda_m=rho[1-cos(2*pi*m/n)]`,

`lambda_1 ~ 2*pi^2 rho/n^2`.

Minimal perfect-reset environment:

`H(E)>=m log 3`, `|E|>=3^m`; equality forces `E ~= Z3^m`.

---

# 12. Current strongest scientific statement

FIN now contains a controlled finite statistical-mechanical example in which:

- a single declared microscopic reversible finite-N process;
- dynamically selected metastable sectors;
- short but quantitatively important memory;
- memory-renormalized effective kinetics;
- exact further lumping;
- global information-preserving dilation;
- and a conditional local reversible transport gate

can all be kept mutually consistent.

The main open fundamental problem has moved from “can dynamics be attached to the FIN potential?” to:

> **Can the information-record factors and reversible transformation structure already present in FIN be promoted, by a typed operational theorem rather than interpretation, into simultaneously existing local subsystems whose incidence and dimensionality are internally sourced?**

That is the correct research frontier after report 255.
