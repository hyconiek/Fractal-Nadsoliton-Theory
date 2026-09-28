# FIN — FULL HANDOFF FROM RESEARCH 327 TO CURRENT STATE

Date: 2026-09-28
Scope: tasks 327–336 plus the unfinished response-law work started after 336 (recorded as 337-WIP).

## 0. Executive state

This handoff begins after the external review of research 300–326 and therefore **does not silently inherit the strongest proof labels from 326**. The first obligation was to replay the global barrier certificate with explicit directed interval decisions. That work is task 327.

The two dominant lanes after 327 are now:

1. **barrier / metastability:** task 327 materially strengthens the 326 separator result to `PASS_DIRECTED_REPLAY`, but repository/external acceptance is not claimed; a complete Arb/MPFR replay of every box remains an optional formal-audit gate;
2. **preparation / observed-process:** tasks 328–336 replace most of the former empirical preparation transfer by exact finite-defect weights, identify the natural control variable `theta=kappa/N`, construct a realizable finite-time preparation mechanism, remove the hard basin oracle in the larger-N regime, freeze an operational control protocol, and precompute a defect-to-post-burn response library.

The current unfinished frontier is **337-WIP: LOW-DEFECT RESPONSE LAW**. At the frozen `theta=2`, only 78 configurations with `D<=2` already carry more than 99.988% of the preparation mass for N=7..10. The intended next step is to derive a cross-N response law for the D=0,1,2 sectors. This cluster/response law is **not yet completed**.

---

# 1. Provenance and review input

The user supplied an external review after research 300–326. Its main recommendations were:

- do not export `Gamma=B4` as an accepted global proof until the full exclusion chain of 326 is independently replayed with rigorous arithmetic;
- prioritize exact finite-N preparation theory rather than another empirical correction or immediate increase of N;
- distinguish fixed `kappa`, fixed per-defect field and fixed information budget;
- distinguish a target preparation distribution from an operational process that actually realizes it;
- retain rejected/unlocalized outcomes explicitly rather than silently conditioning them away.

The original review text is preserved in:

`04_REVIEW_INPUT/EXTERNAL_REVIEW_AFTER_300_326_20260928.txt`

and the response produced after tasks 327–330 is preserved in:

`04_REVIEW_INPUT/REVIEW_RESPONSE_327_330_20260928.md`.

---

# 2. Task-by-task research state

## 327 — AUDITED DIRECTED REPLAY OF 326

Status: **PASS_DIRECTED_REPLAY / external-repository acceptance pending**.

The original 326 used two kinds of global rejection:

- energy lower-bound rejection with an explicit conservative buffer;
- feasibility / fixed-point (`root_possible`) rejection whose full rounding propagation had not been surfaced.

327 replaced the global rejection logic by interval-valued decisions using directed interval arithmetic and independently replayed the tightest cases on a high-precision mathematical kernel.

Main six-dimensional half-wall replay:

- processed boxes: `3,402,039`;
- feasibility/fixed-point exclusions: `431,636`;
- energy exclusions: `1,265,741`;
- local d4 boxes: `3,643`;
- unresolved boxes: `0`;
- maximum depth: `60`;
- smallest positive feasibility margin: `2.236822820353475e-08`;
- smallest positive energy margin: `1.1141814885935524e-09`.

Independent 70-digit interval replay of the two tightest decisions on a separately reconstructed high-precision kernel gave much larger margins:

- feasibility: `>=4.0806488482909208e-4`;
- energy: `>=2.2346782835220542e-7`.

Triple-junction boundary `P0=P1=P2=1/3` was replayed separately in 5D:

- processed boxes: `59,271`;
- unresolved: `0`;
- minimum feasibility margin: `1.1752879963434983e-05`;
- minimum energy margin above d4: `3.6046418241732507e-06`.

The local d4 Krawczyk and Hessian were also replayed at high precision. On the declared radius-0.02 cube:

`lambda_min >= 0.05510885939961514 > 0`.

Preferred high-precision energies:

`V_localized = -0.80531946214230935605557878998...`

`V_d4       = -0.14310032501484691985100395454...`

hence

`B4 = 0.6622191371274624362045748354471982765...`.

The explicit minimum->d4->neighbor path was replayed with interval derivative signs, supporting `Gamma<=B4`; the exhaustive directed separator replay supports `Gamma>=B4` **at this audit standard**.

Current allowed statement:

`Gamma=B4` has **PASS_DIRECTED_REPLAY inside this packet**. Do not write that an external referee or repository maintainer has accepted it. If the project requires the strongest formal standard, replay the entire cover in Arb/MPFR with durable box logs.

Primary report:

`01_RESEARCH_327_330/fin327/AUDITED_DIRECTED_REPLAY_327.md`

## 328 — DISCRETE DEFECT PREPARATION

Status: **PASS for the declared pinning family on tested N=3..12; N=11,12 tails use certified intervals**.

Define

`D=N-n0`, `m=(m1,...,m11)`, `sum m_a=D`.

With

`delta_a=A00-A0a`

and

`B_ab=A_ab-A_a0-A_0b+A00`,

the exact relative pinning weight is

`w_{N,kappa}(m)/w_{N,kappa}(0)`

`= N! / ((N-D)! prod_a m_a!)`

`  * exp[-g sum_a delta_a m_a + (g/(2N)) m^T B m - (kappa/N)D]`.

This is an algebraic rewriting of the same finite-N Gibbs model, not a Gaussian approximation and not a new empirical fit.

The identity was checked against every available J=0 microstate for N=3..10 with log-weight replay residuals of order `1e-14`.

Because

`mu_{N,kappa}/mu_{N,0} proportional exp[-(kappa/N)D]`,

the tail `P_kappa(D>K)` is monotone nonincreasing in `kappa>=0`. Therefore the worst tail for positive pinning is at `kappa=0`.

The central truncation result:

- `D<=4`: 1,365 configurations, worst preparation-TV tail ~`0.9665%` over N=3..12;
- `D<=5`: 4,368 configurations, worst tail ~`0.3859%`;
- `D<=6`: 12,376 configurations, worst tail **<0.1603%**.

For any later linear Markov/readout channel K,

`TV(mu K, mu^(6) K) <= P(D>6)`.

Thus the preparation error is directly controlled without PCA.

Primary report:

`01_RESEARCH_327_330/fin328/DISCRETE_DEFECT_PREPARATION_328.md`

## 329 — CONTROL SCALING PROTOCOLS

Status: **exact identity + finite-N protocol comparison; no physical protocol selected**.

The exact conjugate defect field is

`theta = kappa/N`.

Hence fixed `kappa`, fixed `theta`, and fixed KL information cost are different experiments.

At the already-opened anchor N=6, kappa=12:

- `theta=2`;
- KL = `0.3880089643 bits`;
- `E[D]=0.04221843`.

For N=6..10, fixed kappa=12 visibly weakens with N, whereas fixed theta=2 keeps the defect count and information cost much more stable. A fixed-KL protocol over the same finite range also requires approximately `kappa proportional N`, with theta around 2.2–2.3 by N=8..10.

Do not export `theta=2` as a fundamental FIN constant. It becomes an operational protocol choice only in task333.

Primary report:

`01_RESEARCH_327_330/fin329/CONTROL_SCALING_PROTOCOLS_329.md`

## 330 — PREPARATION GENERATOR REALIZABILITY

Status: **PASS for exact unconditioned biased-Gibbs stationarity; PARTIAL for basin-conditioned realization; PASS for seed-conditioned operational mixing on N=3..10 at theta=2**.

The biased leave-one-out heat-bath

`q_j^prep = softmax_j[(g/N)A(n-e_i) + theta 1_{j=0}]`

satisfies exact detailed balance with the unrestricted biased Gibbs distribution.

Hard reflection at J=0 preserves the conditional weights on each connected component but the basin graph is generally disconnected; numbers of components grow from 1 at N=3 to 112 at N=10. Therefore:

`hard reflection + arbitrary initialization != exact full basin-conditioned target`.

For deep-seed initialization the component containing the seed has essentially all target mass under theta=2. Missing component mass falls from `1.47e-6` at N=4 to `4.86e-14` at N=10.

Worst-state spectral gaps are not the operational preparation time. N=9 has a slow pair with gap ~`0.0372457`, but the deep seed overlaps that pair only at about `-3.5e-8` and `6e-13`.

Direct seed-to-target propagation at theta=2 gives:

- `TV<1%` by about t=1.5–1.75 for all N=3..10;
- `TV<0.1%` by about t=4–4.25.

At N=9:

- t=1: `0.01571`;
- t=2: `0.00596`;
- t=4: `0.000871`;
- t=8: `1.88e-5`.

The remaining issue after 330 was that an exact reflecting wall still uses the declared basin label as an oracle. Task332 later removes that requirement operationally.

Primary report:

`01_RESEARCH_327_330/fin330/PREPARATION_GENERATOR_REALIZABILITY_330.md`

## 331 — DEFECT-BASED OBSERVED PROCESS PREDICTOR

Status: **PASS on already-opened N=7..12; not a new blind holdout and not yet a cheap cross-N post-burn law**.

331 inserted the exact D<=6 preparation into the observed-process pipeline without a new PCA correction.

Preparation contribution comparison:

| N | old empirical preparation component | D<=6 defect preparation component/bound |
|---:|---:|---:|
| 7 | 1.9420% | 0.00273% |
| 8 | 1.4296% | 0.01554% |
| 9 | 2.1691% | 0.05385% |
| 10 | 1.4707% | 0.08493% |
| 11 | 1.4710% | <=0.16028% |
| 12 | 1.8339% | <=0.14310% |

Full microscopic two-time process changes caused only by D<=6 truncation:

- N=7: `2.72e-5 TV`;
- N=8: `1.55e-4 TV`;
- N=9: `5.39e-4 TV`;
- N=10: `8.50e-4 TV`.

With an explicit `U=unlocalized` outcome and the same previously frozen effective generator, maximum two-time joint error versus the microscopic process is:

- N=7: **1.3747%**;
- N=8: **0.8687%**;
- N=9: **0.5915%**;
- N=10: **0.9628%**.

Scientific interpretation: the former 1–2% preparation residual was mostly the cost of an empirical cross-N preparation map, not irreducible hidden memory.

Primary report:

`02_RESEARCH_331_336/fin331/DEFECT_BASED_OBSERVED_PROCESS_331.md`

## 332 — PREPARATION CONTROLLER WITHOUT BASIN ORACLE

Status: **PASS as an operational soft-confinement mechanism on the already-opened finite-N range; low-N escape remains explicit**.

Use the unrestricted biased heat-bath at theta=2. It does not know the basin label.

Equilibrium mass in desired J=0 phase:

- N=3: 95.3170%;
- N=4: 99.1376%;
- N=5: 99.8868%;
- N=6: 99.9808%;
- N=7: 99.99756%;
- N=8: 99.99957%;
- N=9: 99.999946%;
- N=10: 99.999990%.

A killed-process calculation from the deep seed gives the conservative preparation bound at `Tprep=4`:

| N | escape + survivor-TV bound |
|---:|---:|
| 3 | 5.8439% |
| 4 | 1.7500% |
| 5 | 0.3188% |
| 6 | 0.1494% |
| 7 | 0.09489% |
| 8 | 0.08837% |
| 9 | 0.08727% |
| 10 | 0.08871% |

Thus for N>=7 a hard basin oracle is operationally unnecessary if escaped/other/unlocalized outcomes are retained rather than conditioned away.

Primary report:

`02_RESEARCH_331_336/fin332/PREPARATION_CONTROLLER_WITHOUT_BASIN_ORACLE_332.md`

## 333 — FROZEN CONTROL PROTOCOL

Status: **FROZEN BEFORE A GENUINELY NEW HOLDOUT**.

The protocol was frozen from already-opened N<=10 information:

1. initial state: deep seed `n0=N`;
2. unrestricted biased leave-one-out heat-bath at the production g;
3. `theta=2`, equivalently `kappa_N=2N`;
4. preparation duration `Tprep=4`;
5. switch field off after Tprep;
6. continue under the unmodified production generator;
7. no hard J=0 wall;
8. no postselection;
9. retain all localized plus explicit unlocalized/other-phase outcomes.

No future holdout may retune theta or Tprep and still be called a test of this protocol.

### Important checksum note

The JSON contains an internal field

`frozen_sha256 = b7ae70e9...10d3ff`.

The **actual SHA-256 of the final JSON file currently packaged here is**

`bdd761a28fcd0295bac19ce4ba2b3c9e4439aa443e3e27910de33d937ac3d42f`.

The internal value should therefore be treated as a historical freeze identifier / pre-embedding checksum, **not as the checksum of the final serialized file**. For integrity checking use the actual manifest checksum.

Primary artifacts:

`02_RESEARCH_331_336/fin333/FROZEN_CONTROL_PROTOCOL_333.md`

`02_RESEARCH_331_336/fin333/FROZEN_CONTROL_PROTOCOL_333.json`

## 334 — FORMAL MPFR/ARB REPLAY

Status: **NOT EXECUTED**.

This number was reserved for a full MPFR/Arb replay of task327 if repository-level proof acceptance requires a stronger backend and durable complete box logs. Task327 already materially strengthened 326, so 334 remains an optional formal-audit lane rather than a completed research task.

## 335 — FROZEN-CONTROL NEW-g HOLDOUT

Status: **NOT EXECUTED / DELIBERATELY POSTPONED**.

The new-g holdout was not opened because a full cross-g dynamics/response predictor had not yet been frozen. Opening a new g before freezing that law would confound prediction with architecture selection. The frozen task333 control remains untouched and available for a future genuinely new holdout.

## 336 — DEFECT-TO-POST-BURN REDUCED MAP

Status: **PASS for N=7..10 at burn time t=24**.

For each retained defect state x with D<=6, task336 stores

`R_x(a)=P_x[macro outcome a at t_burn=24]`,

where a covers twelve localized labels plus one explicit `unlocalized` outcome.

By linearity,

`p_burn(a)=sum_x w_x R_x(a)`.

Library sizes:

| N | D<=6 J=0 states | full count-state size |
|---:|---:|---:|
| 7 | 2,596 | 31,824 |
| 8 | 5,648 | 75,582 |
| 9 | 9,351 | 167,960 |
| 10 | 11,552 | 352,716 |

At the frozen theta=2, omitted D>6 tail is essentially negligible:

- N=7: `5.13e-11`;
- N=8: `1.89e-10`;
- N=9: `9.06e-10`;
- N=10: `1.36e-9`.

Thus preparation now factorizes as

`control -> exact defect weights -> fixed response library -> post-burn prior`.

The remaining limitation is that the response library is still N-specific and generated from the N-specific microscopic dynamics.

Primary report:

`02_RESEARCH_331_336/fin336/DEFECT_TO_POST_BURN_REDUCED_MAP_336.md`

---

# 3. Current unfinished work — 337-WIP LOW-DEFECT RESPONSE LAW

This work was started after task336 but has **not** been completed into a PASS/FAIL report.

The first exact reduction is already established from the task336 libraries and the frozen task333 protocol.

For theta=2, the entire `D<=2` sector contains exactly

`1 + 11 + C(12,2) = 78`

configuration types (one D=0 state, eleven D=1 states, and 66 weak compositions of D=2 over eleven defect types).

Using the exact task336 weights and adding back the certified D>6 tail, the full preparation mass retained by those 78 states is:

| N | certified mass D<=2 | certified tail D>2 |
|---:|---:|---:|
| 7 | 99.98806636% | 0.01193364% |
| 8 | 99.99094605% | 0.00905395% |
| 9 | 99.99233793% | 0.00766207% |
| 10 | 99.99307156% | 0.00692844% |

Therefore, for the frozen theta=2 protocol on N=7..10, replacing the 2,596–11,552-row D<=6 response libraries by the 78 rows with D<=2 has an immediate process-level TV error bound equal to the corresponding D>2 tail above.

This is a **WIP result**, not yet a new task verdict.

The intended next calculation is a cluster/response decomposition:

- D=0 response `R0`;
- one-defect increments `Delta_a`;
- two-defect interaction residuals `Delta_ab`;
- exploit D12/C12 symmetry to classify defect types;
- test whether appropriately normalized increments obey a transferable law in N (and later g).

At fixed N this decomposition can be made algebraically exact on D<=2. What is still open is the **cross-N/cross-g law** for those increments. No claim of additivity or universal cluster scaling has yet been accepted.

Reproducibility artifacts:

`03_WIP_337/low_defect_theta2_summary.py`

`03_WIP_337/LOW_DEFECT_THETA2_SUMMARY_337_WIP.json`

---

# 4. Current scientific picture

## 4.1 Preparation is no longer primarily an empirical nuisance parameter

The preparation lane now has the chain

`control theta`

`-> exact biased Gibbs / exact defect weights`

`-> finite few-defect representation`

`-> finite-time operational preparation`

`-> explicit localized + unlocalized outcomes`

`-> observed process`.

This is materially stronger than the earlier PCA/initial-slip cross-N preparation map.

## 4.2 The dominant open problem has moved again

After tasks328–336, the principal preparation question is not "how do we fit the prior at each N?". It is now:

**Can the post-burn response of one- and two-defect configurations be derived from a small N/g-transferable mechanical law?**

If yes, the current N-specific response library disappears and the frozen protocol333 can be used for a genuinely new-g holdout without retuning.

## 4.3 Barrier lane remains separate

Task327 strongly supports the 326 global barrier result and repairs the exact proof-status issue raised in the review. Nevertheless, do not let the barrier lane contaminate the preparation lane:

- tasks328–336 do not require `Gamma=B4` to be accepted by an external reviewer;
- a future Eyring-Kramers prefactor does require a clearly stated barrier-proof standard;
- the response-law/new-g preparation test can continue independently.

---

# 5. Epistemic boundaries / nonconclusions

Do **not** claim any of the following from tasks327–337-WIP:

1. FIN derives a physical actuator producing theta=2. `theta=2` is a frozen operational model choice selected from already-opened finite-N data.
2. FIN derives a dimensional physical preparation time. `Tprep=4` is in the declared generator-time units.
3. D<=6 or D<=2 is a universal N->infinity representation. It is a controlled finite-N representation in the presently studied sparse-defect regime.
4. Task331 is a new blind holdout. It reuses already-opened N=7..12 microscopic information.
5. Task332 proves zero escape. It keeps escaped/other outcomes explicit and provides finite bounds.
6. Task336 supplies a universal response law. Its response libraries are N-specific.
7. Task337-WIP has proved cluster additivity. It has only established the extreme D<=2 mass concentration under the frozen protocol.
8. Task327 by itself means a repository maintainer or external referee has accepted `Gamma=B4`. Its internal status is `PASS_DIRECTED_REPLAY`; full Arb/MPFR replay remains available if required.
9. The internal `frozen_sha256` string inside task333 JSON is the current file checksum. It is not; use the manifest checksum.
10. Any of these results derive QM, GR, the Standard Model, spatial incidence, SI units, a physical source of g, or a theory of everything.

---

# 6. Recommended next research program

## P0 — complete 337 LOW-DEFECT RESPONSE LAW

Use only the 78 D<=2 states under frozen theta=2 as the operational preparation basis for N=7..10.

1. construct exact per-N cluster coordinates `R0`, `Delta_a`, `Delta_ab`;
2. quotient them by D12 symmetry;
3. identify the minimal set of inequivalent one- and two-defect response classes;
4. test scaling variables suggested by the exact microscopic generator, not PCA;
5. report direct reconstruction error and retain the certified D>2 tail separately.

Kill-test: if cross-N residuals remain order-percent after accounting for the exact D>2 tail, do not promote a cluster law.

## P0/P1 — response-law holdout

After a response law is frozen, use a state size or gain not used to choose its functional form. Because N=7..10 libraries were already inspected during WIP, they are no longer blind for model selection.

Possible clean options:

- build only the 78 D<=2 responses for N=11 or N=12 and treat them as held-out response-law tests (microscopic N is already opened historically, so this is a response-law holdout, not a new-microscopic holdout);
- preferably freeze a cross-g law and then execute the still-unopened task335 new-g test under protocol333.

## P0 — 335 NEW-g HOLDOUT only after full predictor freeze

Before opening a new gain:

- keep task333 `theta=2`, `Tprep=4`, deep seed, unrestricted preparation, no hard wall, no postselection;
- freeze the response-law dependence on g;
- freeze the production transition law at that g;
- define all outcomes including unlocalized/other phases;
- write hashes before microscopic opening.

No post-hoc change of theta, Tprep, observable definition, response-law basis, or clock law may count as the same held-out test.

## FORMAL AUDIT — optional task334

If repository acceptance of `Gamma=B4` requires it, port the complete 327 half-wall and triple-junction cover to Arb/MPFR with durable complete box logs and independent replay hashes.

## P1 — Eyring-Kramers prefactor

Proceed only with proof status stated explicitly. Do not infer the asymptotic prefactor from a fitted `N^alpha`. Derive it from the actual leave-one-out kinetics, local minimum/saddle geometry, gate multiplicity, valley normalization and discrete count-state lattice.

---

# 7. Minimal continuation instruction for the next agent

Start from task336 and `03_WIP_337`.

**Do not** restart preparation fitting, PCA correction, N=13 generation, or Eyring-Kramers fitting.

First complete the D<=2 cluster-response law. Preserve:

- task333 frozen control `theta=2`, `Tprep=4`;
- no hard basin wall;
- no postselection;
- explicit unlocalized/other outcomes;
- exact defect weights;
- explicit truncation tail;
- separation between already-opened retrospective validation and a genuinely new held-out test.

If the cluster law cannot be made transferable, the correct conclusion is that task336's N-specific response library is still required. Do not hide that failure in another empirical prior correction.
