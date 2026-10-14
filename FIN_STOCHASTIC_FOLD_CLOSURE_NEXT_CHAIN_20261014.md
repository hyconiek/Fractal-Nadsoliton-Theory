# FIN next intrinsic chain: close the microscopic stochastic-fold bridge

Date: 2026-10-14. Repository reference: `ec8baa15` plus the TH intake.
Status: **RESEARCH_PLAN; SF tasks below have not been executed.**
Execution directory: `fin_stochastic_fold_closure_campaign_20261014/`.
Language: English for the execution agent.
Architect decision: [TH-010](fin_intrinsic_theory_review_20261014/TH_010_DECISION.md).

## 1. One scientific objective

Prove a prepared, compact-stopped, process-level critical limit of the already
specified leave-one-out FIN chain; then derive one dimensionless first-passage
law from that same limit. Do not introduce a bath, diffusion coefficient,
clock, graph, fitted memory or new microscopic kinetic rule.

The incoming TH package pays the interacting path-LDP application and gives
microscopic soft coefficients. It explicitly does **not** prove critical
process convergence. This plan replaces its vague "prepared local sequences"
and fixed-original-neighbourhood stopping with a precise **target theorem**.
The restriction is a theorem scope, not an added physical postulate.

There is a promising direct solution: use the full eleven-dimensional
occupation tangent, remove its one known null direction, and prove a stopped
quadratic energy estimate on the ten stable directions. Do not default to a
generic singular-perturbation theorem whose jump/process hypotheses are unpaid.

## 2. Frozen mathematical parent

Use q=12 and A=XX^T with the established exact-decimal kernel parameters and
mathematical circle constant. The exact count generator is

    L_N f(p)=N sum_ij p_i softmax_j[g_N A(p-e_i/N)]
                    [f(p+(e_j-e_i)/N)-f(p)].

Each labelled copy has unit attempt rate. Population normalization remains N.
Do not use empirical refresh for finite-N corrections or auxiliary-Gaussian
Langevin paths for the actual retained process.

Let `(theta_f,g_f,v)` denote the **exact isolated fold root and normalized
null vector** licensed by R7P-031, with source-fixed sign. A decimal midpoint
is a proposal, not an exact fold. Holding its error fixed while N grows changes
`N^(2/3)(g-g_f)` and cannot prove the critical theorem. Use the implicitly
defined exact root in analysis; pay amplified enclosure errors for numerics.

Retain the actual observable `theta_N=g_N X^T p_N`, not a sampled auxiliary
field. Known coefficients are `a=v^T partial_g grad Phi` and `b=D^3 Phi[v,v,v]`.

## 3. Exact input registry

| Alias | Input |
|---|---|
| S01 | AGENTS.md, especially the current TH intake section |
| S02 | fin_intrinsic_theory_review_20261014/INTAKE_20261014.md, CLAIM_REGISTER.json, verification.json |
| S03 | fin_intrinsic_theory_review_20261014/FAST_MODE_PROOF_ROUTE.md |
| S04 | fin_intrinsic_theory_campaign_20261014/TH-002/GENERATOR_TO_ACTION_PROOF.md and LDP_ASSUMPTIONS.md |
| S05 | fin_intrinsic_theory_campaign_20261014/TH-003/OBSERVABLE_SUFFICIENCY.md and PATH_CONTRACTION.md |
| S06 | fin_intrinsic_theory_campaign_20261014/TH-004/SOFT_MODE_MICROSCOPIC_BRIDGE.md and TH-005/STOCHASTIC_FOLD_LIMIT.md |
| S07 | fin_rank7_followup/certificates/R7P-031_simple_fold.json; src/fold_certificate.py, src/coexistence_certificate.py; inputs/fin_handoff_audit/results.json |
| S08 | fin_rank7_mathphysics_next/proofs/MP7-025_simple_fold_normal_form.md and MP7-026_controlled_fold_scaling.md; corresponding result objects |
| S09 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py |
| S10 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/RECURSIVE_LOCAL_CHILD_TEST_131.md and MICROSCOPIC_PROCESS_CONTRACT_138.md |
| S11 | fin_theory_synthesis_20261014/CORE_DERIVED_BRIDGES.md and LITERATURE_REGISTER.json |
| S12 | FIN_POST_F_EVIDENCE_ACTIVATION_CHAIN_20261014.md, for unchanged empirical restrictions |

Hash the actual files you use. The intake replaces only its explicitly
qualified claims, not historical source texts. Source scripts/plots are
diagnostics and input proposals, not proofs of the remaining estimates.

Primary benchmark: [Kraaij's Theorem 2.3](https://arxiv.org/html/1502.05185v3#S2.Thmtheorem3)
licenses the admitted fixed-g path-LDP. It does not replace the growing-horizon
critical weak-limit proof. Cuff's Potts result is a scaling benchmark, not a
FIN mixing theorem. Read any additional theorem fully before applying it.

## 4. Target theorem contract

Fix finite `S>0`, `R>0`, bounded lambda and initial limiting support inside
`(-R,R)`. Define

    g_N=g_f+lambda N^(-2/3),
    Z_N(s)=N^(1/3) v^T[theta_N(N^(1/3)s)-theta_f].

On the occupation tangent define

    Sigma_f=diag(p_f)-p_f p_f^T,
    r_f=Sigma_f Xv,
    ell_f^T=g_f v^T X^T,
    Q_f=I-r_f ell_f^T,
    w_N=Q_f(p_N-p_f).

One sufficiently narrow prepared initial class to prove first is

    Z_N(0) => Z_0,
    E||w_N(0)||^2 <= C/N,
    initial support/interior and tail conditions explicitly frozen.

Start with deterministic bounded z0 and rounded actual type states if that
helps, then state precisely any broader initial law extension. Full Gibbs
equilibrium at g_f is not automatically prepared near this fold; conditioning
or a different initial ensemble must be stated, not silently substituted.

Stop at the soft boundary `|Z_N|>=R`, auxiliary fast boundary
`||w_N||>=N^(-1/4)`, or scaled horizon S. Prove the auxiliary exit probability
vanishes. The desired limiting process is the solution, killed/stopped under
the frozen convention, of

    dZ=-g_f(a lambda+b Z^2/2)ds+sqrt(2g_f)dB

on `(-R,R)`. Neither a globally normalizable cubic-potential equilibrium nor
full escape through distant basins is part of this theorem.

## 5. Task cards

### SF-001 — Exact full-state fast chart and quantitative initial/stop contract

**Scientific question:** does the actual full occupation process have exactly
one soft direction and a usable prepared/stopped scope without new physics?
**Why it matters:** source "prepared" and static retained Hessian stability
alone do not pay stochastic elimination on a growing time horizon.
**Exact inputs:** S01–S11, especially actual root/spectral certificate S07.
**Dependencies:** admitted TH intake; no new experiment or parameter sweep.
**Method:**

1. Bind the exact kernel, root/vector, count generator, clock and finite-N
   convention. Reuse the intake's matching input result where hashes agree.
2. Prove `Jr_f=0`, `ell_f^T J=0`, `ell_f^T r_f=1`, with
   `J=g_f Sigma_f A-I` on the normalized occupation tangent.
3. Prove the ten-dimensional fast restriction is Hurwitz from the triangular
   retained/hidden structure. Declare the chosen basis, norm and gap premise.
4. Construct a derived positive quadratic form P via `B^T P+PB=-I` or an
   explicit equivalent. Numerical P is a proposal if certified bounds are claimed.
5. Freeze g_N, the actual theta/Z coordinate, lambda/R/S domain, initial law
   and stopping convention. Prove the prepared class contains actual type
   sequences, e.g. balanced rounding of `p_f+N^(-1/3)z0 r_f`.
6. Give the coordinate conversion between `ell_f^T(p-p_f)` and actual
   `v^T(theta_N-theta_f)`, including gain-dependent shifts.
7. Distinguish exact symbolic constants, reused interval premises and nominal
   numbers. Do not declare a uniform ten-mode gap from a wrong chart's bound.

**Deliverables:** MODEL_CONTRACT.json, INITIAL_AND_STOP_CONTRACT.json,
FAST_CHART_PROOF.md, GAP_AND_LYAPUNOV_OBJECT.json, INPUTS.sha256.
**Acceptance criterion:** a complete consistent contract, exact chart identities
and stable restriction; a demonstrated nonempty quantitative preparation class.
**Failure/refutation criterion:** an extra nonstable direction, unmatched input,
infeasible preparation or new microscopic assumption is necessary.
**What becomes possible:** generator occupation/containment estimates.
**Must NOT conclude:** a uniform finite-N spectral gap, physical clock or geometry.
**Estimated computational class:** S/M, finite linear algebra and proof.
**Priority:** P0.

### SF-002 — Pay the stopped fast-energy and occupation estimates

**Scientific question:** can actual generator fluctuations remain harmless
over model time `N^(1/3)S` under SF-001's preparation?
**Why it matters:** deleting ten stable means by intuition misses noise,
initial layers and the covariance dependence identified in TH-003.
**Exact inputs:** SF-001, S03/S06/S07/S10 and exact jump/covariance formulas.
**Dependencies:** SF-001 accepted; no fitted macro stochastic process.
**Method:**

1. Derive the local fast expansion with bounded remainders from L_N, retaining
   LOO O(1/N), gain and normalization terms. Leading drift closure does not
   authorize changing finite-N kinetics.
2. Prove a local generator inequality for `V=w^T Pw`, of sufficient form

       L_N V <= -kappa V+C(x^4+|g_N-g_f|^2+1/N),

   where x is the exact chart soft coordinate. All constants must be N-independent
   on the stated frozen domain. Include actual quadratic jump increments.
3. Use stopped Dynkin/nonnegative-energy estimates, with initial and jump
   overshoot payments. Do not assume exponential decay of a process frozen
   at a nonzero exit boundary; distinguish killed occupation and stopped value.
4. Target `E integral_0^(S wedge tau_N) ||w_N(N^(1/3)s)||^2 ds=O(1/N)`.
   Prove or replace it by an explicitly adequate bound with derived exponents.
5. Prove negligible auxiliary tube exit, e.g. a bound of order N^(-1/6)
   from a N^(-1/4) tube and cumulative energy O(N^(-2/3)).
6. Export actual integrated soft-drift residual bounds. Target the mixed term
   O(N^(-1/6)), quadratic-fast term O(N^(-1/3)), and vanishing LOO correction.
7. Establish full occupation convergence to p_f on stopped intervals, sufficient
   for its predictable covariance; leading mean drift alone is inadequate.

**Deliverables:** GENERATOR_ENERGY_PROOF.md, OCCUPATION_AND_EXIT_BOUNDS.md,
REMAINDER_LEDGER.json, constants/domain certificate if numerics are certified.
**Acceptance criterion:** every needed remainder/exit bound follows from the
same microscopic generator, initial law and stopped domain; no unpaid o(1).
**Failure/refutation criterion:** a required integrated term survives, tube
exit is nonnegligible or only a stationary/finite-sample estimate is available.
**What becomes possible:** tightness and identification of the soft process.
**Must NOT conclude:** hidden modes are absent at every N/preparation, a new
information-conservation law or a physical bath has been derived.
**Estimated computational class:** M, analytic proof; small certificate checks only.
**Priority:** P0.

### SF-003 — Prove the compact-stopped martingale-problem convergence

**Scientific question:** is the candidate SDE actually the path limit rather
than a balance of generator characteristics?
**Why it matters:** this is the missing micro-to-mesoscopic bridge, not another
plot of drift versus N.
**Exact inputs:** SF-001/002, accepted a/b/g_f/vector and covariance identities.
**Dependencies:** paid estimates in SF-002. A named lemma alone is insufficient.
**Method:**

1. Write the actual compensated soft martingale and predictable bracket under
   the accelerated clock; prove jump size and third-characteristic bounds.
2. Prove integrated drift convergence to `-g_f(a lambda+b Z^2/2)` and bracket
   convergence to `2g_f`, including hidden occupation and initial-layer scope.
3. Prove stopped tightness in path space using actual characteristics and
   bounded-jump controls. State topology, boundary absorption and any cemetery state.
4. Identify each subsequential limit through tests with controlled Taylor
   errors; prove uniqueness for the bounded-domain one-dimensional problem.
5. Show the auxiliary fast stop disappears in probability and the soft
   hitting/absorption map is licensed by the nondegenerate limiting diffusion.
6. Assemble a self-contained theorem with initial/domain/parameter quantifiers.
   Do not invoke fixed-time LDP as a substitute for weak convergence.
7. Have a small independently runnable checker for any arithmetic certificate;
   retain analytical proof obligations separately from code PASS.

**Deliverables:** STOPPED_CONVERGENCE_THEOREM.md, MARTINGALE_PROOF.md,
BOUNDARY_AND_TOPOLOGY.md, PROOF_DEPENDENCIES.json and claim register updates.
**Acceptance criterion:** a full process convergence theorem under SF-001's
paid contract, or a precisely located failed estimate. No template placeholders.
**Failure/refutation criterion:** only finite-dimensional/generator convergence,
an unpaid fast limit, hidden extra stop or global assertion is supplied.
**What becomes possible:** a limiting first-passage/survival prediction.
**Must NOT conclude:** global metastable mixing/escape, Gamma=B4, stationary
cubic Gibbs density, spatial field propagation or a physical time scale.
**Estimated computational class:** M, primarily proof.
**Priority:** P0.

### SF-004 — One dimensionless first-passage prediction, no fitted law

**Scientific question:** what complete bounded-domain outcome law follows
from the now-proved microscopic limit?
**Why it matters:** a useful effective theory should predict a trajectory
observable, not merely write a canonical equation.
**Exact inputs:** accepted SF-003 theorem and the frozen initial/boundary contract.
**Dependencies:** SF-003 theorem accepted. Otherwise mark NOT_ACTIVATED.
**Method:**

1. Choose before computation one initial point/law, lambda and interval inside
   the theorem domain. Keep both exit sides and survival at finite horizon.
2. Use `c=(4/b)^(1/3)`, `d=c^2/(2g_f)`, `Y=Z/c`, `tau=s/d` to derive
   `dY=(kappa lambda-Y^2)d tau+dW`; transform preparation and boundaries too.
3. Derive the backward killed-generator equation for an exit-side probability
   or Laplace transform. For example `G u=q u` with specified absorbing values.
   Coefficients, time and boundary constants cannot be trajectory-fitted.
4. Provide an analytic integral/ODE representation or a numerical approximation
   with a disclosed/certified deterministic error; do not use Monte Carlo alone
   to define the prediction. Do not report a finite-N error rate unless paid.
5. Use weak convergence to transfer only observables whose continuity/integrability
   has been proved. Bounded finite-horizon laws do not automatically imply
   convergence of an unbounded mean escape time.
6. Optional diagnostics: one deterministic lattice-generator check from existing
   small inputs, at most a few fixed N. No high-N sweep is needed for acceptance.

**Deliverables:** FIRST_PASSAGE_PREDICTION.md, PREDICTION_CONTRACT.json,
DIMENSIONLESS_LAW.json and numerical-error scope where relevant.
**Acceptance criterion:** one falsifiable, source-derived, bounded-domain
mathematical prediction with complete initial/boundary/unit conventions.
**Failure/refutation criterion:** it needs fitted coefficients, unconstrained
global reinjection or unproved hitting/expectation convergence.
**What becomes possible:** later controlled simple-fold universality comparison
or matched-regime theory, after a new architect selection.
**Must NOT conclude:** empirical validation, universal constants of nature,
Potts mixing, physical particles, geometry or quantum/gravity structure.
**Estimated computational class:** S/M.
**Priority:** P1, strictly conditional on SF-003.

### SF-005 — Evidence handoff and single-branch checkpoint

**Scientific question:** which exact micro-to-mesoscopic statement is now
licensed, and what if anything remains to be paid?
**Why it matters:** process convergence must not be buried under completed
file/status counts or expanded into unrelated optional campaigns.
**Exact inputs:** all executed SF artifacts, S01/S02/S12 and unchanged source packages.
**Dependencies:** terminal scoped dispositions for SF-001–003; SF-004 active or not.
**Method:**

1. Build a dependency DAG: source kernel/kinetics -> fold -> preparation ->
   occupation bounds -> process theorem -> exit-law prediction.
2. Separate exact analytic, theorem-application, interval-assisted, finite
   numerical, conditional and OPEN statements. Mark each reused result.
3. Export source hashes, proof constants, failures, replay and checks on a
   temporary copy. Preserve historical files and incoming caches separately.
4. Answer whether ten-mode elimination and the stop were paid, whether a
   uniform finite-N rate was proved, and whether the prediction requires any
   new microscopic premise. Do not let a failed SF-003 activate SF-004.
5. Stop for the architect. Do not choose a graph, source, universality family,
   broader preparation, experiment or TH-006–009 task on your own.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
REPLAY.md, MANIFEST.sha256 and an explicit NEXT_OPEN_ATOM.md if blocked.
**Acceptance criterion:** a small complete replayable package with one honest
terminal theorem/open verdict; every numerical/proof input is bound.
**Failure/refutation criterion:** unchanged open conditions are hidden behind
PASS, historical results relabelled new or new physical premises silently added.
**What becomes possible:** one architect-selected next theoretical step.
**Must NOT conclude:** independent physical evidence, source selection or ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution limits and reporting

- This is a proof campaign. No heavy state enumeration, broad N/g sweep,
  static-digit campaign, new source inventory or external experiment.
- Preserve AGENTS.md, synthesis and prior handoffs. New work only in the
  declared execution directory; return a proposed accepted-claim patch.
- For every executed task export REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
  NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256.
- RESULTS.json must include execution_status, scientific_verdict, supplied/
  derived/fitted objects, proof_scope, arithmetic_status, initial/stop scope,
  dependencies, discrepancies and next_action. No fitted objects are expected.
- An asymptotic proof need not give a numerical finite-N confidence bound;
  do not claim one from schematic O notation or an unvalidated constant.
- Numerical fixed fold centers/vectors require an explicit input-error bridge.
  Existing enclosures license an exact symbolic root, not arbitrary midpoint laws.
- No QW-2191, role transfer, strict/legacy completion, B4, task335, FIN-M0,
  locality/QM/GR/SM, physical calibration or empirical promotion is activated.

## 7. First actions and kill gates

Start **SF-001**, then **SF-002**, then **SF-003**. These are the first three
tasks; no separate repeated TH-001 source inventory is needed.

If the full fast chart has another soft direction, stop the one-dimensional
candidate. If the generator-energy bound is insufficient, return its exact
term/exponent or a rigorously narrower preparation rather than adding noise.
If containment or boundary convergence is unpaid, retain the formal candidate
and path-LDP result, not a stopped theorem. If the proof succeeds, do SF-004
once and return SF-005. No automatic infinite sequence follows.

The expected advance is a derived collective **statistical** law. Its value
for fundamental physics remains conditional on distinct, presently open
kernel/kinetic-source, locality, conservation, units and observational bridges.
