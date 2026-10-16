# FIN after NT: exact auxiliary updates versus a persistent dynamical mediator

Date:2026-10-16. Baseline:3be481a7 plus the NT intake.
Status: **RESEARCH_PLAN; HM-001–HM-005 have not been executed.**
Decision: [mediator architecture](fin_native_transfer_review_20261016/HM_DECISION.md).
Execution directory: fin_auxiliary_mediator_campaign_20261016/.
Language: English for the execution agent.

## 1. One scientific objective

Determine exactly which auxiliary-field constructions preserve the supplied
native FIN dynamics, and which introduce a genuinely different kinetic law.

NT has sufficiently closed the present classical temporal transfer target.
Its Hilbert/critical operator representation is not a physical mediator law.
The important next question is

> Does the exact Hubbard–Stratonovich equilibrium representation license a
> persistent field with feedback, or only an update-level latent variable?

This is executable theory, not another missing-input inventory or sextic-gap
precision campaign. Use the existing MP7 Gaussian joint law and original
LOO heat bath. Diagnostic comparison processes below have separate IDs and
are **not adopted as FIN laws**. Their extra kappa is a supplied comparison
parameter, never a derived constant or physical clock. Any proposal to adopt
one as fundamental requires explicit human approval under a new model ID.

## 2. Exact inputs and reused results

| Alias | Input |
|---|---|
| H01 | Current AGENTS.md and fin_native_transfer_review_20261016/INTAKE_20261016.md |
| H02 | The same review's PROOF_COMPLETIONS_AND_SCOPE.md, MEDIATOR_BRIDGE_SEED.md and verification.json |
| H03 | fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md and MP7-032_gaussian_auxiliary_representation.md |
| H04 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py and accepted exact spectral intervals |
| H05 | fin_generative_memory_review_20261014/INTAKE_20261014.md; existing reversible projection/memory identity |
| H06 | fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md and local-domain contract |
| H07 | fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md §2; UF actual-generator methods |
| H08 | fin_kinetic_invariants_review_20261014/INTAKE_20261014.md and accepted mobility/noise scope |
| H09 | P530 in FIN_Programs_527_536_Auxiliary_Field_Memory_and_Exact_Replay_Report_EN.tex, as prior nonuniqueness context |

Hash actual files. Do not claim source-wrapper custody from actual unpacked
proofs. Reuse generic hidden-dynamics and projection theorems; the new target
is **the exact MP7 joint law and the exact existing native count generator**.
P530's same static Stieltjes response is not already a proof of this full
native path/auxiliary conditional compatibility problem.

## 3. Frozen primary and comparison contracts

### 3.1 Primary FIN remains unchanged

q12, A=XX^T, constant A_jj, parent denominator N, gain g>0, labelled unit-clock
LOO HB. At the critical comparison use the exact g_N=g_c+lambda N^(-1/2),
g_c=12/Lambda6, fixed bounded lambda, N in12N. The old native iid-parent
theorems remain about **that** primary model.

The already admitted equilibrium joint law has

    Pi(n,theta) proportional to Mult(n)
         exp[-N||theta||^2/(2g)+theta.X^T n],
    rho_n=Law(theta|n)=Normal((g/N)X^T n,(g/N)I7).

This is a static representation, not a stochastic trajectory source.

### 3.2 Exact update-only cavity augmentation

For a selected labelled slot, r is the other N-1 counts. Keep denominator N.
Define nu_r proportional to

    exp[-N||theta||^2/(2g)+theta.X^T r] sum_j exp(X_j.theta).

Draw theta from nu_r, then redraw the label with softmax(X theta). Eliminate
the freshly drawn theta. The target is exactly the native leave-one-out law.
This theta is ephemeral per update. It is not a persisting medium or a source
for the already supplied heat-bath law. A circular sampling implementation
that first uses the native label probabilities proves no independent source.

### 3.3 Diagnostic persistent joint sampler: AUX-HB-OU-kappa

For each labelled copy, refresh its label at rate1 using softmax(X theta),
not the native count-only leave-one-out target. Independently evolve

    dtheta=kappa((g/N)X^T n-theta)dt+sqrt(2kappa g/N)dB.

Freeze kappa1 and4 only. This is one **declared diagnostic family**, not a
new law of FIN or a source-free derivation of a physical field. Its joint
invariant Pi is the existing one, while its projected kinetics need not match.

For local analysis suppress outgoing count proposals, with all twelve-label
softmax probabilities unchanged. The stationary comparison is Pi conditioned
on n in D_N; theta remains unbounded. Do not replace reflection by allowed-
label renormalization or claim the old iid trapping/mixing theorem applies
to this new generator. Start the joint diagnostic in its stated stationary
local law; actual new-model iid preparation is outside the initial scope.

## 4. Task cards

### HM-001 — Exact joint law and cavity update bridge

**Title:** What the existing Gaussian representation really generates.
**Scientific question:** does eliminating an update-level auxiliary draw
recover the precise native LOO conditional, including parent normalization?
**Why it matters:** an exact latent update identity is stronger than analogy,
but does not imply a persistent physical mediator.
**Exact inputs:** H01–H04 and section3.2.
**Dependencies:** existing positive-g Gaussian identity and constant diagonal;
no new kinetic premise for the primary model.
**Method:**

1. Write normalized labelled and count joint laws. Distinguish theta|full n
   from theta|the N-1 cavity; they are not the same distribution.
2. Complete the square for every candidate new label. Prove that integrating
   nu_r(theta)softmax_j(X theta) gives native softmax_j[(g/N)Ar].
3. Keep the removed self term and show where constant A_jj cancels. Do not
   renormalize g/N to g/(N-1).
4. For local reflection retain the full proposal then reject departures;
   exporting an allowed-label conditional sampler changes the model.
5. Check tiny cases using completed-square finite weights, not7D Monte Carlo.
   Document that the latent construction does not select the original A7,
   heat bath, attempt clock or distribution of persistent fields.

**Deliverables:** JOINT_AND_CAVITY_CONTRACT.json, EXACT_UPDATE_BRIDGE.md,
SMALL_FINITE_CHECKS.json, NONCONCLUSIONS.md and INPUTS.sha256.
**Acceptance criterion:** an exact finite-N identity with all constants and
conditions, or a preserved counterexample for the actual contract.
**Failure/refutation criterion:** complete-state posterior substituted for
cavity law, changed denominator, dropped multiplicity or boundary renormalization.
**What becomes possible:** a rigorous update-level mediation map with an
explicit distinction from dynamical field evolution.
**Must NOT conclude:** nature's mediator, sourced noise/clock, new force,
physical field propagation, spatial locality or ToE.
**Estimated computational class:** S.
**Priority:** P0.

### HM-002 — Persistent feedback, reversible compression and exactness

**Title:** Equal joint equilibrium is not equal projected process.
**Scientific question:** can AUX-HB-OU-kappa have exactly the native count
semigroup, or even an autonomous stationary projected semigroup?
**Why it matters:** this is the decisive gate before treating an auxiliary
static field as a persistent physical cause of the existing dynamics.
**Exact inputs:** HM-001, H05, section3.3 and the seed projection identity.
**Dependencies:** supplied diagnostic generator and its joint detailed balance;
the family is a comparison, not a changed primary FIN law.
**Method:**

1. Prove invariance/detailed balance of the label and OU parts separately,
   including the count-domain reflection convention when used.
2. Define P=conditional expectation onto the full counts in L2(Pi), Q=I-P.
   Reuse the old positive-memory identity rather than rediscover it.
3. For parity m, compute Q L P m. On a genuine interior type,
   Lm=sigma.softmax(X theta)-m. Use full Gaussian support and the nonconstant
   parity-axis tanh to prove its squared norm is nonzero.
4. Pay continuous-state generator/form domains and second-time differentiability:
   bounded softmax derivatives and Gaussian moments must justify the argument.
   For reflection keep a positive-mass interior type and charge departures.
5. Compare compression's second derivative with the square of its first.
   If unequal, rule out exact autonomous projected semigroup in that scope;
   do not merely claim a mismatch from the first derivative.
6. Export one bounded conditional or correlation witness. No fitted Markov
   replacement and no source-independent physical impossibility claim.

**Deliverables:** DIAGNOSTIC_MODEL_CONTRACT.json, REVERSIBILITY_AND_DOMAINS.md,
COMPRESSION_RIGIDITY_OR_COUNTEREXAMPLE.md, MEMORY_SOURCE_WITNESS.json.
**Acceptance criterion:** exactness is proved or a paid nonzero memory/
intertwining obstruction identifies precisely which equality fails.
**Failure/refutation criterion:** matching equilibrium or averaged rates is
called full dynamic equality, or a formal L² identity lacks its domain payment.
**What becomes possible:** an honest distinction between an exact native
update augmentation and a new feedback model needing its own dynamics.
**Must NOT conclude:** every possible mediator fails, all preparation-specific
approximations are impossible, or physical QM/energy/locality is disproved.
**Estimated computational class:** M, analytical mixed jump/diffusion proof.
**Priority:** P0.

### HM-003 — Derive the diagnostic collective mobility and noise

**Title:** A latent clock can change critical dynamics without changing statics.
**Scientific question:** what does the actual diagnostic generator predict
near the same critical point at fixed kappa1 and4?
**Why it matters:** this turns a missing mediator law into a concrete
kinetic consequence instead of an unmeasurable philosophical ambiguity.
**Exact inputs:** HM-001/002, H04/H06/H08 and the proposed soft chart.
**Dependencies:** correctly defined diagnostic generator; no inherited
primary-heat-bath mean closure is transferred without derivation.
**Method:**

1. Derive the full occupation/field deterministic equations and predictable
   brackets from the mixed generator. Retain label jumps and OU noise.
2. On the parity line set h=sqrt(Lambda6/12)theta6. At criticality derive
   J=[[-1,1],[kappa,-kappa]], its normalized soft left/right vectors and the
   bracket diag(2,2kappa)/N. Prove that mixed jump/continuous brackets vanish.
3. Treat all18 coordinates: prove one soft and17 stable directions for each
   frozen kappa using the exact retained spectral gaps and the four hidden
   count directions. A sampled eigenvalue plot is not a uniform statement.
4. Derive the candidate r_kappa=kappa/(1+kappa), projected cubic/gain drift
   and noise2r_kappa/N. These are generator coefficients, not a process limit.
5. At g0 define the label control directly with theta=0, rather than insert
   g0 into a singular Gaussian density. Match the per-copy label clock.
6. Freeze all coordinate and time conventions before HM-004. Do not fit
   r_kappa or rescale the final clock to make a failed comparison disappear.

**Deliverables:** JOINT_CHARACTERISTICS.md, SOFT_FAST_CHART.json,
DERIVED_MOBILITY_NOISE.json, COEFFICIENT_NEGATIVE_CONTROLS.json.
**Acceptance criterion:** every coefficient/stability assertion comes from
the actual diagnostic generator in the frozen conventions, or a precise failure.
**Failure/refutation criterion:** new macro noise/mobility is inserted,
fast count modes are lost, or equilibrium is used to guess all dynamic coefficients.
**What becomes possible:** a theorem target with a specific, independently
derived clock-dependent prediction at the same static equilibrium.
**Must NOT conclude:** kappa is sourced, the diagnostic is native FIN,
Brownian model time is physical time or a latent field is fundamental.
**Estimated computational class:** S/M, small matrix/symbolic checks and proofs.
**Priority:** P0.

### HM-004 — Pay one fixed-kappa stationary critical process limit

**Title:** A controlled compatibility limit, not a fitted mediator trajectory.
**Scientific question:** does the diagnostic's stationary local slow process
actually converge to the generator-derived candidate?
**Why it matters:** nonzero finite-N memory can coexist with a useful
asymptotic effective law, but only a separate proof determines that law.
**Exact inputs:** HM-001–003 and the actual joint stationary comparison.
**Dependencies:** paid18-coordinate characteristic contract and static
joint moments; not the primary FIN iid trapping theorem for a different generator.
**Method:**

1. Start Pi_D,N(dn,dtheta)=pi_D,N(n)rho_n(dtheta), using the existing MT
   count moments and conditional Gaussian moments for the new field variables.
2. Derive fast containment/occupation estimates or a stationary integrated
   error proof for the actual mixed generator. Account for unbounded theta
   and every positive jump/continuous quadratic variation term.
3. On fixed critical horizons, pay soft drift and bracket replacement,
   jumps, tightness, joint filtration and initialization. Reuse UF methods,
   not its theorem label for a changed generator.
4. Prove or refute

       dZ=r_kappa(a Z-Z^3/3)ds+sqrt(2r_kappa)dB,

   with the same quartic stationary marginal and both fixed kappas. If this
   candidate fails, preserve its original ID and the first failed term.
5. Pay reflected/full joint-process coupling if any unreflected statement
   is made. A native late iid-start theorem cannot be borrowed for AUX-HB-OU.
6. No simultaneous kappa_N limit, huge N sweep, spectral precision or global
   escape/capacity statement. Fixed kappas and fixed scaled horizons suffice.

**Deliverables:** STATIONARY_JOINT_LIMIT_OR_OPEN_LEMMA.md, ERROR_AND_SCOPE.json,
BOUNDARY_AND_FIELD_TAIL_PAYMENT.md and NONCONCLUSIONS.md.
**Acceptance criterion:** a real process theorem under explicit initial/domain/
clock conditions, or the exact unpaid lemma. Generator balance alone is not PASS.
**Failure/refutation criterion:** changed-generator dynamics are inferred
from shared Gibbs measure, an unbounded field/tail is dropped, or mobility is fitted.
**What becomes possible:** quantified separation of static auxiliary
equivalence from actual critical kinetic equivalence.
**Must NOT conclude:** native iid mixing/trapping, physical causal cones,
propagating fields, quantum amplitudes or a unique fundamental mediator.
**Estimated computational class:** M, proof work; no massive state enumeration.
**Priority:** P1 after the three P0 tasks.

### HM-005 — Exact-lift nonidentifiability and architectural stop

**Title:** Which mediator law, if any, was actually selected?
**Scientific question:** do the existing static and native path laws identify
a unique hidden dynamical medium, or exactly what extra premise is missing?
**Why it matters:** it prevents a comparison realization being promoted to
a physical explanation merely because it reproduces some FIN predictions.
**Exact inputs:** all HM results; the explicit exact-marginal reset lifts in
the seed; NT operational/source limits and earlier P530 context.
**Dependencies:** terminal scoped dispositions, not forced success of HM-004.
**Method:**

1. Specialize the known nonunique-lift principle to the exact joint rho_n:
   native count jumps followed by rho_n' refresh, plus arbitrary fixed
   hidden refresh rate. Prove joint detailed balance and exact full count-path
   preservation, not merely the same stationary distribution.
2. Derive one residual-field correlation showing different hidden clocks.
   This is not a new general dilation theorem; it pays the actual MP7 case.
3. Compare three objects: exact ephemeral cavity augmentation, exact-native
   analysis lifts, and persistent feedback diagnostic. Keep their hypotheses
   and path laws separate.
4. State whether HM-004 adds a clock-dependent classical effective law;
   shared statics and a common rescaled shape do not select the hidden clock.
5. Return exactly one architectural conclusion, e.g. STATIC_AUXILIARY_ONLY,
   CONTROLLED_DIAGNOSTIC_MEDIATOR_CLASS or AUXILIARY_KINETIC_SOURCE_UNDERDETERMINED.
   Identify one necessary new physical premise if relevant, then stop.

**Deliverables:** EXACT_LIFT_COMPARISON.md, HIDDEN_CLOCK_WITNESS.json,
HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md, NEXT_OPEN_ATOM.md,
REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** exact primary/diagnostic separation and an evidence-
based missing-law statement, preserving all valid partial results.
**Failure/refutation criterion:** a supplied kappa is named a prediction,
an auxiliary update is called a persistent field, or a new model replaces FIN silently.
**What becomes possible:** a precise physical-mediator source question with
known mathematical alternatives and falsifiable consequences.
**Must NOT conclude:** natural A7/kinetics selection, spatial composition,
physical units/energy, empirical evidence, selector/role transfer or QM/GR/ToE closure.
**Estimated computational class:** S/M.
**Priority:** P0, architect checkpoint after terminal task dispositions.

## 5. Kill tests and stopping policy

- If the cavity identity fails, first audit normalization/self terms; no
  fitted auxiliary density may repair it under the same contract.
- If persistent feedback has nonzero projected memory, stop exact-native
  equivalence for that model. Do not discard a correct asymptotic question.
- If the critical theorem requires a new noise/clock or field-tail premise,
  export that precise open term instead of fitting a trajectory.
- If distinct hidden laws preserve all native data, retain nonidentifiability;
  do not choose one from aesthetic similarity to a physical field.
- A new physical source claim requires independent motivation and a new
  dependency; none of the comparisons supplies units, incidence or nature's A7.

## 6. Execution and handoff

Start HM-001, HM-002, HM-003. Use existing inputs and tiny algebraic fixtures.
No field simulation is required for the first exact identities/no-go. HM-004
is a separate proof gate, not an invitation to a broad kappa/N/g sweep.

Preserve AGENTS, archived campaigns and intakes. New work stays in the new
directory; each task exports REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. Separate supplied/derived/
diagnostic/fitted objects, M/PM/PB/PP/EV, exact identities, formal candidates,
analytic limits and numerical checks. Export failed controls.

No physical adoption of a new model, hardware, spending, messages, outside
experiment, FIN-M0/B4/task335 activation, source/selector inventory replay,
ontology change, legacy-role transfer or quantum/gravity branch is authorized.
