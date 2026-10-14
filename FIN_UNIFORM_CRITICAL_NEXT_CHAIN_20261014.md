# FIN after KI: symmetric-start microscopic critical fluctuations

Date: 2026-10-14. Baseline: `cf1bf9ae` plus current SF/KI intakes.
Status: **RESEARCH_PLAN; UF tasks below have not been executed.**
Execution directory: `fin_uniform_critical_campaign_20261014/`.
Language: English for the execution agent.
Architect decision: [UF_DECISION.md](fin_kinetic_invariants_review_20261014/UF_DECISION.md).

## 1. Why this is a different scientific question

The fold/KI programme derives an honest local statistical law, but starts near
a specifically chosen fold. Ask instead whether a **symmetric reference start**
at an already known FIN instability generates a controlled collective path
law from its own microscopic fluctuations, without selecting a daughter/sign.

Do not claim to discover the k6 pitchfork, its static barrier or mean drift.
Reports131/137 already contain them and reject unrestricted recursive binary
splitting, exact finite-N Markov parity and small-N autonomous-bit promotion.
This campaign must not reopen those exhausted questions. The new object is
a process-level critical limit from the actual generator, on another scaling
window, with initial moments, fast-coupling selection and order of limits paid.

No new microscopic dynamics is introduced. Primary FIN remains HB; MET/BAR
are the already frozen diagnostic controls. An iid-uniform initial distribution
is an explicit mathematical initial-law choice, **not** a derived primordial
state. A mathematical critical gain sequence is not nature's tuning mechanism.
There is no external drive/ramp, added field/noise, graph or physical clock.

## 2. Exact input registry

| Alias | Repository input |
|---|---|
| U01 | AGENTS.md, current KI/SF scope and unchanged ontology/physical gates |
| U02 | fin_kinetic_invariants_review_20261014/INTAKE_20261014.md, CLAIM_REGISTER.json, verification.json |
| U03 | fin_kinetic_invariants_campaign_20261014/KI-002/REVERSIBILITY_AND_MOBILITY.md and KINETIC_RATE_CONTRACTS.json |
| U04 | fin_kinetic_invariants_campaign_20261014/KI-003/KINETIC_CRITICAL_LIMIT.md and METROPOLIS_REGULARITY_GATE.md |
| U05 | fin_stochastic_fold_closure_review_20261014/STOPPED_LIMIT_COMPLETION.md |
| U06 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_generative_operational_campaign_20261014/RG-E04/KINETIC_CONTRACTS.json |
| U07 | fin_rank7_followup/inputs/fin_handoff_audit/results.json, exact active spectral boxes; fin_rank7_followup/certificates/R7P-026_equal_energy_event.json |
| U08 | fin_rank7_mathphysics_next/proofs/MP7-017_first_global_transition.md and MP7-031_finite_N_variational_model.md |
| U09 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/STRENGTHENED_RECURSIVE_LOCAL_CHILD_TEST_131.md and MEMORY_AWARE_METASTABLE_UNIT_VERDICT_137.md |
| U10 | The same report folder: MICROSCOPIC_PROCESS_CONTRACT_138.md and RECURSIVE_LOCAL_CHILD_TEST_131.md |
| U11 | fin_theory_synthesis_20261014/CORE_DERIVED_BRIDGES.md; fin_rank7_mathphysics_next/proofs/MP7-029_conditional_slowing_and_stochastic_laws.md |
| U12 | fin_kinetic_invariants_review_20261014/BENCHMARK_SCOPE.md; its primary Ising Theorem2 is a clock benchmark only, not a FIN theorem |

Bind actual hashes. Reuse exact operator conventions and rate clocks, never
the empirical-refresh finite-N stationary correction or auxiliary Gaussian
Langevin path. Do not replay a large root/saddle campaign. If a relevant earlier
actual-generator critical theorem is found, cite/reuse it and stop its duplicate.

## 3. Dictionary and proposed target

    u_j=1/12, sigma_j=(-1)^j, m=sigma^T p,
    p_m=u+(m/12)sigma,
    chi=Lambda6/12, g_c=1/chi,
    g_N=g_c+lambda N^(-1/2),
    Z_N(s)=N^(1/4) m_N(N^(1/2)s).

The accepted spectrum has Lambda6 as the unique top mode. The uniform
critical point is **not** the localization fold g_f≈3.51564 or first
coexistence g_eq≈3.71834. Use exact g_c from the operator, not a frozen decimal.
Its numerical proposal is about5.12342755; any gain error is amplified by
N^(1/2) in the critical coordinate and must be paid in numerical claims.

The [primary Ising benchmark](https://arxiv.org/html/0712.0790v1#S1.Thmtheorem2)
has a critical single-site-update clock of order n^(3/2). Converting attempts
to the present model clock is a benchmark inference only; FIN's global
equilibrium structure differs, so its global mixing theorem is not transferred.

For HB on p_m the already known drift is `tanh(g chi m)-m`. Its critical
cubic expansion and the actual uniform jump covariance suggest the target

    dZ=[chi lambda Z-Z^3/3]ds+sqrt(2)dB.

This target is **not yet an admitted FIN process theorem**. At uniform all
three g0-matched rules have the same leading equilibrium edge rates and
linear/noise structure, but their nonlinear/finite-N rates remain distinct.
Pay each control separately; a failed MET control need not block primary HB.

The first theorem is compact-stopped, fixed scaled horizon S and bounded
lambda. Only after this, use localization/nonexplosion to remove the compact
soft stop. No global escape or uniform-in-long-time theorem is inferred.

## 4. Key proof hints and traps

### Initialization

For iid uniform labels, occupations have fast second moment O(1/N) and fourth
moment O(1/N^2), while `N^(1/4)m_N(0)->0`. Prove these bounded-label moment
statements. A balanced-uniform deterministic type is a separately tagged
comparison, not the same finite-N initial distribution. No sign-biased seed.

### Dangerous fast coupling

At this scaling a crude `N^(3/4)*|m|*||w||` estimate need not vanish. Do not
copy the fold remainder ledger. For HB, along p_m every nonparity Fourier
feature has zero conditional average in each parity class. Therefore the
parity soft drift derivative in those fast directions vanishes all along the
parity line. Hidden means do not affect HB's leading retained drift, but do
affect finite-N rates/brackets.

For controls, use the subgroup generated by label rotation2: on the normalized
tangent its fixed space is the parity line, and all fast directions are
nontrivial representations. Use actual force/mobility blocks, not symmetry
words alone. A possible bound to prove is

    |b_m(p_m+w,g)-b_m(p_m,g)|
       <= C[||w||^2+||w||(|g-g_c||m|+|m|^3)].

For MET, locally Lipschitz mobility is enough only if the constant and
dangerous linear fast force projections are actually cancelled. Equal-logit
kinks remain present; do not call its microscopic rates globally C3.

### Longer horizon

The physical **model-time** horizon is now N^(1/2)S, not N^(1/3)S. A second
energy inequality with O(1/N) forcing gives cumulative fast energy O(N^-1/2).
The old fast tube N^-1/4 then gives only O(1) by the second-moment exit bound.
Choose, for example, tube N^-1/5 for a vanishing O(N^-1/10) bound, or pay a
stronger fourth-moment maximal estimate. Do not reuse the old N^-1/6 exit rate.
The scaled-time fast occupation target remains O(1/N), so the safe quadratic
fast error is `N^(3/4)*O(1/N)=O(N^-1/4)`.

Microscopic soft jumps are `2 N^(-3/4)` in absolute size; scaled total rate is
O(N^(3/2)); third characteristic is O(N^-3/4). The LOO drift correction O(1/N)
contributes O(N^-1/4), not the fold's old convention-correction exponent.

### Local versus global equilibrium

The candidate quartic limiting stationary density is formally

    rho_lambda(z) proportional to exp[chi lambda z^2/2-z^4/12].

It is integrable, unlike an unrestricted cubic fold potential. But at g_c
the full finite-N Gibbs law is not dominated by u: an existing ordered
coexistence state gives negative free energy there. The local limiting
stationary density must never be substituted for that full equilibrium.
Time/N limit order is a required theorem, not a qualitative caveat.

## 5. Task cards

### UF-001 — Reuse the uniform instability and freeze an unbiased start

**Scientific question:** what genuinely new process theorem remains after
deduplicating the known parity pitchfork and static stability results?
**Why it matters:** this lane must test fluctuation generation from a symmetric
start, not revive the failed recursive-bit/capacity strategy.
**Exact inputs:** U01–U11, especially source spectra and reports131/137.
**Dependencies:** KI intake; no new fold root or microscopic model.
**Method:**

1. Deduplicate existing exact parity-line, g_c, barriers, full stability window
   and finite-N nonlumpability. State the new critical path/initial-law target.
2. Freeze q/X/A, actual generator, E04 clock for controls, exact g_N and m
   normalization. No auxiliary field or a temperature interpretation.
3. Prove one zero/ten stable uniform tangent directions at g_c, with certified
   spectral separation where available; reuse intake's rational propagation.
4. Freeze iid-uniform initial law and its moment/support/tail class. Derive
   initial soft convergence and fast second/fourth moments. Separately label
   balanced-uniform deterministic types, without pooling the two ensembles.
5. Freeze compact R/S/lambda, fast projection/tube and stop semantics. The
   gain is constant along each microscopic trajectory; no ramp/controller.

**Deliverables:** MODEL_AND_INITIAL_CONTRACT.json, DEDUPLICATION.md,
UNIFORM_SPECTRAL_AND_INITIAL_MOMENT_PROOF.md, INPUTS.sha256.
**Acceptance criterion:** a nonduplicate target with one exact model and
non-sign-selecting reference initial law; complete spectral/moment contract.
**Failure/refutation criterion:** an existing theorem already supplies it,
another soft mode is unpaid, or a new drive/physical initial axiom is required.
**What becomes possible:** actual soft-generator coefficient derivation.
**Must NOT conclude:** cosmological initial-state derivation, strict selector
closure, physical binary matter or new spatial incidence.
**Estimated computational class:** S/M.
**Priority:** P0.

### UF-002 — Microscopic quartic drift/noise and fast-coupling selection

**Scientific question:** does the actual scalar parity readout have the proposed
critical characteristics when fast occupation modes are retained honestly?
**Why it matters:** hidden rates invalidate exact finite-N parity Markov closure;
an asymptotic theorem needs a real cancellation/error proof, not that shortcut.
**Exact inputs:** UF-001, U03/U04/U06/U10/U11.
**Dependencies:** fixed actual rates and critical coordinate, not fitted noise.
**Method:**

1. Derive primary HB's parity-line drift and critical cubic term, marking the
   already known deterministic formula as reused. Derive the exact finite-N
   jump drift and bracket, including LOO self exclusion.
2. Prove the dangerous linear fast response vanishes along the entire parity
   line. Give exact Fourier/parity averages, not just a numeric small Jacobian.
3. Derive MET/BAR's actual parity-line rates with the fixed prefactors. Check
   leading linear/cubic/noise coefficients at the critical uniform state;
   do not copy their fold mobility values to this different equilibrium.
4. Pay the subgroup force/mobility projection and Lipschitz remainder for any
   admitted MET control. Return an explicit unpaid term if that proof fails.
5. Export the correct scaled drift, covariance, jump and finite-N error orders.
   Small direct rates/matrices may corroborate, never substitute for bounds.

**Deliverables:** ACTUAL_PARITY_CHARACTERISTICS.md, FAST_SELECTION_RULE.md,
CONTROL_REGULARITY_GATES.md, COEFFICIENT_AND_REMAINDER_LEDGER.json.
**Acceptance criterion:** derived characteristics and no surviving O(1)
fast cross term; every control has its own paid rate/regularity contract.
**Failure/refutation criterion:** a nonvanishing fast influence or different
leading drift/noise survives; no new macro mobility/noise may repair it.
**What becomes possible:** a controlled microscopic critical path limit.
**Must NOT conclude:** exact parity lumpability, arbitrary kinetic universality
or that a Gaussian auxiliary path is the physical microscopic process.
**Estimated computational class:** M, proof and tiny numerical checks.
**Priority:** P0.

### UF-003 — Correct-horizon critical process theorem

**Scientific question:** does the whole observed trajectory converge to the
quartic confining diffusion from the stated symmetric starts?
**Why it matters:** this is the missing actual-law bridge beyond known static
pitchfork data, and its horizon/scaling differ from the completed fold theorem.
**Exact inputs:** UF-001/002, U04/U05, actual count generator.
**Dependencies:** paid primary HB characteristic and fast-coupling gate.
**Method:**

1. Prove the actual-generator fast Lyapunov bound on the N^(1/2) horizon;
   retain the quadratic jump payment and nonnegative stopped exit value.
2. Pay initial random moments/tails and a fast-tube exit bound. Derive the
   occupation error in scaled time, not an assumed stationary fast covariance.
3. Prove integrated drift/bracket replacement, vanishing jumps/third
   characteristic and UI of martingale and square martingale.
4. Use joint stopped-clock/natural-filtration identification and outer-radius
   boundary continuity. Do not infer absorption from interior tests alone.
5. Prove the compact-stopped theorem first. Then pay nonexplosion/tail
   localization for the candidate cubic restoring drift before removing R.
   Fixed-horizon weak convergence is not uniform-in-time convergence.
6. Apply to controls only after their gates pass. A MET failure does not erase
   a primary HB theorem. All initial-law and model-clock scopes remain explicit.

**Deliverables:** UNIFORM_CRITICAL_PROCESS_THEOREM.md, FAST_ENERGY_PROOF.md,
STOPPED_AND_LOCALIZATION_PROOF.md, PROOF_DEPENDENCIES.json.
**Acceptance criterion:** an actual process-level theorem for the primary
readout and specified symmetric initial law, or a precisely located blocker.
**Failure/refutation criterion:** generator balance only, copied fold exponents,
unpaid dangerous cross terms, bad fast tube or circular exit-map use.
**What becomes possible:** controlled non-Gaussian collective response and
local limiting stationary-law analysis without an added noise law.
**Must NOT conclude:** a globally autonomous parity bit, equilibrium global
phase law, capacity/B4, persistent particle, physical time or local field theory.
**Estimated computational class:** M, analytic; no broad trajectories required.
**Priority:** P0.

### UF-004 — Reconcile local dynamical and full Gibbs limits

**Scientific question:** which equilibrium does the quartic limit describe,
and does reversing long-time/N limits change the scientific interpretation?
**Why it matters:** a normalizable local limiting law can still describe a
metastable/local sector rather than typical full-parent equilibrium.
**Exact inputs:** accepted UF-003, U07/U08/U09 and exact finite-N Gibbs law.
**Dependencies:** a proved primary process limit; failed lanes stay inactive.
**Method:**

1. Derive and normalize the candidate limiting stationary density from its
   actual generator. At lambda0 give analytic moment expressions if useful;
   numerical digits are optional and not an acceptance target.
2. Use the already admitted coexistence state and exact relation
   `F_g_c(p_eq)=-(g_c-g_eq)p_eq^T A p_eq/2<0` to show u is not globally dominant.
   Bind the same spectral tuple; do not assume a globally unique localized
   branch at all larger gains or borrow unpaid B4 coverage.
3. Show what the finite-N Gibbs LDP implies for a fixed small neighbourhood
   of u, and contrast with the finite-horizon critical process staying there.
4. State the two limit orders on explicitly named occupation/readout laws.
   Prove the discrepancy for full occupation; a scalar parity discrepancy
   needs its own global phase projection, not an assumption.
5. Distinguish formal limiting invariance, actual convergence as scaled
   time grows, and any full-parent equilibrium or unbounded-time claims.
   Preserve open metastable residence/escape estimates as open.

**Deliverables:** QUARTIC_LIMIT_EQUILIBRIUM.md, LIMIT_ORDER_TABLE.md,
GLOBAL_COMPETITOR_WITNESS.json, NONCONCLUSIONS.md.
**Acceptance criterion:** theorem-level classification of local versus global
equilibrium and stated order of limits; no cubic/global or stationary scope swap.
**Failure/refutation criterion:** local quartic density is simply asserted to
be the full Gibbs marginal, or global phase/moment assumptions remain unpaid.
**What becomes possible:** an honest assessment of spontaneous transient
collective generation versus stable/full-equilibrium structure.
**Must NOT conclude:** arbitrary time-scale separation, global switching
prefactors, natural gain selection, conserved physical energy or cosmic matter.
**Estimated computational class:** S/M.
**Priority:** P0, conditional on UF-003.

### UF-005 — Handoff and architect stop

**Scientific question:** what actual structure can FIN generate from this
symmetric reference start, and what remains initialization/control-dependent?
**Why it matters:** asymptotic symmetry breaking must not be converted into
strict orientation selection, recursively generated particles or a new world.
**Exact inputs:** UF results/failures, U01/U02/U09 and unchanged source packages.
**Dependencies:** terminal scoped dispositions for UF-001–004.
**Method:**

1. Register reused static results separately from new microscopic path/limit
   theorems, nominal computations and open estimates.
2. Export supplied initial/gain/rate assumptions and derived observables,
   coefficients, symmetry selection/cancellation and stationary distinctions.
3. State whether any sign is selected uniquely: a symmetric distribution
   can produce random outcomes without a canonical strict selector.
4. Retain control failures and local/global mismatch. Decide neither a new
   kernel, cosmological preparation, physical graph nor next empirical lane.
5. Export the complete small package, then stop for the architect. No automatic
   quartic precision campaign or reopening of failed recursive-bit research.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NEXT_OPEN_ATOM.md,
NONCONCLUSIONS.md, REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** one honest prepared/reference-start process verdict
and a clear local/global classification, with hashes and independent replay scope.
**Failure/refutation criterion:** reported PASS counts replace proofs or
reference initialization/control/mean-field labels become physical ontology.
**What becomes possible:** one later evidence-based theoretical/source decision.
**Must NOT conclude:** independent source, locality/units/QM/GR/SM/ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution rules and kill tests

- Start UF-001 -> UF-002 -> UF-003. No broad spectrum/gain/N sweep. Do not
  rerun known pitchfork barriers or global saddle/bit capacity campaigns.
- Work only in the new directory. Preserve AGENTS/history and export a
  proposed accepted-claim patch. Use source scripts as proposals, not verdicts.
- Per task: REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md, NONCONCLUSIONS.md,
  REPLAY.md and INPUTS.sha256. Separate execution status and scientific status.
- Exact implicit spectral inputs define g_c; point output is not an interval.
  Existential analytic constants need no fictitious certified N0 or TV rate.
- Primary kill: a surviving scaled m*fast term, unproved stochastic containment
  or another soft mode. Preserve the failure; no noise/drive/controller repair.
- Interpretation kill: full Gibbs dominance disagrees with local limit, or
  random symmetric outcomes are relabelled a unique strict selector. Keep the
  metastable/local statistical description instead of claiming matter/reality.
- No new experiment, physical input, spending, messaging, source/FAR/selector,
  FIN-M0/B4/task335 or quantum/gravity branch is authorized.

The possible advance is a derived collective fluctuation theory from an
unbiased statistical reference ensemble. It does not close the natural
kernel/gain/initial-state, physical-time, locality, conservation or empirical
bridges required for a fundamental theory of reality.
