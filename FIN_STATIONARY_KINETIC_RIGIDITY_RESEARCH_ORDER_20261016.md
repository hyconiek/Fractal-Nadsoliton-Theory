# Research order: stationary critical parity under activity freedom

Date:2026-10-16. Status: **READY_FOR_ANOTHER_EXECUTION_AGENT; NOT EXECUTED HERE.**
Output directory: `fin_stationary_kinetic_rigidity_campaign_20261016/`.
Language: English. No external experiments or new microscopic law required.

## 1. What was actually completed before this handoff

Read [THEOREMS.md](fin_activity_covariance_research_20261016/THEOREMS.md)
and [HANDOFF.md](fin_activity_covariance_research_20261016/HANDOFF.md).
The completed intrinsic preparation derives:

- exact pointwise occupation-covariance inversion into reversible rates,
  with endpoint/positivity compatibility and count multiplicities;
- an exact wheel-stress countermodel preserving Gibbs, even all coordinate
  drifts, and unit attempt capacity while changing microscopic activity/noise;
- the stationary averaged covariance identity and its critical-bracket
  consequence under an explicitly stated concentration condition NC;
- nontrivial same-parity thinning and a localized wheel-spike control.

These are not an independent physical activity source. The original
AS-002/003 and KB-002/003 bridge remains inactive. This order is separately
authorized **intrinsic mathematical research**, not a bypass of that gate.
Its class assumptions are analysis hypotheses, not newly adopted FIN laws.

## 2. One scientific target

Determine the precise scope in which the stationary quartic parity law is
independent of microscopic activity once its exact scalar drift and local
Gibbs law are fixed. The native theorem and KI's three particular kinetics
do not already prove this new class statement.

Target candidate, to be PROVED rather than labelled as an existing result:

    Z_N(s)=N^(1/4)m(p_N(sqrt(N)s))
       => dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB,
          Z0~rho_lambda.

Question beyond the sufficient theorem: is NC implied by the remaining
assumptions, or can admissible microscopic activities produce a different
limiting bracket/history? An invalid sufficient condition is not itself
a different-limit counterexample.

## 3. Inputs and frozen class

Read current AGENTS.md and these actual files, hashing them before work:

1. fin_activity_covariance_research_20261016/THEOREMS.md, RESULTS.json,
   CLAIM_REGISTER.json, verify_research.py and MANIFEST.sha256.
2. fin_rate_coupling_review_20261016/PROOF_AUDIT_AND_ADDITIONAL_CONTROLS.md.
3. fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md.
4. Same campaign: UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md.
5. fin_uniform_critical_review_20261014/INTAKE_20261014.md.
6. fin_critical_metastability_review_20261015/INTAKE_20261015.md,
   PROOF_COMPLETIONS.md and GEOMETRY_CERTIFICATE.json.
7. fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md,
   especially its actual stationary initial-law extension.
8. fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py and accepted exact spectra.

Freeze q12, supplied A7, exact g_N=g_c+lambda/sqrt(N), fixed bounded lambda,
N in12N, actual MT D_N and count weights including multinomial multiplicity.
Let L0_N be native HB reflected ONLY by rejection of outgoing D_N jumps.
Study an explicitly diagnostic count generator Ltilde_N on the same D_N:

- same stationary pi_D,N; retain detailed balance as the licensed class;
- only single-copy replacements and total count-change rate<=CN uniformly;
- EXACT reflected parity drift matching Ltilde_N m=L0_N m at every state;
- stationary initial law pi_D,N;
- first sufficient theorem: NC,
  `E_piD |G_N-E_piD G_N| ->0`, G_N=N Gammatilde_N(m,m).

Do not weaken exact drift matching to equality only at u, on the parity
line, or of the leading ODE. Full coordinate matching is sufficient but
not required. Count Markov dynamics is part of this diagnostic contract;
do not infer it for a projected hidden process from averaged rates.

## 4. Bounded tasks

### KR-001 — Bind the class and nonvacuous examples

Verify the completed stationary covariance lemma and exact support/clock
contracts. Use same-parity thinning eta(g)=1/(1+g^2) as a mandatory non-HB
example: parity drift and Gamma(m) agree exactly while other rates differ.
Keep the local reflected domain and its stationary law; no source fitting.

Export MODEL_CLASS.json, ASSUMPTION_LEDGER.md, NON_HB_EXAMPLES.md and hashes.
Acceptance: precise, nonempty class and correct distinction between NC,
equicontinuity, exact tensor data and averaged covariance. Stop on a failed
inherited input rather than assuming its content.

### KR-002 — Prove the actual stationary critical path theorem under NC

Use exact parity drift matching to transfer the actual UF-002 scalar drift
function and its full selection-rule remainder, not its conditional mean.
Use MT's high stationary moments to pay the integrated errors directly on
sqrt(N) horizons. This stationary route need not assume a new fast mixing
or gap theorem. Pay the reflected boundary difference exactly where native
full-generator expansions are used.

The completed average identity plus NC gives integrated bracket convergence
to2. Jump bound2N^(-3/4) and deterministic bracket cap pay fourth moments
and square-martingale UI. Complete maximum containment using the actual
semimartingale equation, inward cubic drift and vanishing total residual.
Pay tightness, preservation of martingale identities in the JOINT natural
filtration, uniqueness and nonexplosion. Do not replace these by sampled
paths or a generator-balance calculation.

Export STATIONARY_CLASS_PROCESS_THEOREM_OR_OPEN_LEMMA.md,
DRIFT_BRACKET_UI_LEDGER.json and NONCONCLUSIONS.md.
Acceptance: the actual path theorem under the precise class assumptions,
or one exact unpaid lemma with the theorem left unpromoted. No numerical
N0, all-preparation or growing scaled-horizon claim is required.

### KR-003 — Test which noise regularity is genuinely necessary

First retain the balanced-wheel spike as a required negative control for
automatic equicontinuity. Its pointwise covariance differs at u while its
stationary support mass vanishes, so NC and the expected limit can still
hold. Do NOT call it a nonuniversality example.

Then attempt exactly one extension: remove NC but retain stationarity,
reversibility, single-copy support, rate cap and exact parity drift. Either
prove a sufficient integrated-bracket replacement from these hypotheses,
or construct an actual finite-N kinetic family satisfying ALL of them
with a different/random limiting bracket or process. An arbitrary added
noise, a hidden environment not in count state, broken drift, changed
Gibbs weights or an unproved covariance ansatz is not a counterexample.

If no proof/counterexample is obtained in this bounded attempt, export the
precise open bracket-averaging lemma and stop. Do not repair the theory by
adding a fast environment, a graph or a phenomenological clock.

Export REGULARITY_VERDICT.md, PROOF_OR_COUNTERMODEL.md or NEXT_OPEN_LEMMA.md.
A useful sufficient theorem is retained even if this stronger question is open.

### KR-004 — Scope handoff and architectural stop

Separate proved sufficient theorem, constructed kinetic alternatives,
bounded diagnostics and remaining universality questions. State explicitly
whether the independent AS activity-provider gate changed: unless a separate
noncircular source was actually proved, it remains OPEN.

Export REPORT.md, RESULTS.json, CLAIM_REGISTER.json, HANDOFF.md,
NONCONCLUSIONS.md, REPLAY.md, INPUTS.sha256 and MANIFEST.sha256. Include
failed controls. Preserve previous sources, plans, reviews and AGENTS;
return a proposed guardrail patch only. Stop for architectural review.

## 5. Nonnegotiable scope boundaries

No iid-start diagnostic relaxation/trapping, late-iid transfer, arbitrary
initialization, labelled-block record theorem, full operator convergence,
finite-N probability rate or practical sample/clock budget is inherited.
Same count equilibrium is not enough for those statements.

No primary FIN law, ontology, physical A7/kinetic/initial source, units,
locality/incidence/conserved energy, detector/calibration/EV, selector/legacy
roles, FIN-M0/B4/task337/task335 or QM/GR/SM/L_total/ToE gate is changed.
No broad parameter scan, huge count lattice, external experiment, acquisition,
spending, message or closed FAR/source inventory replay is authorized.

An analytic proof, explicit countermodel or exact missing lemma is the
success criterion. A plot, unchanged PASS label or invented provider is not.
