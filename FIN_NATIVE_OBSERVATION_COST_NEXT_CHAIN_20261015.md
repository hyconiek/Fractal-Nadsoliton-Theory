# FIN after DKO: native finite-N observation cost through a smooth corrector

Date:2026-10-15. Baseline:380e9dc0 plus preserved LI/KO/DKO intakes.
Status: **RESEARCH_PLAN; NOC-001–NOC-005 have not been executed.**
Execution directory: `fin_native_observation_cost_campaign_20261015/`.
Language: English for the execution agent.
Decision: [NOC_DECISION.md](fin_dynamic_kernel_observability_review_20261015/NOC_DECISION.md).

## 1. One concrete parent-model theorem

DKO paid fixed-fraction/fixed-fast-horizon observation and a separate Gaussian
resource law. It did not transfer that resource law to a sequence of actual
FIN models with shrinking block fractions and growing records.

Target one native, bounded, **smooth-amplitude** test rather than full path
equivalence or a sign-only triangular theorem. No new microscopic law is
needed: a corrector is an analysis function on existing state variables.

The first target is consistency of a frozen two-model decision based on
M_N block snapshots, where

    k_N=floor(N^(15/16)), alpha_N=k_N/N,
    M_N=floor(N^(1/4)), native spacing h=1,
    W_N(t)=k_N^(-1/2) sum_(a in B_N) f3(xi_a(t_N+t)),
    A_N=(1/M_N) sum_(i=0)^(M_N-1) cos(W_N(i)).

B_N is predetermined, N is in12N and f3(j)=sqrt(2)cos(2*pi*3j/12).
Use the exact integer floors in any finite fixture. These powers are proof
design choices, not physical scaling laws. Eventually k_N<=N/2. No practical
N0 is inherited.

This readout is **not sign-only** and is not invariant under arbitrary
multiplicative amplitude calibration. It uses only the observed block
aggregate at the M grid times; hidden global coordinates below are proof
variables, never additional inputs to the classifier.

## 2. Exact inputs and unchanged premises

| Alias | Repository input |
|---|---|
| N01 | AGENTS.md, current DKO/KO/MT and unchanged physical guardrails |
| N02 | fin_dynamic_kernel_observability_review_20261015/INTAKE_20261015.md, PROOF_COMPLETIONS.md, CLAIM_REGISTER.json and verification.json |
| N03 | fin_dynamic_kernel_observability_campaign_20261015/DKO-001/ACTUAL_CHARACTERISTICS.md and DKO-001/BRACKET_AND_INITIAL_PAYMENT.md |
| N04 | fin_dynamic_kernel_observability_campaign_20261015/DKO-002/JOINT_FAST_PATH_THEOREM_OR_BLOCKER.md and DKO-004/LIMIT_ORDER_AND_TRANSFER_GATES.md |
| N05 | fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md; fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.md and KERNEL_FAMILY_SEED.json |
| N06 | fin_critical_metastability_review_20261015/INTAKE_20261015.md and PROOF_COMPLETIONS.md; fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md |
| N07 | fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md; fin_kernel_observability_campaign_20261015/KO-002/EXACT_FAST_READOUT_KERNEL.md |
| N08 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_rank7_followup/inputs/fin_handoff_audit/results.json |

Hash actual files. Use primary A7 and the frozen diagnostic A7+.01(P3-P4),
with separate model IDs, exact common g_c, original denominator N, labelled
unit-attempt LOO HB and actual iid-uniform preparation. Choose a common
deterministic post-MT origin after both existential local burn-ins. The
required observation horizon is polynomial and must fit both paid residence
windows. No physical source, hardware, graph or kinetic rule is added.

Do not import a growing-horizon diffusion limit or a uniform sampling CLT
that has not been proved. This programme seeks a native generator/decision
comparison that does not require either stronger statement.

## 3. Native state and stable normalized chart

Let b_j count block labels and o_j count outside labels, n=b+o. Their exact
Markov generator has block i->j rate

    b_i softmax_j[(g_c/N) A(n-e_i)]

and outside i->j rate with b_i replaced by o_i. The state update changes
only the corresponding count vector. Self refreshes have zero state increment.
For the auxiliary reflected comparison suppress only jumps whose **total**
count leaves the admitted D_N. Its invariant law is proportional to

    [k!/product b_j!][(N-k)!/product o_j!]
       exp[(g_c/(2N)) n^T A n] 1{n in D_N}.

This is the existing labelled local-Gibbs lift, not a newly assumed initial
ensemble or physical wall. Conditional b|n remains exactly hypergeometric.

For alpha=k/N define hidden proof coordinates

    F_N=N^(-1/2) sum_j f3(j)n_j,
    R_N=W_N-sqrt(alpha)F_N=sqrt(alpha)U_N,
    Phi_N(b,o)=(F_N,R_N).

Normalize private noise this way **before** taking alpha small. Unlike U_N,
R_N has uniformly bounded conditional variance:

    E[R_N|n]=0,
    Var(R_N|n)=(N-k)/(N-1) Var_p(f3).

MT supplies the uniform global-fast second moment. These native moment
facts, not a fixed-alpha Gaussian CLT, are the relevant shrinking-block inputs.

## 4. Analysis generator, corrector and candidate error contract

For each frozen gamma and 0<alpha<=1/2, introduce only an **analysis operator**

    G_gamma,alpha=-gamma f d_f+d_ff-r d_r+(1-alpha)d_rr.

Its stationary Gaussian variances are1/gamma and1-alpha. Let

    g_alpha(f,r)=cos(sqrt(alpha)f+r),
    mu_gamma,alpha=exp[-(alpha/gamma+1-alpha)/2].

For its semigroup P_t, the explicit formula is

    P_t g_alpha(f,r)=exp[-{alpha/gamma*(1-e^(-2gamma t))
                              +(1-alpha)*(1-e^(-2t))}/2]
                      *cos(sqrt(alpha)e^(-gamma t)f+e^(-t)r).

The proposed **discrete** Poisson corrector is

    psi_gamma,alpha=sum_(ell>=0)(P_ell g_alpha-mu_gamma,alpha),
    (I-P_1)psi_gamma,alpha=g_alpha-mu_gamma,alpha.

Prove its linear growth and bounded derivatives through order3 uniformly
for alpha in(0,1/2] and the two fixed positive gammas. No new SDE is adopted
as the finite-N dynamics.

The critical target estimate, **not yet proved by this plan**, is

    E_piD |L_D(v o Phi_N)-(G_gamma,alpha v) o Phi_N|
      <= C ||v||_(derivatives1..3) epsilon_N,k,
    epsilon_N,k=N^(-1/4)+k^(-1/2)+B_N,k,
    B_N,k <= C (N/sqrt(k)) pi_D(departure layer).

The source DKO estimates suggest this route: global drift/FF errors use
MT soft/fast moments; private normalization avoids inverse-alpha blow-up;
block sampling pays RR error; the cross bracket has conditional mean zero;
the third-jump payment is O(k^(-1/2)). Every assertion still needs a proof
uniform in k, not numerical exponent balancing.

## 5. Task cards

### NOC-001 — Native block generator and uniform normalized characteristics

**Title:** Remove the fixed-alpha restriction at the level actually needed.
**Scientific question:** do native stationary moments and generator remainders
remain controlled when k/N decreases?
**Why it matters:** a fixed-alpha OU theorem cannot simply be substituted
into a shrinking-fraction sequence.
**Exact inputs:** N01–N08, especially DKO finite formulas and MT moments.
**Dependencies:** admitted local windows for both frozen models; no new CLT.
**Method:**

1. Freeze the two operators, B_N, exact floors, native clock, readout and
   observation origin. Derive the exact (b,o) count generator and invariant
   reflected law from the labelled parent, including multiplicities.
2. Prove F_N,R_N second moments are bounded uniformly for1<=k<=N/2 in that
   stationary comparison. Use conditional sampling moments for R_N, not a
   Gaussian approximation or unbounded-moment transport through TV.
3. Reexpress all DKO characteristics in (F,R). In particular
   `sum_a |c_a|=2(N-k)` bounds private forcing; do not replace this by a crude
   max|c_a| estimate that creates an unnecessary inverse-alpha factor.
4. Pay FF, RR and FR bracket errors uniformly using conditional block means.
   For general f3 the old-label squared jump is not parity's constant2.
5. Bound the absolute third-jump compensator by C/sqrt(k) plus C/sqrt(N).
   Keep Taylor constants independent of k and both fixed models.
6. Pay the reflection/full-generator discrepancy using the real boundary
   layer, jump size O(1/sqrt(k)) and rate CN. Export the exact function class
   for the candidate section4 L1 comparison, or its first failed term.

**Deliverables:** NATIVE_BLOCK_CONTRACT.json, UNIFORM_MOMENTS_AND_CHARACTERISTICS.md,
GENERATOR_ERROR_THEOREM_OR_BLOCKER.md, BOUNDARY_PAYMENT.md, INPUTS.sha256.
**Acceptance criterion:** a genuine uniform-in-k stationary generator bound
for the required smooth tests, or one explicit unpaid term. No inferred PASS.
**Failure/refutation criterion:** inverse-alpha errors invalidate the proposed
window, a block-count component is inaccessible, or moments are assumed Gaussian.
**What becomes possible:** a direct finite-parent smooth-corrector argument.
**Must NOT conclude:** a uniform full-path Gaussian theorem, actual physical
subsystems, new noise/clock or a finite-N resource theorem already proved.
**Estimated computational class:** M, exact algebra and local moment bounds.
**Priority:** P0.

### NOC-002 — Explicit smooth OU Poisson corrector and ideal mean separation

**Title:** Build the analysis function without fitting the effective law.
**Scientific question:** can the reference corrector be controlled uniformly
in fraction and applied to the native count generator?
**Why it matters:** this gives a quantitative tool beyond qualitative weak limits.
**Exact inputs:** NOC-001 and section4's explicit semigroup.
**Dependencies:** the native function class and two fixed gamma bounds.
**Method:**

1. Verify the semigroup formula and stationary mu directly. Fix the cosine
   frequency at1; do not select it after examining discrimination.
2. Construct the discrete corrector sum. Prove convergence, linear growth
   and derivative1..3 bounds uniform in alpha<=1/2. Exponential contraction
   of the explicit linear semigroup is sufficient; no new stochastic law.
3. Prove the same derivative bounds for P_s psi,0<=s<=1, as required later
   by Duhamel. Track all dependence on the fixed gamma interval.
4. Show that the two ideal means differ by d alpha+O(alpha^2), with d!=0.
   For alpha small prove |mu_A7-mu_control|>=c_* alpha analytically.
5. Record this as an amplitude-sensitive **new** statistic. Do not transfer
   the sign-only channel's arbitrary-gain invariance to cos(W).

**Deliverables:** EXPLICIT_CORRECTOR.md, UNIFORM_DERIVATIVE_BOUNDS.md,
IDEAL_MEAN_SEPARATION.md, FROZEN_COSINE_STATISTIC.json.
**Acceptance criterion:** explicit corrector/function bounds and a strict
shrinking mean separation, not an OU assertion copied without estimates.
**Failure/refutation criterion:** necessary derivatives grow with alpha or
the declared frequency gives no leading contrast. Preserve that outcome.
**What becomes possible:** native fixed-step comparison and a bounded classifier.
**Must NOT conclude:** physical calibration, a sign-only theorem, all-observer
sufficiency or a fitted microscopic law.
**Estimated computational class:** S/M, analytical semigroup calculations.
**Priority:** P0.

### NOC-003 — Fixed-step Duhamel and native record concentration

**Title:** Pay actual temporal dependence without a growing-time diffusion limit.
**Scientific question:** can one control an empirical bounded record in the
real finite-state parent comparison from a one-step generator estimate?
**Why it matters:** Gaussian-family cost does not supply finite-N concentration.
**Exact inputs:** NOC-001–002, native stationary reflection and h=1.
**Dependencies:** uniform corrector bounds and actual invariant moments.
**Method:**

1. Let psi_N=psi o Phi_N. Prove the finite-state/OU Duhamel identity for
   `e_N=P_D,1 psi_N-(P_1 psi) o Phi_N` by differentiating the composed
   semigroups. Phi_N is not assumed to be an exact Markov quotient.
2. Apply native L1 stationarity/contraction and the NOC-001 bound to show
   `pi_D |e_N|<=C epsilon_N,k`. Do not replace pi_D with full global Gibbs.
3. Use the exact identity

       g_N-mu=psi_N-P_D,1 psi_N+e_N,
       sum_(i=0)^(M-1)(g_N(Y_i)-mu)
         =psi_N(Y_0)-psi_N(Y_M)
           +sum_(i=0)^(M-1)[psi_N(Y_(i+1))-P_D,1 psi_N(Y_i)]
           +sum_(i=0)^(M-1)e_N(Y_i).

   Here Y_i is the **full** native block/outside skeleton. The bracketed
   terms are martingale differences, not iid observations.
4. Pay uniform pi_D psi_N^2 from its linear growth and native F/R moments.
   Deduce martingale variance<=CM and endpoint second moment<=C. Bound the
   average defect in L1 by C epsilon_N,k, without a factor M.
5. Derive the actual tail bound sufficient for a mean gap of order alpha:
   `Pr_piD(|A_M-mu|>c alpha)<=C/(M alpha^2)+C epsilon_N,k/alpha`
   with endpoint terms separately shown/absorbed for M>=1.
6. Y_M is a proof terminal variable, not extra observer data: A_M uses only
   M snapshots0,...,M-1. Keep the sample/time convention explicit.

**Deliverables:** FIXED_STEP_DUHAMEL_PROOF.md, NATIVE_MARTINGALE_DECOMPOSITION.md,
DEPENDENT_RECORD_ERROR_BOUND.md, ERROR_LEDGER.json.
**Acceptance criterion:** the native stationary probability bound, or the
exact semigroup/moment term that prevents it. No IID/ESS substitution.
**Failure/refutation criterion:** a full-chain polynomial gap is silently
replaced by a constant, Phi is assumed Markov, or defect errors accumulate unpaid.
**What becomes possible:** a one-sided actual-parent observation-cost theorem.
**Must NOT conclude:** Gaussian experiment equivalence, sign-only cost,
optimal finite-N lower bounds or a practical sample number.
**Estimated computational class:** M, proof gate; no large matrices or trajectories.
**Priority:** P0.

### NOC-004 — The frozen shrinking-block finite-N classifier

**Title:** Transfer a bounded decision event to the unchanged iid-start parent.
**Scientific question:** does this specific native test consistently distinguish
the two supplied operators in one concrete joint asymptotic regime?
**Why it matters:** this is the genuinely missing finite-parent bridge.
**Exact inputs:** NOC-001–003 and the section1 fixed k_N,M_N,h contract.
**Dependencies:** paid uniform generator/concentration bounds, not just balancing powers.
**Method:**

1. Keep alpha_N=k_N/N exactly in both predicted means. Freeze the midpoint
   threshold theta_N=(mu_A7,alphaN+mu_control,alphaN)/2 and all tie treatment
   before computing a test result. No amplitude, frequency or gain retuning.
2. Define the classifier using A_N only. Control predicts the smaller cosine
   mean. The global F_N,R_N, full counts and model label are hidden proof
   variables, not supplied to the decision rule.
3. Apply the stationary native bound with the strict ideal mean gap. Verify
   for the frozen powers

       1/(M_N alpha_N^2)=O(N^(-1/8)),
       N^(-1/4)/alpha_N=O(N^(-3/16)),
       k_N^(-1/2)/alpha_N=O(N^(-13/32)),

   and include the boundary contribution; do not omit integer rounding.
4. Transfer the **whole final bounded error event** to the unreflected
   iid-start parent using the labelled initial coupling and actual departure
   flux over model time M_N. It suffices that initial TV->0; do not multiply
   it by M_N or require it to be o(alpha_N) by comparing per-sample means.
5. Keep initial exceptions and exit probabilities in the final bound. Common
   polynomial preparation plus observation fits the two exponential windows
   asymptotically. No useful N0 or physical runtime is asserted.
6. Prove the provided classifier's error->0 or report the first unpaid term.
   This is a **sufficient upper bound**, not a finite-N necessity/optimality
   theorem, a sign-only result or a transfer of all Gaussian resource laws.

**Deliverables:** FROZEN_NATIVE_TEST.json, FINITE_N_CONSISTENCY_OR_BLOCKER.md,
PARENT_TRANSFER_AND_EXIT_LEDGER.md, NONCONCLUSIONS.md.
**Acceptance criterion:** a genuine unchanged-parent test with paid errors
and declared joint sequence, or an exact remaining obstruction.
**Failure/refutation criterion:** the statistic uses inaccessible global state,
the error is not o(alpha_N), or a Gaussian asymptotic is substituted for a native bound.
**What becomes possible:** theoretical finite-parent observability with a
shrinking fraction and growing dependent record, without new physics inputs.
**Must NOT conclude:** practical feasibility, finite-N minimax cost, calibrated
physical discrimination, source selection, geometry or QM/GR.
**Estimated computational class:** S/M after the proof gates.
**Priority:** P1, dependent on NOC-001–003.

### NOC-005 — Scope handoff and architect checkpoint

**Title:** What crossed from the Gaussian limit back to actual FIN?
**Scientific question:** which finite-parent decision statement is genuinely paid?
**Why it matters:** a useful weaker native theorem is different from unproved
uniform path, sign-only or optimal information equivalence.
**Exact inputs:** all NOC results/failures and N01–N08.
**Dependencies:** terminal scoped dispositions, not forced success.
**Method:**

1. Export the DAG from exact block generator and local moments through
   uniform comparison, corrector, native concentration and final event transfer.
2. Separate reused DKO/KO/MT facts, new native estimates, numerical diagnostics,
   analysis references and actual observer-accessible inputs.
3. State the status of the stronger sign-only/likelihood/optimality questions
   separately. No success here changes their status automatically.
4. Preserve failures and all input/hash/replay objects. No large-N simulation
   campaign is needed to replace a missing inequality.
5. Return one exact next open atom and stop. Physical source, units, locality,
   conservation, apparatus/calibration and EV still need independent bridges.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
NEXT_OPEN_ATOM.md, REPLAY.md, INPUTS.sha256, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** reproducible theorem, counterexample or precise
blocker, with every expanded interpretation unpromoted.
**Failure/refutation criterion:** a smooth-amplitude result is renamed sign-
invariant, hidden state enters a classifier, or a proof scale becomes a physical law.
**What becomes possible:** an evidence-based next architecture decision.
**Must NOT conclude:** natural A7/kinetics/gain/initial-source, physical clock/
units/locality/conservation/EV, ontology/selector/legacy-role, FIN-M0/B4/task335,
QM/GR/SM/L_total or ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution, kill tests and stop rules

Start NOC-001–003. NOC-004 is the dependent native theorem, not an assumed
consequence of the proposed exponent window. If a uniform estimate fails,
export its exact term and stop that branch; any revised window/readout needs
a new contract and architectural review.

Work only in the new directory; preserve AGENTS/history and return a proposed
patch. Per task export REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. RESULTS separates execution/
science verdict, original/reference laws, supplied/reused/derived/fitted
objects, native/limit domains, initial/record/time scopes, arithmetic
assumptions, errors, failures and next action.

No new force/noise/graph/controller, broad N/g/block sweep, heavy count-state
enumeration, physical data acquisition, spending or messaging. No B4,
FIN-M0/task335, source-inventory, selector/legacy-role or quantum/gravity lane
is activated. Any external theorem used must be read and its hypotheses paid.

Success would give one genuine **finite-parent statistical decision bridge**.
It would not derive nature's kernel or clock, create space/matter/conservation,
or validate FIN experimentally. These are independent physical gaps.
