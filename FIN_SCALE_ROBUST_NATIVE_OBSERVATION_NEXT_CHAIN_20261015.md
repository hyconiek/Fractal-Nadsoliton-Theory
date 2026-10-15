# FIN after NOC: a native test invariant to a common readout scale

Date:2026-10-15. Baseline:380e9dc0 plus preserved LI/KO/DKO/NOC intakes.
Status: **RESEARCH_PLAN; SRN-001–SRN-005 have not been executed.**
Execution directory: `fin_scale_robust_native_observation_campaign_20261015/`.
Language: English for the execution agent.
Decision: [SRN_DECISION.md](fin_native_observation_cost_review_20261015/SRN_DECISION.md).

## 1. One operational bridge, not a new microscopic law

NOC paid an actual-parent consistent cosine test along a shrinking block
fraction and growing dependent record. Its amplitude normalization is part
of the observation contract. The next question is whether the same existing
FIN can support a **raw scale-invariant** bounded test.

Do not default to optimizing NOC exponents or proving a minimax lower bound.
Those are non-blocking for this purpose. No extra force, kernel source,
graph, detector dynamics or noise is added. A common-gain transformation is
a declared nuisance/unit class for analysis, not a proved physical instrument.

Freeze simultaneous null and fast readouts from the **same predetermined
block** B_N:

    f1(j)=sqrt(2)cos(2*pi*j/12),
    f3(j)=sqrt(2)cos(2*pi*3j/12),
    W_r,N(i)=k_N^(-1/2)sum_(a in B_N)f_r(xi_a(t_N+i)),
    H_N(i)=1{|W_3,N(i)|>|W_1,N(i)|},
    A_N=(1/M_N)sum_(i=0)^(M_N-1)H_N(i).

Ties are retained with H=0. Observed channels may be transformed into
`a_i W1,a_i W3` for any common finite nonzero multiplier a_i at each sample.
H is unchanged pointwise, even when a_i varies between samples. Independent
channel gains, offsets, erasure a_i=0 and general detector noise are **not**
licensed by this invariance. No actual gain/readout calibration is supplied.

This is a new readout ID, not the NOC cosine or DKO temporal-sign experiment.
The observer uses two block aggregates at M times, never full population
counts, hidden fast variables or an outcome-selected block.

## 2. Exact inputs

| Alias | Repository input |
|---|---|
| S01 | AGENTS.md, current NOC/DKO/KO/MT and unchanged physical guardrails |
| S02 | fin_native_observation_cost_review_20261015/INTAKE_20261015.md, PROOF_COMPLETIONS.md, CLAIM_REGISTER.json and verification.json |
| S03 | fin_native_observation_cost_campaign_20261015/NOC-001/UNIFORM_MOMENTS_AND_CHARACTERISTICS.md and NOC-001/GENERATOR_ERROR_THEOREM_OR_BLOCKER.md |
| S04 | fin_native_observation_cost_campaign_20261015/NOC-002/UNIFORM_DERIVATIVE_BOUNDS.md; NOC-003/FIXED_STEP_DUHAMEL_PROOF.md and NOC-003/NATIVE_MARTINGALE_DECOMPOSITION.md in the same campaign |
| S05 | fin_dynamic_kernel_observability_review_20261015/PROOF_COMPLETIONS.md; fin_dynamic_kernel_observability_campaign_20261015/DKO-001/ACTUAL_CHARACTERISTICS.md |
| S06 | fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md; fin_kernel_observability_campaign_20261015/KO-002/EXACT_FAST_READOUT_KERNEL.md |
| S07 | fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md and INTAKE_20261015.md; fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.json |
| S08 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_rank7_followup/inputs/fin_handoff_audit/results.json |

Hash actual files. Keep primary A7 and diagnostic A7+.01(P3-P4), exact common
g_c, original denominator N, labelled unit-attempt LOO HB, iid-uniform starts,
N in12N and a common deterministic post-MT origin. All previous source,
ontology, units/locality/calibration/EV and selector/role gates remain binding.

The global Gaussian/OU covariance and NOC corrector method are **reused**.
New obligations are the raw angular observable, native anti-concentration,
nonsmooth-to-smooth transfer and its own joint asymptotic window.

## 3. Ideal law and a precise nuisance boundary

At each fixed alpha, the reference snapshot has independent

    W1~Normal(0,1),
    W3~Normal(0,v_gamma(alpha)),
    v_gamma(alpha)=1+alpha(1/gamma-1).

Derive rather than merely quote its angular probability

    p_gamma(alpha)=P(|W3|>|W1|)
                 =(2/pi)arctan(sqrt(v_gamma(alpha))).

For the two fixed gamma values, p_control>p_A7 and

    p_control-p_A7=d_ang alpha+O(alpha^2),
    d_ang=(1/gamma_control-1/gamma_A7)/(2*pi)>0.

Unrestricted separate Gaussian channel gains can match the two variance
ratios exactly. That is a **limiting static Gaussian** nuisance-class overlap,
not a no-go for all native histories. The common-scale invariant theorem
must not be promoted to robustness against independent gains or RG-F's
broader calibrated error classes.

## 4. Native proof route and candidate error contract

The native two-block count state and reflected invariant law remain NOC's.
Use the four-coordinate analysis chart

    F_r=N^(-1/2)sum_(all a)f_r(xi_a),
    R_r=W_r-sqrt(alpha)F_r, r=1,3.

This chart is not assumed to be an exact Markov quotient. The reference OU
analysis operator is the sum of the two NOC operators with gamma1=1 and
gamma3=gamma. All generators, moments and mixed brackets must be paid in
four dimensions, not inferred from separate scalar results.

The raw event boundary is the union of lines W3+W1=0 and W3-W1=0. On the
positive MT cylinder, prove the **native conditional** bound

    P(|W3 +/- W1-z|<=eta | n)<=C(eta+k^(-1/2))

uniformly in shifts z and the allowed total counts. No quantitative Gaussian
CLT is assumed. A useful route selects two label types with different values
of h_plus=f3+f1 or h_minus=f3-f1, conditions on all remaining block counts,
and bounds the maximum mass of the remaining two-type hypergeometric.

Choose a fixed C-infinity compactly supported probability mollifier only
for the proof: H_eta=H convolved with rho_eta. The raw observer never applies
this unit-dependent smoothing or receives artificial noise. Its error is
supported within C eta of the two boundary lines, so stationarity/native
anti-concentration pays E|H-H_eta|<=C(eta+k^(-1/2)).

The intended, **not yet established**, bound is

    P_piD(wrong raw angular classifier)
       <=C[eta^(-2)/(M alpha^2)
             +eta^(-3)epsilon_N,k/alpha
             +(eta+k^(-1/2))/alpha],
    epsilon_N,k=N^(-1/4)+k^(-1/2)+(N/sqrt(k))pi_D(B_D).

Its terms respectively pay corrector martingale variance, native generator
defect and nonsmooth-event replacement. The proof must justify every eta
power; the formula is a target, not a result of this plan.

Freeze a **new**, more conservative joint sequence:

    k_N=floor(N^(63/64)), alpha_N=k_N/N,
    M_N=floor(N^(1/4)), eta_N=N^(-1/32), spacing1.

This observes a larger fraction than NOC's sequence at the same N. It is
not an improved/optimal resource bound and gives no practical N0. Exact
floors, eventual k<=N/2 and boundary/preparation errors remain mandatory.

## 5. Task cards

### SRN-001 — Freeze the raw invariant test and its nuisance scope

**Title:** Two simultaneous channels, one common unknown scale.
**Scientific question:** what invariance and ideal separation does the raw
angular observable actually provide?
**Why it matters:** removing an amplitude premise is an operational bridge;
declaring an unproved detector model would merely add that premise elsewhere.
**Exact inputs:** S01–S08 and sections1–3.
**Dependencies:** accepted NOC/KO mathematics; no new kinetic law.
**Method:**

1. Freeze the two matrices, feature-label maps, block, times, raw event,
   tie convention and new k_N,M_N sequence before producing comparison output.
2. Prove exact invariance under a common nonzero a_i, including finite
   zero/tie cases and time-varying common multipliers. Exclude erasure and
   do not silently add invariance to additive offsets.
3. Derive the Gaussian angular probability by isotropic geometry or a fully
   read primary result. Prove its positive model contrast and linear-in-alpha
   lower separation for a declared small-alpha interval.
4. Construct the exact separate-gain overlap in the limiting static Gaussian
   class. Keep this distinct from native/full-history impossibility and the
   narrower common-scale claim.
5. Export which readout premise a future physical instrument would have to
   license. The mathematics does not prove that a real instrument has it.

**Deliverables:** FROZEN_RAW_READOUT.json, GAUGE_AND_NUISANCE_THEOREM.md,
IDEAL_ANGULAR_SEPARATION.md, OBSERVER_INPUT_CONTRACT.md, INPUTS.sha256.
**Acceptance criterion:** exact invariance and a scoped nonzero ideal contrast,
or a retained failure. No fitted gain/threshold/noise model.
**Failure/refutation criterion:** robustness needs independent gains secretly
restricted, discarded ties, a new clock, or a post-hoc observation choice.
**What becomes possible:** a well-defined raw native target, not just a
smoothed proxy or calibrated-amplitude cosine result.
**Must NOT conclude:** physical calibration, all-instrument robustness,
natural kernel selection or a new microscopic assumption.
**Estimated computational class:** S/M, exact algebra and probability geometry.
**Priority:** P0.

### SRN-002 — Native anti-concentration at the angular boundaries

**Title:** Pay nonsmooth observations without importing a CLT error rate.
**Scientific question:** how much actual conditional block probability lies
near W3+W1=0 or W3-W1=0?
**Why it matters:** bounded weak convergence alone is inadequate when the
required error is o(alpha_N) in a shrinking-fraction experiment.
**Exact inputs:** SRN-001, exact hypergeometric kernel, MT positive coordinate bound.
**Dependencies:** n_j>=cN in the actual comparison domain, k<=N/2.
**Method:**

1. Work conditional on the complete total count; use the exact without-
   replacement law. Keep the shift z arbitrary and all lattice atoms/ties.
2. For h_plus choose labels0,6, whose values differ by4sqrt(2). For h_minus
   choose labels2,8, whose values differ by3sqrt(2). Verify these exact feature
   values and that each selected type has at least cN copies.
3. Condition on how many selected block slots come from that two-type union
   and on remaining block counts. Prove a good-event size between c'k and
   a fixed fraction below the union population. Hypergeometric variance/
   Chebyshev bounds can pay the bad-event probability by C/k.
4. On the good event the two-type hypergeometric maximum atom must be <=C/sqrt(k).
   Prove this with factorial/Stirling bounds or a fully checked finite-
   population result. Do not infer it from a few enumerated cases.
5. An interval of width eta in the normalized sum contains at most
   C(1+eta sqrt(k)) such atoms, regardless of the other-label shift. Deduce
   the uniform C(eta+k^(-1/2)) bound, then the angular-line union bound.
6. Verify tiny support cases independently, as algebra diagnostics only.

**Deliverables:** NATIVE_ANTICONCENTRATION.md, PAIR_CONDITIONING_PROOF.md,
LATTICE_AND_TIE_SCOPE.json, SMALL_EXACT_CHECKS.json.
**Acceptance criterion:** a native bound uniform in total counts, shift and
k in the declared class, or the exact counterexample/unpaid estimate.
**Failure/refutation criterion:** conditioning collapses variance unnoticed,
a lattice atom is dropped, or a Gaussian density bound replaces native control.
**What becomes possible:** honest raw-to-smooth comparison at a shrinking bandwidth.
**Must NOT conclude:** a quantitative full-path CLT, native independence,
physical resolution or an optimally small block requirement.
**Estimated computational class:** M, conditional combinatorics and inequalities.
**Priority:** P0.

### SRN-003 — Four-coordinate corrector with explicit smoothing cost

**Title:** Extend the paid native method to this different nonsmooth target.
**Scientific question:** can the raw-event smoothing error and corrector
derivatives be controlled together in the actual parent?
**Why it matters:** a scale-dependent proof approximation must not become
the observer or conceal an inverse-bandwidth error.
**Exact inputs:** SRN-001–002, NOC generator/Duhamel/martingale results.
**Dependencies:** native anti-concentration and all four normalized moments.
**Method:**

1. Define one compactly supported smooth analysis mollifier and H_eta.
   Prove |H-H_eta| is supported in the C eta boundary strips. Do not add
   random measurement noise to the microscopic dynamics or final data.
2. Extend the native smooth-generator bound to (F1,R1,F3,R3), including every
   mixed bracket. Separate scalar bounds do not prove a joint operator result.
3. Show derivative orders1–3 of H_eta composed with the observed linear map
   cost C eta^(-r). Construct its discrete OU corrector, with derivative
   bound C eta^(-3) and linear growth C eta^(-1)(1+||x||), uniformly in alpha.
4. Pay native corrector second moment <=C eta^(-2) using local/HG moments.
   Repeat the **full-state** fixed-step Duhamel and martingale decomposition;
   Phi is not assumed Markov and no full-chain constant gap is substituted.
5. Obtain the smooth-test concentration bound. Pay both native empirical
   raw/smooth error and the limiting Gaussian mean error separately using
   boundary probabilities, not uniform approximation of the discontinuity.
6. Export the exact eta/alpha/N dependence. If the candidate section4 bound
   fails, give the first failed term rather than a numerical PASS.

**Deliverables:** JOINT_NATIVE_GENERATOR_BOUND.md, MOLLIFIED_CORRECTOR.md,
RAW_EVENT_TRANSFER_BOUND.md, SMOOTHING_ERROR_LEDGER.json.
**Acceptance criterion:** a paid raw-event bound uniform in the stated alpha
and eta ranges, or one explicit estimate blocking it.
**Failure/refutation criterion:** derivative/martingale costs are omitted,
H_eta replaces the raw invariant rule, or finite-N moments are transported through TV.
**What becomes possible:** a genuinely scale-invariant native classifier under
a separately declared joint window.
**Must NOT conclude:** automatic sign-only temporal transfer, optimal native
cost, physical detector noise or an extra informational substrate.
**Estimated computational class:** M, proof gate, no state-lattice enumeration.
**Priority:** P0.

### SRN-004 — Frozen raw native classification and parent transfer

**Title:** Consistency of the actual common-scale invariant decision.
**Scientific question:** does the raw test distinguish the two unchanged
parents under the new shrinking-block/growing-record sequence?
**Why it matters:** this is the native operational result, not a Gaussian-only ratio.
**Exact inputs:** SRN-001–003 and the section4 new exact-floor sequence.
**Dependencies:** paid raw concentration/separation, not exponent balancing alone.
**Method:**

1. Use the exact alpha_N and freeze theta_N=(p_A7(alpha_N)+p_control(alpha_N))/2.
   Declare control when A_N>theta_N and A7 otherwise; keep ties under that rule.
2. Verify the candidate dominant terms tend to zero:

       eta_N^(-2)/(M_N alpha_N^2)=O(N^(-5/32)),
       eta_N^(-3)N^(-1/4)/alpha_N=O(N^(-9/64)),
       eta_N^(-3)k_N^(-1/2)/alpha_N=O(N^(-49/128)),
       eta_N/alpha_N=O(N^(-1/64)),
       k_N^(-1/2)/alpha_N=O(N^(-61/128)).

   Include boundary and exact rounding terms; this is not a practical rate
   or an improvement over the earlier cosine window.
3. Apply the native stationary estimate to the **raw** H observations.
   Unknown common nonzero multipliers leave every recorded H_i unchanged.
   Do not feed hidden F/R, full counts or mollified values to the classifier.
4. Transfer only the final bounded error event by actual initial/exit coupling
   over the polynomial record. Initial TV is added once, not divided by alpha
   or eta and not multiplied by record length.
5. Distinguish M recorded snapshots from any hidden proof endpoint. Preserve
   all outside/exit trajectories as error, never survivor postselection.
6. Prove error->0 or return its first unpaid term. State the nuisance class
   exactly; independent channel gains/offsets/general noise remain unpaid.

**Deliverables:** FROZEN_SCALE_ROBUST_NATIVE_TEST.json,
NATIVE_INVARIANT_CONSISTENCY_OR_BLOCKER.md, PARENT_AND_NUISANCE_TRANSFER.md,
NONCONCLUSIONS.md.
**Acceptance criterion:** a native raw-test consistency theorem for the
declared common-scale class, or a precise failure under the frozen contract.
**Failure/refutation criterion:** only the smoothed/calibrated rule succeeds,
an error exceeds the shrinking mean gap, or observation data are silently enlarged.
**What becomes possible:** removal of one genuine readout premise without
changing the microscopic FIN law.
**Must NOT conclude:** physical instrument licensing, arbitrary gain robustness,
finite-N optimality, full information equivalence or natural A7 selection.
**Estimated computational class:** S/M after the three proof gates.
**Priority:** P1.

### SRN-005 — Scope handoff and next physical-bridge checkpoint

**Title:** What readout restriction was actually removed?
**Scientific question:** which operational claim now follows, and which
physical assumptions are still missing?
**Why it matters:** invariant mathematics is not automatic empirical calibration.
**Exact inputs:** all SRN results/failures and S01–S08.
**Dependencies:** terminal scoped dispositions, not forced success.
**Method:**

1. Export the DAG from the same labelled law through native joint moments,
   anti-concentration, smoothing/corrector, raw decision and parent coupling.
2. Separate reused NOC/KO/DKO facts, new estimates, proof-only variables and
   actual observer data. Preserve negative and open results.
3. Keep cosine consistency, angular common-scale consistency, temporal
   sign-only transfer and minimax necessity as distinct scientific statuses.
4. State missing physical source, instrument/common-gain licensing, units,
   locality/conservation and empirical evidence. Tagging or angular variables
   do not supply spatial incidence or an observer ontology.
5. Stop for architecture review. No automatic frequency/threshold/exponent
   optimization, lower-bound campaign, geometry/quantum branch or experiment.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
NEXT_OPEN_ATOM.md, REPLAY.md, INPUTS.sha256, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** reproducible scoped result or exact blocker,
with every broader physical or mathematical promotion explicitly withheld.
**Failure/refutation criterion:** a proof mollifier is called a detector,
common-scale invariance becomes arbitrary-readout robustness, or a control is adopted as FIN.
**What becomes possible:** one evidence-based physical/operational bridge decision.
**Must NOT conclude:** source/units/locality/conservation/EV, ontology/selector/
legacy-role, FIN-M0/B4/task335 or QM/GR/SM/L_total/ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution and stop rules

Start SRN-001–003. SRN-004 depends on all three; its powers are proposed
proof scales, not a theorem. A failed uniform estimate stops that branch.
Any changed readout/window requires a new identifier and architectural review.

Use only the new directory. Preserve AGENTS/history and return a proposed
patch. Per task export REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. RESULTS separates execution/
science status, original/reference objects, supplied/reused/derived/fitted
inputs, nuisance/readout/initial/time scopes, arithmetic assumptions,
errors, failures and next action.

No broad sweep, heavy lattice enumeration, new microscopic noise/drive/
graph/controller, hardware, experiment, spending or messages. Read any
external theorem actually used and pay its exact hypotheses. Tiny fixtures
verify conventions, not the uniform probability estimates.

Success would strengthen native **statistical** observability under a precise
readout transformation class. It would not derive nature's A7, kinetic law,
physical time/space, conserved energy, detector or quantum/gravity theory.
