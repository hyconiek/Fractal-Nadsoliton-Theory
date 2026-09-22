# FIN: a mathematical-physics research agenda based on the R7O3 handoff

Repository intake update, 2026-09-22: **Target P has now been accepted** in
the declared nonnegative shared-field C4 family. See the
[completed R7O3 audit](fin_r7o3_review/README.md) and
[verification ledger](fin_r7o3_review/verification.json).
The intake obligations MP7-001–006 are now satisfied in that recorded scope;
reuse the verified evidence when its input hashes match instead of repeating
production or the full replay without a reason. The historical candidate/open
wording below describes the state when this plan was written and is superseded
only for Target P by this notice. All proposed phase-alignment, full-X7,
variational, dynamical and finite-N research remains unproved by that intake.

Date: 2026-09-21. Task namespace: **MP7-001–MP7-048**.

Primary working input, selected by the user:
[FIN_R7O3_TARGETP_HANDOFF_20260920](FIN_R7O3_TARGETP_HANDOFF_20260920/HANDOFF.md).
This completed plan supersedes the unfinished R7O2-based draft. R7O2 is retained
as predecessor evidence, not as the current research queue.

This document assigns work to an executing AI assistant. It contains proposed
theorems, proof strategies and conditional model investigations—not newly
accepted research. The present task creates this plan only: no research campaign
is launched and no existing scientific guardrails are changed.

## 1. Assessment from a mathematical-physics perspective

The R7O3 package supplies a **claimed complete covariance-curvature certificate**
for a specified finite model, extending the accepted partial R7O2 result.
The handoff itself calls this a theorem candidate pending supervisory review.
Planning from this latest result does not silently change that acceptance state.
Its scientific value is not the percentage of boxes processed. The important
questions are:

1. What exact invariant statement does the covariance bound establish?
2. Can the model's analytic structure reduce the global minimization problem,
   instead of merely making interval computation faster?
3. Which stationary branches, response functions and critical scalings follow?
4. Which fluctuation predictions can be derived after stating a finite-copy
   statistical model explicitly?
5. Which dynamical or physical conclusions still require additional premises?

The recommended new centerpiece is an **exact phase-alignment argument for the
partition function**, followed by stationary-support classification and a
small, well-chosen global variational target. This could provide a reduction
license that symmetry alone does not supply. It must be proved, including its
zero-amplitude cases, before it is used.

The campaign should produce a small number of meaningful theorems and
dimensionless predictions, with reproducible certificates. Do not replace
mathematical understanding with a new list of hundreds of isolated rectangles.

## 2. Latest working baseline and acceptance boundary

Read the following R7O3 inputs first:

- [Handoff](FIN_R7O3_TARGETP_HANDOFF_20260920/HANDOFF.md).
- [Theorem candidate](FIN_R7O3_TARGETP_HANDOFF_20260920/THEOREM_TARGET_P.md).
- [Report](FIN_R7O3_TARGETP_HANDOFF_20260920/REPORT.md) and
  [state map](FIN_R7O3_TARGETP_HANDOFF_20260920/STATE_MAP.md).
- [Centered-moment proof](FIN_R7O3_TARGETP_HANDOFF_20260920/proofs/R7O3-006_centered_moment.md)
  and [interval-jet contract](FIN_R7O3_TARGETP_HANDOFF_20260920/proofs/R7O3-007_transformed_interval_jets.md).
- [Verification ledger](FIN_R7O3_TARGETP_HANDOFF_20260920/verification.json) and
  [replay instructions](FIN_R7O3_TARGETP_HANDOFF_20260920/REPLAY.md).

The latest handoff reports:

| Quantity | R7O3 reported state |
|---|---:|
| Original residual parents certified closed | 5,432 / 5,432 |
| Active fixed-witness SAFE terminals | 12,425 |
| Remaining unresolved terminals | 0 |
| Producer and clean-directory formula passes | 12,425 each |
| Original compact partition | 18,663 = 13,231 prior SAFE + 5,432 repaired parents |
| Rejected hostile mutations | 12 / 12 |
| PD-test families | 12,411 Sylvester; 14 Gershgorin |
| Smallest reported PD-test lower margin | 13/500000000 |

The last margin is **not** a normalized covariance spectral gap. Its meaning
depends on the tested minor/bound and witness metric. The unbounded-domain join
also uses the already accepted tail theorems; tail hash identity alone does not
re-prove them. Clean-directory replay of the same algorithm establishes
reproducibility, not implementation-independent verification by itself.

Historical audited inputs:

- [R7O2 audit](fin_r7o2_review/README.md).
- [R7O2 verification](fin_r7o2_review/verification.json).
- [R7N accepted results](fin_r7n_review/ACCEPTED_RESULTS.md).
- Current [AGENTS.md](AGENTS.md).

At the historical audited R7O2 baseline:

- 7,340 SAFE leaves have passed complete recertification.
- 3,406 original parents are completely closed; 54 partial parents contain
  171 unresolved terminals; 1,972 original parents remain unprocessed.
- Accepted compact-hull coverage is about 87.6912%, expressed in coordinate
  volume, not probability or confidence.
- Global Target P at `tau0=67/250` remains open in the authoritative register.
- The two fixed-amplitude phase functions already have exactly 60 critical
  points each, with phase-Hessian index counts `(12,24,18,6)`.
- The universal full-X7 index-one statement is false; the stationary-only
  unrestricted version is also false at exact g=5.

Do not use those historical unfinished-parent counts as a new production queue.
R7O3 is the selected basis for the next scientific work. Package A is an intake
gate for its claimed completion, not a request to regenerate the certificates.
No complete fresh mathematical replay was performed while writing this plan.

One concrete editorial issue is already visible: `THEOREM_TARGET_P.md` says
the eigenvalues are ordered increasingly, while its “at most one above tau”
formulation and the min–max proof mean the **second largest** eigenvalue.
Normalize the statement to `lambda1 >= lambda2 >= lambda3 >= lambda4` and
retain the order-free formulation. This is a statement inconsistency to resolve,
not, by itself, a refutation of the certificate method.

Keep two state columns:

| State | Meaning |
|---|---|
| Accepted repository register | What the prior repository audit actually established. |
| Latest working baseline | What R7O3 claims and supplies for independent checking. |

If R7O3 passes independent review, use the resulting Target-P theorem. If it
fails or is incomplete, reopen only the specific missing proof obligations;
the [existing completion plan](FIN_R7O3_TargetP_Completion_Tasks_EN_20260920.md)
remains the fallback. Do not repeat work solely because an older README still
shows the pre-completion counters.

## 3. Model and notation contract

Use the supplied strict finite operator and its rank-seven mediator `A7=X7 X7^T`.
The Fourier feature columns are

```text
sqrt(lambda_k/6) cos(2*pi*k*j/12),
sqrt(lambda_k/6) sin(2*pi*k*j/12),       k=3,4,5,
sqrt(lambda_6/12) (-1)^j,              j=0,...,11.
```

The four-column cosine/alternating restriction is C4. It is not the full X7
space. The classical probability state remains a point of the full 12-label
simplex with an 11-dimensional tangent space.

```text
V_g(p) = D(p||u) - (g/2)(p-u)^T A7(p-u),   u=(1/12,...,1/12),
Phi_g(theta) = ||theta||^2/(2g) - log(mean_j exp((X7 theta)_j)),
g > 0,
p(theta)=softmax(X7 theta),
stationarity: theta = g X7^T p(theta).
```

For the aligned four-amplitude chart, write the unscaled field as

`h_j=J3 cos(pi*j/2)+J4 cos(2*pi*j/3)+J5 cos(5*pi*j/6)+J6(-1)^j`.

Here `Jk=sqrt(lambda_k/6)*s_k` for k=3,4,5, and
`J6=sqrt(lambda_6/12)*s_6`. Do not mix J, s, Fourier complex amplitudes,
probabilities or quantum density matrices.

Write the uniform label vector as `u0` whenever a transformed chart also uses
`u=1-s`. Reserve `s_k` for Cartesian mediator amplitudes; the compact proof
coordinate named `s` is a different scalar. Export a notation table in the
handoff rather than relying on context to resolve these collisions.

Covariance eigenvalues are ordered decreasingly. Negative-index statements
refer to the stated Hessian and coordinates. All coefficients, fields, gains,
copy counts, mobilities and noise laws remain supplied model inputs unless
a separate source theorem is actually proved.

Here D12 means the action generated by label translation `j -> j+1 mod 12`
and reflection `j -> -j mod 12` (24 group elements). A fixed alternating sign
can leave only the 12-element subgroup generated by even translations and
reflections. Do not interchange the group order with the orbit size.

## 4. Proposed analytic centerpiece: phase alignment

### 4.1 Candidate statement to prove, not assume

For fixed nonnegative amplitudes a3,a4,a5 and b, consider

`Z(phi;b)=mean_j exp(sum_{k=3,4,5} a_k cos(2*pi*k*j/12+phi_k)+b(-1)^j)`.

The proposed inequality is

`Z(phi;b) <= Z(0;b)` for b>=0.

If this holds with the appropriate equality classification, it may imply that
every global minimizer of the full dual is D12-equivalent to an aligned
nonnegative C4 representative. The reason is specific to this model: at fixed
Fourier amplitudes the quadratic dual cost is phase-independent.

This is **not** a symmetry-only assertion and does not use arithmetic averaging
of a probability vector with its reflection. The earlier obstruction to those
arguments remains valid. The proposed argument uses the signs of the exact
partition-function coefficients.

### 4.2 Suggested proof route

Derive, or justify through absolutely convergent exponential series,

`exp(a cos x) = sum_{m in Z} I_m(a) exp(i m x)`,

with nonnegative coefficients, strictly positive for every integer m when a>0.
The positivity can be proved directly from the power series for I_m; importing
an external named theorem is not necessary if the finite-series argument is
written fully.

Use

`exp(b(-1)^j)=cosh(b)+sinh(b)(-1)^j`.

Uniform averaging over j imposes a finite character-selection rule. For b>0
and all a_k>0, the supported phase frequencies should form

`L6={m in Z^3 : 3m3+4m4+5m5 = 0 mod 6}`.

For b=0, the proposed support is instead

`L12={m in Z^3 : 3m3+4m4+5m5 = 0 mod 12}`.

If the support and positivity claims are proved, then

`Z(0;b)-Z(phi;b)=sum_m c_m [1-cos(m dot phi)] >= 0`.

Pay absolute convergence before regrouping the Fourier sums. Prove the support
is exact rather than only included in the displayed lattice.

For strictly positive amplitudes and b>0, the frequency vectors
`(2,0,0)`, `(1,1,1)` and `(0,3,0)` are candidates for generators of L6.
Their determinant has absolute value six. For b=0, replace the first vector
by `(4,0,0)` and check the index-twelve statement. Establish the equality
conditions from the entire positive coefficient support.

For b<0, first use an odd label translation to flip the alternating component
and transform the other phases. Do not apply the nonnegative-coefficient
argument directly to an expansion with an unprocessed negative b.

### 4.3 Essential scope restrictions

- Positivity is needed for coefficients of **Z**, not of log Z or an arbitrary
  finite cumulant truncation. The latter can have signed coefficients.
- Equal cosine/sine spectral weights are needed for phase-independent dual
  quadratic cost. Do not transfer the result to a noncirculant perturbation
  without paying the changed metric/cost.
- Zero amplitudes require a separate equality analysis. An arbitrary phase of
  a zero-amplitude coordinate is not a distinct physical or mathematical state.
- Global-minimizer alignment would not align every stationary point or prove
  a unique minimizing orbit. Preserve all existing saddle counterexamples.
- The field/theory, kernel and gain remain supplied; phase alignment would not
  derive a physical selector or an intrinsic source law.

## 5. Work priorities and execution discipline

The eight packages below contain **48 tasks**. Use the order:

```text
A: intake-check the R7O3 completion, without regenerating its production work
B: prove/refute the partition-function alignment mechanism
C: classify supports and attack a sharply defined global variational question
D: formulate coordinate-invariant curvature and response consequences
E: derive local fold and stability laws under explicit dynamics
F: derive controlled finite-N predictions in an explicit statistical model
G: derive the stationary C4-to-X7 bridge and delimit the fixed phase census
H: independent replay, synthesis and handoff
```

B can run independently of a blocked A. C must not use phase alignment until
its proof is accepted within the campaign. D's global conclusions require A;
local D/E/F results can proceed with existing local certificates. G is a
bounded limitation/connection lane, not a new unrestricted atlas campaign.

Create a separate package, suggested name `fin_rank7_mathphysics_next/`.
Keep archives and previous certificates immutable. Use proposed patches for
integration; do not silently rewrite AGENTS.md or mark candidate evidence as
accepted. No research execution is authorized by this document's mere creation.

For each task record: exact claim, domain, input hashes, prerequisites, proof
level, commands, resource limits, complete outputs, rejected approaches and
remaining atom. Distinguish execution DONE from scientific PROVED. Counterexamples
and bounded failures are valid outcomes.

Suggested classes:

- S: at most 5 minutes and 2 GiB per run;
- M: at most 20 minutes and 4 GiB;
- L: at most 60 minutes and 6 GiB, resumable and checkpointed.

Long jobs need checkpoints and explicit worker ownership. Default to one
heavy job; parallelize only under the actual resource/delegation policy.
Do not assume this plan authorizes paid compute, external uploads or a new
unlimited algebra campaign. Stop repeated attempts after two or three failures
with the same cause unless a new mathematical input changes the method.

In particular, do not run the archive's `rm -rf checkpoints/clean_math_shards`
instruction against source evidence. Replay in a fresh output tree and preserve
all original ledgers. Do not inherit the archive's eight-worker setting blindly.
Changing a command to satisfy these safety rules does not change its mathematics.

Priority tiers: A/B first; C and the analytic part of G next; D/E next; F after
the model extension is stated. H is maintained throughout. A partial handoff
is acceptable: mark conditional descendants BLOCKED_BY_DEPENDENCY, and execute
independent tasks instead of pretending every numbered task must succeed.
Do not let margin sharpening or a broad gain atlas consume the budget needed
for the phase-alignment proof and its global variational consequences.

---

## Package A — Admit the latest R7O3 baseline without repeating production

### MP7-001 — Reconcile accepted and claimed states

**Depends:** none. **Class:** S.
Read the accepted R7O2 audit and the newer R7O3 summary. Build two separate
state columns and record the 5,432-parent identity. Note which completion
claims are genuinely new and which are inherited assumptions. Resolve the
increasing/decreasing eigenvalue-order inconsistency explicitly in a proposed
correction; do not edit the frozen source theorem in place.
Validate the complete R7O3 manifest and declared predecessor input identities.
Distinguish archive hashes, extracted-file identities and authoritative
mathematical dependencies; do not recreate deleted archives for this purpose.
**Deliver/accept:** a claim-dependency matrix. Do not rerun the 1,972 old
unprocessed parents merely because the accepted baseline predates the new
candidate; do not accept their claimed completion merely because it is newer.

### MP7-002 — Review the fixed-witness mathematical checker

**Depends:** 001. **Class:** M.
Check saved rational B,c, exact rank and Gram evidence, shared-parameter jets,
normalization, irrational-power bounds, interval rounding and PD criteria.
Confirm that replay performs no hidden eigensolver/proposal generation.
Reconstruct the seven aggregate states from the twelve-label model and confirm
that `inherited/spectral_obs_rounded9.json` encloses the actual supplied
spectral features. A hash protects file identity, not the correctness of the
feature normalization or its derivation.
**Deliver/accept:** an analytic checker audit plus negative controls. A missing
center, rounded sign diagnostic or changed threshold cannot be repaired by
trusting the producer's Boolean flag.

### MP7-003 — Verify the claimed full geometry and global join

**Depends:** 001–002. **Class:** L.
Check all original parent identities, active-version uniqueness, exact trees,
zero unresolved terminals, the R7N compact partition and accepted tails.
Verify the declared 12,425-leaf count from artifacts rather than a report.
Distinguish the 18,663 original compact cells from a fully refined cover:
the latter has 13,231 + 12,425 = 25,656 terminal cells if the two accepted
parts are expanded together. These are different levels of the same cover.
**Deliver/accept:** a complete geometry/join result or an exact list of missing
regions. Equal total volume is not sufficient, and cache/container provenance
must be distinguished from mathematical source identity.

### MP7-004 — Replay all candidate terminal formulas

**Depends:** 002–003. **Class:** L in bounded shards.
Recompute every active SAFE inequality using saved witnesses. Compare another
outward backend on a stratified subset, including the weakest reported margins.
Run deleted-leaf, altered-threshold and corrupted-witness controls.
**Deliver/accept:** a full replay, not deterministic sampling or stored completion
flags. Failed certificates return to an explicit queue; they are not automatically
model counterexamples.

### MP7-005 — Issue the Target-P adjudication

**Depends:** 003–004. **Class:** S.
If all gates pass, record a theorem candidate for `lambda2(M4)<=67/250` on
the declared shared-field domain, with its exact premises. Otherwise identify
the specific missing or invalid atoms and use the earlier completion plan only
for those atoms.
**Deliver/accept:** one unambiguous decision. Target S, full X7, equilibrium
uniqueness and physical gain provenance do not follow from Target P alone.

### MP7-006 — Extract the safe consequences and stop the bookkeeping lane

**Depends:** 005. **Class:** S.
Derive the conditional/global C4 Hessian consequence for `0<g<=250/67`,
preserving endpoint zero-mode caveats. Freeze the accepted proof interface for
later tasks and stop extending masks without a new purpose.
**Deliver/accept:** a short theorem interface and, if needed, a finite residual
queue. Do not present a coordinate-volume percentage as physical confidence.

## Package B — Prove or refute exact phase alignment

### MP7-007 — Derive the exact positive Fourier expansion

**Depends:** model contract only. **Class:** S/M.
Construct the partition-function expansion in Section 4 from the finite Z12
character average. Prove coefficient nonnegativity and absolute convergence
for arbitrary finite nonnegative amplitudes and b>=0.
**Deliver/accept:** a self-contained lemma and exact small-order checks. Do not
use positivity of log-Z coefficients, and do not confuse this all-orders
argument with the already-known quartic resonance approximation.

### MP7-008 — Prove the inequality before discussing equality

**Depends:** 007. **Class:** S/M.
Establish `Z(phi;b)<=Z(0;b)` using the real Fourier-pair sum or an equivalent
triangle argument. Handle b<0 by an exact label translation. Check consistency
with direct finite sums at adversarial phase choices and both alternating signs.
**Deliver/accept:** a proved inequality or an admissible counterexample.
If it fails, retire the reduction route immediately rather than repairing it
by a numerical maximization assumption.

### MP7-009 — Classify equality in the strictly positive interior

**Depends:** 008. **Class:** M.
Prove the support lattices L6/L12 and their annihilators. Verify the proposed
generators and indices exactly. Relate the resulting phase locks to label
translations and reflection symmetry, including the alternating-sign convention.
**Deliver/accept:** a complete equality statement for positive active amplitudes.
Numerical agreement with six or twelve maxima is not the equality proof.

### MP7-010 — Classify all zero-amplitude strata

**Depends:** 009. **Class:** M.
For each active subset of {3,4,5} and b=0 or b>0, derive the restricted support
lattice. Use the image of `m -> sum k m_k` modulo 6 or 12 to classify equality.
Discard inactive phases as redundant coordinates, not as distinct states.
**Deliver/accept:** a finite stratum table including the all-zero field and the
pure alternating field. No boundary case may be omitted from a global-minimum
reduction.

### MP7-011 — Prove the global-minimizer reduction in the exact dual

**Depends:** 008–010 and accepted primal/dual equivalence. **Class:** M.
Show that the quadratic dual cost is unchanged by phase alignment and sign-
correcting label translations. Determine whether this proves only existence
of an aligned minimizer or also that every minimizer is D12-equivalent to one.
Use the equality classification for the stronger claim.
**Deliver/accept:** a precise reduction license from full X7 minimization to
the nonnegative C4 chart, or the weaker statement actually established.
Do not infer that all stationary points are aligned or that the orbit is unique.

### MP7-012 — Challenge the reduction and document why it is new

**Depends:** 011. **Class:** M.
Check the earlier symmetry-only and reflection-averaging obstructions against
this different nonlinear dual argument. Test changed quadratic metrics,
noncirculant perturbations and broken conjugate-mode degeneracy as negative
controls. State exactly which model features the proof uses.
**Deliver/accept:** a robust scope statement. If an old obstruction genuinely
contradicts the new proof, resolve the contradiction before using the reduction;
do not declare the old result wrong solely because the new route is attractive.

## Package C — Stationary supports and a bounded global thermodynamic question

### MP7-013 — Classify self-consistent boundary supports

**Depends:** 007–011. **Class:** M.
Prove when a Fourier mean is strictly positive under an aligned nonnegative
field, using the subgroup generated by the active frequencies. Test the candidate
self-consistent supports `empty`, `{6}`, `{4}`, `{4,6}`, `{3,6}` and
`{3,4,5,6}`. Check constrained-minimum KKT conditions at zero amplitudes.
**Deliver/accept:** an exact support classification or corrected list. A missing
amplitude with strictly negative inward energy derivative cannot support a
constrained minimum.

### MP7-014 — Exclude or certify boundary branches in a declared gain window

**Depends:** 013. **Class:** M/L.
Use explicit one- or two-harmonic stationary equations on each boundary
support. A natural first window is `0<g<=250/67`, but a smaller paid interval
is acceptable. On monotone subsystems, try finite supersolution iteration
followed by a local contraction box rather than an unbounded root scan.
**Deliver/accept:** exact boundary exclusions or isolated branch certificates.
The known nonzero two-harmonic g=5 witness is a required negative control;
no exclusion may accidentally claim all gains.

### MP7-015 — Define one full-domain amplitude stationary problem

**Depends:** 011–014. **Class:** M.
Choose exact g=37/10 as a first useful target, not an arbitrary broad gain atlas.
Use `s=g E_s[C4]` to derive a compact stationary box and the accepted isotone
map to build sound box contractors. Keep uniform and boundary roots explicit.
**Deliver/accept:** a complete four-amplitude root-domain specification and
candidate list. A three-root numerical pattern is not yet exhaustive.

### MP7-016 — Prove a first global variational improvement if feasible

**Depends:** 015. **Class:** L.
Attempt a complete stationary complement exclusion at the chosen gain or an
independent global objective lower bound. Certify energies of every surviving
minimum candidate. If uniform is uniquely global at g=3.7, combine monotonicity
in g with the already-certified upper witness at 3.71835.
**Deliver/accept:** a valid sharper bracket or a bounded remaining domain.
Do not assume the expected uniform/localized/saddle catalog is complete.

### MP7-017 — Only then test whether the local crossing is globally first

**Depends:** a successful 016 and the relevant local event certificates.
**Class:** L.
Attempt parameter-uniform amplitude exclusion on a small gain strip containing
the local equal-energy event. Prove root-tube joins, boundary exclusions,
energy ordering and crossing transversality throughout the required strip.
**Deliver/accept:** identification of the first attaining orbit only if all
global obligations pass. Otherwise retain the improved bracket and the local
event as distinct results.

### MP7-018 — Publish the minimizer classification with its unresolved remainder

**Depends:** outcomes of 011–017, including bounded failures. **Class:** S.
Separate phase alignment, boundary-support exclusion, stationary exhaustion,
global energy ordering, orbit multiplicity and selection of one orbit member.
State which combination was actually proved.
**Deliver/accept:** a concise variational theorem register. Even a unique
minimizing orbit does not choose a physical label or supply QW-2191 closure.

## Package D — Invariant curvature, robustness and response

### MP7-019 — Identify what is invariant in the covariance theorem

**Depends:** model contract; 006 for a global Target-P application. **Class:** M.

1. Start from feature matrix C and the Euclidean quadratic dual cost.
2. For an invertible change `theta=T eta`, derive the transformed feature
   matrix `C'=C T`, quadratic metric `G'=T^T T`, and covariance
   `M'=T^T M T`.
3. Express the curvature statement using generalized eigenvalues of `(M',G')`,
   not ordinary eigenvalues of M' after discarding G'. Prove congruence
   invariance of Hessian inertia.
4. Distinguish linear coordinate changes from nonlinear charts: away from a
   stationary point, a nonlinear Hessian has additional gradient terms.
5. Produce one exact diagonal-rescaling example in which raw covariance
   eigenvalues change but the generalized statement and inertia do not.

**Deliver/accept:** `proofs/MP7-019_metric_invariance.md` and an exact regression
test. A numerical threshold expressed in a supplied feature normalization must
not be described as a normalization-free physical constant.

### MP7-020 — Extract a genuine spectral margin, if the certificates support it

**Depends:** 005, 019. **Class:** M/L, lower priority than B/C.

For each fixed witness set `G=B^T B` and
`K=tau0 G-E[(B^T F-c)(B^T F-c)^T]`. Obtain a verified bound `K >= m G`
with `m>0`. This gives a meaningful bound `lambda2(M4)<=tau0-m` on that cell.
One conservative route is `m >= lower(lambda_min(K))/upper(lambda_max(G))`;
derive the implication and bound both quantities rigorously. Determinant or
leading-minor positivity alone is not a numerical lower eigenvalue bound.

Aggregate all relevant cells, including the previously accepted compact SAFE
part. Audit the tails separately: strict compact bounds do not supply a uniform
strict margin on an unbounded domain. Return a compact-only margin if that is
all that is paid. Any gain extension requires a global normalized margin.

**Deliver/accept:** a margin ledger with the kind of bound and its domain.
Do not reinterpret the R7O3 value `13/500000000` as a uniform spectral gap.
No automatic Target-S claim: its actual sharper threshold must be checked.

### MP7-021 — Separate scaling equivalence from genuine model robustness

**Depends:** 019; 020 only for a quantitative global neighborhood. **Class:** M.

Prove the equivalence of `A7 -> c A7`, `g -> g/c` for c>0 at the level of V,
with the corresponding mediator coordinate change. This tests identifiability
of gain and operator normalization, not their physical provenance.

For nonuniform perturbations of positive retained spectral weights, first hold
the unscaled fields J fixed. The probability law then stays fixed, while
`C' = C S` and `M'=S M S` for a diagonal feature scaling S. Derive a valid
ordered-eigenvalue perturbation/comparison bound; do not assume ordinary
eigenvalues are invariant under this congruence. Then state how stationary
solutions move when J is not held fixed, using an implicit-function bound.

**Deliver/accept:** separate equivalence and robustness propositions. A
coefficient rectangle inside this finite model does not derive the strict
kernel, establish a legacy completion map, or transfer legacy physical roles.

### MP7-022 — Derive static response with an explicit source convention

**Depends:** local nondegenerate stationary certificates; 019. **Class:** M.

For an aligned stationary branch `s=g mu(s)`, `H=I/g-M`, differentiate to
check the candidate identity `ds/dg=H^{-1}s/g^2`. Certify its components or
selected invariant contractions on one existing local branch box.

Distinguish two source experiments:

- Adding `-f^T s` to the dual potential gives `ds/df=H^{-1}` at a stable root.
- Adding a microscopic feature field h inside the exponential family gives
  `s=g mu(s+h)`. The total feature response is
  `d mu/dh=(I-g M)^{-1} M`, where nonsingularity is required.

Derive both formulas rather than transferring a susceptibility from one source
definition to the other. Check symmetry, signs where justified, and the g=0
limit of the microscopic response. Covariance is the Fisher information metric
of this supplied exponential family, not automatically a spacetime metric.

**Deliver/accept:** a source/observable/response table, exact derivation and
verified local intervals. Do not invert through a fold singularity.

### MP7-023 — Identify the collective soft direction and its limitations

**Depends:** 006, accepted C4 cooperativity, 022. **Class:** M.

Use the accepted entrywise nonnegative C4 covariance. Pay irreducibility or
strict positivity on each region where a unique positive Perron vector is
claimed. At a positive stationary branch, compare that vector with the
certified fold null vector and response direction, using interval angles or
projector norms in the declared metric.

Target P gives at least three C4 Hessian directions with curvature
`>=1/g-tau0` for `g<1/tau0`. It does not determine the sign of the remaining
direction. A smooth leading spectral projector requires a paid spectral
separation; it is not automatic at the uniform state or a degenerate boundary.

**Deliver/accept:** a collective-mode explanation and quantified scope.
Positivity in amplitude coordinates does not choose a label, a spatial
orientation, a particle species, or a physical time direction.

### MP7-024 — Export a reusable response theorem interface

**Depends:** available results of 019–023. **Class:** S.

Collect formulas with their exact variable conventions, metric, gain domain,
source type and nondegeneracy requirements. Provide a callable routine or
small exact test fixture for each numerical use. Report both the mathematical
response and any additional premise required to interpret it physically.

**Deliver/accept:** one compact response report and machine-readable records.
Prefer a few certified observables, such as energy slope and susceptibility
along a declared direction, over an unlabelled matrix of floating numbers.
Missing global acceptance must remain visible on every global corollary.

## Package E — Local bifurcation laws and explicitly conditional dynamics

### MP7-025 — Reconstruct the simple-fold normal form from accepted data

**Depends:** accepted R7P-026–031 local event certificates. **Class:** M.

Read the original fold equations and root boxes; do not use the approximate
gain `3.51564471684` as an exact input. Choose and record a normalized null
vector v and a transverse complement. With `F=grad Phi`, verify one simple
zero mode and an invertible positive transverse Hessian in the chosen scope.
For a full-X7 statement, use the full local certificate or Package G, not
Target P alone.

Carry out Lyapunov–Schmidt elimination of the transverse variables. Derive
and certify the coefficients in

```text
epsilon = g-g_fold,
Psi(xi,epsilon) = Psi0(epsilon) + a epsilon xi + (b/6) xi^3
                 + controlled higher terms,
a = v^T partial_g F,       b = D^3 Phi[v,v,v]
```

at the exact fold. Fix the sign of v so the proposed convention is `a<0,b>0`,
if the certified signs allow it. These coefficients depend on the declared
normalization of xi; the final observable prediction must transform correctly.

**Deliver/accept:** a normal-form certificate with coefficient intervals,
invertibility bounds and an explicit neighborhood, not a cubic fit.

### MP7-026 — Pay the remainder and derive fold scaling

**Depends:** 025. **Class:** M/L.

Bound the reduced remainder and its first two xi derivatives uniformly.
Prove, rather than presume, a two-branch neighborhood for epsilon>0.
Derive and verify the leading forms

```text
xi_± ~ ±sqrt(-2 a epsilon/b),
Delta Psi_local ~ (2^(5/2)/3) |a|^(3/2) epsilon^(3/2)/sqrt(b),
soft curvature ~ sqrt(2 |a| b epsilon).
```

Check the constants against the chosen xi normalization and transverse
elimination; include explicit error bounds over a stated epsilon interval.
Use interval continuation at several points only as a cross-check of the
uniform argument. A logarithmic regression slope is not a proof of an exponent.

**Deliver/accept:** controlled branch separation, local energy difference and
soft-curvature laws. The local saddle-to-minimum energy difference is not a
proved global escape barrier without a separating-path/mountain-pass theorem.

### MP7-027 — Derive branch thermodynamics without assuming global dominance

**Depends:** existing stable and saddle branch certificates; 022, 025–026
where applicable. **Class:** M.

At a stationary point check `Phi_g(theta)=V_g(p(theta))` and derive
`d Phi_branch/dg=-||theta||^2/(2g^2)=-||mu||^2/2`. For the localized-minus-
uniform energy difference, bound its derivative near the certified local
equal-energy event near `3.71834489812038`. Pay branch existence over the
interval on which the derivative is used.

Quantify the response divergence near the fold under the source convention
of MP7-022. Separate a local equal-energy crossing, a metastable spinodal,
and a global equilibrium transition. Use “latent heat” only after introducing
and justifying an actual temperature/energy convention; otherwise report a
dimensionless slope discontinuity or order-parameter jump.

**Deliver/accept:** a branch event table with local/global status on each row.
No first-transition claim unless MP7-017 has paid the global complement.

### MP7-028 — Compare two declared gradient dynamics

**Depends:** the potential and local stability certificates. **Class:** M.

Study the explicitly added dynamical assumption `dot(s)=-L(s) grad Phi(s)`.
Compare L=I with one stated positive diagonal mobility. Prove the Lyapunov
identity and check existence of solutions in the chosen domain. For the
nonnegative C4 cone, check its boundary: at `s_i=0`, the inward flow for
positive diagonal L follows from `partial_i Phi=-mu_i<=0` where the accepted
nonnegative-mean result applies. A general SPD off-diagonal mobility need not
preserve this cone and needs its own boundary treatment.

At a critical point show that the linearization is `-L H`; derivatives of L
multiply a zero gradient there. Prove similarity to
`-L^(1/2) H L^(1/2)` for SPD L. Thus instability count is determined by Hessian
inertia, while rates and paths can change with mobility.

**Deliver/accept:** two dynamics with the same equilibria and controlled
stability comparison. Do not derive a unique evolution law from Phi alone.

### MP7-029 — Quantify conditional slowing down and distinguish stochastic laws

**Depends:** 025–028; 031–032 for a finite-copy comparison. **Class:** M.

For the declared mobility, derive the relaxation rate near a stable fold
branch, including its mobility-dependent prefactor and remainder. Explain
which part of the square-root law comes from the potential and which part
requires a regular, nonvanishing mobility.

If adding a stochastic extension, state it completely. A constant-mobility
dual Langevin process with noise covariance `2L/N` has a candidate stationary
density proportional to `exp(-N Phi)` on the full mediator space. Verify
normalization and stationarity, with boundary conditions if a restricted cone
is used. MP7-032 can relate the full-space equilibrium density to an auxiliary
field of a finite-copy model; it does not equate their dynamical trajectories.

As an alternative, give two finite-occupation reversible Markov chains with
different symmetric attempt rates but the same Gibbs law. Derive detailed
balance. Do not claim a Kramers prefactor, physical switching time or dominant
escape saddle without the required global/dynamical hypotheses.

**Deliver/accept:** a conditional kinetic proposition and an explicit
nonuniqueness-of-rates example. No physical clock is sourced.

### MP7-030 — Synthesize the local mechanism and its nonconclusions

**Depends:** available results of 025–029. **Class:** S.

Give one clear account of the fold, stable/metastable branches, local crossing,
susceptibility and mobility dependence. Label every curve segment by whether
it is certified, asymptotic with controlled remainder, or exploratory only.
Optional plots must display the validated domain, not extrapolate a theorem
to a broad gain range.

**Deliver/accept:** a local bifurcation report usable without inspecting logs.
Explicitly distinguish static landscape, chosen kinetics, noise law and
physical interpretation. Neither a Lyapunov function nor a fold supplies the
missing kernel/gain source, selector or laboratory realization.

## Package F — An explicit finite-copy statistical extension

This package adds a declared dimensionless model; it does not claim that the
repository has already derived a physical population size or a temperature.
Use the full rank-seven operator, not the four-column restriction, for the
unrestricted 12-label ensemble. N is a positive integer copy count supplied
by this extension. It is not automatically the nadsoliton's bit count.

### MP7-031 — Define the finite-N ensemble and derive its variational limit

**Depends:** model contract. **Class:** M.

For labels `x_a in {0,...,11}`, define the normalized partition sum of weights

```text
12^(-N) exp[(g/(2N)) sum_{a,b=1}^N A7[x_a,x_b]].
```

For occupations `n_j>=0`, `sum n_j=N`, and `p_j=n_j/N`, derive the exact
weight

```text
(N! / product_j n_j!) 12^(-N) exp[(Ng/2) p^T A7 p].
```

Show why `A7 u0=0` identifies the leading variational rate with V_g. Prove
the finite-state large-N variational limit with a uniform entropy bound,
including zero occupations via `0 log 0=0`; do not assume Stirling is uniform
at the boundary. Check the all-pairs diagonal convention: dropping self-pairs
changes the normalization by a computable constant if A7 has constant
diagonal, and must not be done silently.

**Deliver/accept:** an exact finite model, its partition normalization and a
proved variational limit. N and g remain supplied parameters.

### MP7-032 — Derive the exact auxiliary Gaussian-field representation

**Depends:** 031. **Class:** M.

Use completion of a seven-dimensional Gaussian square to prove or correct

```text
Z_N(g) = (N/(2*pi*g))^(7/2)
         integral_{R^7} exp[-N Phi_g(theta)] dtheta,       g>0.
```

Account for every constant and prove integrability using the quadratic cost
and bounded finite features. Construct the joint label/auxiliary-field law.
Its proposed conditionals are independent labels with distribution
`p(theta)=softmax(X7 theta)`, and
`theta | p ~ Normal(g X7^T p, (g/N) I7)`.

Derive exact consequences such as
`Cov(theta)=(g/N)I7+g^2 Cov(X7^T p)` in the full ensemble. Test N=1, g tending
to zero, and at least one very small N by direct enumeration. Handle g=0
as a limit or directly in the label model, not by substituting into a singular
Gaussian prefactor.

**Deliver/accept:** a self-contained exact bridge between two equilibrium
representations. The auxiliary field is not a newly derived physical field,
and equilibrium equivalence does not determine a microscopic dynamics.

### MP7-033 — Reconcile 11D probability and 7D mediator fluctuations

**Depends:** 031–032; a stable interior stationary certificate. **Class:** M.

Let `Sigma=diag(p)-p p^T`, `M7=X7^T Sigma X7` and `G7=I7-g M7`. In the
simplex chart with tangent columns `e_i-e_11`, i=0,...,10, derive

```text
H_chart = B0^T (diag(1/p)-g A7) B0,
det(H_chart) product_{j=0}^{11} p_j = det(G7).
```

Prove the identity by the matrix determinant lemma and check the g=0
normalization. Derive the leading stable-phase empirical covariance

```text
Cov(p) ~ [Sigma + g Sigma X7 G7^(-1) X7^T Sigma]/N,
Cov(X7^T p) ~ M7 (I7-g M7)^(-1)/N.
```

Independently derive the exact finite-N fluctuation-response identity for an
external microscopic feature source f: `d E[mu]/df=N Cov(mu)` under tilt
`exp(N f^T mu)`. Distinguish this exact finite-N relation from the single-phase
Gaussian approximation and from response to a source conjugate to theta.

**Deliver/accept:** compatible determinant, covariance and source formulas.
Using M4 here would omit the three sine fluctuations and is inadmissible for
the full ensemble. A nonpositive G7 invalidates the stable Gaussian formula.

### MP7-034 — Compute local phase weights and symmetry multiplicity

**Depends:** 032–033 and certified stable local minima near coexistence.
**Class:** M/L.

For each isolated stable minimum of the full mediator potential, derive its
Laplace prefactor `det(I7-g M7)^(-1/2)`. Cross-check using the multinomial
simplex calculation. Prove the size of the localized label orbit from its
actual stabilizer; a unique peak and reflection symmetry may imply 12 copies,
but do not import that count without checking the certified root.

Use disjoint caps around the uniform root and the translated localized roots.
If the localized orbit has 12 members, the proposed local-cap weight ratio is

```text
R_N ~ 12 exp[-N Delta V(g)] sqrt(det G7_uniform / det G7_localized),
Delta V = V_localized - V_uniform.
```

Prove a local cap gap and specify whether the caps live in mediator or
empirical-probability space. Their exact finite-N conditioning is not identical;
use MP7-032 to justify any asymptotic comparison.

**Deliver/accept:** certified local prefactor intervals and a conditional
weight-ratio formula. Without a global complement bound, this is not a claim
that these caps contain essentially all equilibrium mass.

### MP7-035 — Derive finite-size rounding with explicit error control

**Depends:** 027, 034; 017 only for global-phase claims. **Class:** M/L.

If both branches persist and the local crossing is transverse, derive the
equal-cap-mass shift suggested by

```text
Delta V(g_N) = [log(12) + 0.5 log(det G7_uniform/det G7_localized)]/N
              + controlled error.
```

Compute the sign from the certified derivative of Delta V, not intuition.
Evaluate whether a window `g=g_eq+c/N` yields a controlled logistic cap-weight
law. Replace 12 if the orbit calculation gives a different multiplicity.

Pay the remainder: positive minimum probability for a simplex approach;
bounded third/fourth derivatives; positive Hessian; cap separation; and a
bound outside the Gaussian core. A core radius `N^(-alpha)` with
`1/3<alpha<1/2` is a possible starting route, not a substitute for bounds.
Report an explicit N range if obtained; very large sufficient N is an honest
outcome. Local fold scaling and coexistence rounding use different limiting
regimes and must not be conflated.

**Deliver/accept:** a controlled finite-N prediction or an explicit missing
remainder/global-gap atom. No unqualified laboratory or universe-size claim.

### MP7-036 — Validate the extension and mark empirical checks correctly

**Depends:** available results of 031–035. **Class:** M.

Use exact small-N label/occupation sums for normalization, symmetry and
response identities. For larger N, optional sampling must record seeds,
initialization, autocorrelation/mixing diagnostics and uncertainty. Near
coexistence, failure to switch between caps is not evidence of absent mass.
Do not enumerate `12^N` states when occupation counts or the exact Gaussian
representation give a smaller validation problem.

**Deliver/accept:** analytic tests and a reproducible diagnostic notebook or
script, with sampled statements labelled NUMERICAL_EVIDENCE. The extension
must remain explicitly supplied, even if its internal predictions are exact.

## Package G — The stationary C4-to-X7 bridge and phase-census scope

### MP7-037 — Prove angular stability at interior aligned stationary points

**Depends:** 007–010; no global Target-P assumption needed. **Class:** M.

This is a high-priority analytic task, not merely another numerical Hessian
scan. Differentiate the positive Fourier expansion at aligned phases. The
candidate identity is

```text
Q := Hess_phi Phi|_{phi=0}
   = [sum_m c_m m m^T] / Z(0;b).
```

Pay termwise differentiation. If all three paired amplitudes are positive,
prove that the positive support spans R^3 and hence Q is positive definite,
including b=0 with the correct support lattice.

For Cartesian radial amplitudes `s3,s4,s5>0`, set D=diag(s3,s4,s5).
Derive the nonlinear-chart identity

```text
Q = D H_odd D - diag(s_k partial_{s_k} Phi).
```

At a full stationary point the gradient term vanishes, so
`H_odd = D^(-1) Q D^(-1) > 0`. Verify this carefully, including the sign
convention for sine columns. At a nonstationary field the gradient term cannot
be discarded.

**Deliver/accept:** an analytic positive-odd-block theorem for the stated
interior stationary family, or a precise refutation. Reconcile it with the
accepted nonstationary index-two example and the boundary two-harmonic g=5
stationary counterexample; neither has the full required interior premise.

### MP7-038 — Combine curvature and angular stability without overextension

**Depends:** 006, 037; 014 for any boundary-exhaustive claim. **Class:** M.

At a reflection-even root, prove the Cartesian Hessian splitting
`H7=H_even4 direct_sum H_odd3`. Combine positive odd curvature with Target P
to bound the full negative index at interior aligned stationary points for
`0<g<=250/67`. Keep endpoint zeros and the leading even mode explicit.

Then list what would be required to cover boundary aligned roots: the support
classification, actual branch exclusions, and any transverse blocks that
survive at zero amplitudes. MP7-011 aligns global minimizers, not arbitrary
stationary points; it does not remove nonaligned stationary families from an
all-stationary theorem.

For one existing branch near fold/coexistence, export an explicit positive
odd-block margin on a certified parameter tube if feasible. A pointwise
strict analytic theorem and a uniform numerical margin are distinct outputs.

**Deliver/accept:** a scoped 4D-to-7D stationary bridge and an obligation list,
not an unrestricted restoration of an already-refuted index-one conjecture.

### MP7-039 — Investigate the known boundary counterexample structurally

**Depends:** existing exact g=5 stationary witness and support classification.
**Class:** M/L, optional after the main analytic/global targets.

Continue only its two-harmonic invariant family over a declared bounded gain
interval. Track the relevant transverse eigenvalue blocks and, if feasible,
isolate a zero crossing with a validated augmented system. If no crossing is
isolated within the budget, report the certified tube and the remaining box.

Determine the isotropy action on the critical subspace before naming a
bifurcation. A simultaneous two-dimensional crossing cannot be treated as a
generic one-dimensional pitchfork. For a period-four pattern, investigate
the residual translation/reflection representation and whether a cubic
invariant such as `Re(z^3)` is symmetry-allowed; derive this for the actual
representation rather than assuming a universal normal form.

**Deliver/accept:** one symmetry-resolved local mechanism or a bounded open
atom. Do not turn this task into an unrestricted full-X7 gain atlas.

### MP7-040 — Test which fixed-phase critical points can be full equilibria

**Depends:** accepted fixed-fixture exactly-60 census. **Class:** M.

Keep its decimal amplitudes exact as declared and do not replace them by
coexistence amplitudes. A phase critical point solves only angular equations.
For nonzero theta, define `mu=X7^T p(theta)` and test whether

```text
(I7-theta theta^T/||theta||^2) mu = 0,
theta dot mu > 0,
g_candidate = ||theta||^2/(theta dot mu).
```

These are the conditions for a positive g with `theta=g mu`. Derive the
equivalence. Use exact root boxes/intervals to exclude full stationarity where
the projected residual is bounded away from zero. A residual interval merely
containing zero is inconclusive; if pursuing a root, use a separately specified
coupled amplitude/phase problem. Treat theta=0 separately.

**Deliver/accept:** a classification of what the 60 phase roots do and do not
imply. Do not call them 60 physical states, full equilibria or global minima.

### MP7-041 — Build a local amplitude-to-phase continuation interface

**Depends:** accepted nondegenerate phase root boxes; 040. **Class:** M/L.

For a selected nonzero-amplitude fixture and a nonsingular phase Hessian,
apply a validated implicit-function theorem to obtain `phi=phi(a)` on an
explicit amplitude box. Compute the effective amplitude Hessian by the Schur
complement

```text
H_eff = H_aa - H_a_phi H_phi_phi^(-1) H_phi_a.
```

Keep separate the Hessian of log Z and that of Phi; their angular signs are
opposite. A phase block with a negative direction for Phi excludes a full
local minimum on that branch. The complement identity requires the exact
coordinate metric and a nonsingular eliminated block.

Start with one branch that answers a live question. Expand to other branches
only if it changes a theorem, not simply to produce a larger catalog.

**Deliver/accept:** one reusable local continuation certificate, or an
explicit failure domain. Do not infer a global amplitude homotopy or reuse
the fixed-fixture census outside its paid amplitude box.

### MP7-042 — Reconcile global minima, aligned saddles and the phase catalog

**Depends:** available results of B/C and 037–041. **Class:** S.

Produce a quantifier table for: all fields; all stationary points; interior
aligned stationary points; boundary aligned stationary points; all global
minimizers; and fixed-amplitude phase critical points. Record independently
whether each statement is 4D, 7D or on the probability simplex.

Explain why a global-minimizer alignment theorem can be powerful even while
a complete 7D stationary census remains open. Conversely, explain why a
positive local Hessian or a list of phase roots cannot establish globality.

**Deliver/accept:** a logically consistent scope map preserving both the new
positive results and all previously accepted counterexamples.

## Package H — Independent challenge, synthesis and portable handoff

### MP7-043 — Audit the proof dependency graph for circularity

**Depends:** maintained from 001 onward; finalize after available packages.
**Class:** S/M.

Create one graph whose nodes are propositions/certificates, not just scripts.
Mark inherited, newly proved, candidate, numerical and added-model nodes.
Verify in particular that phase alignment does not assume global minimum
classification; global classification does not assume the expected three-root
picture; and finite-N asymptotics do not assume an unproved global gap.

**Deliver/accept:** an acyclic dependency graph with exact artifact references.
Every new global statement must have a visible path to full-domain coverage
or an analytic global theorem.

### MP7-044 — Run targeted mathematical and implementation negative controls

**Depends:** relevant produced proofs/checkers. **Class:** M.

At minimum test: swapped eigenvalue ordering; omitted metric; rank-deficient B;
altered c; a missing or duplicated leaf; compensated overlap/gap; a changed
threshold; a stale input hash; direct use of negative b in the positivity
argument; omitted zero-amplitude strata; nonstationary removal of the polar
gradient term; and substituting M4 for M7 in finite-N fluctuations.

Also check failed interpretation gates: local energy equality is not global
coexistence; a phase root is not a full stationary point; an auxiliary Gibbs
law does not source a physical dynamics; and equal mass among translated
phases does not select one label.

**Deliver/accept:** expected rejection reasons and observed outcomes. Semantic
mathematical counterchecks may be documented exact examples rather than code
tests. Do not count inherited tests again as newly designed controls.

### MP7-045 — Replay the new results from a clean environment

**Depends:** available artifacts and 043–044. **Class:** bounded L shards.

Copy the executable proof inputs into a new directory, disable dependence on
producer caches and use relative paths. Verify source identities and enumerate
all records before checking them. Report versions, arithmetic backend, hardware,
per-task resources and complete failure logs. Never delete the source evidence
to force a clean run. No unapproved package/network installation is implied.

Separate exact algebra, interval-assisted proof, numerical exploration and
proof-assistant formalization. If the same checker is replayed twice, say so;
only genuinely different implementation/derivation earns that independence
claim. Reuse the complete R7O3 intake replay if its hashes match; do not spend
another full pass just to duplicate its count under a new task ID.

**Deliver/accept:** a fresh replay ledger or exact reproducibility blockers.
No concealed producer-specific absolute paths or success-by-cache shortcuts.

### MP7-046 — Write the scientific report around conclusions, not job counts

**Depends:** all available scientific outputs. **Class:** S/M.

Lead with the strongest proved statement and its domain. Then give the
physical interpretation permitted by the supplied model, limitations and
counterexamples. Report unsuccessful approaches if they eliminate a plausible
route. Place exact theorem statements and proof sketches before large tables.

Use separate columns for execution state and scientific state. Useful latter
values are `PROVED_ANALYTIC`, `PROVED_INTERVAL_ASSISTED`, `REFUTED_IN_SCOPE`,
`NUMERICAL_EVIDENCE`, `OPEN`, and `CONDITIONAL_ON_ADDED_MODEL` with a separate
proof-level field for the last case. A resource stop is not a refutation.

**Deliver/accept:** `REPORT.md`, `CLAIM_REGISTER.json`, `NONCONCLUSIONS.md`
and a finite `NEXT_ATOMS.md`. Do not replace a specific missing lemma with
another generic “derive everything from the kernel” campaign.

### MP7-047 — Prepare conservative integration proposals only

**Depends:** 043–046. **Class:** S.

Draft `AGENTS_PROPOSED_PATCH.md` and `ACCEPTED_CLAIM_CANDIDATES.md` with exact
domains, prerequisite hashes, proof levels and explicit nonconclusions. Any
R7O3 acceptance proposal must correct the eigenvalue-order wording and preserve
the Target-P/Target-S distinction. Any new minimizer or stationary bridge
proposal must state its separate analytic reduction premises.

**Deliver/accept:** review-ready text, not a silent update of AGENTS.md or an
authoritative merge. Supervisory/user review remains the acceptance boundary.
No selector, physical clock, gain source, dimensional unit, laboratory record,
legacy role transfer, SM/GR, L_total or ToE closure may be smuggled into a
downstream finite-model theorem.

### MP7-048 — Deliver the portable handoff and stop cleanly

**Depends:** all attempted tasks, including blocked/skipped records. **Class:** S.

Package the exact inputs, new proofs, executable checkers, certificates,
negative controls, manifest and a minimal clean replay command list. Include
every task ID with its outcome and a pointer to the smallest remaining atom.
Record what was not attempted and why. State whether any worker remains
running; normally terminate the campaign cleanly with no unattended jobs.

**Deliver/accept:** a self-contained handoff following the template below.
It must allow a supervisor to distinguish genuine new mathematics, repeated
replay, conditional statistical extensions and unresolved conjectures without
reconstructing the research conversation.

## 6. Concrete mathematical recipes for the executing assistant

These are starting derivations to check, not extra accepted results or
additional task IDs. Correct them explicitly if an exact model comparison
reveals a mismatch. Never force the artifacts to fit the plan.

### 6.1 Freeze the model before applying any formula

Build `MODEL_CONTRACT.json` containing the exact upstream coefficient source,
its hashes, the retained modes, column order, Fourier signs, label origin,
spectral enclosures, feature metric, gain convention and probability-to-dual
maps. The lambdas in this plan are the supplied positive mediator spectral
weights; do not substitute adjacency eigenvalues for Laplacian-sector weights.

For the R7O3 proof chart, verify the proposed correspondence

```text
r = exp(-2 J3),     A = sqrt(r) = exp(-J3),
s = exp(-3 J4/2),   t = exp(-J5/2),
y = exp(-2 J6),     z = t^sqrt(3).
```

Match the seven aggregate weights to actual label multiplicities and feature
values, rather than trusting their names. Their physical-coupling meaning is
shared dependence on these supplied fields, not a new physical source.

Consult current repository guardrails before opening a different research
subtree. In particular, the mandatory prerequisite notes for work in
`fundamental_action_reconstruction` remain binding; this campaign does not
authorize reopening its closed source/selector lanes or changing its files.

### 6.2 A self-contained route to the Fourier positivity lemma

For n>=0, derive

```text
I_n(a) = sum_{r>=0} (a/2)^(2r+n) / [r! (r+n)!],
I_(-n)(a) = I_n(a).
```

This follows from multiplying the two exponential series in
`exp[(a/2)e^(ix)] exp[(a/2)e^(-ix)]`. Sum absolute coefficients at x=0
to control convergence. For differentiated phase series, pay the additional
polynomial factors in the integer frequencies, for example through exponential
majorants or direct factorial estimates.

After the Z12 average, a product coefficient survives with weight cosh(b)
when `3m3+4m4+5m5=0 mod 12`, and with weight sinh(b) when the residue is six.
The two residue classes are disjoint. When b=0 the latter weights vanish;
when a paired amplitude is zero only its zero Fourier index survives.

For an active index set S, define the homomorphism from `Z^S` into Z6 or Z12.
Its kernel is the frequency-support lattice. Equality phases are characters
annihilating that kernel. Show that the resulting character of its finite
cyclic image extends to the ambient cyclic group. This is a useful route to
proving that the equality fields are label translates even on boundary strata.
Count distinct fields, not phase coordinates of absent amplitudes.

### 6.3 Derive the primal/dual bridge without a minimax shortcut

Consider the joint function

```text
L(p,theta) = D(p||u0) - theta^T X7^T p + ||theta||^2/(2g).
```

Minimizing over theta gives V_g(p), since `X7^T u0=0`. Minimizing over p
gives Phi_g(theta) by the finite Gibbs variational identity. This is an
infimum-over-two-variables identity, not an unjustified interchange of a
minimum and a maximum. Prove existence/coercivity and the correspondence of
minimizers, including the entropy convention at zero probability.

This argument must use the supplied positive rank-seven factorization. It
does not license applying the same dual formula to an arbitrary indefinite
kernel. Once checked, it explains why phase alignment of theta can settle
a full-simplex energy question despite not being arithmetic symmetrization
of p.

### 6.4 Boundary supports: use exact low-dimensional equations

Prove that a positive aligned field has a strictly positive Fourier mean
precisely at retained frequencies in the subgroup generated by its active
frequencies. Classify the possible generated subgroups of Z12 by their gcd.
This gives the candidate saturated support table:

| Generator gcd with 12 | Retained stationary support |
|---:|---|
| 12 | empty (uniform) |
| 6 | {6} |
| 4 | {4} |
| 3 | {3,6} |
| 2 | {4,6} |
| 1 | {3,4,5,6} |

Then verify the following equations using the actual normalization:

- Pure alternating support: `J6=(g lambda6/12) tanh(J6)`.
  The elementary `tanh(x)<x` for x>0 gives a useful low-gain exclusion.
- Pure k=4 support: with `mu4=(exp(J4)-exp(-J4/2)) /
  (exp(J4)+2 exp(-J4/2))`, solve `J4=(g lambda4/6) mu4`.
  A nonzero branch can appear through a fold before loss of uniform linear
  stability; do not use only the derivative at zero to exclude it.
- Support {4,6}: verify factorization under the order-three and order-two
  character coordinates. If it holds, the preceding two stationary equations
  decouple exactly.
- Support {3,6}: writing J=J3, K=J6, derive

```text
x = sinh(J)/(cosh(J)+exp(-2K)),
y = (cosh(J)-exp(-2K))/(cosh(J)+exp(-2K)),
J = (g lambda3/6) x,       K = (g lambda6/12) y.
```

For a declared maximum gain, an upper starting rectangle follows from x,y<=1.
If the map is order-preserving, iterating its upper corner gives progressively
smaller upper bounds on every fixed point. Terminate with an actual local
contraction/uniqueness proof around zero, or retain the surviving rectangle.
Small floating iterates alone do not prove that the only root is zero.

All coefficient intervals must come from the exact supplied spectral source.
If this route finds a boundary branch in the intended window, certify and
include it; do not discard it because the expected catalog omitted it.

### 6.5 Complete-domain root exclusion at one gain

For `gamma=(sqrt(lambda3/6),sqrt(lambda4/6),sqrt(lambda5/6),
sqrt(lambda6/12))`, every nonnegative stationary amplitude satisfies
`0<=s_i<=g gamma_i`. In J coordinates the corresponding upper bounds are
`g lambda_k/6` and `g lambda6/12`. Verify these directly from bounded features.

Cover this exact stationary hull using a deterministic queue. Each box needs
one of: an equation component excluding zero; a sound fixed-point contractor;
an interval/Krawczyk exclusion; a verified unique-root inclusion with an
accounted collar; or a recorded subdivision. Known roots do not make their
unverified surrounding collars safe. Handle shared faces consistently, and
do not discard a boundary box by a test requiring strictly positive amplitudes.

Keep the objective route available: to prove uniform globality, an exact
global lower bound for V or Phi can replace full stationary exhaustion.
Whichever route succeeds must cover the entire domain. An assumed count of
three roots, random restarts and positive C4 cooperativity are not substitutes
for that coverage.

If V is nonnegative at g=3.7, monotonicity follows from A7 being positive
semidefinite. Explicitly derive the lower-gain conclusion. If no proof closes
within the bounded campaign, return the exact residual boxes and best paid
bracket, without identifying the local crossing as globally first.

### 6.6 Minimal successful outcomes and fallback order

A strong campaign need not prove every ambitious target. Prefer these
outcomes in order:

1. Reliable R7O3 intake, with its precise statement and no unnecessary
   regeneration of completed production work.
2. An analytic phase-alignment/equality theorem, or a decisive counterexample
   exposing its exact failed premise.
3. An analytic stationary angular-stability bridge with all coordinate terms
   paid, together with the boundary-support classification.
4. One useful global variational improvement, if the analytic reduction is
   available and the bounded complement proof closes.
5. Controlled local response/fold predictions and an exact finite-N
   auxiliary-field representation, explicitly conditional on the model.

If B fails, stop all tasks using its reduction and continue independent
response/fold/statistical tasks. If A is delayed, prove analytic results that
do not need Target P and keep their Target-P corollaries conditional. If C
does not close, retain local finite-N cap ratios without global mass claims.
Do not respond to a blocked global question by launching an unlimited search.

## 7. Required package structure and task record

Suggested layout (future artifacts, not files already produced by this plan):

```text
fin_rank7_mathphysics_next/
  HANDOFF.md
  REPORT.md
  MODEL_CONTRACT.json
  INPUTS.json
  CLAIM_REGISTER.json
  TASK_LEDGER.json
  NONCONCLUSIONS.md
  NEXT_ATOMS.md
  REPLAY.md
  AGENTS_PROPOSED_PATCH.md
  ACCEPTED_CLAIM_CANDIDATES.md
  MANIFEST.sha256
  proofs/
  src/
  tests/
  certificates/
  results/
  logs/
```

Each task record must contain at least:

```json
{
  "task_id": "MP7-001",
  "execution_state": "NOT_STARTED",
  "scientific_state": "OPEN",
  "claim": "Exact proposition or explicitly bounded question",
  "quantifiers_and_domain": "All relevant variables, faces, gains and metrics",
  "model_premises": [],
  "input_hashes": {},
  "dependencies": [],
  "proof_level": "NONE",
  "commands": [],
  "resource_usage": {},
  "artifacts": [],
  "negative_controls": [],
  "remaining_atom": "Smallest missing condition, if any"
}
```

Use a structured rational encoding for exact bounds, not a float-only JSON
field. Logs and summaries must agree on active proof versions. Counts refer
to unique mathematical objects, not the number of times a command was run.
Hash the final package contents except the manifest itself and verify the
manifest after copying. Do not claim a missing historical archive was hashed
when only its extracted files were available.

## 8. Handoff template for supervisory review

The returned `HANDOFF.md` should answer these questions in this order:

1. **What is genuinely new?** State the strongest proved result and its exact
   domain in a few sentences. If no new theorem was proved, say that directly.
2. **What happened to the R7O3 input?** Give intake verdict, exact threshold,
   descending eigenvalue convention, dependency status and replay coverage.
   Distinguish repository acceptance from campaign-level verification.
3. **Does phase alignment hold?** Supply the inequality, equality strata,
   convergence argument, full-X7 reduction scope and any failed premises.
4. **What is globally known?** Separate boundary exclusion, stationary
   exhaustion, global energy ordering, first-transition identification,
   orbit multiplicity and label selection.
5. **What is locally known?** Give fold, branch, response and Hessian results
   with certified gain/coordinate boxes and asymptotic-error bounds.
6. **What was added as a model assumption?** List finite N, ensemble, source,
   mobility, noise and units. Keep an auxiliary field distinct from a sourced
   physical field, and static equivalence distinct from kinetic equivalence.
7. **What failed or remains open?** Include exact counterexamples, invalid
   methods, stopped jobs and finite residual queues; do not hide them behind
   DONE counts or global PASS labels.
8. **How is it replayed?** Provide a minimal safe command sequence, dependency
   versions, input identities, expected counts and resource requirements.
9. **What should be integrated?** Point to proposed text only, with scope
   restrictions and a list of claims that must not be promoted.
10. **What is the next smallest research atom?** Recommend at most three
    genuinely new questions justified by the resulting proof frontier.

Final acceptance rule: a complete certificate plus a valid analytic implication
can support a theorem in its stated finite model. It does not by itself supply
physical provenance, a selector, dimensional units or evidence that the model
describes nature. Preserve that distinction even when the mathematics succeeds.
