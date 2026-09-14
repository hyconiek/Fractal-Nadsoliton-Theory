# FIN post-handoff research master plan for an executing AI agent

Version: 1.0 — 2026-09-13.

Prepared from the accepted local handoff audit, not from the original chat's
optimistic summary. This document is a research execution specification. Its
proposed calculations, lemmas, and milestones are **not already completed
research**. The present writing task does not launch the campaign.

## 1. Mission and definition of completion

Execute a broad, falsification-first mathematical research campaign around
the supplied strict rank-seven model. Convert recoverable numerical claims
into reproducible results; certify selected statements when possible; find
and preserve counterexamples when statements fail; identify precisely what
remains open. Return a complete, independently auditable English handoff.

The plan contains **128 individually identified tasks in 16 work packages**.
Some tasks are conditional on earlier discoveries. “Complete the plan” means
that every task has an honest terminal disposition and associated evidence,
not that every conjecture must be proved. A refutation, an explicitly bounded
unsuccessful computation, or a documented missing prerequisite can be the
correct result. Never fill a dependency gap with an assumption and later
describe the conclusion as unconditional.

The principal scientific deliverables are:

1. A strengthened, independently replayable exact-certificate foundation.
2. Certified local stationary/coexistence/fold information where feasible.
3. A separate resolution or accurately bounded status of the **four-amplitude**
   Ising, intraparity, and off-face curvature problems.
4. A reconstruction of the phase calculations with complete root coordinates,
   derivative conventions, and, where feasible, proof-grade root/cover data.
5. A full-seven-coordinate stationary atlas with rigorously scoped indices,
   without resurrecting the already-refuted everywhere curvature bound.
6. A carefully limited global rank-seven variational result, if the preceding
   information actually supplies the necessary hypotheses.
7. A next handoff containing source code, proof objects, failures, and exact
   replay commands—not just summaries, plots, or counts.

Do not confuse an interesting conditional landscape with a strict-derived
physical law. No task below authorizes claiming a source of active gain,
selector closure, physical units, apparatus, SM/GR, or a theory of everything.

## 2. Mandatory reading and source precedence

Read in this order before implementing mathematical extensions:

1. The applicable current [AGENTS.md](AGENTS.md), including the latest
   “Local pre/post-Discord handoff intake and rank-seven scope correction”
   section and the relevant earlier ST293, ST344–ST361, ST389–ST410,
   ST458/ST459, and quantum-package guardrails.
2. [Intake report](fin_handoff_audit/REPORT.md), completely.
3. [Accepted corrected proofs](fin_handoff_audit/PROOF.md), completely.
4. [Audit implementation](fin_handoff_audit/research.py),
   [tests](fin_handoff_audit/test_research.py), and
   [verifier](fin_handoff_audit/verify.py).
5. The recorded `fin_handoff_audit/results.json` and `verification.json`.
6. The entire supplied
   [handoff](FIN_research_artifacts_pre_and_post_Discord/FIN_full_chat_research_handoff_pre_and_post_Discord.md)
   and its [manifest](FIN_research_artifacts_pre_and_post_Discord/MANIFEST.txt).
7. Only the additional historical source files needed by the selected task.
   Discover their real paths with `rg`; do not invent a filename from a program ID.

Precedence: current repository guardrails and accepted corrected proofs outrank
the imported handoff's informal status labels. The imported bundle is primary
evidence of what was reported, not automatic evidence that it was proved.
If sources genuinely conflict, record the conflict and inspect the underlying
argument. Do not resolve mathematical disagreements solely by file date.

Do not browse or re-scrape the shared ChatGPT conversation. The local package
is the authorized research input. If later literature is needed, use the
environment's applicable browsing skill and primary sources only. Literature
access is optional for most core tasks; never pretend to have read an
unavailable source.

This campaign should normally use a new root-level package. If work actually
enters `fundamental_action_reconstruction/`, first read all six mandatory
documents listed at the top of AGENTS.md. This plan does not reopen frozen
upstream FAR lanes by implication.

## 3. Frozen state: accepted, numerical, open, and false

| Item | Starting status | Permitted use |
|---|---|---|
| Strict W, A, sector multiplicities and outward spectral provider | Existing accepted inputs | Reuse with provenance and replay. |
| Weighted edge decomposition 11+55 | Analytic identity, numerical reconstruction available | Structural input, not a new information substrate. |
| Passive Schur positivity | Accepted in the specified passive model | No negative static stiffness from this elimination. |
| Rank-seven vertex/uniform budget threshold at g=4 | Exact interval-certified | Not the first global transition or an information-dimension theorem. |
| Exact seven-coordinate Gibbs duality | Accepted | Do not constrain the probability state to the active linear subspace. |
| CRT/Ising representation and probability identities | Accepted in the specified sector | Not a global phase-reduction or general correlation-inequality theorem. |
| Local crossing near g=3.71834489812038 | Independently reconstructed numerical root | Seed for a new certificate, not yet a first-global-transition theorem. |
| Four-amplitude saddle near the crossing | Numerical reconstruction | Seed; inspect the full seven-coordinate Hessian separately. |
| Corrected one-dimensional resolvent positivity | Exact strict-interval face certificate | Preserve r=sech, not exp. |
| Extreme face s4=s5=0, s3,s6>=0 | Exact strict-interval second-curvature certificate | Equality only on its stated compactified boundary. |
| Global boundary-Ising curvature bound | Open: no replayable complete certificate supplied | Must rebuild, prove, or refute. |
| Dominant-mass/intraparity-W global closure | Open and dependency-sensitive | Do not promote the imported summary table. |
| Positive-orthant four-amplitude global ceiling | Open | Includes off-face and limit-domain obligations. |
| Universal cooperativity/isotonicity | Not accepted by the intake audit | Prove exact hypotheses or provide a counterexample. |
| Sixty-point phase census | Imported numerical claim | Reconstruct roots, multiplicities, derivatives, and complement coverage. |
| Everywhere full-seven-coordinate Hessian index <=1 | **Refuted** | Keep the explicit regression counterexample; do not try to prove it. |
| Index <=1 at all stationary points | Separate open question | The nonstationary counterexample alone does not decide it. |
| Global rank-seven minimizer orbit at the crossing or g=4 | Not closed by this intake | No transfer from the full-rank positive-edge result. |
| Quantum separability/discord results | Existing separate conditional packages | Preserve distinctions between marginals, correlations, and access models. |

The accepted counterexample is especially important: at
`h_j=2 cos(pi*j/2)` two orthogonal full-seven-coordinate covariance directions
are each greater than `313/960`. Thus the seven-coordinate dual Hessian has
at least two negative directions for every supplied `g>=3.7`. At the reported
crossing the two eigenvalues are numerically about `-0.06125315`. This field
has not been asserted stationary. The four-coordinate restriction can omit
one of these directions.

## 4. Notation contract: do not overload symbols

Use the following names in new code and reports even when older files use
shorter ambiguous symbols:

| Name | Meaning |
|---|---|
| `W_kernel` | The real symmetric 12-by-12 strict weighted adjacency matrix. |
| `A_full` | `diag(W_kernel 1)-W_kernel`. |
| `lambda[k]` | Fourier eigenvalue of A_full; lambda[0]=0. |
| `A7` | Top-sector rank-seven positive-semidefinite mediator, not a positive-edge graph Laplacian. |
| `X7` | 12-by-7 normalized Fourier feature matrix with A7=X7 X7^T. |
| `C4` | Four cosine/alternating columns of X7, in order 3,4,5,6. |
| `theta7` | Full Cartesian dual coordinate. |
| `s4` or `s_amp` | Four-dimensional amplitude coordinate; do not confuse the vector name with component s_4. |
| `J3,J4,J5,J6` | Coefficients of the unscaled trigonometric observables in h. |
| `p` | Twelve-label classical probability distribution, not a density matrix. |
| `rho_C` | Quantum density `I/12+W_kernel/20`, when that separate lane is used. |
| `W_par` | Intraparity covariance part; never call it W_kernel. |
| `q_even` | Probability of the even parity class. |
| `M4`,`M7` | Covariances of C4 and X7 under the declared p. |
| `H4`,`H7` | Cartesian dual Hessians `I/g-M4`, `I/g-M7`. |
| `sigma` | Candidate global four-coordinate second-curvature ceiling; accepted on the stated face only. |
| `index(H)` | Number of strictly negative eigenvalues. Call it Morse index of a critical point only when stationarity is established. |
| `g_eq` | Gain at the reconstructed/certified local equal-energy branch event. |
| `g_global` | First gain with a global nonuniform competitor, only if separately defined/proved. |

For j=0,...,11, the X7 column order is

```text
(3cos, 3sin, 4cos, 4sin, 5cos, 5sin, 6alt).
X7[j,kcos] = sqrt(lambda_k/6) cos(2*pi*k*j/12), k=3,4,5.
X7[j,ksin] = sqrt(lambda_k/6) sin(2*pi*k*j/12), k=3,4,5.
X7[j,6alt] = sqrt(lambda_6/12) (-1)^j.
C4 = X7[:, (0,2,4,6)].
```

Consequently `Jk=sqrt(lambda_k/6)*s_k` for k=3,4,5, and
`J6=sqrt(lambda_6/12)*s_6`. The field norm, the theta norm, the s norm,
and amplitudes of an alternative complex Fourier convention are not
interchangeable.

Eigenvalue notation in curvature tasks is decreasing:
`lambda_1(M)>=lambda_2(M)>=...`. Spectral sector labels `lambda_k(A_full)`
are labels, not that ordering convention. Write the matrix argument where
there is any ambiguity.

## 5. Execution contract for a lower-reasoning-budget agent

### 5.1 Standard loop for every task

1. Read the task card and the named prerequisites.
2. State the exact domain, assumptions, target claim, and proof level in a
   machine-readable task record **before** implementing the search.
3. Produce the smallest diagnostic or symbolic identity that can falsify the
   proposed route. Use fixed seeds for stochastic exploration.
4. Separate discovery code from the checker. Numerical optimizer output is
   a candidate for certification, never the certificate itself.
5. Save the entire candidate or proof object. Preserve failed candidates too,
   tagged with the reason for rejection.
6. Run positive tests, negative controls, and the appropriate regression subset.
7. Write a short result note: what changed, what did not change, and the exact
   next missing premise. Update dependencies, then choose the next ready task.

If three attempts fail for the same mathematical reason, do not merely change
seeds, increase precision, or extend a grid. Name the obstruction and take the
task's alternative branch. A genuinely new coordinate chart, inequality,
provider, or exact factor can justify another attempt; repetition alone cannot.

### 5.2 Status vocabulary

Every task has one execution status and a separate scientific claim status.

Execution: `NOT_STARTED`, `RUNNING`, `DONE`, `RESOURCE_STOP`,
`BLOCKED_DEPENDENCY`, `SUPERSEDED_BY_COUNTEREXAMPLE`, `OUT_OF_SCOPE`.

Scientific: `EXACT_PROVED`, `INTERVAL_CERTIFIED`, `NUMERICAL_REPRODUCED`,
`NUMERICAL_NEW`, `COUNTEREXAMPLE_CERTIFIED`, `COUNTEREXAMPLE_NUMERICAL`,
`CONDITIONAL_LEMMA`, `UNRESOLVED`, `NOT_A_SCIENTIFIC_CLAIM`.

`DONE + UNRESOLVED` is allowed when the prescribed bounded investigation
was completed and its limitations are honestly recorded. `DONE` alone must
never be read as “the conjecture is true.”

### 5.3 Resource classes and safe stopping

These are proposed conservative defaults for the later execution campaign,
not permission to run anything now or exceed environment limits:

| Class | Typical task | Default maximum per bounded run |
|---|---|---|
| S | Algebraic identity, unit test, low-dimensional diagnostic | 5 minutes, 2 GiB. |
| M | Root isolation, moderate symbolic reduction, structured sample | 20 minutes, 4 GiB. |
| L | Adaptive cover, interval atlas, complete finite census | 60 minutes, 6 GiB, resumable checkpoints. |

Use one heavy process at a time on the recorded 16-GB host unless the current
environment and user explicitly authorize another arrangement. No unlimited
Groebner campaign, duplicate worker, unattended installation, network upload,
external compute purchase, or automatically restarted stopped campaign.
Checkpoint long runs every 30–60 seconds; capture exit status and last complete
checkpoint. A parent command timeout must not leave a heavy child alive.

For initial covers use a 10,000-leaf pilot. A later bounded pass may use up to
200,000 leaves if the pilot justifies it. Larger jobs require an explicit
resource decision in the task record and compliance with user/environment
authority. A cap is a method/resource stop, not a counterexample.

The lanes can be developed independently where dependencies permit, but the
default execution is serial. This document is not an instruction to spawn
subagents; obey the executing environment's delegation policy.

### 5.4 Keep a small working memory without losing scientific context

At the start of each execution turn, load only the state map, current task
card, required proof dependencies, and the current candidate/checkpoint after
the mandatory initial reading is complete. Maintain this short working note:

```text
Current task:
Exact target and quantifiers:
Coordinate/state space:
Already paid premises:
Unpaid premises (must not be assumed):
Current smallest check:
Result file/checkpoint:
Next action if the check passes / fails:
```

Read a prerequisite's actual proof or certificate interface before using it;
do not rely on a previous agent's one-line paraphrase. A completed prerequisite
task with scientific status UNRESOLVED does not supply a proved lemma.
However, it can supply a useful counterexample, partial domain, or failure
reason. Continue an affected task only as an explicitly conditional or
dependency-independent route.

Use this selection rule:

```text
Choose a NOT_STARTED task whose required inputs actually exist.
Prefer: new exact witness -> small local certificate -> missing boundary
lemma -> bounded global cover -> optional interpretation.
If a required theorem is false, retire every argument that needs it.
If a required theorem is unresolved, either take an independent route
already allowed by the task or mark BLOCKED_DEPENDENCY with the exact atom.
Never repeatedly select the same failed task without a new mathematical input.
```

### 5.5 Useful numerical seeds, explicitly not certificates

These are orientation-specific seeds from the accepted audit implementation.
Use its column convention and verify the residual; a D12-transformed
representative has different coordinates but the same invariant content.

```text
Local equal-energy gain candidate:
  g_eq ~ 3.7183448981203795
Localized C4 coordinate:
  (1.8199035812800786, 1.9139895546687176,
   1.9145691325468461, 1.3672032801954980)
Barrier C4 coordinate at that numerical gain:
  (0.9409570673214581, 1.0014394023775526,
   0.9621088536282435, 0.6864149839950574)
Localized maximum probability:
  ~ 0.8363652266510411
Imported fold gain candidate, not interval-certified by the intake:
  ~ 3.515644716839593
Corrected resolvent face seed:
  r_sech ~ 0.6608385948, s3 ~ 1.703819057,
  S ~ 0.057549460989
Extreme-face saturation constants:
  sigma ~ 0.267443244229, t_star ~ 0.4264779295
```

For a point-root proof, convert the proposed center and preconditioner to
explicit rationals and then certify. For a spectral inequality, obtain
constants from the exact interval provider, not this seed block. The exact
uniform spinodal formula is `12/lambda_6`; its decimal display is not an
independent theorem about dynamics.

If the executing agent receives only this plan and not the repository or
referenced input package, first request the missing local files. Do not
manufacture the prior audit, proof objects, or baseline test output from the
seed values printed here.

## 6. File layout and evidence contract

Use a new package, proposed as `fin_rank7_followup/`, to avoid contaminating
the immutable source bundle or the accepted intake baseline:

```text
fin_rank7_followup/
  README.md
  STATE_MAP.md
  CLAIMS.json
  TASKS.json
  WORKLOG.md
  config.json
  environment.json
  src/                  # Small modules, shared conventions.
  tests/                # Unit, property, mutation, and regression tests.
  proofs/               # Human-readable derivations with exact assumptions.
  candidates/           # Full numerical candidates, including rejected ones.
  certificates/         # Rational endpoints, root boxes, cover trees, witnesses.
  logs/                 # Commands, stdout/stderr, exit codes, elapsed resources.
  results/              # Derived result tables, never the sole proof objects.
  figures/              # Optional visual aids generated from saved data.
  verify.py             # Read-only verification by default.
  build_results.py      # Explicit regeneration, separate from verification.
  HANDOFF.md
  MANIFEST.sha256
```

Task evidence should be named `R7P-NNN_*`. Exact rational numbers should be
serialized as numerator/denominator pairs or strings, never silently converted
to JSON binary floats. Keep “display decimal” separate from “certified bound.”

Example task record:

```json
{
  "id": "R7P-025",
  "execution_status": "NOT_STARTED",
  "claim_status": "UNRESOLVED",
  "domain": "C4 Cartesian stationary equation near a supplied local seed",
  "assumptions": ["frozen strict spectrum", "supplied gain convention"],
  "dependencies": ["R7P-012", "R7P-024"],
  "inputs_sha256": {},
  "commands": [],
  "random_seeds": [],
  "outputs": [],
  "acceptance_checks": [],
  "rejected_alternatives": [],
  "next_missing_atom": null
}
```

A cover certificate must contain the original domain, exact split rules,
every leaf or a losslessly replayable tree, the accepted inequality on each
leaf, and every unresolved leaf. “200,000 boxes passed” without the boxes
and inequalities is insufficient. A root certificate must contain the exact
system, center, box, preconditioner, residual/Jacobian enclosures, strict
inclusion margins, and uniqueness scope. An eigenvalue claim must name the
matrix, dimension, ordering, domain, and treatment of multiplicities.

## 7. Dependency architecture and recommended scheduling

```text
A intake and reproducibility -> B coordinates and analytic infrastructure
                                      |
                  +-------------------+--------------------+
                  v                   v                    v
             C face upgrades     D local events       E full-7D tests
                  |                   |                    |
                  +--> F Ising algebra -> G boundary proof |
                                      |                   |
                              H intraparity bounds        |
                                      |                   |
                              I off-face 4D proof         |
                                                          |
B -> J cooperativity test                                 |
B -> K phase/cumulants -> L phase census -----------------+
                                                          v
                                               M full-7D global analysis

B -> N passive/finite-N interpretation     B -> O quantum comparison
                         all lanes -> P independent review and handoff
```

An arrow means mathematical dependency, not compulsory chronological blocking
of unrelated lanes. In particular, K/L need not wait for a successful I.
J is a hypothesis-validation lane, not an assumed input to every other lane.
N/O are secondary interpretive lanes: keep them bounded and do not let them
displace the core mathematical tasks merely because they are easier.

Recommended batches:

1. A and B; obtain a clean baseline and correct coordinate/derivative machinery.
2. C plus the cheapest D/E tasks; prioritize exact local wins and immediate
   stationary/full-dimension falsification.
3. F/G and J; independently reconstruct the missing boundary proof and test
   whether the claimed cooperative structure is available at all.
4. H/I, conditional on their actual dependencies; in parallel in the schedule,
   not necessarily in processes, pursue K/L phase reconstruction.
5. M only after the relevant local, phase, and domain tools are ready.
6. N/O in bounded secondary slots; P checkpoints throughout and final assembly.

Start the final handoff skeleton early. Do not wait until the end to discover
that root coordinates, seeds, rejected boxes, or source hashes were not saved.

---

## Work package A — Intake, provenance, and reproducibility

### R7P-001 — Freeze the actual starting state

**Depends:** none. **Class:** S.
Record `git status --short`, current commit if available, Python/platform and
installed numerical package versions. Read the mandatory sources. Hash the
accepted audit inputs and imported manifest. Separate pre-existing dirty files
from proposed new outputs; do not reset or clean them.
**Accept/output:** `environment.json` and a baseline provenance note, with no
unexplained source mutation. If the recorded host differs, revise resource
assumptions rather than blindly using the historical 16-GB number.

### R7P-002 — Replay the accepted audit without overwriting it

**Depends:** 001. **Class:** M.
Run `python3 fin_handoff_audit/verify.py`, not `--record`, first. Save stdout,
stderr, exit code, duration, and the exact command. The baseline expects 19
intake tests and 56 inherited regressions. Investigate failures as environment,
implementation, or mathematical discrepancies; do not regenerate expected
results to make a failure disappear.
**Accept/output:** a PASS replay or a precise discrepancy record blocking only
the affected dependencies. No new scientific result is counted for replay.

### R7P-003 — Convert the handoff ledger into explicit claims

**Depends:** 001. **Class:** S.
Create one claim record per substantive audited statement. Include quantifiers,
coordinate dimension, finite versus limiting domain, dependencies, and original
handoff section. Split compound claims such as “boundary proved, therefore
global” into separately testable atoms. Tag already-refuted full-7D claims.
**Accept/output:** `CLAIMS.json`; every claim has a status and a cited local
source. No optimistic status is inherited merely from a CSV filename.

### R7P-004 — Inventory missing proof objects

**Depends:** 003. **Class:** S.
List what is absent for the 60-root phase census, Ising cover, dominant-mass
bound, numerical fold, and global resolvent minimum. Distinguish missing source
code, missing raw output, missing exact inputs, and missing mathematical lemma.
Map each missing item to a rebuilding task below.
**Accept/output:** a missing-evidence matrix. It must explicitly say that a
summary count, plotted curve, or decimal coefficient table is not a replayable
global proof object.

### R7P-005 — Establish schemas and immutable/raw boundaries

**Depends:** 003–004. **Class:** S.
Create the proposed package, task schema, candidate schema, interval encoding,
and proof-certificate schema. Make raw imported artifacts read-only by policy;
new corrected results go elsewhere. Define an output-to-input provenance link
for every derived file.
**Accept/output:** schema validators and fixture tests rejecting missing domain,
missing quantifiers, or float-only “exact” endpoints. Do not require a complex
framework when a small standard-library validator suffices.

### R7P-006 — Separate verifier and result generator

**Depends:** 005. **Class:** S.
Implement a read-only checker entry point and a separately named regeneration
entry point. Verification must fail on altered proof inputs, invalid interval
endpoints, or unresolved cover leaves when a global PASS is requested.
**Accept/output:** command-line contracts and tests showing that verification
does not silently rewrite a baseline. Add a deliberate one-coefficient
mutation and confirm rejection before trusting later certificate workflows.

### R7P-007 — Make every expensive computation resumable

**Depends:** 005–006. **Class:** S.
Implement task IDs, fixed seeds, bounded batch sizes, periodic checkpoints,
lock/PID handling, and safe interruption. Simulate an interruption during a
small dummy cover and resume without dropping or duplicating leaves.
**Accept/output:** resumed and uninterrupted runs produce the same scientific
state. A stale lock has an explicit recovery procedure; a live lock prevents
duplicate heavy execution. Do not restart any historical stopped campaign.

### R7P-008 — Select the first live frontier

**Depends:** 001–007. **Class:** S.
Write `STATE_MAP.md` with ready tasks, blocked tasks, and why each proposed
next step adds a new witness, lemma, or certificate. Choose small C/D/E wins
before a global cover. Create all 128 task records, including conditional and
secondary ones, so later skipped work remains visible.
**Accept/output:** first batch and explicit stop/pivot rules; no generic
bridge-source, selector, or closed-lane replay is scheduled.

## Work package B — Coordinates, derivatives, domains, and proof primitives

### R7P-009 — Build and cross-check both feature spaces

**Depends:** 002,005. **Class:** S.
Construct W_kernel, A_full, X7, C4, and A7 independently from the frozen
definition. Check A7=X7 X7^T, zero column means, sector norms, and the exact
embedding C4=X7[:,(0,2,4,6)]. Test random theta and s fields under both routes.
**Accept/output:** a coordinate module and a mapping table. Deliberately swap
a sine-column sign and require a test to detect the convention mismatch.

### R7P-010 — Implement the D12 action explicitly

**Depends:** 009. **Class:** S.
Generate all 24 label permutations, induced real feature transformations,
and their actions on fields and probabilities. Verify group closure and
objective equivariance. Distinguish a reflection-fixed subspace from a
fundamental domain; the former does not cover generic states.
**Accept/output:** exact/trigonometric identities plus numerical round-trip
tests. Stabilizer detection initially remains numerical; certified stabilizers
require later exact invariance and separation arguments.

### R7P-011 — Re-derive primal/dual correspondence

**Depends:** 009. **Class:** S.
Derive joint minimization of `D(p||u)+||theta||^2/(2g)-theta^T X7^T p`.
Prove the softmax formula, stationarity correspondence, equality of global
infima, and the conditions under which boundary probability minimizers are
excluded. Keep g>0 explicit. Record what happens at g=0 separately.
**Accept/output:** a short theorem with an implementation identity test.
Do not call this a minimax theorem or restrict p-u to Range(A7).

### R7P-012 — Implement derivatives through the required orders

**Depends:** 009,011. **Class:** M.
Compute gradients and covariance Hessians using stable log-sum-exp. Implement
third and, where required, fourth directional derivatives as centered moments
or cumulants. Compare against symbolic small examples and finite-difference
diagnostics at moderate scales. Keep ambient, spherical, and phase derivatives
as separately named functions.
**Accept/output:** derivative tests and documented normalization factors.
Finite differences are diagnostics, not the proof-grade interval derivative.

### R7P-013 — Harden transcendental interval primitives

**Depends:** 006,009. **Class:** M.
Reuse the existing strict rational spectral enclosures. Add only the missing
sqrt, exp, log, tanh, and sech enclosures needed by later tasks, with rational
range reduction and explicit remainder bounds or an approved verified library.
Test negative arguments, near-zero values, large fields, and zero denominators.
**Accept/output:** exact endpoint tests and independent high-precision spot
checks. Higher precision alone is not outward rounding or a proof.

### R7P-014 — Prove the correct Hessian-index transfer at stationary points

**Depends:** 011–013. **Class:** M.
Use the joint primal/dual Hessian and Schur complements on the probability
tangent space to derive how negative directions of the 11D primal and 7D dual
correspond at a common interior stationary state. Identify any additional
strictly positive directions and all nondegeneracy assumptions.
**Accept/output:** an exact inertia theorem and tests. Do not transfer Hessian
signatures through nonlinear coordinates at noncritical points without the
extra gradient terms.

### R7P-015 — Derive compact domains before global searches

**Depends:** 011. **Class:** S.
Let `R=max_j ||X7_j||`. Derive `||theta||<=g R` at stationary points from
theta=g X7^T p. Derive a possibly larger sublevel enclosure using
`log(mean exp(X7 theta))<=R||theta||`. State which bound is for roots and
which is for global negative sublevels; account for zero/equality cases.
**Accept/output:** exact/coarsely interval bounds and normalized domain files.
Do not use the 4D amplitude box as a cover of the 7D stationary domain.

### R7P-016 — Set up asymptotic/domain charts honestly

**Depends:** 009,013,015. **Class:** M.
Enumerate relevant supports in positive-field limits, and separate that problem
from arbitrary signed full-7D limits. Check whether each proposed exponential
coordinate change is algebraic. In particular, odd-sector cos(5-angle) values
involve sqrt(3)/2: independently varying exp(J5/2) and exp(sqrt(3)J5/2)
relaxes the model unless their relation is enforced.
**Accept/output:** chart definitions, physical constraints, and forbidden
relaxations. Keep interval-transcendental charts when polynomialization is invalid.

## Work package C — Finish low-dimensional exact face results

### R7P-017 — Independently replay the two accepted face certificates

**Depends:** 006,013. **Class:** S.
Reconstruct their polynomials from the formulas rather than reading saved
Bernstein coefficients. Verify the exact endpoint factor 1-q, denominator
positivity, and every interval lower bound. Cross-check the corrected sech
parameter against a direct covariance calculation.
**Accept/output:** a second finite certificate implementation or a clearly
labelled dependent replay. Do not count identical copied code as independent
verification; preserve the exp-versus-sech rejection test.

### R7P-018 — Derive the resolvent derivative polynomial

**Depends:** 017. **Class:** S.
Write S=N/D for the corrected one-dimensional formula. Derive N'D-ND', factor
only factors known nonzero on the open interval, and record the reduced degree.
Do not assume the handoff's “quintic” before checking cancellations such as
the factor 1+r.
**Accept/output:** exact symbolic polynomial identities, denominator proof,
and coefficient intervals from the strict spectrum. Compare derivative signs
with direct differentiated formulas at fixed rational test points.

### R7P-019 — Certify uniqueness of the one-dimensional minimum

**Depends:** 018. **Class:** M.
Isolate the candidate r near 0.6608385948. Prove derivative sign on the
complement and a unique root in its box using interval monotonicity/Newton
or a coefficient-robust exact root-count method. A Sturm sequence computed
for rounded coefficients does not certify the strict transcendental polynomial.
**Accept/output:** one isolated minimum plus complete interval coverage, or a
precise unresolved interval list. Endpoint limits must be included explicitly.

### R7P-020 — Enclose the minimum value and original coordinates

**Depends:** 019. **Class:** S.
Propagate the isolated r interval through S(r), q=1/(1+r),
t=sqrt(1-r^2), and s3=arcosh(1/r)/sqrt(lambda3/6), using validated primitives.
Record interval widths and a conservative positive safety margin.
**Accept/output:** certified face minimum and location intervals if 019 passes;
otherwise only conditional bounds on the candidate box. No assertion that
this is the minimum over all four fields is allowed.

### R7P-021 — Audit every equality case of the extreme-face proof

**Depends:** 017. **Class:** M.
Track equality in the scalar-channel estimate, the 45 block, the physical
constraint x^2<=2q-1, and the Bernstein steps. Treat q=1, x=t_star as a
compactified point, not a finite s6 state. Check endpoints x=0 and q=1/2.
**Accept/output:** an equality-locus theorem or a corrected weaker statement.
Never infer uniqueness merely because one plotted curve touches the bound once.

### R7P-022 — Quantify a face gap away from the equality point

**Depends:** 021. **Class:** M.
On a declared compact subset excluding a rational neighborhood of the equality
locus, prove a positive gap below sigma. Near the equality locus, derive a
validated local expansion with remainders in physical coordinates. Allow a
piecewise bound instead of forcing an unnecessarily sharp closed form.
**Accept/output:** explicit domains and gap constants reusable by an off-face
cover. No eigenvalue differentiation across a double root without a matrix
or spectral-subspace treatment.

### R7P-023 — Verify local parity-mixing asymptotics

**Depends:** 012,021. **Class:** M.
At the double root, perform degenerate perturbation in epsilon=1-q. Derive
the two reported first-order shifts and certify their signs. Then distinguish
fixed-t shifts from the reoptimized-envelope slope near 0.1312828584;
the latter requires a separate local optimization argument.
**Accept/output:** exact first-order formulas, interval remainder bounds if
available, and a clear status for the reoptimized slope. Do not conflate them.

### R7P-024 — Export a reusable face theorem interface

**Depends:** 017–023, allowing unresolved subresults to remain tagged. **Class:** S.
Define a small checker API accepting certified spectral intervals and returning
only the conclusions actually verified. Include parameter-domain validation,
endpoint handling, and mutation tests that break one hypothesis at a time.
**Accept/output:** `face_certificate.json`, proof note, and checker. If uniqueness
or quantitative gaps remain open, export basic positivity/curvature conclusions
without pretending the stronger API fields were proved.

## Work package D — Certify local roots, coexistence, and folds

### R7P-025 — Reconstruct localized and saddle seeds

**Depends:** 009,012,024. **Class:** S.
Solve the four-amplitude stationarity equations near the recorded localized
and saddle seeds. Store full coordinates, p, gain, objective, gradient, H4,
H7, and primal tangent Hessian diagnostics. Repeat with an independently
implemented residual before any interval proof.
**Accept/output:** complete candidate records with residuals and conventions.
Solver success flags without residual verification are rejected.

### R7P-026 — Certify the local equal-energy event

**Depends:** 013,025. **Class:** M.
Use the five-equation system `(s/g-E[C4], Phi_g(s))=0` in `(s,g)`.
Choose a rational box around the recorded crossing and apply a validated
Krawczyk/interval Newton inclusion. Verify positivity of g and interior p.
**Accept/output:** a unique root in the stated five-dimensional box or an
honest failed inclusion. This proves a local equal-energy event, not the first
global transition or absence of another branch with lower energy.

### R7P-027 — Certify the localized root's full stability type

**Depends:** 014,026. **Class:** M.
At the entire root box enclose H4, H7, and the primal tangent form. Use
verified LDL/congruence or rigorous eigenvalue bounds. A positive H4 alone
does not prove positive H7. If the full form is positive, state local
minimality and its neighborhood scope.
**Accept/output:** an inertia certificate across all required directions, or
the specific uncertain eigenvalue interval. Do not accept floating stationarity
or positive sample eigenvalues as a theorem.

### R7P-028 — Prove local crossing transversality

**Depends:** 026–027. **Class:** M.
Derive the branch energy derivative using stationarity, in either primal or
dual coordinates, and certify a nonzero interval sign at the event. Establish
a locally unique crossing and, if affordable, a small connected parameter
tube with explicit overlaps.
**Accept/output:** local crossing direction, gain interval, and uniqueness
scope. Disconnected individual root boxes are not a continuous branch tube.

### R7P-029 — Certify the barrier saddle at a declared gain

**Depends:** 013,025. **Class:** M.
Choose either a rational gain near coexistence or a parametric treatment tied
to the certified event. Isolate the saddle and certify H4 and H7 signatures.
If using a rational gain, do not describe it as the exact event gain.
**Accept/output:** a local saddle certificate and the energy barrier relative
to clearly named states. It is not automatically the globally lowest mountain
pass or a dynamical transition rate.

### R7P-030 — Rebuild the numerical fold with the correct augmented system

**Depends:** 012,025. **Class:** M.
For the 4D stationary equation F, solve `F=0`, `D_s F v=0`, and a fixed
normalization of v, in nine unknowns `(s,g,v)`. Use the reported gain near
3.51564471684 only as a seed. Check the other eigenvalues and numerical
transversality before attempting an interval proof.
**Accept/output:** full augmented candidate, normalization, Jacobian, and
residual—not only p_max and a near-zero Hessian eigenvalue.

### R7P-031 — Certify a simple local saddle-node if possible

**Depends:** 013,030. **Class:** L.
Isolate the augmented root, prove the stationary Jacobian has exactly one
kernel direction, and enclose the two nonzero fold coefficients involving a
left nullvector, F_g, and F_ss[v,v]. Check full-7D transverse directions
separately before any full-model interpretation.
**Accept/output:** a scoped simple-fold theorem or separate statuses for root
existence, kernel dimension, and transversality. Failure of one stage cannot
be hidden behind success of another.

### R7P-032 — Assemble a local branch diagram with certified joins

**Depends:** 026–031 as available. **Class:** L.
Connect certified local roots by validated continuation, using overlapping
root/tube certificates and explicit joins. Stop at the declared resource cap.
Distinguish the local fold, local energy crossing, and uniform spinodal
12/lambda6. Preserve unbridged intervals rather than drawing a continuous
certified curve across them.
**Accept/output:** a local diagram with each segment labelled certified,
numerical, or missing. No claim of exhaustive hysteresis or physical time.

## Work package E — Test the stationary-only full-seven-coordinate problem

### R7P-033 — Permanently register the known full-7D counterexample

**Depends:** 009,013. **Class:** S.
Reproduce the two exact trial directions at `h=2 cos(pi*j/2)` and the lower
bound 313/960. Compute the stationarity residual to demonstrate why this
does not itself settle a critical-point-only theorem.
**Accept/output:** a compact analytic witness and regression test. Any later
claim of an everywhere full-7D index-one theorem must fail automatically against
this witness, even if a restricted four-dimensional test suite passes.

### R7P-034 — Derive the invariant two-harmonic stationary system

**Depends:** 010–012. **Class:** S.
Study `h=J cos(pi*j/2)+K(-1)^j`. Derive, rather than merely assume,
`x=sinh(J)/(cosh(J)+exp(-2K))`,
`y=(cosh(J)-exp(-2K))/(cosh(J)+exp(-2K))`.
Verify the full stationarity reduction
`J=g(lambda3/6)x`, `K=g(lambda6/12)y`, and all omitted residuals.
**Accept/output:** a two-variable exact reduction with parameter domain and
symmetry justification. This is a particularly cheap route to testing the
stationary-only question, not a promised counterexample.

### R7P-035 — Search the two-harmonic stationary branches

**Depends:** 034. **Class:** M.
At g=37/10, the local crossing neighborhood, g=4, and selected nearby gains,
find all candidates in the compact bounds inherited from stationarity. Include
zero, signed J, and any K branches allowed by the declared domain. Compute
full H7 signatures, especially sine directions omitted by C4.
**Accept/output:** reproducible roots and possible index-two-or-higher witnesses.
A finite scan is not an exhaustive theorem unless its complement is covered.

### R7P-036 — Certify or bound the two-harmonic outcome

**Depends:** 013,035. **Class:** M/L.
If a nondegenerate higher-index stationary candidate exists, isolate it and
certify two negative H7 directions. If none is found, prove a bounded
two-variable exclusion or record the unresolved boxes. Use the full residual
identity from 034 to lift a certified reduced root.
**Accept/output:** a genuine stationary counterexample, or an explicitly scoped
restricted-family result. Never declare the full stationary-only conjecture
true merely because this family yielded no counterexample.

### R7P-037 — Build a reproducible full-7D stationary discovery atlas

**Depends:** 010,012,015. **Class:** L.
Use a fixed mixture of symmetry-generated seeds, low-mode seeds, Sobol points,
and seeded random starts within the stationary ball. Store failed solves and
deduplicate under all 24 D12 actions using a stated tolerance.
**Accept/output:** full root vectors, orbit candidates, objective values,
H7 spectra, and coverage statistics. “No new roots in another batch” is
numerical saturation only, not stationary exhaustion.

### R7P-038 — Certify representative roots and orbit separation

**Depends:** 013,037. **Class:** L.
Prioritize roots with unexpected index, trivial stabilizer, or competitive
energy. Certify root boxes and Hessian inertia, then prove boxes from distinct
candidate orbits cannot overlap under D12. Exact stabilizers need invariance
and exclusion of additional group elements.
**Accept/output:** a partial certified atlas with exact orbit counts only for
the isolated components. A representative list is not a global root count.

### R7P-039 — Classify restricted/full stability mismatches

**Depends:** 027,029,038 as available. **Class:** S.
For every certified or numerical reflection-fixed root, compare H4 with H7
and decompose the latter into even/odd reflection blocks. Explain whether
extra negative modes are phase-like, amplitude-like, or mixed in the chosen
Cartesian basis.
**Accept/output:** a mismatch table with proof level per row. Do not describe
a restricted minimum as a full minimum when the odd block is negative.

### R7P-040 — Decide the next stationary-only theorem target

**Depends:** 033–039. **Class:** S.
If a stationary index>=2 witness is certified, retire the universal
stationary-only index-one conjecture and formulate narrower valid questions.
Otherwise select a precise parameter/domain restriction for further proof,
with an explicit exclusion obligation.
**Accept/output:** a claim-status update and a non-repetition decision. This
decision does not affect the already-valid four-amplitude face theorem.

## Work package F — Reconstruct the missing boundary-Ising algebra

### R7P-041 — Rebuild the exact four-state boundary model

**Depends:** 009,013,016. **Class:** S.
Enumerate the six even labels, aggregate their multiplicities into four
states `(Aspin,Y)`, and verify the observables and Ising weights directly.
Preserve the degeneracy term `H_Y=3J4/4-(log 2)/2`.
**Accept/output:** exact four-state features, weights, and forward/inverse
parameter maps for strictly positive probabilities. A two-spin model with
uniform degeneracies is a different model and must fail a comparison test.

### R7P-042 — Characterize the admissible probability domain

**Depends:** 041. **Class:** M.
Derive the three polynomial inequalities from positive probability ratios,
normalization, and nonnegative J3,J4,J5. Prove the positive-interior converse
using log ratios. Investigate zero-probability boundary points separately:
which are actual limits, and which merely satisfy necessary inequalities?
**Accept/output:** an exact domain theorem or a labelled outer relaxation.
Any counterexample found only in an unproved relaxation is not yet a
counterexample to the original exponential family.

### R7P-043 — Derive the covariance polynomial invariants

**Depends:** 041–042. **Class:** M.
Symbolically construct the 3-by-3 covariance and its characteristic polynomial
in three independent probabilities. Derive trace, second elementary symmetric
coefficient, and determinant. Verify the determinant factor
`3 lambda3 lambda4 lambda5 p1 p2 p3 p4/8` independently.
**Accept/output:** factored polynomial expressions and coefficient provenance.
Square-root scalings may be removed by a justified congruence where useful;
remember that congruence preserves inertia, not eigenvalues themselves.

### R7P-044 — Design a logically sufficient eigenvalue-count test

**Depends:** 043. **Class:** M.
Translate `lambda2(M)<=sigma` into a sound inertia or shifted-characteristic
criterion. One candidate is a sufficient Descartes sign-variation bound for
`det((sigma+z)I-M)`; derive its exact implication and equality handling before
coding. Alternatively use verified LDL with pivot charts.
**Accept/output:** a standalone lemma tested on matrices with zero, one, two,
and three eigenvalues above sigma, including threshold multiplicities. Do not
use the sign of the determinant alone to count supercritical eigenvalues.

### R7P-045 — Prove the special double-root identities exactly

**Depends:** 021,043–044. **Class:** S.
Substitute the reported double-root probabilities expressed in t_star into
the covariance invariants. Prove the two threshold eigenvalues and identify
the remaining one. Verify all active physical-domain constraints.
**Accept/output:** an exact equality certificate retaining the dependence of
sigma and t_star on the same spectral parameters. Replacing these dependent
quantities by unrelated decimal constants is not an exact identity test.

### R7P-046 — Enumerate boundary strata of the probability domain

**Depends:** 042–045. **Class:** M.
List p_i=0 faces and active-constraint intersections, their attainable support
ranks, and which require separate parameter limits. Analyze each reduced
covariance before launching a bulk interior cover.
**Accept/output:** a finite stratum table with proofs or unresolved obligations.
If a lower-dimensional stratum has a simple analytic bound, export it as a
boundary certificate reusable by the global checker.

### R7P-047 — Recover and certify relaxed-domain counterexamples

**Depends:** 042–044. **Class:** M.
Independently relax each physical constraint and search for lambda2>sigma.
Save full probabilities and reconstruct fields whenever possible. Rationalize
a witness and certify the eigenvalue violation plus the exact violated
physical constraint.
**Accept/output:** negative controls demonstrating why domain constraints
matter. Failure of an enlarged domain is a method warning, not a failure of
the properly constrained model.

### R7P-048 — Freeze a boundary proof specification

**Depends:** 041–047. **Class:** S.
Choose the exact domain, polynomial criterion, equality neighborhood, boundary
lemmas, and cover representation for G. Define the expected leaf outcomes
and the treatment of dependency uncertainty. Estimate complexity from a
small pilot rather than planning an unbounded elimination.
**Accept/output:** a versioned proof-specification JSON and a mathematical note.
If the target was refuted in the admissible domain, retire G's positive-proof
route and redirect it to counterexample verification.

## Work package G — Complete or refute the boundary-Ising certificate

### R7P-049 — Implement a small auditable Bernstein engine

**Depends:** 006,013,048. **Class:** M.
Support the required low-degree multivariate polynomials on boxes or simplex
charts. Verify power-to-Bernstein conversion against exact point evaluations
and symbolic examples. Use exact rational interval coefficients for uncertain
strict inputs; save chart transformations.
**Accept/output:** conversion and subdivision tests, including a polynomial
with a known interior negative region. A tool that certifies only positive
examples is insufficiently tested.

### R7P-050 — Prove and test the leaf classifier

**Depends:** 044,046,049. **Class:** M.
Each leaf must be classified by an explicit implication: outside the physical
domain, covered by a boundary lemma, satisfies a sufficient inertia condition,
belongs to the equality neighborhood, or remains unresolved. Check interval
sign directions carefully.
**Accept/output:** leaf reason codes, inequalities, margins, and adversarial
fixtures. Never accept a leaf merely because an optimizer found no violation
inside it or its midpoint has the desired sign.

### R7P-051 — Reproduce a bounded global boundary cover

**Depends:** 048–050. **Class:** L.
Run the 10,000-leaf pilot, then at most the documented next bounded pass.
Record the complete split tree and every unresolved leaf. Compare the actual
unresolved hull with the imported “collapse onto double root” claim without
requiring the old numerical counts to match.
**Accept/output:** replayable coverage statistics and unresolved geometry.
Persistent remote leaves indicate a missing inequality or a possible witness;
they may not be deleted as presumed numerical noise.

### R7P-052 — Derive a local equality-neighborhood certificate

**Depends:** 045,051. **Class:** M/L.
Use coordinates adapted to active constraints near the double root. Derive
the physical tangent cone and the leading polynomial terms. Seek an exact
factorization, copositive form, or interval remainder bound that controls the
whole neighborhood. Do not differentiate an individual ordered eigenvalue
through the double eigenvalue as if it were simple.
**Accept/output:** a local neighborhood theorem with explicit radius and
remainder bounds, or a list of uncontrolled cone directions.

### R7P-053 — Resolve remote ambiguous components

**Depends:** 051–052. **Class:** L.
For unresolved components away from the equality neighborhood, try alternative
inertia charts, lower-dimensional boundary reduction, or validated local
optimization/KKT isolation. Treat imported numerical distances and derivative
values only as hints.
**Accept/output:** certified resolution of each component or retained exact
unresolved boxes. A KKT point does not exhaust a constraint component unless
the component and its boundary are independently covered.

### R7P-054 — Check the assembled certificate independently

**Depends:** 051–053. **Class:** M.
Build a checker that reads only the frozen proof specification and the saved
certificate tree, recomputes every inequality, and verifies coverage. Corrupt
a split, remove a leaf, flip a sign, and alter one spectral interval to test
rejection.
**Accept/output:** independent replay logs and mutation failures. If any
unresolved leaf remains, the result must be a partial certificate, not a
global PASS with a caveat hidden in prose.

### R7P-055 — Export the boundary theorem or certified refutation

**Depends:** 054. **Class:** S.
State exactly whether the bound holds on the original attainable domain,
an explicitly larger semialgebraic domain, or only covered subdomains. Include
all equality cases and whether they occur at finite fields or only limits.
**Accept/output:** one theorem/refutation record with links to every proof
dependency. If incomplete, export a precise partial result and identify the
smallest remaining mathematical atom rather than saying “essentially proved.”

### R7P-056 — Propagate the boundary outcome through dependencies

**Depends:** 055. **Class:** S.
Update H/I and the overall state map. A certified boundary counterexample may
refute the proposed global four-amplitude ceiling; verify attainability before
making that inference. An unresolved boundary proof blocks only conclusions
that depend on it, not the already-accepted face theorem or independent phases.
**Accept/output:** a dependency impact report and a ready-task queue reflecting
the actual result, not the original expected proof path.

## Work package H — Intraparity complementarity and dominant-mass bounds

### R7P-057 — Derive conditional sector distributions from shared fields

**Depends:** 009,016,041. **Class:** M.
Compute C_plus, C_minus, their means, and q_even directly from the same
J3,J4,J5,J6. Prove which quantities are independent of J6. Derive the proposed
q0>=1/2 at J6=0 or find a counterexample; do not silently import it from the
handoff's cooperative interpretation.
**Accept/output:** exact formulas and a proof/status for the parity-weight
bound, with all shared-field constraints preserved.

### R7P-058 — Reconstruct the C_minus covariance domain

**Depends:** 057. **Class:** M.
Derive the three-probability `(u,d)` parameterization and its actual domain
from nonnegative shared fields. Verify every factor in the displayed 2-by-2
covariance. Test whether the claimed supremum
`(3 lambda4+lambda5)/32` is attainable or only an upper envelope.
**Accept/output:** a scoped exact supremum/bound theorem or a numerical
counterexample. Distinguish arbitrary three-probability distributions from
the physically constrained parameter subset.

### R7P-059 — Derive the dangerous-sector condition exactly

**Depends:** 058,013. **Class:** M.
Prove the smaller C_minus eigenvalue is below sigma in the stated domain,
then derive the determinant inequality equivalent to its larger eigenvalue
being at least sigma. Solve for the admissible d^2 threshold and u interval,
with denominator and endpoint signs certified.
**Accept/output:** an exact dangerous-set description. Squaring or clearing
denominators without sign checks invalidates the equivalence and must be
caught by tests against direct covariance evaluations.

### R7P-060 — Reconstruct the dominant C_plus mass from the same fields

**Depends:** 057–059. **Class:** M.
Derive the lower bound on the largest conditional even-sector probability
implied by a dangerous C_minus. Recover the one-dimensional optimization
leading to the reported candidate near u=0.5280356754, if the reduction is
valid. Save the exact expression rather than only its minimum.
**Accept/output:** a proof of the reduction or the missing premise that prevents
it. Optimizing two parity classes independently is not an admissible substitute.

### R7P-061 — Certify a sufficient coarse dominant-mass bound

**Depends:** 059–060. **Class:** M/L.
Prefer a conservative bound with visible margin over chasing the exact
minimum. Test the imported bounds d/u>=0.8515 and the ratio bound >=11.45
with strict intervals on the full dangerous set. Recompute their implication
for the dominant probability; do not simply copy 0.7112098557.
**Accept/output:** a complete one-dimensional/constraint cover and a rational
lower bound, or exact unresolved intervals and a weaker valid substitute.

### R7P-062 — Prove the covariance envelope under a dominant mass

**Depends:** 043,061. **Class:** M.
Derive the claimed upper envelope for the second elementary symmetric
covariance invariant with one probability at least alpha. Verify any
monotonicity only on its stated domain. Use `lambda2<=sqrt(e2)` only after
PSD and eigenvalue ordering justify it.
**Accept/output:** an exact constrained envelope or a validated bound with
all extremizer branches considered. A numerical maximizer of e2 is not
automatically its global envelope.

### R7P-063 — Assemble the intraparity Weyl argument

**Depends:** 055,057–062. **Class:** M.
Split on whether lambda1(C_minus)<=sigma. In the dangerous case combine the
proved dominant-mass result with the C_minus bound and the actual admissible
q_even interval. Check both endpoints of the affine Weyl combination and
every condition needed for the proposed q>=1/2 shortcut.
**Accept/output:** `lambda2(W_par)<=sigma` in an explicitly declared domain,
or a conditional lemma naming each unpaid prerequisite.

### R7P-064 — Falsify the decoupling shortcuts and export the result

**Depends:** 057–063. **Class:** M.
Search for violations when shared-field constraints or the lower q bound are
removed, and certify simple witnesses where possible. Package the valid
intraparity theorem separately from the invalid relaxed models.
**Accept/output:** a regression set showing that physical coupling constraints
are actually checked. If H remains unresolved, do not label I's input as a
proved intraparity bound merely because its numerical margin looks large.

## Work package I — The genuinely open off-face four-amplitude problem

### R7P-065 — Define the exact off-face target and stable representation

**Depends:** 016,024,056,064. **Class:** S.
State whether the target is lambda2(M4)<=sigma or an equivalent Mtilde inertia
bound, on `s3,s4,s5,s6>=0`. Use the Schur denominator eta>0, not a full
resolvent inverse near singularities. If H/G are unresolved, mark their use
conditional or choose a direct method that does not need them.
**Accept/output:** one precise target specification and a dependency-clean
choice between the direct and factored proof routes.

### R7P-066 — Perform a bounded adversarial off-face search

**Depends:** 012,065. **Class:** M.
Use fixed-seed multistart/DE and structured slices near the known face,
double-root neighborhood, dangerous C_minus region, and large-field limits.
Save full candidate fields and physical constraints. Directly recompute the
covariance at every putative violation in an independent implementation.
**Accept/output:** certified-candidate requests or a labelled numerical search
ledger. Search saturation does not establish a global ceiling.

### R7P-067 — Validate compactification and asymptotic boundary coverage

**Depends:** 016,046,065–066. **Class:** M/L.
Construct a finite domain-plus-tail decomposition or constrained compactified
charts. Preserve irrational exponential relations and shared fields. For
each tail chart, bound discarded probabilities and covariance perturbations;
include simultaneous and multiscale large-field limits.
**Accept/output:** a proved covering of the intended unbounded domain or a
finite-domain theorem specification. “All s_i<=12” is not the whole orthant.

### R7P-068 — Certify a local physical-cone neighborhood of the extreme face

**Depends:** 022–023,065,067. **Class:** L.
Near the saturation point use matrix blocks, spectral projectors, or a
polynomial inertia criterion with validated remainders. Test one-sided
s4,s5 directions and their mixed terms on the nonnegative cone. An indefinite
ambient Hessian can still be positive on this cone, but that needs a proof.
**Accept/output:** an explicit neighborhood certificate or the unresolved cone
direction. Do not infer global transverse monotonicity from this local result.

### R7P-069 — Build a finite off-face complement cover

**Depends:** 065,067–068. **Class:** L.
Cover the remaining compact domain using verified inertia charts or
Bernstein bounds where genuinely polynomial. Use previously proved face
gaps to prune cells. Save complete coverage and unresolved leaves, with
adaptive split criteria tied to actual interval uncertainty.
**Accept/output:** a bounded partial or complete certificate. Avoid blindly
uniform 4D grids when dependency-aware subdivision is available.

### R7P-070 — Certify any real violation before revising the ceiling

**Depends:** 066 or unresolved candidates from 069. **Class:** M.
If lambda2>sigma appears, enclose the field and covariance and certify two
supercritical directions, or a suitable exact inertia witness. Verify
nonnegative fields and the original shared-field model. Determine whether
the violation is finite or only a limiting construction.
**Accept/output:** an admissible counterexample or a rejected numerical
artifact. Do not move sigma upward and continue as though the original
conjecture had been proved.

### R7P-071 — Assemble the strongest justified 4D theorem

**Depends:** 067–070 and all actually used lemmas. **Class:** S.
Possible outcomes include the original global ceiling, a weaker ceiling
sufficient on a certified gain interval, a compact-domain theorem, or a
counterexample. Compare any proved ceiling with the reciprocal of the actual
certified gain interval, not a floating g_eq.
**Accept/output:** a precisely scoped four-amplitude Hessian consequence.
It does not determine full H7 indices, global minimizers, or a physical selector.

### R7P-072 — Record equality, optimality, and nontransfer boundaries

**Depends:** 071. **Class:** S.
If a bound is proved, identify whether it is sharp, attained, or only approached
at infinity; identify all proved equality cases. List explicitly which phase
directions are omitted. If incomplete, give the exact residual domain and
best certified subresults.
**Accept/output:** the final I-lane status and a machine-readable prohibition
on promoting it to the already-refuted everywhere full-seven-coordinate claim.

## Work package J — Is the claimed cooperative structure actually valid?

### R7P-073 — Define the precise covariance-sign conjecture

**Depends:** 009,012,041. **Class:** S.
State which observables, which positive coupling domain, which normalization,
and which covariance entries are claimed nonnegative. Separate nonnegative
means, pair covariances, and the Jacobian of the particular fixed-point map.
They are distinct propositions.
**Accept/output:** finite formulas for every required sign condition. The word
“ferromagnetic” and a general theorem name do not themselves pay the hypotheses.

### R7P-074 — Search systematically for sign violations

**Depends:** 073. **Class:** M.
Test small fields, coordinate faces, large-field limits, and fixed-seed
interior samples. Use higher-precision recomputation only to distinguish
roundoff from a candidate violation. Save exact inputs for every negative
entry and its relation to the proposed physical domain.
**Accept/output:** a reproducible sign table and candidate counterexamples.
If a counterexample exists, stop invoking universal isotonicity downstream.

### R7P-075 — Attempt a direct finite-group proof

**Depends:** 073–074. **Class:** M.
Expand the finite partition function or duplicated covariance expression
in a basis that might expose nonnegative coefficients. Track finite-group
character identities and any parity obstruction. If invoking an external
correlation inequality, check its exact hypotheses against this discrete
model rather than transferring an XY/Ising statement by analogy.
**Accept/output:** a direct proof, a correctly sourced applicable lemma, or
a precise reason the named general theorem does not apply.

### R7P-076 — Prove a restricted replacement if the universal claim fails

**Depends:** 074–075. **Class:** M.
Possible useful replacements are local sign positivity near certified roots,
selected Jacobian entries only, a face-restricted monotone map, or a weaker
invariant cone. Derive the smallest replacement actually needed by a later
task; do not invent a stronger claim to preserve the original narrative.
**Accept/output:** a scoped lemma/counterexample and a revised dependency map.
Purely numerical local positivity remains numerical until enclosed.

### R7P-077 — Certify fixed-point Jacobians at known roots

**Depends:** 026,029,073. **Class:** M.
At certified local roots enclose `DT=g Cov(C4)` and its spectrum. A positive
Perron vector may be inferred only under the required positivity/irreducibility
conditions. Distinguish an eigenvalue above one from convergence properties
of a separately specified iterative algorithm.
**Accept/output:** local Jacobian spectra, verified signs, and the conditions
under which the “collective escape direction” interpretation is legitimate.

### R7P-078 — Analyze monotone iteration only with paid premises

**Depends:** 075 or 076, plus 077. **Class:** M.
If a valid isotone map exists, construct explicit subsolutions/supersolutions
and a trapping order interval. Prove invariance and convergence to a fixed
point under the appropriate hypotheses. If isotonicity is unavailable, compare
ordinary fixed-point iteration numerically without calling it a monotone theorem.
**Accept/output:** a theorem for a declared algorithm or a bounded numerical
convergence study, not a sourced physical dynamical law.

### R7P-079 — Quantify alignment without overstating it

**Depends:** 025,077. **Class:** S.
Recompute the alignment between a saddle's leading direction and the
localized-minus-saddle displacement. Give interval enclosures if root/eigenvector
certificates permit them. Explain that near alignment does not prove an exact
one-dimensional invariant path or minimum-action transition route.
**Accept/output:** a reproducible scalar diagnostic and explicit nonconclusions.
Avoid treating an almost-unit cosine as an exact reduction theorem.

### R7P-080 — Close the cooperativity claim ledger

**Depends:** 073–079. **Class:** S.
Separate universally proved signs, locally certified signs, numerical signs,
and counterexamples. Remove unpaid uses of universal cooperativity from
all downstream notes and code assumptions.
**Accept/output:** a compact cooperative-structure report and tests enforcing
its domain. If nothing universal is proved, preserve that honest outcome;
the exact CRT identity and other independent results remain valid.

## Work package K — Rebuild the phase and cumulant calculations

### R7P-081 — Identify exactly which amplitudes the phase census used

**Depends:** 003–004,009. **Class:** S.
Recover the reported fixed amplitudes and the convention relating them to
h, theta, or complex Fourier coefficients. Distinguish angular coexistence
amplitudes from radial g_eq amplitudes. If only rounded amplitudes are
available, define that decimal fixture explicitly and create a separate
interval-parameter extension task rather than pretending it is an exact root.
**Accept/output:** an unambiguous phase fixture and provenance note. Missing
amplitude definitions block comparison with the reported 60-root count.

### R7P-082 — Derive cumulants two through four exactly

**Depends:** 009,012,081. **Class:** M.
Expand the uniform finite average of powers of the field and derive
`K4=kappa2/2+kappa3/6+kappa4/24` for a mean-zero field. Record every amplitude,
phase, sign, and normalization factor. Compare direct finite sums with the
resonance expansion.
**Accept/output:** an exact symbolic cumulant implementation and tests. The
fourth cumulant is not simply the fourth raw moment; subtract its disconnected
variance contribution correctly.

### R7P-083 — Enumerate cubic phase locks and their scope

**Depends:** 010,082. **Class:** M.
Derive the resonance conditions from frequency sums mod 12. For strictly
positive amplitudes and each sign of the alternating component, solve the
simultaneous cubic phase-maximization conditions modulo 2pi. Compute D12
orbits and stabilizers in the declared convention.
**Accept/output:** an exact count if justified, including the degeneration
when an amplitude vanishes. A fixed-positive-amplitude cubic maximum is not
yet a maximum of full log-mgf or of the unconstrained radial landscape.

### R7P-084 — Test the quartic compatibility claim

**Depends:** 082–083. **Class:** M.
List all quartic phase-sensitive terms with their actual signs after the
cumulant subtraction. Test whether the cubic locks simultaneously maximize
the relevant quartic contribution. Search deliberately for amplitude choices
that invalidate a universal compatibility claim.
**Accept/output:** an exact theorem with an amplitude domain, or a counterexample
and narrower valid statement. Do not infer compatibility from one coexistence
fixture or a list of resonance frequencies alone.

### R7P-085 — Construct full phase derivatives and remainder diagnostics

**Depends:** 012,081–084. **Class:** M.
Implement full log-mgf and K4 values, phase gradients, and phase Hessians in
the same coordinates. Reproduce a fixed-seed Sobol diagnostic and save every
sample used to estimate remainder norms. Check chain-rule terms from the
nonlinear phase parametrization.
**Accept/output:** complete diagnostic arrays and maximum locations. Sampled
maxima of a remainder are lower estimates of its supremum, not certified
uniform upper bounds.

### R7P-086 — Prove usable uniform C2 remainder bounds

**Depends:** 013,082,085. **Class:** L.
Derive a remainder estimate for values, gradients, and Hessians on the stated
phase torus and amplitude box. Use rigorous finite-sum/log bounds, analytic
series tails, or an adaptive interval cover. Keep dependency inflation visible
and subdivide only where it buys a proved margin.
**Accept/output:** explicit uniform bounds, or a bounded failed attempt with
the exact factor preventing the desired structural-stability conclusion.
Do not substitute the imported 65,536-sample extrema for these bounds.

### R7P-087 — Choose a complete real phase-equation representation

**Depends:** 082,085. **Class:** M.
Consider periodic angle boxes, sine/cosine polynomial equations with circle
constraints, or tangent-half-angle charts. If using algebraic equations,
exclude nonphysical complex solutions and track chart singularities. Define
how roots crossing the 0/2pi seam are represented once.
**Accept/output:** a root/cover specification on the full real three-torus.
Never treat one half-angle chart as a global parameterization without its
missing infinity faces.

### R7P-088 — Run a pilot and select the phase proof route

**Depends:** 081–087. **Class:** M.
Estimate root separation, weakest Hessian margins, and complement-gradient
behavior with a bounded pilot. Decide whether direct full-function isolation,
quartic-to-full continuation, or a combination is most realistic. State the
required C2 margin and the available certified bound.
**Accept/output:** a phase proof plan with explicit success/failure thresholds,
not an unconditional promise that all 60 roots will persist.

## Work package L — Phase census, continuation, and angular events

### R7P-089 — Reconstruct the complete numerical quartic candidate list

**Depends:** 087–088. **Class:** L.
Run reproducible root solves from structured and seeded starts. Save every
distinct candidate's phase triple, residual, Hessian spectrum, and symmetry
relations. Deduplicate on the torus using a documented tolerance.
**Accept/output:** the actual root catalog, whether it contains 60 points or
another number. A target count must not influence merging tolerances to force
agreement with the imported summary.

### R7P-090 — Interval-isolate quartic critical points

**Depends:** 013,089. **Class:** L.
Certify existence and local uniqueness for each representative candidate.
Enclose the phase Hessian and establish its index. Prove distinctness of
all translated/chart-related boxes and reconstruct symmetry copies only
under a justified action.
**Accept/output:** a partial or complete locally certified catalog. Report
unresolved nearly singular candidates individually instead of dropping them
from the count or assigning an index from midpoint eigenvalues.

### R7P-091 — Prove or explicitly fail quartic phase exhaustion

**Depends:** 087,090. **Class:** L.
Cover the complement of certified root neighborhoods and prove a nonzero
gradient component or a gradient-norm lower bound on each cell. Treat torus
seams and singular charts. Use exact algebraic root counts only if they
count the correct real constrained domain.
**Accept/output:** an exhaustive theorem or a map of unresolved complement
cells. Euler-characteristic agreement is a consistency check, never a
substitute for complement exclusion.

### R7P-092 — Continue isolated quartic roots to the full function

**Depends:** 086,090. **Class:** L.
Use a homotopy `K_t=K4+t(K_full-K4)` with validated local continuation,
or isolate full roots directly near the candidates. Prove that each tracked
branch remains nondegenerate if claiming preservation of index.
**Accept/output:** root-to-root correspondence with certified tubes or explicit
local boxes. A small numerical displacement and matching count do not exclude
birth/death of additional roots elsewhere during the homotopy.

### R7P-093 — Exclude additional full phase roots

**Depends:** 086,091–092. **Class:** L.
Either transfer a certified quartic complement-gradient gap using a strictly
smaller uniform gradient remainder, or perform a direct full-function cover.
Check the inequality on the exact domain outside the actual certified
neighborhoods, not outside arbitrary fixed-radius balls.
**Accept/output:** a full-function census theorem if all cells pass; otherwise
a partial correspondence plus an explicit remaining complement problem.

### R7P-094 — Test robustness to amplitude uncertainty

**Depends:** 081,090–093. **Class:** L.
Replace a fixed decimal amplitude fixture by an interval box justified by an
upstream angular/event certificate if available. Certify root persistence
and index preservation on that box, with overlap checks if subdivided.
**Accept/output:** a parametric phase theorem or an exact statement that the
result applies only to a supplied fixed fixture. Do not transfer between
angular and radial coexistence amplitudes without the actual map.

### R7P-095 — Re-derive the first pure-alternating angular instability

**Depends:** 010,012,013. **Class:** M.
At fixed Cartesian dual radius compute the **constrained spherical Hessian**
of log-mgf around the pure k6 direction. Include the Lagrange multiplier
term and compare every transverse sector. Recover or correct the proposed
k3 resonance equation and certify which sector crosses first.
**Accept/output:** a local angular-instability theorem and radius interval,
or a corrected competing-sector result. Ambient covariance eigenvalues alone
are not the constrained spherical Hessian.

### R7P-096 — Reconstruct angular fold/coexistence and compare truncations

**Depends:** 088,094–095 as available. **Class:** L.
Define the constrained angular stationary systems for the full and quartic
models; reconstruct the reported fold, equal-value event, and saddle barrier.
Certify the most accessible local event using bordered systems and correct
radius normalization.
**Accept/output:** a scoped angular diagram with proof level per landmark.
The angular 2-to-12 transition is not the radial localization transition in g,
and neither automatically proves global rank-seven minimality.

## Work package M — Full-seven-coordinate global variational questions

### R7P-097 — Define the global transition through an exact ratio problem

**Depends:** 011,015. **Class:** M.
For nonuniform p with `Q=(p-u)^T A7(p-u)>0`, study
`2D(p||u)/Q` as the candidate characterization of the first energetic
transition. Treat Q=0 and the limit p->u separately. Prove precisely when
an infimum is attained or approached and how it relates to local g_eq.
**Accept/output:** a correct variational definition and valid initial lower/
upper bounds, not an identification of g_global with the numerical local crossing.

### R7P-098 — Build sign-aware global lower bounds

**Depends:** 011,015,097. **Class:** M.
Develop interval box-simplex bounds or full-dual box bounds that respect A7's
signed real-space entries. Use PSD ordering `A7<=A_full` only for consequences
it actually licenses: for g>=0, `V7>=V_full` pointwise is valid, whereas
transfer of full-rank minimizers or positive-edge cap geometry is not.
**Accept/output:** tested lower-bound primitives and their tightness gaps on
known candidates. Reject bounds that are not rigorous on entire cells.

### R7P-099 — Certify an initial global gain bracket

**Depends:** 026,028,097–098 as available. **Class:** L.
Combine a proved uniform-global lower-gain region with an explicit normalized
competitor having strictly negative energy at an upper gain. Use rational
probability witnesses and outward objective evaluation where practical.
An equal-energy candidate alone is not a strict upper-side witness.
**Accept/output:** a rigorous bracket for a first energetic change, with
unresolved width stated. Do not declare the crossing unique or identify its
attaining orbit unless separately proved.

### R7P-100 — Combine the stationary atlas with certified exclusion tools

**Depends:** 038,091/093 where relevant,098. **Class:** L.
Use certified local minima and saddles to focus, but not replace, a global
complement cover. Bound energies on regions without isolated roots or exclude
stationarity there. Account for simplex boundary behavior through the analytic
entropy argument instead of ignoring small probabilities.
**Accept/output:** a partial or exhaustive stationary/minimum atlas with
explicit uncovered regions. A Morse polynomial with the correct alternating
sum is not proof that no additional cancelling pairs exist.

### R7P-101 — Test fixed-gain global orbit uniqueness independently

**Depends:** 027,038,098–100. **Class:** L.
Choose one precise supplied gain, preferably g=4 or a rational near g_eq,
and attempt to prove exactly one minimizing D12 orbit using disjoint certified
neighborhoods plus strict complement separation. Do not assume a peak cap
threshold imported from the full positive-edge model.
**Accept/output:** a scoped fixed-gain global theorem, a competing certified
orbit, or a remaining objective gap. Local convexity and a positive unresolved
global gap cannot be reported as exact attainment/uniqueness.

### R7P-102 — Audit all symmetry-reduction licenses used by M

**Depends:** 010,039,093,100–101. **Class:** M.
For every reduction, name the theorem that forces a minimizer into that
subspace or the complete fundamental-domain representation. Include zero
amplitudes, enhanced stabilizers, generic orbit-24 states, and phase-chart
singularities. Use the full objective, not a truncated phase model alone.
**Accept/output:** a reduction-license table. Any unsupported reduction
downgrades the affected global result to its actually searched subspace.

### R7P-103 — Study connections only after specifying a dynamical law

**Depends:** 029,038,100, plus a declared law. **Class:** M/L.
If studying saddle connections, choose and label an explicit gradient,
replicator, or other mobility law. Numerically integrate stable/unstable
manifolds and, only if feasible, use validated trajectories or isolating
blocks to certify selected connections.
**Accept/output:** trajectories and connection claims with proof levels.
The energy alone does not determine time, mobility, transition rates, basin
weights, or physical nucleation probabilities.

### R7P-104 — Publish the global-frontier status without overclaiming

**Depends:** 097–103. **Class:** S.
Distinguish local event, global energy bracket, fixed-gain minimizer orbit,
stationary exhaustion, Morse indices, and dynamical connections. State which
are proved, numerical, refuted, or unresolved, and which conclusions depend
on a particular supplied gain or amplitude class.
**Accept/output:** a global state-map update. If no global theorem is achieved,
the successful local, boundary, and counterexample results remain valid and
must not be hidden by an all-or-nothing verdict.

## Work package N — Secondary passive-memory and finite-N checks

These tasks are useful interpretation safeguards and finite mathematical
extensions. They are not a license to reopen generic gain-source speculation.

### R7P-105 — Independently verify the weighted Hodge decomposition

**Depends:** 009,013. **Class:** S.
Re-derive the rank-11 image and rank-55 cycle projection, orientation covariance,
and contraction identities. Add exact small-graph fixtures and weighted
relabeling tests. Distinguish an algebraic cycle space from crossings in a
particular planar drawing.
**Accept/output:** a reusable structural theorem/checker with declared graph
assumptions. Do not count an established decomposition as a newly sourced
physical interaction or an informational layer under the nadsoliton.

### R7P-106 — Investigate the reported 28 positive cycle eigenvalues

**Depends:** 105. **Class:** M.
Reproduce the numerical count with stated tolerances. Use D12 representation
blocks to predict multiplicities, then seek interval separation of the
distinct positive values if the block structure makes this affordable.
**Accept/output:** an exact multiplicity result, certified separated clusters,
or a numerical-only count. Matrix rank 55 by itself does not prove that there
are exactly 28 different positive eigenvalues.

### R7P-107 — Derive the exact finite-N empirical generator

**Depends:** 009,105. **Class:** S/M.
For independent walkers write transitions of occupation counts and derive
drift and state-dependent quadratic variation. Verify total-probability
conservation, the uniform equilibrium, and the exact multinomial covariance.
Use small N for exact finite-state tests before simulation.
**Accept/output:** an exact generator specification and regression showing
that A/(6N) is only the uniform-state noise covariance, not an all-state law.

### R7P-108 — Compare the jump process with the OU approximation

**Depends:** 107. **Class:** M.
At several explicitly supplied N values, compare equilibrium mode variances,
short-time covariance, and non-Gaussian diagnostics using reproducible
simulations or exact small-N transitions. State the scaling limit and any
error theorem separately from empirical convergence plots.
**Accept/output:** a finite-N/OU comparison with sample uncertainty. Noise
seeding all modes is not evidence that passive noise supplies active gain
or selects a unique localized phase.

### R7P-109 — Reconstruct the even/odd memory realization

**Depends:** 009,013. **Class:** M.
Block A_full, compute the exact transfer/self-energy form, and identify
hidden pole groups and visible residue ranks. Seek a proof of the reported
minimal visible dimension five using observability/controllability or residue
rank arguments, not only a numerical rank threshold.
**Accept/output:** a minimality theorem or a clearly numerical realization
diagnostic. A hidden dimension of six is not automatically a minimal visible
memory dimension of six.

### R7P-110 — Test the memory passivity statement over its actual domain

**Depends:** 109. **Class:** M.
Prove the relevant positive-real/Stieltjes or static Schur property under
the supplied symmetric passive assumptions. Distinguish frequency-dependent
memory, effective static stiffness, and growth-rate statements. Add a
negative-loading example only as an explicitly changed premise.
**Accept/output:** a scope-safe passivity result. A counterexample outside
the passive class does not evade ST293 without a new strictly sourced law.

### R7P-111 — Specify a conditional finite-copy Gibbs realization of V_g

**Depends:** 011,107. **Class:** M.
If useful, define an explicit finite-copy energy, normalization, and reference
measure whose large-N rate function has the supplied entropy-minus-quadratic
form. Derive all factors, including pair-count and diagonal/self terms.
State beta, J, and their product g as supplied parameters.
**Accept/output:** a conditional model with finite-N corrections or a bounded
derivation attempt. It does not establish the FIN provenance of beta, J,
the pump, the clock, or the chosen mediator rank.

### R7P-112 — Produce a passive/active interpretation firewall

**Depends:** 105–111. **Class:** S.
List exactly which passive identities are proved, which active terms are
added, and which choices are necessary for a dynamics or equilibrium law.
Check all new narrative sections for accidental source claims.
**Accept/output:** a short firewall note and terminology tests/review checklist.
If no new source object was introduced, the gain-source problem remains open;
this is not a universal impossibility theorem for every future completion.

## Work package O — Secondary quantum/spectral comparison

The quantum packages are a distinct conditional branch. Use them to clarify
operational claims and common spectral inputs, not to invent a causal bridge.

### R7P-113 — Freeze the quantum theorem scopes relevant to comparison

**Depends:** 001–003. **Class:** S.
Read the applicable proofs and current guards in the Hartree-equivalence,
separable-stationarity, and discord-robustness packages. Separate canonical V,
the mixed-flow family, and the larger pure-flow-equivalent family. Record
individual versus average marginal assumptions.
**Accept/output:** a scope matrix. The 22-product rank-132 construction belongs
to the separable-stationarity report; it is not a complete summary of the
later discord-robustness report.

### R7P-114 — Recompute the spectral comparison with interval semantics

**Depends:** 013,113. **Class:** S.
Express exact density spectral gaps through lambda6 and lambda6-lambda5.
Keep the previously certified rational lower enclosures Delta_L and delta_L
distinct. Recompute their product bound without replacing inequalities by
equalities with rounded transcendental quantities.
**Accept/output:** a table of exact expressions, certified lower intervals,
and numerical displays. Preserve the established discord theorem unless a
separate, fully verified improvement is actually proved.

### R7P-115 — Audit every proposed causal arrow between discord and localization

**Depends:** 073,113–114. **Class:** S.
List what the two research branches share (for example spectral numbers)
and what they do not share (state space, generator, preparation, operational
access, gain source). Attempt to write each causal claim as a theorem with
an explicit coupling map and premises.
**Accept/output:** a causal-claim obligation matrix. Missing maps remain
missing; numerical correlation or a common gap factor is not their substitute.

### R7P-116 — Verify local-indistinguishability versus joint distinguishability

**Depends:** 113. **Class:** M.
Replay the stationary channel mixture examples with fixed supplied programs.
Verify identical local channels, separability at balanced mixing, and joint
Swap discrimination using the existing exact formulas. Distinguish this
comparison from the separate passive-instrument equivalence example.
**Accept/output:** a reproducible operational comparison, marked as existing
theorem replay unless a genuinely new result is derived. No claim that local
FIN geometry identifies the joint preparation law.

### R7P-117 — Study controlled perturbations without losing assumptions

**Depends:** 013,113–116. **Class:** M.
If exploring small changes in loading or spectral parameters, certify density
positivity, simplicity/gap assumptions, and relevant block invertibility before
reusing a discord bound. Keep uniform-zero-loading and degeneracy limits as
adversarial tests.
**Accept/output:** a scoped robustness interval or a precise failure boundary.
An estimate proved for a fixed marginal cannot silently become an unrestricted
comparison over a changing marginal class.

### R7P-118 — Audit preparation resources and accessible observations

**Depends:** 113,116. **Class:** S.
For each discussed channel identify supplied program copies, classical shared
randomness, heralding, quantum routing, coherent control, and joint measurements.
Separate fixed-known-target preparation from a universal unknown-input channel
acting on states possibly entangled with a reference.
**Accept/output:** an operational resource table. Separable outputs do not by
themselves prove universal LOCC implementability, and symmetry does not itself
exclude coherent control of a known physical Hamiltonian.

### R7P-119 — Keep any legacy comparison genuinely separate

**Depends:** 113,117–118. **Class:** S/M.
If the comparison uses a legacy program, read its actual definition and
independently verify positivity and the required qualitative hypotheses.
Do not copy strict numerical constants or reinterpret signed legacy weights
as classical jump rates.
**Accept/output:** a separate scoped check or `OUT_OF_SCOPE` if no new comparison
is needed. No completion bridge or physical-role transfer follows from similar
downstream algebra.

### R7P-120 — Decide whether the spectral bridge adds any new theorem

**Depends:** 113–119. **Class:** S.
Classify the outcome as a new quantitative robustness result, a structural
analogy, an operational nonidentifiability clarification, or no new theorem.
Count existing replayed results separately from new research.
**Accept/output:** a concise comparison report preserving all physical
nonclosure gates. Do not manufacture a discovery merely because this secondary
lane was included in a long campaign.

## Work package P — Independent verification and next handoff

### R7P-121 — Audit claim-to-evidence traceability

**Depends:** all completed scientific tasks; run incrementally. **Class:** S.
For every accepted theorem or numerical finding, trace the exact source,
domain, proof file, checker, and raw output. Reject orphan conclusions with
only a table/figure citation. Verify that every task, including blocked and
superseded ones, has a terminal disposition before final completion.
**Accept/output:** a complete claim dependency graph and an orphan-claim report
with zero unaddressed accepted claims.

### R7P-122 — Perform independent and mutation-based checker review

**Depends:** 006 and all exported certificates. **Class:** M/L.
Recheck the most consequential root, inertia, and cover certificates with
independent routines where feasible. Mutate signs, denominators, endpoints,
spectral intervals, coverage leaves, and coordinate embeddings.
**Accept/output:** rejection logs and a clearly stated independence level.
Two wrappers around the same underlying unchecked calculation are not two
independent proofs. Any discovered checker bug triggers review of all dependent
results, not just the failing fixture.

### R7P-123 — Run the entire regression suite and compare the baseline

**Depends:** 121–122. **Class:** M/L.
Run baseline and new tests with read-only verification. Compare imported
artifact hashes with the starting manifest. Compare exact regenerated proof
layers exactly and numerical layers with predefined tolerances; do not choose
tolerances after seeing discrepancies.
**Accept/output:** verification.json with versions, commands, counts, durations,
source hashes, and explicit non-replayed historical material. Record tests
that are skipped and why; never include them in the passed count.

### R7P-124 — Assemble result tables and optional figures from saved evidence

**Depends:** 121,123. **Class:** S/M.
Generate concise tables for local events, stationary orbits, phase counts,
curvature domains, and global gain brackets. Give each row a proof-status
column. Plot only data saved in the package, and visually distinguish
certified regions from numerical interpolations.
**Accept/output:** reproducible result artifacts with source commands. A smooth
plot is not a certified continuation tube; large displayed precision must not
exceed the justified enclosure or numerical reliability.

### R7P-125 — Write the main English research report

**Depends:** 104,112,120,121–124. **Class:** S.
Lead with what is actually new relative to the intake audit. Present precise
theorems, counterexamples, numerical findings, failures, and remaining gaps
in separate sections. Include the distinction between four-amplitude,
phase-restricted, full-seven-coordinate, and 11D primal statements.
**Accept/output:** an English report with no unsupported global or physical
claims. Produce source-only Markdown by default; do not generate a PDF unless
the user later asks for one.

### R7P-126 — Prepare scoped repository guardrail updates

**Depends:** 125. **Class:** S.
Draft concise additions to AGENTS.md stating accepted new results, corrected
claims, exact nontransfer limitations, and the next genuinely live frontier.
If the later execution authorization includes repository integration, append
them carefully; otherwise deliver the proposed patch without applying it.
**Accept/output:** new guardrails or a reviewable patch. Never replace the
existing audit guardrails with an optimistic campaign summary.

### R7P-127 — Build a portable replay bundle and machine-readable manifest

**Depends:** 123–126. **Class:** M.
Package the new code, tests, proof objects, raw candidates, logs, task/claim
ledgers, report, handoff, and source-input references. Hash all included files
except the manifest itself; record how the manifest is generated. Test a clean
directory replay using only bundled files and declared dependencies.
**Accept/output:** a portable directory or archive plus its SHA-256 manifest.
Do not require hidden local files, live notebooks, private services, or an
unrecorded random seed to reproduce a promoted result.

### R7P-128 — Deliver the final handoff and stop the campaign cleanly

**Depends:** 121–127. **Class:** S.
Use the template below. Include every unresolved obligation, counterexample,
bounded no-go, stopped process, and exact next recommended atom. Verify that
no worker or child process is accidentally left running. Do not launch a new
research campaign merely because the present report lists future questions.
**Accept/output:** HANDOFF.md, final verification summary, package location,
and a clean stopping state ready for review by the next higher-reasoning agent.

---

## 8. Detailed algorithm recipes and common failure modes

These recipes reduce the amount of mathematical improvisation required during
execution. They do not waive the need to verify their hypotheses in each task.

### 8.1 Stable softmax and covariance calculation

For a field vector h, subtract `max(h)` before exponentiation. Compute means
and the centered second moment, rather than subtracting two nearly equal large
raw moments when possible. For a conditional parity class, normalize its six
weights directly using its own log-sum-exp. Do not obtain a tiny odd-class
probability by subtracting a rounded q_even from one and then divide by it.

This is important at large positive J6, where ordinary binary floating point
can produce q_even=1 long before the actual odd class disappears. Such a
rounding artifact can create false infinities, false zero covariance, or a
spurious resolvent violation. Exact limiting distributions should be handled
by their own boundary formulas, not by overlarge finite-field exponentials.

For intervals, use a proved common shift or interval-safe normalization.
Simply applying floating log-sum-exp to interval midpoints does not enclose
the true probabilities.

### 8.2 Correct root-certificate structure

For a square system F and rational box B centered at x0:

```text
1. Propose x0 and an approximate inverse Y numerically.
2. Freeze rational x0, rational Y, and rational B.
3. Enclose F(x0) and J(B) with verified arithmetic.
4. Form K(B) = x0 - Y F(x0) + (I - Y J(B))(B - x0).
5. Check the exact inclusion and the hypotheses of the chosen
   Krawczyk/interval-Newton theorem.
6. Record a positive strict inclusion margin in every coordinate.
7. Prove any claimed uniqueness scope explicitly; do not extend it
   beyond the theorem's certified domain.
```

Use a standard precise theorem formulation in the proof note; do not rely
on this schematic recipe as a substitute for the existence/uniqueness
hypotheses. If the inclusion fails, inspect which coordinate or dependency
causes the failure. Scaling, a better preconditioner, or a different box may
help. A smaller box can fail just as a larger one can; do not assume monotonic
improvement as the radius decreases.

Positive/negative definiteness should be certified across the root box, not
only at x0. For a parameter tube, separately certify overlaps or a common
unique branch argument. A chain of disjoint successful boxes is not a tube.

### 8.3 Covariance eigenvalues and multiple roots

When only the number of eigenvalues above sigma matters, use inertia of
`sigma I-M`; do not insist on enclosing each ordered eigenvalue independently.
Options include:

- Exact/interval LDL decomposition with verified nonzero pivots and pivot charts.
- Congruence to a simpler block form, with invertibility of the transformation.
- A verified sufficient sign-variation criterion for the shifted characteristic
  polynomial, including zero coefficients and equality cases.
- A rigorously chosen subspace on which the quadratic form has a definite sign.

A congruence does not preserve eigenvalues, only inertia. Ordered eigenvalues
are generally nonsmooth at multiplicities. At the double root use matrix
perturbation restricted to the degenerate eigenspace, or a polynomial/inertia
criterion. Numerical eigenvector derivatives near a crossing are unstable
and are not a safe foundation for a global proof.

For a refutation of lambda2(M)<=sigma, a certified two-dimensional subspace
with `v^T M v>sigma ||v||^2` for every nonzero v is enough. Save the subspace
basis and the positive-definiteness certificate for its compressed form.
Do not accept two unrelated positive Rayleigh quotients if their span contains
an uncontrolled cross term.

### 8.4 The parity Schur reduction and its exact limitations

In C4 coordinates, write

```text
M4 = W_par + b b^T,
W_par[6,:] = W_par[:,6] = 0,
eta = 1 - b6^2/sigma > 0,
Mtilde = W_par[345,345] + b345 b345^T/eta.
```

The accepted statement is equality of the negative inertia counts of
`sigma I4-M4` and `sigma I3-Mtilde`.

The scalar `S=1-b^T(sigma I-W_par)^(-1)b` additionally requires an invertible
denominator matrix. Its “one direction becomes two” sign interpretation
requires exactly one W_par eigenvalue already above sigma and none equal.
Record which condition is established on each cell. Near singular W_par,
prefer the valid Schur representation or another inertia chart.

The corrected one-dimensional variable is `r=sech(sqrt(lambda3/6)*s3)`.
Its endpoints are limits in s3 where appropriate. Any code that reintroduces
`r=exp(-sqrt(lambda3/6)*s3)` must be rejected by a regression test.

### 8.5 Adaptive cover discipline

Use a queue of exact cells, each with a chart, parent, and split history.
Classify a cell using only verified inequalities on the **whole** cell.
Suggested prioritization: uncertain physical constraints first, then near-zero
inertia indicators, then large interval overestimation. Split a coordinate
that materially reduces the dominant uncertainty; do not blindly split the
longest coordinate if its effect on the expression is negligible.

Every final leaf has exactly one status:

```text
EXCLUDED_FROM_DOMAIN
CERTIFIED_BY_INEQUALITY
CERTIFIED_BY_BOUNDARY_LEMMA
CERTIFIED_BY_LOCAL_EQUALITY_LEMMA
UNRESOLVED_RESOURCE_STOP
COUNTEREXAMPLE_CANDIDATE
```

Global PASS requires no unresolved or unchecked candidate leaves and a
verified complete covering of the original domain. Do not discard measure-zero
boundaries when the claimed theorem includes them. Keep overlap rules and
shared boundaries explicit. A true inequality on a larger outer domain proves
it on the model; a violation in that larger domain need not refute the model.

### 8.6 Parameter dependencies and exact zeros

The constants sigma, t_star, q_min, and q_crit are functions of the same
strict spectral values. Independent wide interval evaluation can destroy
exact cancellations. Prefer symbolic reduction of known identities before
interval substitution. For example, use

`t_star^2 = 1 - sigma/(lambda3/6)`

and factor the exact endpoint zero `(1-q)` in the accepted face polynomial.
Never set a small computed coefficient to zero because it is “obviously”
roundoff. Either derive its exact vanishing or enclose it as nonzero/uncertain.

Conversely, do not use a dependency identity beyond its domain: independently
perturbed coefficients may no longer satisfy it. A robustness theorem must
state whether perturbations preserve the defining spectral relations.

### 8.7 Correct phase and spherical derivatives

For a nonlinear phase parametrization theta(phi), the phase Hessian contains
both the pulled-back Cartesian Hessian and the gradient contracted with the
second derivative of theta(phi). Omitting the second term changes critical
point classification. At a sphere-constrained stationary point, the tangent
Hessian also includes the Lagrange multiplier correction.

Use three separate APIs:

```text
cartesian_dual_hessian(theta7, g)
fixed_amplitude_phase_hessian(phi, amplitudes, sign6)
fixed_radius_spherical_hessian(theta7, radius)
```

The phase torus has seams; amplitude-zero strata collapse phase coordinates;
the D12 action is finite, not a continuous gauge direction. Do not remove
a zero mode as “symmetry” unless there really is a continuous symmetry or
a redundant coordinate proven in the chosen formulation.

### 8.8 Local events versus global physics

Maintain this implication discipline:

```text
small residual
  -> numerical candidate only
validated root box
  -> local existence/uniqueness in that box
positive full tangent Hessian
  -> strict local minimum
equal energy with uniform
  -> local branch crossing
strict global complement exclusion
  -> global minimizer statement in the declared domain
```

None of these arrows supplies g, mediator rank, a clock, an apparatus, a
particular orbit member, or empirical evidence. A D12 orbit of minimizers
is not a selector. A variational landscape is not a unique dynamical law.
Energy equality does not prove a unique first transition.

## 9. Acceptance gates for major campaign claims

| Proposed claim | Minimum evidence required |
|---|---|
| Unique corrected 1D resolvent minimum | Complete derivative-root count/exclusion, isolated root, endpoint comparison, strict spectral intervals. |
| Local rank-seven coexistence | Validated augmented root and explicit local scope; no global claim. |
| Localized state stable in full model | Stationary root plus positive H7 or equivalent full primal tangent certificate. |
| Simple fold | Validated augmented root, one-dimensional kernel, both paid transversality coefficients. |
| Stationary-only index-one conjecture refuted | A certified full stationary root with two certified negative directions. |
| Boundary-Ising ceiling | Exact admissible domain, complete cover/analytic proof, local double-root treatment, all boundary cases. |
| Intraparity-W ceiling | Every actual Weyl, parity-weight, dominant-mass, and boundary prerequisite proved. |
| Global positive-orthant 4D ceiling | Domain/tail coverage, off-face control, equality treatment, and no unpaid lemmas. |
| Universal cooperative map | Exact required covariance signs and map/domain conditions, not only CRT terminology. |
| Exactly 60 phase critical points | Correct fixed-amplitude fixture, all root boxes, separation, Hessian indices, complete complement exclusion. |
| Same full/quartic phase topology | Certified continuation or direct roots plus no new roots on the complement; actual amplitude-domain scope. |
| Full rank-seven global orbit uniqueness | Full-domain reduction license or cover, local minima, strict outside separation, boundary exclusion, orbit/stabilizer proof. |
| Physical active-gain/source mechanism | A genuinely new explicitly typed source law and provenance evidence; none of the landscape tasks alone suffices. |

If the required evidence is absent, export a weaker result explicitly. The
word “certificate” must always identify the checked finite object and its
analytic connection to the claimed theorem.

## 10. Ready-to-use handoff template

The executing agent should fill the following structure in
`fin_rank7_followup/HANDOFF.md`. Do not replace it with a short success summary.

```markdown
# FIN rank-seven follow-up research handoff

## 0. Identity and execution state
- Campaign version and date:
- Repository baseline commit and dirty-worktree note:
- Environment and dependency versions:
- Executed task range and terminal task count:
- Package path and manifest hash:
- Running/stopped process inventory:

## 1. Executive result
- Most important genuinely new theorem:
- Most important certified counterexample:
- Most important reproduced numerical finding:
- Most important unresolved gap:
- What was already known before this campaign:

## 2. Source and coordinate contract
- Input hashes and baseline verification command:
- X7/C4 column conventions and parameter definitions:
- Exact spectrum provider and interval provenance:
- Imported fixtures versus derived exact states:

## 3. Complete task ledger
| Task ID | Execution status | Claim status | Main output | Remaining atom |
|---|---|---|---|---|

## 4. Accepted new theorems
For each theorem:
1. Exact statement, quantifiers, and domain.
2. Assumptions and dependencies.
3. Proof sketch and full proof file.
4. Exact certificate location and checker command.
5. Certified constants, not only display decimals.
6. Equality, degeneracy, and boundary cases.
7. Explicit nonconclusions.

## 5. Counterexamples and rejected conjectures
For each witness:
- Original claim.
- Full coordinates/parameters and admissibility.
- Exact or numerical verification level.
- Which claim is refuted and which nearby claims remain open.
- Regression test preventing reintroduction.

## 6. Numerical findings not yet certified
- Full candidate files, seeds, residuals, and tolerances.
- Solver and stopping settings.
- Search domain and coverage limitations.
- Candidate-to-theorem obligations.

## 7. Root and continuation atlas
- Root boxes, orbit/stabilizer information, and Hessian dimensions.
- Explicit certified joins and missing intervals.
- Local crossing/fold versus global transition distinctions.

## 8. Cover certificates and unresolved domains
- Domain/chart definitions and exact constraints.
- Leaf counts by status and complete tree files.
- Equality-neighborhood treatment.
- Unresolved boxes with coordinates and failure reasons.
- Tail/compactification coverage and remaining asymptotic strata.

## 9. Phase results
- Fixed-amplitude fixture and sign conventions.
- Quartic and full root coordinates and index counts.
- Uniform remainder bounds versus sampled estimates.
- Exhaustion status, chart seams, amplitude-zero cases.

## 10. Full-seven-coordinate/global results
- Known nonstationary index-two witness preserved.
- Status of stationary-only index claims.
- Local versus global minimizer results.
- Certified gain brackets and unresolved objective gaps.

## 11. Conditional interpretation and resource assumptions
- Passive/active law distinction.
- Quantum marginal/correlation/access distinctions.
- Which clocks, gain coefficients, and preparation resources are supplied.
- No selector/legacy-role/physical closure unless separately proved.

## 12. Failures, resource stops, and do-not-repeat list
- Exact failed method and mathematical/resource reason.
- Last valid checkpoint and resumability.
- New information required before retrying.

## 13. Verification and portability
- Clean replay command and expected output.
- Tests passed, failed, skipped, and not attempted.
- Mutation-test outcomes and independence level.
- Files required but not included, if any.

## 14. Proposed AGENTS.md update
- Only accepted new scope statements.
- Corrections and explicit nontransfer rules.

## 15. Ranked next frontier
- At most five concrete next tasks.
- Each must name one new witness, lemma, chart, source atom, or certificate.
- State what would count as completion or falsification.
```

## 11. Final instructions to the executing agent

Work persistently through ready tasks, but preserve epistemic honesty over a
preferred narrative. Do not ask for clarification when the present plan already
gives a safe, reasonable default. Do request direction when execution would
require new authority, missing private data, external paid compute, or a
material expansion beyond the user's authorized campaign.

Small exact results are valuable. A well-certified counterexample can be more
important than another large numerical scan. A precisely localized remaining
gap is more useful to the next reviewer than the phrase “almost closed.”

When prerequisites fail, continue independent tasks and record conditional
ones as blocked or superseded. Do not silently shrink a theorem's domain
while keeping its global title. Do not silently expand a theorem's domain
because nearby numerical results look similar.

Finish by handing over the actual research machinery and evidence. The next
reviewer should be able to determine, without reading the original chat,
exactly what was proved, what was computed, what was falsified, and what remains
to be done.
