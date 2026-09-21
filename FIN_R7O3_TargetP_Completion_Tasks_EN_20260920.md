# FIN R7O3: controlled completion of the physical-coupled Target-P cover

Date: 2026-09-20. Task namespace: **R7O3-001–R7O3-040**.

Prepared after reading the R7O2 handoff, its state ledgers, the physical-coupled
checker and the three generations of workers. **This is a task plan, not an
integration or acceptance of R7O2 research.** No R7O2 formula-level replay was
performed in preparing it, no worker was launched, and AGENTS.md was not changed.

The assistant executing this plan should return a new handoff for supervisory
review. It must not merge scientific claims into the authoritative repository
or edit AGENTS.md without a later explicit instruction.

## 1. Supervisory decision

The next priority is no longer to invent another broad enclosure architecture.
R7O2 supplies a concrete candidate improvement: a **shared-parameter interval
Taylor enclosure of a centered second moment**. Its correctness and replay
interface must be checked first; if those checks pass, finish the pending
finite computation using that method.

Recommended order:

1. Establish trustworthy provenance, a canonical parent registry, and a
   replayable certificate format.
2. Validate the new enclosure and the already-produced proof objects.
3. Finish the first-pass work on the currently unprocessed original parents.
4. Freeze the resulting repair queue and refine **only its unresolved leaves**.
5. Independently replay every accepted terminal formula and the complete
   geometry before attempting a global Target-P theorem.

Do not restart the obsolete R7N weight-box campaign, run the dynamic worker
blindly over all parents, chase the exploratory 2.26% figure, reopen the already
closed fixed-fixture phase census, or move to Target S/full X7 while this
finite completion question is unresolved.

This is not a promise that Target P is true or that every pending parent will
close. The permitted outcomes are a complete proof, a certified admissible
counterexample, or an explicit bounded residual with a genuinely new next atom.

## 2. Read these sources and preserve their status

Read in order:

1. Current [AGENTS.md](AGENTS.md), especially the accepted September 19 and
   September 20 rank-seven audit sections.
2. [Accepted R7N results](fin_r7n_review/ACCEPTED_RESULTS.md) and
   [R7N audit README](fin_r7n_review/README.md).
3. [R7O2 handoff](FIN_R7O2_CONTINUATION_HANDOFF_20260920/HANDOFF.md) completely.
4. [R7O2 current state](FIN_R7O2_CONTINUATION_HANDOFF_20260920/CURRENT_STATE.json).
5. [Physical-coupled checker](FIN_R7O2_CONTINUATION_HANDOFF_20260920/fin_rank7_targetp_physical_full_20260920/src/physical_centered_moment.py).
6. The depth-five, dynamic, chunk and four-parent repair source files named in
   the artifact map below before executing or adapting them.

The accepted R7N baseline remains:

- Exactly 60 critical points for each of the two declared fixed phase
  functions; no need to repeat that research here.
- Target P, `lambda2(M4)<=67/250`, globally open.
- 13,231 accepted compact-hull leaves and 5,432 unresolved leaves in the
  independently audited R7N partition.
- The unresolved R7N fraction is approximately 0.3631049472406795 of the
  compact coordinate hull.
- The sharp sigma target, full-X7 globality, source/selector and physical-role
  claims remain separate and open as previously stated.

R7O2's additional mathematical claims are **pending acceptance**, even where
their stored metadata says `fully_closed` or `INTERVAL_CERTIFIED`.

## 3. Starting state: what was checked for planning

Read-only structural inspection confirmed the following, without evaluating
the new mathematical leaf inequalities:

| Quantity | Structural result |
|---|---:|
| Frozen residual entries | 5,432 |
| Initial chunk-layer original indices | 540 |
| Dynamic-v2 stream original indices | 160 |
| Depth-five stream original indices | 2,760 |
| Unique recorded processed indices | 3,460 |
| Cross-layer index duplicates | 0 |
| Unprocessed original indices | 1,972 |
| Largest processed index | 3,461 |
| Unprocessed holes below that maximum | 3,454 and 3,458 |

The entire frozen `refined_failed` list and root hull are identical to those
in the available extracted R7N predecessor. This is stronger than matching
only a count. It is not a replay of the R7O2 certificates.

The handoff reports 3,406 closed parents and 54 depth-five repair parents.
The four earlier difficult parents 332, 338, 357 and 363 have a designated
replacement proof, not an additional contribution to the parent count.

### 3.1 Two different residual-volume measures

The handoff's reported fraction 0.12318558366584603 uses conservative
**whole-parent accounting**:

```text
not-yet-processed parent volume       ~ 0.1230813421248372
whole volume of 54 repair parents     ~ 0.00010424154100882935
conservative parent residual total   ~ 0.12318558366584603
```

Summing only the saved unresolved terminal boxes inside those 54 repair
parents instead gives approximately `0.000006331724266361647` of the hull.
This smaller number is a **structural sum from unreviewed proof records**, not
a newly accepted scientific bound. It illustrates why parent and terminal
accounting must be kept separate.

After validation, report both measures explicitly. Never call a whole-parent
upper estimate the exact remaining terminal residual. Never count an accepted
parent and its own accepted leaves twice. The large remaining work is primarily
the unprocessed set, not the current 54-parent repair set.

## 4. Concrete issues found in the implementation

These are implementation/replay obligations, not a declaration that the new
enclosure is mathematically wrong.

### 4.1 Missing center in saved certificates

`physical_centered_moment.certify` computes and returns a rational `center_c`,
but the worker `pack` routines do not preserve it in terminal records. A fixed
center is part of the centered-moment proof object. Do not rely on regenerating
the same numerical eigensystem and mean in a different NumPy environment.

Add a checker accepting **saved fixed B and c**. For older records, either
recover and justify the original center or choose a new rational center and
issue an explicitly labelled replacement certificate for the same cell.
Successful recertification with a different center is legitimate; calling it
byte-for-byte replay of the old certificate is not.

### 4.2 Floating diagnostic endpoints

The sign decisions inside the checker use rational intervals, but the workers
serialize d1/d2/d3 and the Gershgorin bound as floats. These are display data,
not portable exact proof endpoints. Serialize exact rational or explicitly
outward-rounded endpoints, including the method that paid positivity.

### 4.3 Volume equality is not coverage

The current workers check that summed terminal volume equals parent volume.
Equal volume alone does not rule out overlapping boxes compensated by gaps.
Dynamic paths record L/R while the split axis is not retained at every
internal node. A new verifier must reconstruct or record the actual split
tree and prove exact geometric coverage, not just compare scalar volumes.

### 4.4 Parent selection and silent malformed-file handling

The two stream workers have different exclusion sets. The older dynamic-v2
selection predates the depth-five layer and would duplicate work if launched
unmodified. Some scans use `except: pass`, and the depth-five worker treats
filenames in the v2 directory as done without validating every file's content.

Replace this with explicit validation and quarantine. A corrupt, truncated,
unknown-version or wrong-parent file is neither accepted proof nor silently
completed work.

### 4.5 Hard-coded predecessor location and absent local ZIP

Workers expect `/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920`.
That path is not portable. During this planning inspection the predecessor
ZIP was absent locally, but the extracted `FIN_R7N_HANDOFF_20260920/` directory
was available. Do not require the user to restore a deleted ZIP merely to
read already-available verified scientific inputs.

If the exact ZIP exists, compare the handoff's required SHA-256. Otherwise
use the extracted predecessor manifest and explicit dependency hashes, recording
that the archive-container hash itself was not rechecked. Do not regenerate a
ZIP and claim it has the identity of the missing original container.

The bundled `predecessor_minimal/` is not automatically a runnable replacement:
its root `results/` directory is absent, while inherited imports can read
`results/safe_union_v2_audited.json` at import time. Resolve actual dependencies
or refactor a pure mathematical dependency slice; do not create dummy JSON
files to silence the loader.

### 4.6 Runtime budgets are currently soft

The stream budget is tested after a parent returns from `imap_unordered`.
A difficult parent can exceed the nominal budget before any return occurs.
Implement explicit parent-level and batch-level supervision, cancellation
records and resumable ownership. Retain atomic writes, but do not confuse
atomic file replacement with mathematical validity or collision-free scheduling.

## 5. Exact mathematical contract for the new primitive

Let F_j be the four normalized C4 feature vectors and let M4 be their covariance
under the stated shared nonnegative fields. Let B be an exact rational 4-by-3
matrix of full column rank, and set `Z_j=B^T F_j`. Let c be any fixed rational
3-vector. Then

`Cov(Z) <= E[(Z-c)(Z-c)^T]`

because their difference is `(E[Z]-c)(E[Z]-c)^T`, which is PSD.

Therefore proving

`K = (67/250) B^T B - E[(Z-c)(Z-c)^T] > 0`

on an entire cell is sufficient to prove that M4 has at most one eigenvalue
above 67/250 there. B need not be exactly orthonormal, but its rank and the
Gram matrix must be handled exactly. Numerical eigenvectors are proposals,
not proof inputs until rationalized and checked.

The internal chart is

```text
A = sqrt(r), u = 1-s, v = 1-t, y = exp(-2 J6),
r = A^2, s = 1-u, t = 1-v, z = t^sqrt(3).
```

The seven aggregate weights, in the current representative-state order
`[0,4,6,2,3,5,1]`, are

```text
1,
2 s t^3,
A^2 t^4,
2 A^2 s t,
2 A t^2 y,
2 A s t^2 y / z,
2 A s t^2 z y.
```

All share the same variables and normalization. The proof engine must retain
that dependence through its interval jets. Independently intervalizing every
weight is an outer relaxation and was the earlier method's bottleneck.

For each matrix entry f in E, the intended enclosure is midpoint Taylor:

`f(X) subset f(m) + grad f(m) dot (X-m)
              + (1/2) sum_ij H_ij(X) (X_i-m_i)(X_j-m_j)`.

A conservative remainder using absolute Hessian bounds is acceptable. Validate
all first/second derivative rules, normalization inverses and domain hypotheses.
B and c must be constant during these derivatives. They can be proposed at
the original cell midpoint even if the midpoint of the square-root chart maps
to a different point; this is not an error because c is arbitrary and fixed.

The worker currently uses `bounded_rationals(9)`: outward rounding on a
10^-9 grid, not a floating tolerance. Prove that the exact spectral/coordinate
inputs remain enclosed. Coarse rounding can cause false rejection, but may
not be replaced with rounding to nearest or a positive-sign tolerance.

## 6. Execution protocol and output structure

Create a separate workspace such as `fin_rank7_targetp_completion/`.
Do not append new products inside the supplied handoff archive.

```text
fin_rank7_targetp_completion/
  README.md
  STATE_MAP.md
  TASKS.json
  INPUTS.json
  config.json
  src/
  tests/
  inherited/                 # references or verified copied dependencies
  parent_registry.json
  fast_pass/parents/
  repairs/parents/
  certificates/
  checkpoints/
  logs/
  results/
  verify.py                  # read-only mathematical/geometry replay
  HANDOFF.md
  AGENTS_PROPOSED_PATCH.md    # proposal only; do not apply
  MANIFEST.sha256
```

Each original parent keeps its immutable original index, original path and
cell hash. Each proof version has its own hash and an explicit `supersedes`
link. There must be one active proof version per parent, although historical
versions remain available.

Maintain separate states:

- `UNPROCESSED`;
- `PRODUCED_UNCHECKED`;
- `CERTIFIED_CLOSED`;
- `PARTIAL_REPAIR_REQUIRED`;
- `RESOURCE_STOP`;
- `INVALID_INPUT_OR_PROOF`;
- `COUNTEREXAMPLE_CERTIFIED`.

Do not map the supplied `fully_closed=true` directly to `CERTIFIED_CLOSED`
before verifying its formulas and geometry.

Suggested default limits:

- S: 5 minutes, 2 GiB;
- M: 20 minutes, 4 GiB;
- L: 60 minutes, 6 GiB, durable checkpoints;
- one worker until profiling, integrity and cancellation tests pass;
- later worker count must be set explicitly for the actual host and must not
  create duplicate parent ownership;
- each repair parent has both a wall-time limit and a certificate-call/depth
  limit. A timed-out parent stays visible in the residual registry.

These are safeguards for later authorized execution, not an instruction to
start computation now or to purchase external resources. Follow the actual
environment's permission and delegation rules.

---

## Package A — Provenance and canonical state

### R7O3-001 — Freeze inputs without merging them

**Depends:** none. **Class:** S.
Read the mandatory sources, record git state and versions, and hash the supplied
handoff, current state, source checker and frozen 5,432-parent input. Preserve
user deletions and unrelated dirty files. Record which R7O2 claims are pending
review versus already accepted R7N prerequisites.
**Deliver/accept:** INPUTS.json and a baseline note. No AGENTS.md edits, no
scientific merge and no restatement of imported PASS flags as new evidence.

### R7O3-002 — Verify predecessor identity with an honest fallback

**Depends:** 001. **Class:** S/M.
If the predecessor ZIP is available, verify the stated hash. Otherwise verify
the extracted predecessor manifest and required scientific dependency hashes.
Compare all residual rows, their order, original paths and the hull against
the accepted R7N checkpoint, not only the count 5,432.
**Deliver/accept:** an exact parent-list identity check and an explicit archive-
hash status. Do not manufacture a replacement container identity.

### R7O3-003 — Make dependency resolution portable and fail-closed

**Depends:** 002. **Class:** S/M.
Replace `/mnt/data/...` assumptions in the new working copy with one configured
predecessor root. Enumerate import-time data dependencies, including the inherited
safe-union registry. Use the existing full predecessor if available; treat the
minimal tree as incomplete until its needed dependencies are verified.
**Deliver/accept:** a pure import smoke test in a clean directory, with no search
or write triggered by import. Missing files cause explicit errors, not dummy inputs.

### R7O3-004 — Build the canonical original-index registry

**Depends:** 002–003. **Class:** S/M.
Parse all three indexed layers and the four-parent replacement. Verify file hashes,
content indices, cell identities, duplicate/conflicting records and supersession.
Recompute processed and missing index sets rather than using maximum filename.
Explicitly retain holes 3454 and 3458 unless new validated records fill them.
**Deliver/accept:** a registry for all 5,432 parents with one active record each.
Malformed or unindexed files are quarantined and reported, never silently skipped.

### R7O3-005 — Reconcile both volume ledgers

**Depends:** 004. **Class:** S/M.
Compute exact rational parent-level and terminal-level volumes separately.
Account for the four repaired parents as replacements, not extra accepted volume.
Compare the reconstructed totals with CURRENT_STATE.json and explain rounding.
**Deliver/accept:** exact conservation equations and a discrepancy report.
The reported 12.3186% parent residual is not silently equated to the smaller
terminal residual, and neither figure is treated as physical probability.

## Package B — Validate and harden the certificate primitive

### R7O3-006 — Write the precise centered-moment proof

**Depends:** 003. **Class:** S.
Derive the Loewner inequality, the role of fixed rational c, the full-column-rank
condition on B, and the min-max implication for lambda2(M4). State the threshold,
feature order and shared-field domain. Show why floating orthonormality is not
required and why the exact Gram matrix is.
**Deliver/accept:** a short standalone proof and rank-deficient/incorrect-Gram
negative controls. This remains conditional mathematics, not a source of gain.

### R7O3-007 — Audit the transformed interval jets

**Depends:** 006. **Class:** M.
Verify the map from a rational `(r,s,t,y)` cell to the outer `(sqrt(r),u,v,y)`
box, all seven aggregate weights, the z=t^sqrt(3) value/derivative enclosures,
and positivity of every divided denominator. Check multiplication/inverse
second derivatives and the midpoint-Taylor remainder.
**Deliver/accept:** analytic derivations plus independent diagnostic fixtures.
Finite-difference agreement alone is not a proof that a Hessian interval encloses
the full box.

### R7O3-008 — Separate proposal generation from proof checking

**Depends:** 006–007. **Class:** M.
Provide `propose(cell)` returning rational B and c, and
`certify_fixed(cell,B,c,policy)` performing no eigensolver or optimizer call.
Freeze the spectral input and arithmetic policy. Verify rank using exact minors;
check the matrix signs from interval entries.
**Deliver/accept:** a checker usable on another platform without regenerating
the proposal. Different valid proposals are allowed, but their certificates
must be recorded as different proof versions.

### R7O3-009 — Export complete leaf certificates

**Depends:** 008. **Class:** S/M.
Add missing `center_c`, exact basis, chart bounds, threshold, source hashes,
precision policy, exact/outward sign endpoints and proof method to each SAFE
leaf. Preserve float diagnostics only as display fields. Record axes and exact
split points at every internal node.
**Deliver/accept:** schema tests rejecting missing c, float-only proof endpoints,
wrong threshold, incomplete rank evidence and unknown backend versions.
Do not relabel incomplete legacy records as complete without recertification.

### R7O3-010 — Validate arithmetic and negative controls before mass execution

**Depends:** 007–009. **Class:** M.
Compare unrounded rational and outward-rounded backends on a fixed small corpus.
Test one-time versus repeated initialization, preserved input enclosures,
zero-denominator rejection and false-positive resistance. Mutate z coupling,
Gram normalization, c handling and one derivative coefficient.
**Deliver/accept:** a paid method-validation gate. If the enclosure proof or tests
fail, stop the large run and fix this atom first; do not generate thousands
of certificates from an unvalidated primitive.

## Package C — Verify and normalize inherited R7O2 proof objects

### R7O3-011 — Implement an independent geometric tree checker

**Depends:** 004,009. **Class:** M.
Prove that each child pair exactly covers its parent and that all terminals
are accounted for once. Use explicit split data where available; otherwise
reconstruct a consistent partition from exact boxes, refusing ambiguity.
**Deliver/accept:** rejection of missing leaves, overlapping leaves with equal
total volume, shifted boundaries, wrong-parent cells and inconsistent L/R paths.
Scalar volume equality is an additional check, not the coverage theorem.

### R7O3-012 — Validate the four historical replacement parents

**Depends:** 008–011. **Class:** M/L.
Check the complete repairs for original parents 332, 338, 357 and 363, including
all their leaves and coverage. Mark the raw earlier parent proofs superseded.
If c is absent, recertify the same leaves with an explicitly saved rational
center and record that this is a replacement certificate.
**Deliver/accept:** four independently replayable parent proofs or an explicit
repair queue. Never add their volume on top of the raw parent volume.

### R7O3-013 — Replay a stratified inherited corpus

**Depends:** 008–012. **Class:** M.
Test whole-cell closures, shallow and deep trees, both PD methods, near-zero
reported margins, boundary cells, and the 54 partial parents. Replay fixed
saved B wherever available and migrate centers transparently.
**Deliver/accept:** benchmark and discrepancy list. A sampled PASS licenses
confidence in the implementation workflow, not acceptance of every inherited
leaf. Any mathematical discrepancy blocks the affected bulk method.

### R7O3-014 — Replay all inherited SAFE leaves in bounded shards

**Depends:** 013. **Class:** L, repeat as bounded checkpoints.
Run the fixed-witness checker on every inherited active SAFE terminal and the
geometry checker on every processed parent. Store cumulative exact counts and
failed leaf identities. This may run in scheduled slots beside the new first
pass, but no inherited parent becomes certified merely from its old flag.
**Deliver/accept:** complete inherited replay or a precisely enumerated pending
replay set. Global completion requires this set to be empty.

### R7O3-015 — Freeze the verified/pending baseline registry

**Depends:** 012–013 and a completed checkpoint of 014. **Class:** S.
Distinguish valid produced records, independently certified closed parents,
partial parents, invalid records and genuinely unprocessed indices. Publish
which scientific checks remain pending if work is being staged.
**Deliver/accept:** versioned parent_registry.json, exact two-level volumes,
and a derived ready queue. Do not use an unqualified single “done” counter.
The full inherited replay may remain active while the first pass is scheduled;
it must be complete before global acceptance in Package G.

## Package D — Finish the fast first pass

### R7O3-016 — Derive the work queue by set difference

**Depends:** 004,010,015. **Class:** S.
Start from all original indices 0–5431 and subtract only valid produced records
under the active registry. Requeue invalid records explicitly. The initial
unprocessed set should contain 1,972 entries if the structural starting state
has not changed, including 3454 and 3458.
**Deliver/accept:** an immutable queued-index list with cell hashes. A contiguous
range starting at the highest processed number is forbidden.

### R7O3-017 — Implement safe worker ownership and cancellation

**Depends:** 009,016. **Class:** M.
Use per-parent ownership, atomic artifact publication and explicit state
transitions. Add parent and batch timeouts rather than relying only on the
post-result budget check. Preserve or requeue interrupted parents, and reject
truncated JSON instead of silently ignoring it.
**Deliver/accept:** interruption/resume and duplicate-worker tests. Indexes must
be reconstructible from validated immutable parent artifacts.

### R7O3-018 — Run a small first-pass pilot under the fixed policy

**Depends:** 010,017. **Class:** M.
Use the whole-cell test followed by `t -> r -> s -> t -> r`, depth at most five.
Run one worker first; record certificate-call cost and peak resources. The split
can be proposed numerically, but its exact rational endpoint is the geometry.
**Deliver/accept:** a pilot with independent formula/partition replay. Increase
worker count only after validating resource use and collision-free ownership.

### R7O3-019 — Complete every remaining first-pass parent

**Depends:** 018. **Class:** L, repeat with durable checkpoints.
Process the entire verified missing set in bounded chunks. Do not revisit closed
parents or discard slow ones. Each original parent finishes either with a complete
safe tree, a complete safe-plus-unresolved partition, or an explicit resource-stop
record. Track actual IDs, not just file counts.
**Deliver/accept:** all 5,432 originals have a produced, well-formed active state.
This is first-pass completion, not yet global Target-P proof.

### R7O3-020 — Freeze the actual repair queue

**Depends:** 019 and geometry checks. **Class:** S.
Combine the initial 54 partial parents with newly discovered failures and any
invalidated inherited proofs. Save each unresolved terminal, its original
parent ID, exact box, prefix proof and source hash. Do not assume the final
queue still has 54 entries.
**Deliver/accept:** an immutable repair queue and separate whole-parent and
terminal residual volumes. All already-safe prefixes remain credited only once.

## Package E — Repair the frozen unresolved terminals only

### R7O3-021 — Build a repair-only worker

**Depends:** 020. **Class:** M.
Accept an explicit parent/terminal queue rather than scanning all original IDs.
Import each verified safe prefix and replace only its unresolved leaves. Preserve
the original parent identity, prior proof hash and replacement relation.
**Deliver/accept:** a worker that cannot accidentally reschedule unrelated closed
parents. Its output is a new active proof version, not an append-only mixture
of overlapping old and new leaves.

### R7O3-022 — Execute a bounded dynamic local fallback

**Depends:** 010,021. **Class:** M/L per bounded batch.
Initially consider split axes r,s,t. Use immediate child successes and diagnostics
only to propose an axis; every accepted child needs the same fixed-witness proof.
Start with at most seven additional levels and an explicit call/time cap per
parent. Save exact split points and both children.
**Deliver/accept:** closed repairs or complete residual subtrees. A favorable
float d3 score is navigation, not evidence for accepting an unchecked child.

### R7O3-023 — Independently validate every repaired parent

**Depends:** 021–022. **Class:** M/L.
Replay the formulas with saved B,c and verify the replacement geometry, including
the inherited safe prefix. Check that each old unresolved leaf is wholly covered
by its replacement subtree and that no sibling was lost.
**Deliver/accept:** a per-parent VERIFIED result with counts and hashes.
Update the registry atomically only after this verification succeeds; leave
the previous proof version available for diagnosis.

### R7O3-024 — Separate arithmetic failure from geometric failure

**Depends:** surviving leaves from 022–023. **Class:** M.
On a small frozen difficult corpus compare outward precisions, e.g. 9, 12 and
18 decimal places, with the same saved witnesses. Separately compare a new
proposed witness. Classify coarse-rounding inflation, wide-cell Taylor remainder,
poor basis/center and possible genuine violation.
**Deliver/accept:** a measured diagnosis. Higher precision is useful only if it
changes a paid bound; it is not automatically stronger scientific evidence.
Never alter a threshold or introduce a sign tolerance to manufacture PASS.

### R7O3-025 — Freeze the post-repair state and decide the branch

**Depends:** 023–024. **Class:** S.
If every original parent is verified closed, proceed to the final global audit.
If not, stop broad computation and publish the surviving microscopic terminals
with exact domains and failure reasons. Do not silently increase depth for all
5,432 parents or resurrect a noncanonical exploratory checkpoint.
**Deliver/accept:** either a candidate-complete global proof set or one immutable
residual queue. The next package applies only to actual surviving terminals.

## Package F — Analyze survivors, not the whole hull again

### R7O3-026 — Recompute the actual M4 at each surviving locator

**Depends:** nonempty queue from 025. **Class:** M.
Use the original shared-field model and fixed threshold, with stable normalization.
Compare center values and selected rational nearby points with interval enclosures.
Preserve the r,s,t,y domain and all physical coupling relations.
**Deliver/accept:** a distinction between negative-gap method failures and
candidate positive-gap points. Midpoint safety does not prove the whole cell;
a failed certificate does not prove a violating point.

### R7O3-027 — Test one targeted enclosure improvement

**Depends:** 024,026. **Class:** M.
Choose only the dominant identified issue: a tighter z derivative enclosure,
a smaller Taylor chart, a better fixed center, a different rational subspace,
or a validated subdivision adapted to the remainder. Keep a before/after
difficult-cell benchmark.
**Deliver/accept:** one genuinely new improvement or a bounded no-improvement
result. Do not launch unrelated polynomial elimination or another unconstrained
global optimizer as a substitute for this diagnosis.

### R7O3-028 — Permit a new split direction only with evidence

**Depends:** 027 if r/s/t fallback stalls. **Class:** M.
If the common normalization or y direction dominates the remainder, run a
bounded comparison including a y split or a justified alternative chart.
Record why this is a new blocker-cut rather than more of the same subdivision.
**Deliver/accept:** a local validated repair or a precise unsuccessful result.
The initial r/s/t preference is a cost policy, not a theorem forbidding y.
No uncontrolled change of variables may relax the physical model silently.

### R7O3-029 — Certify any genuine counterexample

**Depends:** a stable candidate from 026–028. **Class:** M/L.
Prove admissibility and certify two supercritical covariance directions, for
example using an exact rank-two basis U with
`U^T(M4-(67/250)I)U > 0` on a declared parameter enclosure. Recompute
independently from the original model.
**Deliver/accept:** a certified counterexample or a rejected numerical artifact.
Only an admissible violation decides Target P. A violation in a relaxed
independent-weight or independent-parity model does not.

### R7O3-030 — Enforce the research stop rule

**Depends:** 025–029. **Class:** S.
If survivors remain after the bounded local alternatives, return their exact
subdomains and one missing analytic/algorithmic atom. If no survivors remain,
mark unnecessary survivor tasks OUT_OF_SCOPE with that reason and proceed.
**Deliver/accept:** no invisible work queue and no unlimited campaign. Resource
exhaustion is not mathematical nonexistence; successful local repair is not
global closure until the complete proof set is checked.

## Package G — Independent global Target-P audit

### R7O3-031 — Check all 5,432 original parents exactly once

**Depends:** candidate completeness from 025/030, or partial audit if incomplete.
**Class:** L in shards.
Rebuild the active proof set from validated artifacts, not from summary counters.
Verify original index, path and cell hash; exact child coverage; active-version
uniqueness; and zero unresolved leaves for a requested global PASS.
**Deliver/accept:** a complete parent-coverage theorem or a precise failure.
Checking 5,432 files is insufficient if indices repeat or a parent is omitted.

### R7O3-032 — Replay every terminal mathematical certificate

**Depends:** 008–010,031. **Class:** L in resumable shards.
Recompute the coupled weights, midpoint/Taylor enclosure, fixed B,c moment
matrix, exact rank and positive-definiteness test for every active SAFE leaf.
Use no stored Boolean or float eigenvalue as the proof decision.
**Deliver/accept:** a full all-leaf replay with exact failed identities if any.
Sampling is useful for debugging only. A global theorem requires the complete
scientific replay and its documented analytic premises.

### R7O3-033 — Join the new proof to the accepted R7N cover and tails

**Depends:** 031–032. **Class:** M.
Show that the 5,432 certified parents are exactly the previously unresolved
R7N part of the compact hull. Include the 13,231 already accepted compact
leaves and the accepted unbounded tails. Check common boundaries explicitly.
**Deliver/accept:** a global domain-union argument. Proving the old residual
alone is not the whole nonnegative-field theorem without this join.
Do not add W_par safety as if it were automatically M4 safety.

### R7O3-034 — Run hostile proof and scheduling mutations

**Depends:** 031–033. **Class:** M.
Delete a leaf; overlap two leaves while preserving total volume; omit one
low-index hole; corrupt a center, Gram matrix, rank minor, split or threshold;
duplicate a parent across layers; truncate a file; and alter a dependency hash.
**Deliver/accept:** all invalid constructions are rejected with precise reasons.
The verifier must also reject a requested global PASS when even one residual
or unchecked terminal remains.

### R7O3-035 — State the theorem or partial result with its gain consequence

**Depends:** 031–034. **Class:** S.
Only after all gates pass may the report state `lambda2(M4)<=67/250` globally
in the declared nonnegative shared-field family. Then derive the C4 Hessian
consequence: at most one strictly negative eigenvalue for supplied
`0<g<=250/67`. Preserve zero-eigenvalue caveats unless separately excluded.
**Deliver/accept:** a theorem candidate for supervisory review, or a clearly
partial result. No transfer to sigma, X7, physical gain provenance or selector
closure is allowed.

## Package H — Handoff, portability and no automatic integration

### R7O3-036 — Produce exact final accounting and evidence tables

**Depends:** 025,030,035. **Class:** S.
Report processed, produced-unchecked, certified-closed, partial, invalid and
unprocessed original indices. Give exact parent and terminal volume fractions
and proof-leaf counts, separately from optimizer statistics and wall time.
**Deliver/accept:** a table derivable entirely from active verified artifacts.
Never mix the exploratory R7O 2.26% estimate into the canonical R7O2/R7O3 state.

### R7O3-037 — Test clean-directory mathematical replay

**Depends:** 032–036. **Class:** M/L.
Copy only declared inputs and code into a clean directory, resolve all paths
without the producer's `/mnt/data` layout, and replay geometry and mathematics.
Separate full checks from smoke tests and stored-record checks. Exclude cache
from scientific manifests, or explicitly classify any inherited cache entries.
**Deliver/accept:** reproducible commands, dependencies and actual test counts.
Do not claim a clean replay from a manifest check plus sampled inequalities.

### R7O3-038 — Write the scientific report and proposed guardrails

**Depends:** 035–037. **Class:** S.
Lead with what changed from accepted R7N and from the unreviewed R7O2 handoff.
Describe the new enclosure, completion/refutation/partial status, analytic
dependencies, quantitative bounds and nonconclusions.
**Deliver/accept:** REPORT.md, accepted-claim candidates and AGENTS_PROPOSED_PATCH.md.
Do not apply the patch or alter the authoritative accepted-result register:
the user requested a handoff before integration.

### R7O3-039 — Package every required proof object

**Depends:** 037–038. **Class:** M.
Include the exact frozen input list, active parent registry, complete tree
objects, B and c witnesses, interval bounds, source versions, checkers,
mutations, raw failures, logs and dependency manifest. Preserve superseded
proofs with explicit lineage, or package their hashes and archival locations.
**Deliver/accept:** a portable package with no missing centers, hidden global
variables, missing predecessor data or presumed original ZIP identity.

### R7O3-040 — Return the final handoff and stop

**Depends:** 036–039. **Class:** S.
Use the template below. State whether Target P is a fully verified theorem
candidate, a certified refutation, or still partial. Confirm that no worker is
left running unintentionally. Provide at most three ranked next atoms.
**Deliver/accept:** HANDOFF.md and the package path for supervisory review.
Do not launch Target S, full X7, another phase campaign or repository integration
automatically after finishing this computation.

---

## 7. Required certificate and registry fields

### 7.1 Complete SAFE-leaf record

```json
{
  "original_index": 0,
  "original_parent_sha256": "...",
  "proof_version": "R7O3-fixed-witness-v1",
  "cell": [["lo", "hi"], ["lo", "hi"], ["lo", "hi"], ["lo", "hi"]],
  "chart": "sqrt(r),1-s,1-t,y",
  "threshold": "67/250",
  "basis_num": [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
  "basis_den": 100000,
  "center_c": ["rational", "rational", "rational"],
  "rank_minor_rows": [0, 1, 2],
  "rank_minor_exact": "nonzero rational",
  "spectral_source_sha256": "...",
  "checker_source_sha256": "...",
  "arithmetic_policy": {"kind": "outward rational intervals", "grid_digits": 9},
  "internal_chart_bounds": [],
  "positive_denominator_bounds": [],
  "moment_entry_enclosures": [],
  "gram_matrix_exact": [],
  "pd_method": "Sylvester or Gershgorin",
  "pd_bounds_exact": [],
  "status": "PRODUCED_UNCHECKED"
}
```

The zero matrix above is a schema placeholder, **not** a valid example
certificate. A real record must pass a nonzero exact rank-minor check.
The independent checker, not the producer's assignment, upgrades its status.

Every internal tree node records:

```text
node ID, exact parent cell, split axis, exact split coordinate,
left/right child IDs and boxes, parent index, and source/proof version.
```

For a reconstructed legacy tree, identify it as a geometric reconstruction.
It need not reproduce the original search order, but it must prove the same
complete, non-overlapping-up-to-boundaries parent cover.

### 7.2 Active parent record

```text
original index / original path / original cell hash
reported source state
independently verified state
active proof hash
superseded proof hashes
geometry-check result
formula-replay result and progress
SAFE / unresolved terminal counts
whole-parent volume
SAFE-terminal volume
unresolved-terminal volume
owner / lock / checkpoint / resource stop
```

Record a source-file name and its internal parent ID independently and check
their agreement. Do not infer validity from filename existence.

## 8. Worker-selection rules that prevent duplicate work

The first-pass set is computed from the canonical registry:

```text
ALL = {0,...,5431}
VALID_PRODUCED = IDs with a well-formed active parent record
FIRST_PASS_QUEUE = ALL - VALID_PRODUCED
INVALID_REPLAY_QUEUE = explicit IDs whose existing records were rejected
```

The repair queue is not another call to the first-pass selector:

```text
REPAIR_QUEUE = unresolved terminal IDs in valid partial-parent partitions
```

The independent formula-replay queue is also separate:

```text
REPLAY_QUEUE = active SAFE terminal certificates not independently verified
```

Do not equate these queues. A parent can have a produced complete tree while
its formula replay is pending. It is not yet eligible for the final global
theorem. A parent with a few unresolved terminals may already contain useful
safe subdomains, which must be preserved during repair.

Use one writer per parent proof version. A cancellation must release or expire
ownership without marking the parent safe. On restart, load only valid complete
JSON records with expected hashes/schema; list every other file as an error
or abandoned temporary artifact.

## 9. Independent verification recipe

For each active parent:

```text
1. Match original index, original path and exact box to the frozen R7N list.
2. Resolve the active proof version and all replacement links.
3. Verify the complete geometric partition, not only total volume.
4. For every SAFE leaf:
   a. validate its physical compact domain and threshold;
   b. reconstruct the coupled seven-state model from fixed input data;
   c. verify the exact rank of saved B;
   d. hold B and c fixed and enclose the normalized second moment;
   e. verify the PD inequality with exact/outward endpoints;
   f. record the actual positive margin and source hashes.
5. Refuse global closure if any unresolved or unchecked leaf remains.
6. Recompute exact parent and terminal volumes as consistency diagnostics.
```

A valid geometric partition plus valid SAFE inequalities proves a parent.
A sum of positive margins without coverage does not. A coverage tree whose
leaf inequalities are only stored assertions does not.

## 10. Logical gates for the final mathematical statement

The chain required for global Target P is:

```text
accepted R7N safe compact leaves
  + all 5,432 exact original residual parents, fully verified
  = complete compact hull

complete compact hull
  + accepted unbounded-tail coverage
  = all declared shared nonnegative fields

all declared shared nonnegative fields
  -> lambda2(M4)<=67/250
  -> index_negative(I4/g-M4)<=1 for 0<g<=250/67.
```

The last line concerns the four-amplitude Cartesian Hessian only. At the gain
endpoint, a non-strict ceiling does not itself exclude zero eigenvalues.
If a stronger uniform physical gap is desired, it requires a separate direct
matrix argument and quantitative normalization; do not import the withdrawn
Schur-gap shortcut.

Target P is weaker than Target S. A proof at 67/250 is not a proof at sigma.
A certified violation above 67/250 in the same admissible model would also
violate the smaller sigma ceiling, but a violation in a relaxed model would
not have that implication.

There is no inference here to a sourced gain, physical clock, selector,
QW-2191 discharge, apparatus, laboratory evidence, legacy-role transfer,
SM/GR, L_total or a theory of everything.

## 11. Required next-handoff template

```markdown
# FIN R7O3 Target-P completion handoff

## 0. Identity and provenance
- Package version, date, source hashes and environment.
- ZIP identity checked / unavailable; extracted-tree identity checked separately.
- Exact original 5,432-parent input and its hash.
- No source archive modified; no unauthorized merge performed.

## 1. Executive result
- Target P: theorem candidate / certified refutation / partial.
- What changed after R7O2, and what remains only imported evidence.
- Strongest accepted new mathematical statement and its exact scope.

## 2. Primitive validation
- Centered-moment and min-max proof.
- Transformed jets, denominator bounds and irrational-power coupling.
- Backend policy and comparison tests.
- Separation of proposal generation from fixed-witness checking.

## 3. Canonical parent registry
| State | Parent count | Whole-parent fraction | Terminal residual fraction |
|---|---|---|---|
- All original IDs accounted for, including noncontiguous holes.
- Replacement lineage for 332,338,357,363 and later repairs.
- No duplicate active proof or silently ignored malformed file.

## 4. First pass
- Starting queue, exact processed set, interruption/resume behavior.
- Fixed split policy and actual runtime/call statistics.
- Frozen final repair queue, not assumed equal to the original 54.

## 5. Repairs
- Original unresolved leaf -> replacement subtree map.
- Saved B,c witnesses, exact bounds and split geometry.
- Surviving leaves, if any, with precise failure classification.

## 6. Independent replay
- Full geometry coverage, all-leaf formula checks and mutation outcomes.
- Source/dependency hashes and actual verification commands.
- Tests passed, failed, skipped, historical-only or not executed.

## 7. Global join and theorem statement
- Join to accepted R7N compact safe leaves and tails.
- Exact threshold, domains, equality/zero-mode caveats and gain consequence.
- No transfer to sigma or X7 unless separately proved.

## 8. Counterexamples / bounded failures
- Actual admissibility and interval witness for any claimed violation.
- Method failure versus mathematical failure versus resource stop.
- No unsupported replacement by the exploratory 2.26% state.

## 9. Full R7O3-001--040 task ledger
| ID | Execution state | Scientific state | Artifacts | Remaining atom |
|---|---|---|---|---|

## 10. Portability
- Complete file manifest and clean replay.
- No hidden /mnt/data paths, missing centers or unstated input files.
- Cache/temporary-file policy and nonreplayed dependencies.

## 11. Integration proposal only
- Proposed AGENTS.md patch and accepted-result index.
- Await supervisory review; do not apply automatically.

## 12. Next atoms
- At most three specific tasks justified by the actual residual or theorem.
- Confirm no background worker remains active unintentionally.
```

## 12. Final direction

Keep the campaign narrow. The finite original residual is frozen, a promising
coupled enclosure has been supplied, and the main challenge is now controlled
completion plus independently checkable evidence. Improve bookkeeping and proof
serialization before producing more bulk results.

If all parents close, return a complete global Target-P proof candidate for
review. If some survive, return only those exact cells and the smallest new
mathematical obligation. Do not change the target, hide unfinished work,
or start another branch automatically. **Do not merge this continuation yet.**
