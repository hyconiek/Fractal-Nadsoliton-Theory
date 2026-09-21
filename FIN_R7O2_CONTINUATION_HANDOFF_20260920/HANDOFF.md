# FIN R7O2 continuation handoff — physical-coupled Target P

## 0. Predecessor identity

This handoff is the **direct successor** to the final R7N bundle:

- predecessor archive: `FIN_R7N_HANDOFF_20260920.zip`
- required predecessor SHA-256: `ba56cc0ecf69ac970e1d6fe719e0770753437213a4ef345949ffcf1411941d9c`
- predecessor scientific state: exact fixed-fixture phase census (quartic and full log-mgf both exactly 60 critical points); Target P globally open with 5,432 compact-hull residual cells occupying `0.3631049472406795` of compact-hull volume.

The predecessor remains immutable. Work after that bundle is kept in separate continuation directories.

## 1. Executive scientific delta after R7N

The main new development is a **physical-coupled centered-moment enclosure** for Target P at the exact practical threshold

`tau0 = 67/250 = 0.268`.

The old Target-P interval machinery lost too much dependence by enclosing state weights separately. The new checker instead keeps the common physical dependence on `(r,s,t,y)` through a second-order interval jet, including the coupled terms

- `t^3`,
- `t^4`,
- `t^(2+sqrt(3))`,
- `t^(2-sqrt(3))`,
- and the common normalization denominator.

For a rational 3-dimensional subspace basis `B` and a fixed center `c`, it uses the rigorous one-sided inequality

`Cov(Z) <= E[(Z-c)(Z-c)^T]`

and proves positive definiteness of

`tau0 * (B^T B) - E[(Z-c)(Z-c)^T]`

by interval Sylvester minors or, secondarily, Gershgorin. This avoids the unstable interval subtraction in direct covariance enclosure while preserving the physically coupled weights.

This was the missing primitive identified in the R7N handoff.

## 2. Exploratory R7O campaign — discovery only

Directory: `fin_rank7_physical_coupled_20260920/`.

This earlier continuation established that the new physical-coupled enclosure is much stronger than the R7N weight-box relaxations.

A rigorously consolidated experiment on 300 original R7N residual parents closed all 300:

- 140 as whole cells,
- 160 through local exact-cover refinement,
- incoming residual: `36.31049472406795%` of compact hull,
- after those 300 parents: `25.81476446699962%` remained,
- exact reduction: `10.49573025706833` percentage points of compact hull.

Evidence:

- `fin_rank7_physical_coupled_20260920/results/physical_coupled_complete_local_tree.json`
- local repair example: `results/repair_4694_depth5_r.json`

A later exploratory adaptive run suggested residuals near ~2.26% were reachable. **Do not promote that number to the canonical proof state.** Some intermediate old proof chunks were not retained completely, so this exploratory path is evidence for method selection, not the preferred global proof object.

That is why R7O2 below was restarted from the original frozen set of all 5,432 R7N residual cells.

## 3. Canonical fresh R7O2 campaign

Directory: `fin_rank7_targetp_physical_full_20260920/`.

Input is frozen in:

`inputs/target_p_residual_5432.json`

It contains exactly the 5,432 residual parents exported by the R7N Target-P cover. R7O2 is independent of the partially lost exploratory adaptive proof objects.

### Base proof policy

For each original residual parent:

1. run the physical-coupled certificate on the whole cell;
2. if needed use the fixed local split sequence
   `t -> r -> s -> t -> r`;
3. every terminal SAFE leaf must independently satisfy the physical-coupled interval certificate;
4. children must exactly cover the parent by rational endpoints;
5. checker failure is `UNRESOLVED`, never a counterexample.

A small number of difficult parents can receive a separately stored deeper local repair tree. The first such set was parents `332, 338, 357, 363`; all four now have a complete exact repair proof.

Evidence:

- `results/repair_four_parents_complete.json`
- `checkpoints/full_cover_chunks/index.json`

## 4. Exact current R7O2 state

This is the state that the next campaign should continue from.

Original R7N residual:

- parent count: **5,432**
- volume fraction of compact hull: **0.3631049472406795**

Fresh R7O2 state:

- unique original parents processed: **3,460 / 5,432**
- rigorously fully closed parents: **3,406**
- fast-depth-5 parents on repair queue: **54**
- not yet processed by the fresh cover: **1,972**

Exact volume accounting relative to the full compact hull:

- newly certified by fresh R7O2: **0.23991936357483346**
- current repair queue: **0.00010424154100882935**
- not yet processed: **0.1230813421248372**
- total fresh unresolved: **0.12318558366584603**

Thus the fresh portable continuation has already certified about **66.0744% of the incoming R7N residual volume**. Equivalently, the unresolved part of the full compact hull has dropped rigorously from `36.3105%` to **12.3186%** in this fresh run.

This is still a **partial-domain certificate**. Target P is **not yet globally proved**.

Machine-readable summary: `CURRENT_STATE.json`.

## 5. Current proof-object layout

### A. Initial full-cover chunks

`checkpoints/full_cover_chunks/index.json`

- first 540 parents accounted for,
- four raw difficult parents `332,338,357,363` are replaced by the exact repair proof,
- effective closed count = 540.

### B. Dynamic-v2 stream

`checkpoints/stream_v2/index.json`

- 160 additional parents,
- all 160 closed,
- cumulative base after this layer = 700 closed parents.

### C. Fast depth-5 stream

`checkpoints/depth5_stream/index.json`

- 2,760 additional parents processed,
- 2,706 closed,
- 54 unresolved parent indices frozen for later repair.

Current 54-parent repair queue is stored explicitly in that index. No parent on this queue is a mathematical counterexample.

The combined processed set contains 3,460 unique original indices. Two indices below the current numerical frontier (`3454`, `3458`) are simply unprocessed; this is an execution-order artifact, not a geometry gap.

## 6. Important bottleneck lessons

### What worked

- preserving common physical weight dependence is the major improvement;
- centered second moment is more interval-stable than direct covariance subtraction;
- local refinement is extremely effective after physical coupling is preserved;
- the old conclusion that `t` is always the privileged split axis no longer holds universally: after physical coupling, many difficult parents close under several axes;
- deterministic fixed depth-5 closes the overwhelming majority cheaply; deeper fallback should be reserved for the small repair queue.

### What did not work well

- direct interval covariance `E[ZZ^T]-E[Z]E[Z]^T` reintroduced dependency inflation;
- matrix preconditioning gave only marginal additional closure;
- global scheduling/ML was useful only for navigation and is not part of any proof;
- repeating global grid refinements from R7N is obsolete and should not be resumed.

## 7. Exact continuation protocol

### Step 1 — verify predecessor dependency

The current source scripts expect access to the R7N predecessor sources. Verify the old ZIP hash first:

`sha256sum FIN_R7N_HANDOFF_20260920.zip`

It must equal:

`ba56cc0ecf69ac970e1d6fe719e0770753437213a4ef345949ffcf1411941d9c`.

The run that produced this handoff restored it under:

`/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920`

### Step 2 — finish the fast depth-5 first pass

From `fin_rank7_targetp_physical_full_20260920/` run the depth-5 worker in bounded chunks, e.g.:

`python src/stream_depth5.py 5 20 10000`

The worker is append-safe/atomic at the parent-file level and rebuilds its index. Continue until all 5,432 original parents are either closed or present on the frozen repair queue.

Do **not** reinterpret an interval failure as a counterexample.

### Step 3 — freeze the complete repair queue

After the fast pass is complete, make a single immutable list of all unresolved original parent indices. Record:

- count,
- exact volume fraction,
- hashes of their depth-5 proof objects.

### Step 4 — repair only that queue

Use the dynamic local fallback logic from `src/stream_cover_v2.py` as a **repair-only worker**. Do not rerun it blindly over all parents, because its current selection logic predates the depth-5 stream and would duplicate work.

Recommended repair policy:

- retain fixed depth-5 result as the parent proof prefix;
- on each unresolved terminal, choose among `r,s,t` only;
- axis selection is navigation only: each child must pass the same physical interval certificate;
- exact rational child coverage is mandatory;
- use a bounded local extra depth first (current prototype used up to 7 extra levels);
- if a terminal survives, isolate it and analyze its center margin / one-axis splits before increasing depth globally.

### Step 5 — global audit only after residual = 0

Target P may be promoted to global only if every one of the 5,432 frozen R7N residual parents has a complete exact-cover proof tree and the pre-existing R7N tails/direct-safe regions are included.

Required final checks:

- 5,432/5,432 original residual parents accounted for exactly once;
- no unresolved terminal leaf;
- exact child-volume conservation at every repair tree;
- formula-level replay of every terminal physical-coupled SAFE certificate;
- mutation tests: deleted leaf, shifted boundary, altered threshold, altered fixture/model dependency;
- combine with the already accepted R7N tails / non-residual compact cover;
- only then state `lambda_2(M4) <= 67/250` globally.

## 8. Scientific nonconclusions

At this handoff:

- **Target P remains open globally.**
- No certified Target-P counterexample has been found.
- Target S at the sharp `sigma` is still open.
- Nothing in R7O/R7O2 upgrades tau0 results to sigma.
- No new full-X7/global-energy theorem is implied.
- The exploratory ~2.26% residual is not the canonical portable proof state.
- The exact phase results from R7N remain accepted and unchanged.

## 9. Strongest next atom

The next atom is no longer “invent a better enclosure.” That part succeeded.

**Finish the remaining 1,972 fast-depth-5 parents, freeze the resulting repair queue, and run the physical-coupled dynamic fallback only on that queue.**

The key stop condition is binary:

- if all repair parents close: build the full 5,432-parent audit and attempt global Target P;
- if some survive bounded local fallback: stop broad computation and analyze only those surviving microscopic terminal cells.

Do not start Target S or full X7 before this decision is resolved.
