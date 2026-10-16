# FIN after BA: one native count-partition feasibility problem

Date:2026-10-16. Status: **RESEARCH_PLAN; FP-001–FP-004 NOT EXECUTED.**
Execution directory: `fin_count_partition_feasibility_campaign_20261016/`.
Language: English. Decision: [partition decision](fin_reversible_bracket_averaging_review_20261016/PARTITION_DECISION.md).

## 1. New target, not another equivalent open averaging label

BA has paid exact finite parity-cut/Dirichlet/conditional-noise rigidity.
KR under NC/IBA and the finite Poisson criterion are reused, not rerun.
Universal q12 IBA without NC is still open.

Test the single candidate in
[COUNT_PARTITION_FEASIBILITY_SEED.md](fin_reversible_bracket_averaging_review_20261016/COUNT_PARTITION_FEASIBILITY_SEED.md).
It is an intrinsic function of existing count variables, not a new hidden
state or physical controller. Its limiting masses and IF-feasible IBA
obstruction have been derived. What remains is **actual admissible count
conductances**, with exact statewise drift and one N-independent rate cap.

One valid uniformly capped closed-partition family would give

    liminf E sup_(t<=S)|integral_0^t (G_N-EG_N)ds|>=S M>0,
    M=E_rho0[Z^2 exp(-Z^2)]>0.

No such family is claimed here. Infeasibility kills this frozen candidate,
not every countermodel or the universal IBA proposition. Do not tune a
sequence of thresholds until one appears to work.

## 2. Freeze the object and inputs

Read and hash current AGENTS.md and:

| Alias | Input |
|---|---|
| P01 | fin_reversible_bracket_averaging_review_20261016/INTAKE_20261016.md, PROOF_AUDIT_AND_SCOPE.md, CLAIM_REGISTER.json |
| P02 | Same review: COUNT_PARTITION_FEASIBILITY_SEED.md and verification.json |
| P03 | BA source BA-001/COMPONENT_LEDGER.md and CLASS_AND_FLUX_CONTRACT.json |
| P04 | KR review PROOF_COMPLETIONS_AND_SCOPE.md; sufficient process theorem |
| P05 | Actual MT review/proofs, exact D_N definition, conditional Gibbs law and moment bounds |
| P06 | Actual UF-002 native scalar drift expansion and versioned safe_seed/spectral inputs |
| P07 | AC full-coordinate wheel and KR scalar spike negative controls |

Fix lambda0, g_N=g_c=12/Lambda6 exactly, N in12N, the actual MT reflected
domain and count weights W_N including multinomial multiplicity. Freeze
epsilon=1/8 and the existing integer k3 Fourier sums C(n),S(n):

    Z_N=N^(1/4)m, R_N=2(C(n)^2+S(n)^2)/N,
    p(z)=1/2+(1/8)exp(-z^2),
    T(z)=-(2/gamma3)log(1-p(z)), gamma3=1-Lambda3/Lambda6,
    A_N={R_N<=T(Z_N)}, B_N=D_N\A_N.

The candidate diagnostic may set rates between the two phases to zero.
Every other hypothesis stays unchanged: reversibility with the same pi_D,N,
single-copy support, exact native reflected parity drift at EVERY state
and total exit rate<=CN with one finite C independent of N. There is no
extra state variable or modification of primary FIN. The diagnostic class
does not require full coordinate drift matching or D12 generator covariance;
do not silently add or claim either. The partition itself is D12 invariant.

## 3. Task cards

### FP-001 — Verify the frozen phase contract against the actual domain

**Question:** are the phase geometry, static limits and exact necessary
balance identities correctly bound to the current MT/q12 objects?
**Inputs/dependencies:** P01–P07; current analytical seed, not an assumed counterfamily.
**Method:** bind exact Fourier normalizations, gamma3 and the domain's D12
invariance. Reuse the seed's radius CDF and moment UI rather than recount
them as discoveries. Retain ties and exact count classification; a floating
boundary cannot decide scientific feasibility. Export the actual statewise
source b_x=W_x d_native(x), not only the leading cubic drift. Check exact
phase balance and record any finer internal-component condition not paid.
**Deliverables:** FROZEN_PARTITION.json, INPUTS.sha256, EXACT_DOMAIN_AND_SOURCE.md,
STATIC_SCOPE_LEDGER.json and REPORT.md.
**Acceptance:** one correctly typed fixed partition and exact graph/source.
**Failure/stop:** using a large artificial cylinder, rounded g_c, modified
Gibbs law, appended phase state or unavailable conditional equilibrium.
**Must not conclude:** existence of a closed-phase generator or IBA failure.
**Priority/class:** P0/S–M.

### FP-002 — Prove feasible native flows or an exact dual obstruction

**Question:** can within-phase count conductances satisfy every microscopic
drift and capacity constraint with a uniform finite C?
**Inputs/dependencies:** FP-001 and the seed's exact primal/dual alternative.
**Method:** seek c_e>=0 with Ac=b and Bc<=CN W on legal within-phase edges.
Same-parity edges may be omitted for nonnegative feasibility. Either give
an analytic construction with paid nonnegativity and a fixed cap, or give
a valid dual certificate with unrestricted real phi and beta>=0 satisfying
every edge inequality but negative total pairing.
Prove the statement on the true all-N family, not by solving a tiny unrelated
LP. A uniform unbounded subsequence plus native generators elsewhere is
allowed only if it defines a complete admissible sequence with one cap.
No cap growing with N, near-drift equality or missing boundary constraints.
**Deliverables:** CONDUCTANCE_CONSTRUCTION_OR_DUAL.md, POSITIVITY_AND_CAP_PROOF.md,
FAMILY_OR_OBSTRUCTION.json, RESULTS.json and REPLAY.md.
**Acceptance:** a genuine asymptotic primal construction or a correctly
scoped obstruction/open feasibility lemma. An abstract phase-total balance
identity alone is insufficient.
**Failure/stop:** only interface-edge rejection or vanishing-mass repair
is proposed. The seed proves that successful rates must change Gamma in
stationary L1 by at least M>0. A local spike cannot achieve it.
**Must not conclude:** a numerical infeasibility at one small N settles
the asymptotic problem, or rejecting this partition proves universal IBA.
**Priority/class:** P0/M, one fixed proof attempt; no large lattice enumeration.

### FP-003 — Derive the actual temporal consequence only after feasibility

**Question:** does the constructed family yield the promised nonvanishing
time-integrated bracket, not merely a large sufficient Poisson bound?
**Inputs/dependencies:** a paid FP-002 construction; P02/P04/P05.
**Method:** apply conditional stationarity to each closed phase, transfer
the weighted native drift error and fourth-moment UI, and obtain the two
nonzero residual means. Use Jensen to prove the positive S M lower bound.
No separate limiting SDE is required just to refute IBA. If additionally
claiming a different process limit, pay its tightness/joint-martingale
identification separately. Optionally use the already derived native mixture
eta_N=N^-2 to obtain strictly positive irreducible rates, paying its full
count-path coupling error O(eta_N N^(3/2)). Do not claim a uniform lower
rate or native diagnostic mixing from that positivity.
**Deliverables:** IBA_OBSTRUCTION_OR_FAILURE.md, CONDITIONAL_MOMENT_LEDGER.json,
STRICT_POSITIVITY_SCOPE.md and NONCONCLUSIONS.md.
**Acceptance:** a genuine actual-q12 temporal counterfamily, if the primal
exists. If FP-002 is negative/open, mark this task NOT_ACTIVATED.
**Failure/stop:** assuming phases are closed before constructing rates,
or replacing the actual covariance functional by a new noise ansatz.
**Must not conclude:** native HB is refuted, a physical memory device exists,
or every stationary quartic approximation fails.
**Priority/class:** P1/S–M, strictly after construction.

### FP-004 — Exact architectural disposition and stop

**Question:** what did this single frozen candidate actually establish?
**Inputs/dependencies:** terminal FP-001–003 dispositions and current claim registers.
**Method:** choose NATIVE_IBA_COUNTERFAMILY, FROZEN_PARTITION_INFEASIBLE or
OPEN_UNIFORM_FLOW_LEMMA. Retain KR under NC/IBA in every case. A failed
candidate does not prove the universal theorem. Preserve all failures and
return for architecture; no automatic threshold changes or new graph/clock.
**Deliverables:** REPORT.md, HANDOFF.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch only.
**Acceptance:** reproducible scope with no forced success.
**Failure/stop:** simulation or an unchanged open statement is renamed
a construction, or a comparison generator is adopted as physical FIN.
**Must not conclude:** independent activity/A7 source, units/locality,
conservation/EV, selector/legacy-role, QM/GR/SM/L_total/ToE closure.
**Priority/class:** P0/S, final checkpoint.

## 4. Execution boundaries

Start FP-001/002 only. The present intake has derived the static candidate
and conditional consequences, not these missing feasibility rates. Per task
save REPORT/RESULTS, proof/open lemma, NONCONCLUSIONS, REPLAY and input hashes.
Work only in the new directory; preserve AGENTS, old plans/reviews/sources.

No external experiments, hardware, spending or messages. No broad source
inventory, classifier sweep, huge count enumeration, new physical noise,
graph, controller, bath or ontology premise. All physical source, units,
incidence/locality/conserved energy/calibration/EV, FIN-M0/B4/task337/task335,
selector/legacy-role and quantum/gravity/ToE gates stay unchanged.
