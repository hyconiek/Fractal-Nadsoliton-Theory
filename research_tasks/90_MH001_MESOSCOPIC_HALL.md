# MH-001 — One frozen growing complement: lattice shell and adjacent cancellation

Execution prompt for GPT-6.1 Sol. **UNEXECUTED; OPTIONAL P2**.
Cwd: `/media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison`.
Write only `fin_mesoscopic_hall_20261019/`.

## 1. Problem

Prove or isolate the SHELL and ADJACENT estimates for the single radius
L=1/10 fast-coordinate ball in the current DH seed. Then, only if both are
paid, prove Hall for its dense-complement cut family.

## 2. Significance

This reaches a positive-mass growing cluster instead of another finite
orbit of holes. It is a technical advance toward the frozen carry ansatz,
not an intrinsic activity source or full IBA counterexample.

## 3. Inputs

- `AGENTS.md`, `research_tasks/00_EXECUTION_RULES.md`.
- `FIN_MESOSCOPIC_HALL_AND_INTERMEDIATE_SOURCE_NEXT_CHAIN_20261017.md`.
- `fin_dense_orbit_fibre_transport_review_20261017/CLUSTERED_COMPLEMENT_SEED.md`.
- `fin_dense_orbit_fibre_transport_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md`.
- `fin_dense_orbit_fibre_transport_campaign_20261017/DH-001/DENSE_COMPLEMENT_HALL_THEOREM.md`.
- `fin_balanced_phase_cut_research_20261016/COMPACT_PHASE_PREFIX_LOCAL_LIMIT.md`.
- `fin_balanced_phase_cut_research_20261016/EXACT_COMPACT_CURRENT_MARGINALS_AND_GLUING.md`.
- `fin_balanced_phase_current_campaign_20261016/MARGINAL_AND_SUPPORT_CONTRACT.json`.
- `fin_balanced_phase_current_review_20261017/INTAKE_20261017.md`.
- `fin_rbcut_weighted_degree_review_20261016/PROOF_AUDIT_AND_CORRECTIONS.md`.

## 4. History and duplicate check

KR/BA/FP/RB/RBCUT/SC/PC/DH are predecessors, not unfinished tasks to restart.
Separated holes and one bounded 36-predecessor cluster are already paid.
The strong shielding-weight route has a real counterexample. The current
MH order explicitly leaves SHELL and ADJACENT unproved.
Check `rg -n 'A_ball|ADJACENT|SHELL|shrinking.annul'` in the named DH/PC
packages and any newer mesoscopic-Hall output before proceeding.

## 5. Paid geometry

Keep N in12N, exact g_c, lambda=0, actual D_N/E_N, phases, full H support
and corrected I',O'. On adjacent compact fibres |z|<=1 put

    p_m=u+m sigma/12, Y_N=sqrt(N)(p-p_m),
    A_ball={||Y_N||<=1/10}.

The 36 moves v_ab=e_b-e_a-sigma/6 have zero mean, second moment P_fast/6
and norm²=5/3. The seed already proves predecessor erosion by
1/(8sqrt(N)); reuse it. Charge true E_N holes separately.

## 6. Exact missing lemmas

For B_right the un-eroded right ball, uniformly over the compact band,

    I'_(l+1)({L-1/(8sqrt(N))<||Y||<=L})
       >=c_* N^(-1/2) J_A,l, c_*>0;                          SHELL
    |O'_l(A_ball)-I'_(l+1)(B_right)|
       =o(N^(-1/2) J_A,l).                                 ADJACENT

Do not replace these by weak convergence or unrelated o(1) statements.

## 7. Proposed improved route

First verify that even/odd count-sum constraints reduce the relevant fast
lattice translations to finitely many A5 direct-sum A5 residue classes
(candidate indexing: even count modulo6 when N in12N). This is a proof
route to establish, not a new assumed lattice law.

Then seek a uniform weighted lattice-annulus estimate with an error below
its width scale, and a **differenced** adjacent-level Stirling/Edgeworth
expansion. Exploit exact carry/source cancellation before bounding errors.
Separate the bulk expansion, residue-class discrepancy and exponential E
correction. Do not use a smooth Gaussian Riemann sum with an unproved sharp-
boundary error of the same order as the desired shell mass.

## 8. Actions

1. Bind exact lattice cosets, membership and the frozen carry formulas.
2. Attempt one shell count estimate with uniform lower margin.
3. Attempt the adjacent cancellation at the required smaller scale.
4. Charge E terms and apply O'(A)-I'(R) only if both estimates pass.
5. Report precisely which larger cut families remain outside the theorem.

## 9. Premises

The ball is an analysis norm, not physical space. Keep exact count
multiplicities, phase ties, source and allowed departure-layer interphase
edges. Do not replace the full H by a stronger within-phase graph. Unknown
asymptotic N0 is allowed if the theorem genuinely proves eventual uniformity.

## 10. Results

One mesoscopic-family Hall theorem, a valid counterexample to the frozen
estimate, or a specific open lattice/cancellation lemma. No full Hall or
actual globally glued rates follow automatically.

## 11. Verification

Independently verify the lattice translation and residue calculation on
small exact count fixtures. Reconcile prefactors/covolume with PC's actual
slice theorem. Numerical lattice counts may expose a false estimate but
cannot prove all-N uniformity. A failed sufficient estimate is not proof
that all H-supported flows are impossible.

## 12. Dependencies

No ID/RP result required. This is an optional mathematical lane and does
not supersede the fundamental source priority. Full coupling and global
gluing remain later, separately gated tasks.

## 13. Failure

Retain the exact order/sign of the failing remainder. Do not tune L,
threshold, phase partition or marginal ansatz. Return to Astra before a
different growing family is considered.

## 14. Stop

One residue/lattice-expansion route plus one justified repair. No whole
count-simplex enumeration, new bounded-hole variants or cap-digit campaign.

## 15. Outputs

`LATTICE_CONTRACT.json`, `SHELL_PROOF_OR_OPEN_LEMMA.md`,
`ADJACENT_PROOF_OR_OPEN_LEMMA.md`, `MESOSCOPIC_HALL_DISPOSITION.md` and the
common handoff files. Verdict: `ONE_MESOSCOPIC_FAMILY_PAID`,
`PARTIAL_LEMMA_ONLY`, `COUNTEREXAMPLE` or `OPEN_LEMMA`. Stop.
