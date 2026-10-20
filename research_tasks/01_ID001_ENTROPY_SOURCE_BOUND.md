# ID-001 — Native time resolution versus finite-entropy deterministic source

Execution prompt for GPT-6.1 Sol. Date: 2026-10-19. **UNEXECUTED; P0**.
Work from `/media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison`.
Write only `fin_intrinsic_entropy_source_bound_20261019/`.

## 1. Mathematical problem

Derive a quantitative necessary source-entropy/error budget for a stationary
deterministic process to reproduce the actual native FIN N2 difference
process at time resolution delta. Test a theorem, not a new microscopic law.

Target: a measure-preserving flow T_t on a standard probability space,
H=h_mu(T_1)<infinity, and deterministic measurable q12 readout. At
delta=1/m, let Y_k be its stationary readout and X_k the native N2 skeleton.
For a stationary joining with mismatch e_delta=P(X_0!=Y_0)<=11/12, seek

    h_delta <= delta H + h2(e_delta) + e_delta log(11),             (* )

where h2 is binary entropy in nats, not a static spectral histogram.
Allow a different readout/joining for each delta if the proof supports it;
keep the same finite H bound. Mark the actual supported quantifiers.

## 2. Why it matters

This tests whether finite-dimensional smooth deterministic informational
dynamics could be an exact all-resolution source of the native jump law.
It is upstream of identifying that law with physical time. It does not
derive the source or exclude infinite-dimensional/singular-limit alternatives.

This is a conditional test of an effective coding of intrinsic relational
histories after a comparison clock chart has been declared. It does not
assume that the primordial nadsoliton has an external time, finite dimension
or finite KS entropy. Do not infer those properties from the twelve labels.

## 3. Exact inputs

- `AGENTS.md`, `research_tasks/00_EXECUTION_RULES.md`.
- `fin_strategy_audit_20261019/INTRINSIC_EVOLUTION_ASSESSMENT.md`.
- `fin_strategy_audit_20261019/QW206_INSPECTION.json`.
- `fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md` §§1–2.
- `fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py` (read formulas; do not run its writer).
- `FIN_ST383_Asymptotic_Markov_Rate_Distortion_Bounds.json`.
- `FIN_ST372_ST386_Global_Collision_Envelope_Exact_Invariant_Bases_and_Continuum_IR_Report_EN.tex`, ST383 subsection.
- `fin_unified_strict_evolution_audit_20261017/AUTONOMOUS_COMPOSITION_THEOREM_OR_OBSTRUCTION.md` §§1–2.
- `fin_kinetic_composition_research_20261016/THEORETICAL_RESULTS.md` §3.
- `FIN_Ten_Research_Programs_Monograph.md`, Program 15.
- `fin_intermediate_clock_bridge_research_20261017/TRACE_CLOCK_AND_INSTRUMENT_THEOREM.md` §§1–4.
- `fin_strategy_audit_20261019/references/lecomte_2007_sections_2_2_to_3_2.txt`
  and the sibling original PDF. Read §§2.2–2.3 and the distinction at the end
  of §3.2; do not use its regularized jump-sequence entropy as metric entropy.

## 4. Prior campaigns and duplicate result

Program 15 already proves a finite-dimensional unitary-bath obstruction.
KCR gives a regular-observable compression bound. ST383 already uses
per-symbol Fano for a supplied stationary strict-heat source at a fixed
sampling convention. The intermediate model already distinguishes trace
time, last-visible records and raw paths. Do not reproduce those as the
main discovery. The general infinite metric entropy of nontrivial CTMC
histories is established literature, not a new FIN theorem.

The architect's final search found these equivalents, not the full native
resolution-dependent source/joining inequality. Before execution run a
targeted `rg -n -i 'stationary.joining|entropy.*resolution|ID-001|finite.entropy'`
over the listed predecessor directories and any newly present
`fin_intrinsic_entropy_source_bound_*`. If an equivalent accepted result
now exists, reuse it and return DUPLICATE_TARGET, without recomputation.

## 5. What is already paid

For any fixed finite g and the supplied circulant A7, p_d is proportional to
exp[(g/2)A_0d], p_d>0, and the exact difference generator is

    B=2(1 p^T-I),
    P_delta=e^(-2 delta)I+(1-e^(-2 delta))1 p^T.

This quotient and its clock are model inputs. Use the exact p(g), not the
rounded matrix or an assumed uniform p at positive g.

## 6. Genuinely open part of this assignment

Bind the entropy-rate comparison to this actual process and a deterministic
source class, including approximation and resolution quantifiers. Determine
whether it yields a nontrivial lower error bound as delta tends to zero.
The static QW206 histogram is already identified as insufficient; another
histogram run is not needed.

## 7. Proposed proof route

First establish the block-entropy/Fano joining bound, then pass to entropy
rates. Bound the readout entropy by h_mu(T_delta). For delta=1/m, prove or
properly cite the transformation power rule
h_mu((T_delta)^m)=m h_mu(T_delta), giving h_mu(T_delta)=H/m.
This avoids an unproved arbitrary clock-change theorem.

Use the exact P_delta to expand

    h_delta=-sum_i p_i sum_j P_delta(i,j)log P_delta(i,j),
    rbar=2(1-sum_i p_i^2)>0.

Target h_delta/delta=rbar log(1/delta)+c(p)+O(delta|log delta|)
with a proved remainder for fixed p, or a declared compact positive-p class.
Do not silently impose p_min uniformly over unbounded g.

## 8. Concrete actions

1. Freeze source space, invariant measure, deterministic observation, shared
   model-time convention, joining and entropy definitions.
2. Prove (*) or provide a valid counterexample to its exact stated scope.
3. Derive the native small-delta expansion from the semigroup, including
   diagonal/no-jump terms and a remainder.
4. Test the proposed consequence `liminf e_delta/delta >= rbar` for a common
   finite H bound. It is a TARGET, not a premise. Also retain the finite-delta
   inequality as the useful source budget if a sharper corollary fails.
5. Apply the result conditionally to the QW R26 polynomial Hamiltonian family
   on a supplied compact invariant/coercive energy set. Prove finite metric
   entropy from appropriate regularity, or leave that application conditional
   on H<infinity. Do not apply it automatically to the NS/RP PDE.

## 9. Assumptions and limits

The stationary joining of whole sampled sequences is stronger than finite-
window weak convergence or a single two-time covariance match. No implication
from the latter is allowed without proof. An external random readout/clock
adds information and leaves the deterministic-source class. Constant clock
rescaling does not source seconds. If H grows with resolution, report the
required growth instead of claiming a universal no-go.

## 10. Expected results

One scoped entropy-budget theorem and native application; or a precise open
lemma/counterexample. Record exact realization and approximation separately.
No new law or empirical success is requested.

## 11. Acceptance and independent checks

Check the entropy bound by a separate elementary block-chain-rule argument.
Use the exact N2 reset matrix to verify the leading coefficient independently
of a generic generator expansion. Small numerical evaluations are diagnostics
only. Mandatory controls: fixed positive sampling interval; discrete-time
Bernoulli dynamics with finite entropy; infinite-entropy/path-space source;
same static p with different path laws; stochastic readout; and unchanged
two-time correlations without the required full joining. Positive KS must
not be identified with thermodynamic irreversibility.

## 12. Dependencies

No RP/Hall/B4 completion is required. This task uses existing model laws as
targets, not as sources of themselves. No independent physical hypothesis
is adopted; H<infinity specifies the candidate class being tested.

## 13. If the proof fails

Identify the first failed entropy/regularity/joining implication. Preserve
any valid finite-delta inequality. Do not insert hidden randomness or change
the metric to rescue the advertised conclusion. Return to the architect.

## 14. Stop rules

One route plus one explicit correction. No chaos simulation, broad g/N scan,
new microscopic principle or full NS PDE entropy campaign. Stop after the
scoped theorem/counterexample/open lemma and its controls. A known general
theorem alone, without the native budget/application, is not a new result.

## 15. Output and return

Export `SOURCE_CLASS.json`, `ENTROPY_BUDGET_PROOF.md`,
`NATIVE_RESOLUTION_BOUND.md`, `COUNTERCONTROLS.md`, bounded verifier/results
and the common handoff files. Verdict: `SCOPED_SOURCE_BUDGET_PROVED`,
`PARTIAL_BUDGET_ONLY`, `COUNTEREXAMPLE` or `OPEN_LEMMA`.
Return to Astra. Do not automatically propose a bath, continuum or physical
generator after excluding this one source class.
