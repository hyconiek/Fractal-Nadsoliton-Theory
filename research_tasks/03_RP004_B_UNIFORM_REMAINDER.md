# RP004-B — Trace-admissible nonlinear remainder on the frozen log window

Execution prompt for GPT-6.1 Sol. **UNEXECUTED; P1, optional field track**.
Cwd: `/media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison`.
Write only `fin_rp004_flux_continuation_20261019/RP004-B/`.

## 1. Problem

Control the actual prepared nonlinear solution uniformly up to
T_epsilon=12 log(1/epsilon), in a norm that pays boundary derivative traces.
No positive flux is assumed or required at this stage.

## 2. Importance

The current O_T(epsilon^j) expansions have uncontrolled T dependence.
An estimate on the target window is needed before a small radiative term
can dominate the exact flux errors.

## 3. Exact inputs

- `AGENTS.md`, `research_tasks/00_EXECUTION_RULES.md`.
- `research_tasks/RP004_SHARED_CONTRACT.json`.
- `fin_coherent_mode_radiation_review_20261019/RETARDED_FLUX_PROOF_ROUTE.md`.
- `fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md` §§1–2,6–7.
- `fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md`.
- `fin_nonlinear_stability_transfer_review_20261017/INTAKE_20261017.md`.
- `fin_coherent_mode_radiation_campaign_20261017/MODE_AND_FLUX_CONTRACT.json`.
- `fin_coherent_mode_radiation_campaign_20261017/ERROR_SCOPE.json`.
- `FIN_RP004_RETARDED_FLUX_CONTINUATION_20261019.md`.

## 4. History and duplicates

NS already provides full-field orbital persistence, moving-centre continuation
and local transfer. RP provides fixed-time smooth expansions. None supplies
the growing-time error in the latest reviewed files. Do not redo NS or infer
this estimate from its PASS. Targeted freshness check:
`rg -n -i 'logarithmic|uniform.*remainder|trace.admissible'` in the listed
RP/NS reviews and any actual RP004 continuation.

## 5. Paid prepared reduction

The actual action/data preserve the real e1 species sector, phi odd and
psi even. This is an invariant Cauchy class, not a truncation for arbitrary
data. Set K=tanh x, eta=phi-K, L=-d_x²+4-6sech²x, L_a=L+a+2,
v=sqrt(3/4)sech²x cos(sqrt(a+2)t), and use the zero-data retarded u from A.

The exact equations are

    eta_tt+L eta=-6K psi²-6eta psi²-6K eta²-2eta³,
    psi_tt+L_a psi=-12K eta psi-6eta²psi-g4 psi³-g6 psi⁵.

For r=psi-epsilon v and w=eta-epsilon²u, use exactly the remainder equations
in RETARDED_FLUX_PROOF_ROUTE. All independent r,w Cauchy data are zero.

## 6. Missing estimate

Seek uniform bounds for r,w in H2 x H1 (or a proved equivalent graph norm).
Possible sufficient targets are O(epsilon³(1+t)^p) and
O(epsilon⁴(1+t)^q), with finite explicit powers and constants. Values p=2,
q=3 are proposed bootstrap guides, not premises. A different proved bound
is acceptable only if its actual errors still vanish at the frozen log time.

## 7. Method

Absorb the full variable-coefficient linear operator into its wave group.
Use its positive energy on the actual parity classes (odd scalar gap and
positive transverse gap), graph-norm equivalence and 1D product estimates.
Estimate u from its forced linear equation. Bootstrap the small nonlinear
coefficients rather than replacing the linear dynamics by a coarse
Gronwall factor exp(Ct) without checking its exponent at t=12log(1/epsilon).

## 8. Actions

1. Prove the relevant domain/norm equivalence and parity-preserved gap.
2. Set explicit bootstrap inequalities for eta,psi,r,w, with fixed finite
   signed g4,g6 and constants depending on declared finite bounds M4,M6.
3. Close them on the proposed interval uniformly for a in[-3/4,-1/4], or
   show precisely which bound prevents it.
4. Translate the result to endpoint derivative traces and time-integrated
   products required by the original flux.
5. State epsilon0 and constant dependence honestly; numerical values are not
   required if explicit finite analytical bounds suffice.

## 9. Scope

Do not claim uniformity over unbounded nonlinear coefficients. Bounds M4,M6
specify a conditional class, not natural values or fitted parameters. Do not
pin the wall: odd eta removes the even translation direction for this datum.
Do not remove the odd scalar shape mode. All stronger norm claims need their
own domain/product proof. No Phi6 or additional forcing is allowed.

## 10. Output types

Uniform comparison theorem; explicit failure of the specified bound; or an
open product/domain/bootstrap inequality. Fixed-T control alone is partial.

## 11. Independent checks

Re-derive remainder signs by substituting phi=K+epsilon²u+w and
psi=epsilon v+r into the original equations. Verify zero initial data and
epsilon parity independently. Check the endpoint norm embedding explicitly.
Small symbolic checks are support, not the infinite-dimensional proof.

## 12. Dependencies

Independent of A's eventual flux lower bound, but bound to the same u,
normalization and time/radius contract. It may run beside A. C needs both.

## 13. Failure handling

Keep any valid norm estimate and exhibit why it does not pay the target
window or flux trace. No silent smaller horizon, unknown derivative trace,
or numerical fitting of the error exponent.

## 14. Stop

One graph-energy bootstrap and one directly justified repair. Stop at a
proved bound or exact open lemma; no large lattice or long-time simulation.

## 15. Deliverables

`UNIFORM_REMAINDER_OR_FAILURE.md`, `STRONG_NORM_AND_TIME_SCOPE.json`,
`TRACE_ERROR_BUDGET.json`, common proof/report/handoff files.
Verdict: `LOG_WINDOW_COMPARISON_PAID`, `PARTIAL_FIXED_TIME_ONLY`,
`TARGET_ESTIMATE_FAILED` or `OPEN_LEMMA`. No flux promotion at this stage.
