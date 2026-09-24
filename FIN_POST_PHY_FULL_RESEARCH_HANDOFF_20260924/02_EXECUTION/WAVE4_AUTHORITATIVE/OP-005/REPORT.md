# OP-005 — One joint error ledger for carrier, channel and memory

Status: **PASS_CONDITIONAL_END_TO_END_LEDGER_WITH_SPLIT_DATA_AND_PREDICTABLE_REUSE_RULES**.

## Why composition was needed

The repository contains individually valid conditional pieces: a calibrated/no-fit channel protocol, checkpoint-valid matrix rank testing, carrier/readout identifiability conditions, and explicit detector-memory aliases. Their guarantees do not automatically survive adaptive reuse of the same data or post-hoc protocol choice.

## Frozen protocol graph

The safe default workflow is:

1. **Preregistration.** Freeze hypotheses, null classes, transformations, allowed controls, all outcomes (including return/no-escape/no-click), stopping rules and budgets before scoring data.
2. **Calibration block `D_cal`.** Independently establish the nuisance set for control ratio, efficiency/confusion support, detector reset/memory bound and preparation map. C13's channel design already reserves `delta_cal=0.025` for its external calibration gate. If the set is empty or too wide, output `NO_DECISION`.
3. **Memory block `D_mem`.** Use a typed broad/null class and either a fixed/checkpoint matrix certificate or another theorem-valid procedure. If detector memory is not independently bounded, the only allowed attribution is to the combined system+detector process, not to the system alone.
4. **Channel/carrier block `D_chan`.** Use a frozen no-fit statistic (the C13 two-control classifier or the clock-map-free escape profile). The choice between candidate channel tests may depend on earlier blocks, but the scoring block must be fresh or covered by a simultaneous/predictable-reuse theorem.
5. **Joint claim.** Emit a conjunction only if every prerequisite gate passes; otherwise abstain at the first failed gate. Raw attempted records are never dropped by postselection.

A concrete 5% ledger can reserve

`delta_cal=0.025`, `delta_mem=0.0125`, `delta_chan=0.0125`.

Existing 5% subtests cannot simply be inserted unchanged into this 5% joint claim; their thresholds/sample counts must be recomputed at the allocated stage level. If one insists on keeping the existing 5% rank test together with the 2.5% calibration gate and another 5% channel claim, the resulting familywise bound is at best 12.5% by this ledger, not 5%.

## Composition theorem

Let `F_k` be the information available before stage `k`. If each reached stage satisfies, uniformly over the declared composite null/nuisance set,

`P(E_k | F_(k-1), stage k reached) <= delta_k`,

then for any adaptive ordering/stopping rule based only on past information,

`P(any false promoted claim) <= sum_k delta_k`.

No independence between stages is required. Abstention only reduces the chance of reaching later stages.

For e-values, multiplication is licensed only when the next factor is conditionally safe,

`E[E_k | F_(k-1)] <= 1`,

so that the running product is a nonnegative supermartingale/e-process. Marginal validity of dependent e-values is insufficient: if `E1=E2=2` with probability `1/2` and `0` otherwise, each has expectation `1` but `E[E1 E2]=2`.

## Shared nuisance variables

Calibration uncertainty must be carried as a **joint confidence set**, not as separately optimized best cases for topology, readout and memory. Every downstream guarantee is required to hold uniformly over the surviving nuisance set. Calibration drift is handled by an explicit drift envelope or triggers abstention/recalibration; it cannot be silently absorbed into the physical model.

## Adversarial ordering / stopping

Because stage validity is conditional on the full past, an analyst may stop after a failed gate, choose among preregistered downstream branches using prior blocks, or alter sample allocation predictably. What is not allowed is choosing a feature/test after seeing the same scoring outcomes unless a simultaneous selection correction has been proved.

## Operational consequence

This supplies a coherent **conditional** workflow, but does not make the currently enormous memory-rank cost practical. MEM-009's broad-null bracket remains roughly `2.64e5` necessary-design to `4.29e9` sufficient all-time resets, and OP-002's clock-free channel route is also expensive. The ledger prevents these costs from being hidden by adaptive reuse.

No empirical execution, calibration custody or physical success is claimed.
