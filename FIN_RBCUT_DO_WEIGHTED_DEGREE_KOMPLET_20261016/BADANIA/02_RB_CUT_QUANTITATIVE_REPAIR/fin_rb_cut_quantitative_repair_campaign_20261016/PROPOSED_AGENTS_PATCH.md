# Proposed guardrail patch ONLY — not applied

In the RB-CUT / rare-boundary flow register, consider adding:

- A new finite capacitated construction lemma: `K_N^* = max(0,sup_S sum_S s / c0_H(out S))`, with positive/zero interpreted as infinity. Finite K_N^* is equivalent to old lower-ideal feasibility; uniform K_N^* is a **sufficient** exact-drift and O(N)-rate construction, not yet paid for actual q12.
- Signed allowed-edge repairs must be permitted. A nonnegative-only augmentation can falsely rule out legitimate final conductances.
- Single deleted upward parity edges admit no upward H-path replacement with the same endpoints (one-step parity grading). Collective cancellation is still possible.
- An unbounded native-relative multiplier does not refute the original uniform absolute vertex-rate constraint.
- No q12 all-N K bound, flow family, negative ideal, IBA counterfamily, independent physical kinetic source or ToE status is changed by this package.
