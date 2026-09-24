# EMP-002 — Preregister a discriminating pilot, not a full rank campaign by default

Status: **NO_FEASIBLE_REAL_PILOT_WITHOUT_PLATFORM__PARAMETERIZED_LOWEST_COST_CARRIER_PILOT_FROZEN**.

EMP-001 admits no real platform, so this task cannot honestly declare a feasible physical pilot. It can, however, freeze the lowest-cost **conditional** pilot that should be considered first if a future platform clears the readiness gate.

## Cost/value comparison

The current programme strongly disfavors beginning with full memory rank or the clock-map-free wave/unitary profile:

- broad memory rank: `~2.64e5` necessary-design to `4.29e9` sufficient all-time resets;
- clock-map-free wave/unitary profile: `~4.71e8` expected attempts at the displayed strict-model point;
- CHANNEL-05: `31,299,102` main attempts plus unknown calibration shots;
- embedded-jump carrier relation: hundreds to thousands of **recorded exits**, with a particularly strong frozen shell-swap negative control having TV `0.3348080497`, KL `0.2996465192` nats and Chernoff information `0.07654211034`.

The carrier observation therefore removes a large class of wrong relational carriers at much smaller ideal information cost and is exactly invariant to global clock rescaling.

## Parameterized pilot `EMP-002-CARRIER-01`

**Prerequisites:** EMP-001 READY; independent detector-memory/confusion calibration; a preparation that identifies one starting state without using the scoring exits; all exits/no-clicks retained.

**Primary no-fit endpoint:** from each prepared state, the two empirically dominant destinations must assemble into one connected 2-regular graph on 12 states. Graph recovery is label-free; only after recovery are shell distances assigned.

**Frozen scoring options:**

- minimal one-state negative-control discrimination: maximum 31 recorded exits for the *specific* strict-vs-shell-1/2-swap pair under the ideal iid model and equal-prior Chernoff calculation;
- structural carrier recovery: 435 recorded exits for a single prepared state top-two set at family error <=5%; for full simultaneous destination-frequency protection, 612 exits per state (7,344 recorded exits total for 12 states).

These are mathematical record counts, not apparatus attempt counts: a real efficiency/no-click model must convert attempts to exits before feasibility can be claimed.

**Abstention:** any calibration failure, detector-memory ambiguity, missing outcome, insufficient exits by the maximum budget, or graph that is not connected 2-regular returns `NO_DECISION/FAIL_CARRIER`, not a refit.

**Held-out prediction:** after the graph is recovered, shell-group destination frequencies on fresh scoring data must match the frozen strict profile within a preregistered simultaneous confidence region. No cyclic labels or absolute clock are imported.

## Verdict

This is the preferred **future first pilot by current certified information scale**, but because no platform has passed EMP-001 the present empirical verdict remains `NO_FEASIBLE_REAL_PILOT`. No experimental evidence is generated here.
