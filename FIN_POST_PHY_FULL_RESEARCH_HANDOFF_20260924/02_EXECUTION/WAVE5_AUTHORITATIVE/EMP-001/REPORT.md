# EMP-001 — Platform readiness and information-budget gate

Status: **NOT_READY_NO_REAL_PLATFORM_OR_INDEPENDENT_CALIBRATION_RECORD_SUPPLIED**.

This task is an empirical gate, but external data collection/contact/purchases are not authorized and no named apparatus record is present in the repository handoff. The correct output is therefore a precise readiness failure, not a simulated experiment.

## Required records that are currently missing

A platform cannot enter OP-005 until it supplies, independently of scoring data:

- a preparation/reset specification and measured reset fidelity for the exact null/process class being tested;
- a detector-memory/reset bound sufficient to separate system memory from detector memory (OP-004); otherwise only the composite system+detector rank is identifiable;
- a complete outcome map retaining no-click, return and no-escape attempts;
- a control-map calibration supporting the CHANNEL-05 conditions `eta<=0.02` and first dimensionless time in `[0.23,0.24]`, with calibration failure budget `0.025`;
- an efficiency/confusion model. If the C13 support-preserving row-stochastic route is used, its primitive safe radius is `0.00941745678649` at relative efficiency uncertainty `r=0.02`; the support pattern itself must be independently justified;
- calibration drift bounds over the scoring run;
- per-reset/per-attempt time, failure and monetary costs, so the information budget can be compared with certified sample requirements;
- raw-data custody/provenance sufficient to prove that calibration, training and scoring blocks were not silently reused contrary to OP-005.

## Current certified cost scales to be matched by a real platform

- clock-free carrier recovery: 435 recorded exits for one prepared state to recover the top-two set, or 612 exits per state for simultaneous `12x11` frequency protection at family error <=5% in the ideal embedded-jump model;
- specific strict-vs-shell-swap negative control: 31 exits from one known preparation reaches the conservative ideal Chernoff bound for Bayes error <=5%;
- CHANNEL-05: 31,299,102 total main attempts, plus an **unknown** calibration-shot count because no apparatus calibration noise law is supplied;
- OP-002 clock-map-free unitary-vs-wave profile: about `4.71e8` expected attempts at the displayed strict-model point;
- broad memory-rank exclusion: current necessary-design to sufficient all-time bracket is roughly `2.64e5` to `4.29e9` resets.

Without a platform rate/cost/fidelity record, none of these counts is a feasible experimental recommendation.

## Verdict

`NOT_READY` is forced by missing nature-side inputs, not by failure of the mathematical protocols. No empirical validation, equipment recommendation or purchase is made.
