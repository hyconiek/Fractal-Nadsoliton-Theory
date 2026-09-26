# MEMORY-AWARE-METASTABLE-UNIT-VERDICT-137
## Updated physical verdict after reports 131–136

Date: 2026-09-26

The strengthened programme changes the status of the proposed binary
refinement substantially.

## 1. What failed

### Recursive local pitchfork

The representative main localized minimum has no detected Z2-odd soft mode
from g=3.516 to g=50.

So the simple rule

    every localized state -> two stable children

is not supported.

### Exact Markov parity bit

The exact leave-one-out process is not strongly lumpable onto parity count.

Hidden microstate information changes transition activities by O(1).

So a projected parity bit is not exactly Markov.

### Small-N metastable autonomy

Exact D12-equivariant committors for N=3..6 overwhelmingly escape into the
localized sector before completing a +k6 -> -k6 switch.

Thus local mean-field minima alone do not make a usable finite-N bit.

## 2. What survived

### Endogenous binary incidence

The uniform -> ±k6 pitchfork still gives:
- two stable mean-field daughters in a finite g window;
- a common parent saddle;
- an exact one-dimensional barrier formula.

This remains the cleanest endogenous binary-incidence prototype in FIN.

### Decaying projection memory

Although exact lumpability fails, the Mori-Zwanzig memory kernel for parity
projection decays rapidly at N=3..6.

At N=6 its normalized norm is roughly

    8.64e-2 at t=1,
    2.80e-2 at t=2,
    5.87e-3 at t=4,
    3.43e-4 at t=8,
    1.32e-6 at t=16.

The observed tail decay rate is about 0.69.

Slow hidden eigenmodes exist, but the slowest ones are symmetry-decoupled from
the parity memory source.

This is the key new memory result.

### Approximate Markov closure is plausible, not exact

For N=6 the exact projected semigroup differs from the instantaneous Markov
closure by about:
- 0.13% at t=0.1;
- 1.72% at t=1;
- 2.49% at t=2.

So memory is measurable but short-lived in the tested small-N regime.

## 3. Correct next architecture

The next theory should not demand exact self-similar binary refinement.

It should seek:

    microscopic Markov process
      ->
    metastable basin partition
      ->
    short memory / fast intrabasin relaxation
      ->
    controlled effective generator
      ->
    next coarse level.

A coarse state may have 2, 3, 4, or another number of children.

The branching number must come from capacity/time-scale separation.

## 4. Required four-time-scale ledger

Every proposed effective unit must now report:

    tau_mem
      memory-decay time;

    tau_relax
      intrabasin equilibration time;

    tau_switch
      desired transition time between coarse states;

    tau_escape
      leakage time to states outside the chosen coarse unit.

A useful approximately Markovian metastable unit needs an ordering such as

    tau_mem, tau_relax
      <<
    tau_switch
      <<
    tau_escape.

If instead tau_mem is comparable to tau_switch, the memory kernel must remain
in the effective theory.

## 5. Current k6 status

At the saddle-level optimum:

    tau_mem (small-N parity projection)
      is O(1);

    tau_relax (mean-field k6 slow mode)
      is ~118;

    the metastable exponent is only
      B_switch~1.35e-5.

So a parametrically clean switching regime would require extremely large N.

This makes k6 a useful mathematical prototype but currently a poor candidate
for the first physically economical emergent unit.

## 6. Recommended next system

The previously derived three four-state sectors are now more promising than an
imposed binary tree.

They already arise from the lowest barrier hierarchy.

The next P0 target should therefore compare:
- parity-bit coarse graining;
- 3x4 metastable-sector coarse graining,

using the SAME leave-one-out generator and the SAME tests:
capacity, memory decay, semigroup error and parameter count.

Whichever has the cleaner time-scale separation should define the next
effective level.
