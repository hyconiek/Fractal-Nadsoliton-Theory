# CONSERVATION-SOURCE-265
## Existing FIN naturality principles do not select composition conservation; a stronger record-content principle would

Date: 2026-09-27

Status:
exact comparison of the three admissible Z3 reset-dilation gates.

The family

    F_c(x,y)=(y+c,x+c)

contains:
- c=0: SWAP;
- c=+1;
- c=-1.

Reports 240 and 257 establish:
- all three are bijective;
- all produce the same full reset with a fresh uniform ancilla;
- F_+ and F_- are mutual inverses;
- a reversible process can mix them symmetrically;
- only SWAP protects nontrivial additive composition charges.

The question here is whether already accepted FIN principles force c=0.

## 1. Principles that do NOT distinguish c

### Global information preservation

All F_c are bijections.

No selection.

### Minimum environment dimension

All act on the same 3 x 3 pair state space.

No selection.

### Correct one-unit MaxEnt reset

All produce the same uniform subsystem output with a fresh uniform environment trit.

No selection.

### Reversibility of the stochastic process

F_+ and F_- can appear with equal rates and satisfy detailed balance jointly.

No selection.

### Exchange and Z3 covariance

The three-gate family survives these symmetries.

No selection.

Therefore current principles do NOT derive the composition conservation law.

## 2. A stronger candidate principle: transport records without rewriting their content

Treat the pair as two information records carrying values x and y.

Allow their LOCATIONS to be exchanged, but penalize rewriting their record values.

Define content-mutation cost as:
the minimum number of value edits needed after optimally matching output records to input records.

Exact enumeration gives:

    SWAP:
      average cost = 0
      worst cost = 0

    F_+:
      average cost = 4/3
      worst cost = 2

    F_-:
      average cost = 4/3
      worst cost = 2.

Equivalently:

    boxed:
    SWAP uniquely preserves the unordered input value multiset {x,y} for every pair.

Thus the principle

    "elementary transport may move records but not rewrite their contents"

selects SWAP uniquely.

It simultaneously implies conservation of every empirical label count.

## 3. Epistemic boundary

This is a good SOURCE CANDIDATE, not yet an accepted FIN law.

Abstract bijectivity permits content recoding.
To forbid F_± one needs a persistent notion of record content/identity that is stronger than Shannon information preservation.

The next proof target is therefore:

    derive stable record identity from FIN's relational continuity,
    or state composition conservation explicitly as an additional microscopic law.

## Verdict

The conservation problem is now sharply localized.

No tuning parameter is needed.

One discrete law distinguishes the two universality classes:

    content-preserving transport
      ->
    SWAP
      ->
    conserved densities
      ->
    gapless diffusion;

    content-rewriting reversible gates
      ->
    reaction gap
      ->
    no protected hydrodynamic density mode.
