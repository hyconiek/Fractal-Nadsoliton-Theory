# PHYS-011 — Q12-ARCHITECTURE-AND-N2-MAPPING

Date: 2026-09-29

Status: **EXACT MODEL MAPPING + SOURCE-CODE ARCHITECTURE AUDIT**

## Result

The VRSPAD-144 RTL/host architecture is structurally capable of two 12-state
Potts nodes without changing the FPGA topology.

The key hardware facts are:

- the latch chain contains 144 state latches;
- `cuts` define the boundaries between Potts nodes;
- within a region with no cuts, pulses mutually inhibit the other latches;
- the host model accepts arbitrary `node_sizes`;
- the stored weight matrix is 144 x 144;
- each stored weight byte is interpreted as signed two's-complement by the
  energy pipeline, so positive and negative pair couplings are available.

The 16-channel demonstrator is **not q=12**: its published Potts conversion
uses four q=4 nodes. It is therefore only a stochastic-rate/system-fidelity
reference, not proof of a q=12 run.

## Exact N=2 FIN mapping

For two labelled FIN copies at N=2, the diagonal A7 terms are constant and

    P(i,j) ∝ exp[(g/2) A7[i,j]].

Because A7 is circulant, the marginal distribution of

    d = (j-i) mod 12

is

    P(d) ∝ exp[(g/2) A7[0,d]],

which is exactly the frozen PHYS-007 12-bin conditional fingerprint.

Therefore the stronger hardware realization is:

- node 1: q=12;
- node 2: q=12;
- 24 physical stochastic channels total;
- pair energy E(i,j)=-(g/2)A7[i,j];
- measure only d=(j-i) mod 12.

No anchor, postselection or second observable is required.

## Why this is preferable to the single-node shortcut

A single q=12 node with twelve programmed biases can reproduce the same
histogram, but it bypasses the physical pair coupling. The two-node construction
actually uses the A7 pair interaction while preserving the already-frozen
observable and countermodel predictions.

## Source boundary

This proves architectural representability and exact mathematical mapping.
It does not show that the physical device actually attains the target
distribution at TV<=0.003.
