# Architect decision after MT

Date: 2026-10-15.
Decision: **INTRINSIC_OBSERVATION_AND_CLASSICAL_RECORDS**.
Status: selected research programme; observation tasks have not been executed.

The local critical regime is now sufficiently controlled for an effective-
theory next layer. Do not make optimal polynomial-gap constants or another
trapping-exponent refinement the default prerequisite.

## Next scientific question

> What observable information about the collective mode can a small existing
> subset, or a history of an existing microscopic copy, actually contain?

The point is to derive a readout channel/noise from the **same microscopic
law**, not to insert a detector, observer, storage rule or extra stochastic
field. Copy tags are observation choices, not physical cells or space.

## Two different candidate bridges

At a time within the MT window the conditional labelled law is exchangeable.
Given total parity count M, a fixed k-copy block has exactly a hypergeometric
parity count. Together with the critical global fluctuation scale this suggests
an instantaneous readout threshold k~sqrt(N). At k/sqrt(N)->alpha>0, the target
scaled channel is `Y=Z+alpha^(-1/2)G`. This is a proposed observation theorem,
not a new result already supplied by MT.

An instantaneous noisy record may have no continuous critical-time path law:
its private fluctuations can decorrelate in O(1) model time while the signal
evolves on sqrt(N). Testing that obstruction comes before writing a measured
signal SDE by analogy.

For a **fixed** k-copy block, use a different observable: the existing path
functional

    Y_N(s)=N^(-1/4) integral_(t_N)^(t_N+sqrt(N)s) block_parity_mean(t) dt.

For primary HB, a promising target is

    dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB,
    dY=Z ds+sqrt(2/k)dW,

with independently generated limiting Brownian noises. The new noise would
come from tagged-copy refresh martingales; it is not supplied by hand. The
integral is a mathematical observation, not a derived physical memory device.

## The key premise and proof gate

Counts alone do not determine tagged-copy histories. First bind the actual
labelled single-copy refresh generator and forbid silent count-invisible
permutations. A diagnostic permutation process can show the ambiguity without
being admitted as a new FIN law.

For HB, the candidate proof has a concrete starting identity: for fixed block
contrast `v=block_mean-global_mean`, the labelled generator gives
`L v=-v+O(1/N)`. The integrated martingale identity pays endpoint/remainder
errors. Joint bracket calculations, not assertions of independent noise,
must identify the effective channel. MET/BAR are conditional sensitivity
checks; their tagged-generator errors cannot simply be bounded by copying
the HB identity.

## Why this is a physics-facing step

It connects microscopic states -> collective law -> observable records with
explicit statistical errors and information limits. It can distinguish
recoverability from information conservation and a static noisy readout
from a continuous observation process. This adds operational content without
adding a microscopic physical postulate.

It still does not derive a physical apparatus, units, spatial location,
measurement backreaction, natural kernel/initialization or QM. A shared
classical latent variable is not entanglement, Born rule or emergent spacetime.

Execute [the OBS-C01–C05 English plan](../FIN_INTRINSIC_OBSERVATION_NEXT_CHAIN_20261015.md),
then stop for the architect. No experiment or further automatic campaign.
