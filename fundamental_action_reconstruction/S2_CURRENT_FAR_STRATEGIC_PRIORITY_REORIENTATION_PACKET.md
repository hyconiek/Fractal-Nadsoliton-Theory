# S2 Current FAR Strategic Priority Reorientation Packet

Status: `S2_UPDATED_LEGACY_KERNEL_RESTORED_AS_INTERMEDIATE_BRIDGE_KERNEL_NO_FALSE_PASS`
As of: `2026-05-31`

## Goal

Record the current strategic reorientation forced by the author's latest
correction: the legacy kernel is **not** a discarded dead end.  It is restored
as an intermediate kernel on the path toward identifying the strict kernel.

This packet supersedes the older retirement-only reading.  It does not license
silent identity, silent physical-role transfer, QW-2191 discharge, or ToE
closure.

## Critical update: legacy kernel restored as intermediate bridge kernel

The current bridge discipline is:

```text
K_legacy_ont(d) = alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)
    -> intermediate / incomplete bridge kernel

K_strict_gate(d) = cos(omega*d+phi)/(1+beta*d^eta)
    -> completed / enriched strict kernel
```

The intended bridge direction is therefore not “retire legacy and forget it”,
but:

```text
legacy kernel + missing strict-side characteristics -> strict kernel
```

The current repository already contains finite certificates showing that the
strict completion adds structure that the legacy kernel did not fully encode,
including:

1. the unique A/P/D completion ansatz on the audited finite domain,
2. GF(2) phase-sign reconstruction and rank-11 uniqueness,
3. carrier/edge/node commutative-diagram bookkeeping,
4. path-cohomology and cycle-closure parity constraints,
5. strict damping/compression behavior, including the nonlinear `d^eta`
   compression absent from the legacy linear `beta_tors*d` denominator.

Therefore the legacy kernel is restored as an input/intermediate object, but it
remains incomplete relative to the strict kernel.

## Priority 1: legacy -> strict completion bridge

The highest theoretical priority is now:

```text
construct and audit the explicit completion bridge
K_legacy_ont -> K_strict_gate
```

A valid bridge must specify at least:

1. the amplitude/normalization passage,
2. the phase/frequency/topological bit passage,
3. the damping/compression passage from legacy linear torsion damping to strict
   nonlinear compression,
4. the selector/source premise responsible for the certified phase/topological
   data,
5. which strict-side characteristics are additions rather than legacy contents.

## Priority 2: role-transfer audit after bridge completion

When the bridge is fully specified, the next mandatory step is a separate
legacy-role transfer audit.

That audit must decide, claim by claim, whether each legacy physical role:

1. survives unchanged under strict completion,
2. survives only in modified/compressed form,
3. becomes a strict-side successor statement with different semantics, or
4. is rejected.

This especially applies to legacy claims involving:

- `sin^2(theta_W)=alpha_geo/12`,
- `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)`,
- `beta^N` gravity hierarchy,
- any proposed `beta_tors -> chi_11` orientation/torsion role.

## Priority 3: selector / source obstruction

`QW-2191` remains a real strict-core selector obstruction.  The GF(2), path
cohomology, and cycle-closure certificates identify where the bit lives and how
it is constrained, but they do not by themselves export its strict source.

Any proposed source such as `beta_tors -> chi_11`, spontaneous symmetry
breaking, or observer-readout must be stated as a bridge/source theorem and must
not be smuggled in as an already proven strict-core selector closure.

## Consequence for current FAR work

The current priority order is now:

1. explicit `legacy -> strict` completion bridge,
2. role-transfer audit after the bridge is fully specified,
3. selector/source theorem for the remaining bit/orientation obstruction,
4. only then auxiliary local direct-route mass decomposition.

## What S2 does not claim

`S2` does not claim:

- that `K_legacy_ont == K_strict_gate` as a raw identity,
- that the legacy kernel already contained all strict nadsoliton
  characteristics,
- that legacy physical roles automatically transfer to the strict kernel,
- that `beta_tors -> chi_11` is already proven,
- that `QW-2191` is discharged,
- that ToE is closed.

But it DOES claim:

- that the legacy kernel is restored as an intermediate bridge kernel,
- that the strict kernel should be audited as a completed/enriched legacy
  continuation rather than an unrelated replacement,
- that bridge completion must be followed by explicit role-transfer auditing,
- that compression and other strict-side additions must remain visible.

## Product

This packet is a guardrail update, not a closure theorem.
