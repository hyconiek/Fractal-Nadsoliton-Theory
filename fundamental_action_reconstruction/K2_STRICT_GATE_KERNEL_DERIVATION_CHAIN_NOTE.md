# K2 Strict Gate Kernel Derivation Chain Note

Status: `K2_CURRENT_REPO_STATE_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

Record, in one place, how the later strict gate kernel was actually obtained in
the current repo lineage.

This note exists to prevent one recurring misreading:

```text
the later strict gate kernel was not directly derived from the old 4.4
ontological kernel by a finished bridge theorem
```

Instead, the current repo shows an operational chain:

```text
refreeze scan -> refrozen baseline -> micro support repair ->
micro/Stage-C intersection -> freeze bundle -> renormalization interpretation
```

## Step-by-step chain

### 1. QW-2038: derivation-compatible refreeze scan

[RAPORT_QW2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.md)

What it does:

1. scans a family of later-form kernels,
2. selects the best row under mass/flavor/GW gate constraints,
3. outputs a candidate tuple.

Best row:

```text
omega = 0.18575
phi   = 0.16250
beta  = 1.00000
eta   = 1.80000
```

This is a candidate-selection result, not yet a legacy-to-strict bridge.

### 2. QW-2039: derivation-compatible refrozen kernel gate

[RAPORT_QW2039_DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2039_DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE.md)

What it does:

1. promotes a refrozen baseline,
2. keeps `beta` inside derivational CI95,
3. freezes the baseline for downstream testing.

Selected refrozen kernel:

```text
omega = 0.18575
phi   = 0.16250
beta  = 0.92000
eta   = 1.80000
```

So even here the final `beta = 1.0` working gate has not yet appeared.

### 3. QW-2041: canonical vs refrozen reparameterization audit

[RAPORT_QW2041_CANONICAL_REFROZEN_REPARAMETERIZATION_AUDIT.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2041_CANONICAL_REFROZEN_REPARAMETERIZATION_AUDIT.md)

This is the key semantic warning.

It compares:

```text
canonical TeX:  omega = pi/4, phi = pi/6, beta_tors = 0.01, eta = 1
refrozen 2039: omega = 0.18575, phi = 0.16250, beta = 0.92, eta = 1.8
```

Verdict:

```text
CANONICAL_REFROZEN_REPARAMETERIZATION_FAIL
CANONICAL_SEMANTIC_DRIFT_CONFIRMED_INTERNAL
```

Required next step:

```text
DEFINE_EXPLICIT_BRIDGE_OPERATOR_BETWEEN_CANONICAL_AND_EFFECTIVE_SEMANTICS
```

So the repo itself already says the bridge is missing.

### 4. QW-2045: phase-conditioned pointwise micro derivation

[RAPORT_QW2045_PHASE_CONDITIONED_POINTWISE_MICRO_DERIVATION.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2045_PHASE_CONDITIONED_POINTWISE_MICRO_DERIVATION.md)

This stage targets the refrozen `QW-2039` kernel and gives partial micro
support, but still fails strict pointwise identifiability:

- too few bins,
- phase condition too weak.

### 5. QW-2048: spectral phase-locked pointwise derivation

[RAPORT_QW2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.md)

This is the methodological repair step.

The old pointwise phase estimate is replaced by a phase-lock from the
signed-dynamic operator spectrum.

Result:

- pointwise identifiability is repaired,
- the target is still the refrozen `QW-2039` tuple
  with `beta_target = 0.92`, `eta_target = 1.8`.

### 6. QW-2046 and QW-2049: micro / Stage-C intersection

- [RAPORT_QW2046_MICRO_STAGEC_INTERSECTION_GATE.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2046_MICRO_STAGEC_INTERSECTION_GATE.md)
- [RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md)

This is where the later strict working tuple is selected.

Selected kernel:

```text
omega = 0.18575
phi   = 0.16250
beta  = 1.00000
eta   = 1.80000
```

This is an intersection-gate result:

```text
micro support + Stage-C pass + external hard flags
```

It is not yet a canonical-ontology bridge theorem.

### 7. QW-2050: freeze bundle

[RAPORT_QW2050_SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2050_SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE.md)

This freezes the selected strict working kernel for independent execution.

### 8. QW-2064: micro-derived renormalization constants

[RAPORT_QW2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.md)

This is a later interpretation layer:

- `z_beta_target = 100`,
- `delta_eta_target = 0.8`,
- micro medians are compared against those targets.

So `QW-2064` does not generate the strict kernel from scratch.
It interprets part of the already frozen strict kernel in renormalization
language.

## Bottom line

The strict gate kernel currently comes from:

```text
operational refreeze and gate selection chain
```

not from:

```text
completed theorem-level derivation from the old 4.4 ontological kernel
```

## Consequence for FAR

This means:

1. the strict gate kernel is operationally usable,
2. but its provenance is later-pipeline and gate-operational,
3. therefore it cannot be silently substituted for the old ontological kernel
   inside action-first ontological reconstruction,
4. unless an explicit bridge is added.

## What K2 does not claim

`K2` does not claim:

- that the strict gate kernel is false,
- that the refreeze chain is invalid,
- that the legacy kernel is false,
- that no bridge can ever exist,
- that ToE is closed.
