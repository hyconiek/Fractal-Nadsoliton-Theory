# P62 Legacy Physical Role Transfer To Strict Gate Kernel Probe

Status: `P62_EXECUTED_LEGACY_PHYSICAL_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F4`, the repo now explicitly classifies three legacy physical roles as
still carried only by the legacy kernel package:

1. the legacy Weinberg-angle role,
2. the legacy fine-structure role,
3. the legacy gravity-hierarchy role.

`P62` asks the next honest question:

```text
does the current repo already export a rigorous transfer of any of those
legacy physical roles onto K_strict_gate?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_LEGACY_PHYSICAL_ROLE_CLAIMS_AND_STRICT_KERNEL_PIPELINE_BUT_NO_RIGOROUS_LEGACY_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_AFTER_P62
```

## Why it fails

`P62` confirms all of the following at once:

1. `F4` classifies the three legacy physical roles as legacy-side only,
2. `P47/N50` still leave the kernel bridge absent,
3. `F2` still classifies the strict gate kernel as a later-pipeline
   operational object,
4. `QW-2005` still downgrades the legacy roles rather than re-establishing
   them as strict proofs.

Therefore the current repo still lacks:

1. an explicit legacy-to-strict kernel bridge operator or theorem,
2. an explicit claim-specific role-transfer theorem for the legacy
   Weinberg-angle semantics onto `K_strict_gate`,
3. an explicit claim-specific role-transfer theorem for the legacy
   fine-structure semantics onto `K_strict_gate`,
4. an explicit claim-specific role-transfer theorem for the legacy
   gravity-hierarchy semantics onto `K_strict_gate`,
5. an explicit retained-vs-replaced partition theorem for legacy physical
   roles on the strict side.

## Real reduction after `P62`

`P62` does not say that the strict kernel is false.

It says only this stronger current-repo-state result:

```text
the repo exports legacy physical-role claims and a later strict kernel
pipeline, but it does not yet export a rigorous transfer of those
legacy physical roles onto the strict gate kernel
```

So the frontier is no longer:

```text
maybe those old physical roles came along automatically with the new kernel
```

It is now:

```text
either export claim-specific role-transfer theorems, or keep those legacy
roles explicitly outside the strict-kernel inheritance package
```

## What P62 does not claim

`P62` does not claim:

- theorem-level proof that the strict gate kernel is wrong,
- theorem-level proof that the legacy formulas are wrong,
- theorem-level proof that no future role transfer can ever exist,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. either add one explicit role-transfer theorem or bridge packet,
2. or formalize the current non-transfer result theorem-level,
3. while keeping strict operational outputs separate from legacy role
   inheritance.
