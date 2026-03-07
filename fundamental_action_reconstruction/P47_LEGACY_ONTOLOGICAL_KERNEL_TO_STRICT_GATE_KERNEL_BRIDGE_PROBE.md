# P47 Legacy Ontological Kernel To Strict Gate Kernel Bridge Probe

Status: `P47_EXECUTED_LEGACY_TO_STRICT_KERNEL_BRIDGE_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `K1`, the repo now explicitly names a kernel split:

```text
K_legacy_ont(d) := alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)
K_strict_gate(d) := cos(omega*d+phi)/(1+beta*d^eta)
```

`P47` asks the next honest question:

```text
does the current repo already export a rigorous bridge identifying
K_legacy_ont(d) with K_strict_gate(d) or with a strict renormalized image of it?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_LEGACY_AND_STRICT_KERNELS_BUT_NO_RIGOROUS_LEGACY_TO_STRICT_KERNEL_BRIDGE_AFTER_P47
```

## Why it fails

`P47` confirms all of the following at once:

1. the legacy kernel is explicitly present in canonical legacy sources,
2. the strict gate kernel is explicitly present in Release `4.9` / `QW-2049`
   sources,
3. `F1` restores the old ontological parameter layer
   `D_f / alpha_geo / beta_tors`,
4. `QW-2005` already downgrades several old legacy formulas from strict proof
   to heuristic/model-level status,
5. but the current repo still does **not** export a rigorous bridge object
   between the two kernels.

In particular, the current repo still lacks:

1. an explicit amplitude normalization or absorption map explaining the loss
   of visible `alpha_geo`,
2. an explicit renormalized damping map
   `beta_tors -> (beta, eta)` or an equivalent strict translation rule,
3. an explicit phase/frequency bridge from
   `(pi/4, pi/6)` to `(0.18575, 0.16250)`,
4. an explicit kernel-specific retained-vs-replaced partition theorem,
5. an explicit kernel-level ontological anchoring witness showing that the
   strict gate kernel is an internal route of one informational nadsoliton,
   rather than a law of a separate informational layer.

## Real reduction after `P47`

`P47` does not say that the new strict kernel is false.

It says only this stronger current-repo-state result:

```text
the repo exports two kernel objects and a partial revalidation downgrade,
but it does not yet export a rigorous bridge identifying them
```

So the frontier is no longer:

```text
maybe the new kernel replaced the old one somehow
```

It is now:

```text
either export a bridge, or keep the kernel split explicit
```

## What `P47` does not claim

`P47` does not claim:

- theorem-level proof that the strict gate kernel is wrong,
- theorem-level proof that the legacy kernel is wrong,
- theorem-level proof that no bridge can ever exist,
- full ToE closure,
- discharge of `QW-2191`,
- closure of the selector route,
- automatic restoration of old `alpha_geo/beta_tors` formulas on the strict
  side.

## Recommended next move

The correct next move is now:

1. either build one explicit bridge packet
   `legacy kernel -> strict gate kernel`,
2. or formalize the current nonidentification result theorem-level for the
   current repo state,
3. while keeping the ontological correction from `AX9/K1` explicit.
