# F2 Strict Gate Kernel Provenance And FAR Input Classification Packet

Status: `F2_EXECUTED_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `K1/P47/N50`, the repo already knows that the legacy ontological kernel
and the later strict gate kernel are not rigorously identified.

After `K2`, the derivation chain of the later strict gate kernel is also
explicit.

`F2` converts that knowledge into one current-FAR input classification rule.

## Result

`F2` establishes the following current-repo-state classification:

1. the later strict gate kernel is a
   `later_pipeline_operational_gate_selected_kernel`,
2. the old `D_f / alpha_geo / beta_tors` layer remains the only currently
   restored ontological parameter layer for action-first FAR,
3. therefore `A1/A4/A8` must not silently substitute the strict gate kernel
   for the old ontological kernel or for its canonical structural parameter
   layer,
4. the strict gate kernel may enter FAR only as a later operational control,
   benchmark, or downstream strict-pipeline import unless a bridge is added.

## Classification rule

For current FAR:

```text
K_legacy_ont / (D_f, alpha_geo, beta_tors)
  = canonical-ontology-supported action-first source layer

K_strict_gate
  = later-pipeline operational strict working kernel
```

So the honest current rule is:

```text
do not silently inherit full ontological roles of K_legacy_ont into K_strict_gate
```

## A1 consequence

`A1` may use:

- the nadsoliton ontology from `AX9`,
- the canonical parameter layer from `F1`.

`A1` may **not** silently replace that ontological source layer with the later
strict gate kernel.

## A4 consequence

`A4` may use the later strict gate kernel only as:

- operational comparison,
- later-pipeline control,
- downstream consistency target.

It may **not** treat `QW-2049` as if it already supplied the ontological shell
measure or the canonical fractal substrate semantics.

## A8 consequence

`A8` may use the later strict gate kernel only as:

- operational gravity-side control or benchmark,
- later strict-pipeline input.

It may **not** silently convert that into proof that the old
`alpha_geo/(2 beta_tors)` hierarchy semantics have already been inherited by
the strict gate kernel.

## What F2 does not claim

`F2` does not claim:

- theorem-level proof that the strict gate kernel is wrong,
- theorem-level proof that the strict gate kernel can never be bridged back to
  the ontological kernel,
- full selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is:

1. keep FAR action-first stages under this kernel-input classification,
2. only later add an explicit bridge or non-bridge theorem,
3. keep light-before-observer ordering and the corrected `AX9` ontology
   explicit.
