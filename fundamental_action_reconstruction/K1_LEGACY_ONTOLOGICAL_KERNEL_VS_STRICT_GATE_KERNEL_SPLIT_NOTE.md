# K1 Legacy Ontological Kernel vs Strict Gate Kernel Split Note

Status: `K1_CURRENT_REPO_STATE_KERNEL_SPLIT_NOTE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

Make explicit one structural ambiguity that is currently easy to miss in the
repo:

```text
the symbol K(d) is used for two materially different kernel objects
belonging to two different layers of the FIN program
```

This note does **not** claim a bridge theorem.

It records the strongest honest current reading:

1. the legacy FIN program contains an ontological/effective kernel,
2. the later strict pipeline contains a gate-selected strict working kernel,
3. the repo does not yet export a rigorous derivation identifying one with the
   other.

## Kernel object 1: legacy ontological/effective kernel

Canonical legacy sources:

- [DIAGRAMS_KERNEL_TRANSFORMATION.md](../DIAGRAMS_KERNEL_TRANSFORMATION.md)
- [TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex](../TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex)
- [RAPORT_QW1729_NADSOLITON_KERNEL_CHARACTERISTICS_MAP.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW1729_NADSOLITON_KERNEL_CHARACTERISTICS_MAP.md)

The legacy ontological/effective kernel is:

```text
K_legacy_ont(d) := alpha_geo * cos(omega*d + phi) / (1 + beta_tors*d)
```

with the legacy canonical interpretation:

- `alpha_geo = 4 ln 2` as info-geometry amplitude,
- `omega = pi/4`,
- `phi = pi/6`,
- `beta_tors = 0.01` as hierarchy/torsion damping datum.

This object is not just an arbitrary fit in the legacy program.

Inside the legacy program it is presented as:

1. an effective compression of the more detailed legacy `K_total` route,
2. a carrier of ontological structure tied to the informational nadsoliton,
3. a source used downstream in legacy formulas for hierarchy, EW angle, and
   related model-level identities.

## Kernel object 2: strict gate working kernel

Canonical strict sources:

- [RELEASE_4_9_TEXTBOOK_EN_PL.md](../RELEASE_4_9_TEXTBOOK_EN_PL.md)
- [RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md)

The strict gate working kernel is:

```text
K_strict_gate(d) := cos(omega*d + phi) / (1 + beta*d^eta)
```

with the selected `QW-2049` working tuple:

- `omega = 0.18575`,
- `phi = 0.16250`,
- `beta = 1.0`,
- `eta = 1.8`.

This object is not documented as a simple restatement of the legacy kernel.

Inside the current strict program it is presented as:

1. the selected common structural function for the Release `4.9` strict gate,
2. a kernel supported by the repaired micro/Stage-C intersection route,
3. a route later tied to micro-supported renormalization constants
   `Z_beta`, `delta_eta`.

## Why this is a real split

The split is structural, not cosmetic:

1. `alpha_geo` is explicit in `K_legacy_ont` and absent from
   `K_strict_gate`,
2. damping changes from linear `d` to nonlinear `d^eta`,
3. the legacy tuple `(pi/4, pi/6, 0.01)` is not numerically close to the
   strict tuple `(0.18575, 0.16250, 1.0, 1.8)`,
4. the legacy route reads the kernel as ontological/effective structure,
5. the strict route reads the kernel as selected gate-level structural
   function after methodological repair and micro support.

Therefore the current strongest honest statement is:

```text
K_legacy_ont(d) and K_strict_gate(d) are not yet rigorously identified in the repo
```

## What current repo evidence already says about this split

The repo already contains partial revalidation language supporting this split.

In
[RAPORT_QW2005_PRE1700_TEX_REVALIDATION_MATRIX_EN_PL.md](../material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2005_PRE1700_TEX_REVALIDATION_MATRIX_EN_PL.md):

- the old frozen kernel `(omega=pi/4, phi=pi/6, beta_tors=0.01)` is marked
  `PARTIAL / CONFLICTED`,
- `sin^2(theta_W)=alpha_geo/12` is marked
  `HEURISTIC / NOT STRICTLY DERIVED`,
- `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)` is marked
  `PARTIAL / MODEL-LEVEL`,
- exact gravity hierarchy from `beta^N` is marked
  `MODEL-CONSISTENT, NOT FULL INDEPENDENT PROOF`.

So the repo already does **not** treat the entire old legacy kernel edifice as
current strict proof.

## Relation to the ontological correction from AX9

`AX9` restored the canonical ontology:

```text
the nadsoliton is itself the primordial information of the universe in a
solitonic state; there is no separate informational layer underneath it
```

That correction changes how the legacy kernel should be read.

After `AX9`, the honest reading is:

1. `K_legacy_ont(d)` should **not** be described as a coupling law sitting on
   top of an independent informational substrate,
2. instead, `K_legacy_ont(d)` is best read as an internal effective coupling
   pattern of one informational nadsoliton,
3. this strengthens the obligation to keep `D_f`, `alpha_geo`, and
   `beta_tors` visible in `A1/A4/A8`,
4. this also preserves the route discipline
   `nadsoliton -> light -> matter -> emergent observer`,
   with observer downstream.

However, `AX9` does **not** eliminate the kernel split.

It fixes ontology, but it does not produce a theorem of the form:

```text
K_legacy_ont(d) == K_strict_gate(d)
```

## Ontological misreading risk on the strict side

There is a specific documentation risk on the strict side that should now be
named explicitly.

In the later strict presentation, the kernel is introduced as a common
structural function of a deeper informational-dynamical structure, but that
presentation does not always re-state the corrected `AX9` ontology with enough
force.

Because of that, the strict-side documentation can be misread as if it assumed
the following background picture:

```text
informational substrate/layer -> nadsoliton or sectors -> kernel law
```

rather than the corrected canonical reading:

```text
informational nadsoliton -> internal light/matter/observer projections -> local kernels
```

This is a real risk of interpretation because:

1. the later strict kernel is presented operationally as a selected working
   gate function,
2. its explicit `alpha_geo / beta_tors` ontological markers disappear,
3. the strict text can then look as if the kernel belongs to an underlying
   informational medium rather than to one nadsoliton interpreted
   informationally.

After `AX9`, the safest correction is:

1. `K_strict_gate(d)` should also be read, whenever possible, as a later
   working-gate projection of one informational nadsoliton route,
2. not as evidence for a separate informational layer outside the nadsoliton,
3. and not as a theorem that the strict kernel has already inherited all
   ontological roles of the legacy kernel.

So the current honest repo state is:

```text
the strict kernel is operationally usable,
but ontologically under-anchored unless it is read through AX9/K1 discipline
```

## Consequence for FAR

For the current `fundamental_action_reconstruction` lane, the strongest honest
reading is:

1. `A1/A4/A8` may and should track the legacy ontological parameter layer
   (`D_f`, `alpha_geo`, `beta_tors`) as
   `canonical-ontology-supported`,
2. the strict gate kernel from `QW-2049` may still be used on the later strict
   working pipeline,
3. but the repo should not silently pretend that these are already one and the
   same kernel object.

## Open bridge frontier

The missing bridge is now sharp.

The repo still lacks an explicit packet/theorem establishing some combination
of:

1. amplitude normalization or absorption map explaining the loss of explicit
   `alpha_geo`,
2. renormalized damping map
   `beta_tors -> (beta, eta)` or an equivalent strict translation rule,
3. phase/frequency bridge
   `(pi/4, pi/6) -> (0.18575, 0.16250)`,
4. a theorem stating which legacy formulas are retained only as
   ontology-supported/model-level relations and which are replaced by the
   strict gate pipeline.

Until such a bridge exists, the safest notation is:

```text
K_legacy_ont(d)   !=?   K_strict_gate(d)
```

with the equality question still open.

## What K1 does not claim

`K1` does not claim:

- that the legacy kernel is invalid,
- that the strict gate kernel is invalid,
- that `AX9` solves the split,
- that the old `alpha_geo / beta_tors` formulas are strict-derived,
- that the strict gate kernel fully supersedes every ontological role of the
  legacy kernel,
- that ToE is closed.

## Recommended next move

The next honest move is to add one explicit bridge note or theorem-spec:

```text
legacy ontological kernel -> strict gate kernel bridge or non-bridge result
```

so the repo no longer relies on implicit layer-switching under one shared
symbol `K(d)`.
