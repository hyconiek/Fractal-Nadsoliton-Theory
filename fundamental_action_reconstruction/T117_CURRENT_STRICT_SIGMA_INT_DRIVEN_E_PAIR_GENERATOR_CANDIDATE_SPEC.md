# T117 Current Strict Sigma-Int Driven `E_pair` Generator Candidate Spec

Status: `T117_CURRENT_STRICT_SIGMA_INT_DRIVEN_E_PAIR_GENERATOR_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `T115/F268/N380` the strict-side lane exports one explicit
fractal-light pair-map-rule **candidate form** which reduces a finite
path-ensemble carrier field `E_pair` to pair-indexed candidate phases.

After `T116/F269/N381` the repo also exports one explicit **template**
carrier-field instance `E_pair_nad12_uniform_template_v1`.

The next honest gap is therefore no longer “invent a map form”.
It is:

```text
export one explicit source-side generator (still candidate-level)
for producing E_pair from an internal strict-core datum
without theta inputs and without populated-instance inputs
```

`T117` proposes one admissible **candidate** generator that uses:

1. an internal `Z2` sigma-int datum value `sigma_int_input ∈ {+1,-1}`,
2. a finite `12`-slot recurrence scaffold (nad12 motif),

to generate a finite, normalized, pair-indexed carrier field `E_pair`.

This spec does **not** assume any identification between candidate and strict
sigma-int sources. The same `Z2` input can be instantiated by:

1. `sigma_int_candidate` (`B4`) (candidate object; hybrid-supported), or
2. `sigma_int_strict_derived_v1` (`F307/N418`) (strict-side source upgrade; premise provenance).

No equality theorem between these two objects is used or implied here.

This is intentionally weaker than any claim of:

- strict derivation of `E_pair`,
- bridge/export-map discharge,
- theta export,
- pair population,
- selector closure,
- ToE closure.

## Candidate generator form

### Inputs

1. `sigma_int_input ∈ {+1,-1}` (internal `Z2` sigma-int datum value),
2. a fixed amplitude parameter `eps ∈ [0,1]` (candidate),
3. a fixed pair-indexed sign mask `b_{i,k} ∈ {+1,-1}` for:
   - pair slot `i ∈ {1,2}`,
   - octave/path index `k ∈ {0,1,...,11}`.

Notation: in the formulas below, `\sigma_{int}^{in}` denotes the chosen input
value `sigma_int_input`.

### Output

The generator outputs a pair-indexed path-ensemble carrier field:

```text
E_pair := { pair1: {paths: [...]}, pair2: {paths: [...]} }
```

with the `T116` contracts:

1. distances `d >= 0`,
2. weights `w >= 0`,
3. `sum_k w_k = 1` for each pair slot,
4. noncyclic: no `theta` input and no populated-instance input,
5. observer-free: no `K_obs` primary selection.

### Concrete admissible mask choice (minimal `Z2` grading on the 12-ring)

One admissible fixed sign mask is:

```math
b_{1,k} := (-1)^k,
\qquad
b_{2,k} := (-1)^{k+1}.
```

This mask is:

1. pair-indexed (pair2 is a fixed shifted mask, not an observer choice),
2. finite and fully explicit,
3. `Z2`-valued.

### Explicit weight law

For each pair slot `i ∈ {1,2}` and each `k ∈ {0,...,11}` define:

```math
d_{i,k} := k,
\qquad
w_{i,k} := \frac{1 + \sigma_{int}^{in}\, \varepsilon\, b_{i,k}}{12}.
```

Then define:

```text
E_i(pair) := { (d_{i,k}, w_{i,k}) }_{k=0..11}.
```

Nonnegativity is automatic for `eps ∈ [0,1]`.
Normalization holds by construction.

## Intended use (composition only; no upgrade)

The only intended use of this generator is to supply a noncyclic input
`E_pair` to the already exported reduction form `T115`, i.e.:

```text
sigma_int_input
  -> E_pair (T117 candidate generator)
  -> (theta_1^cand, theta_2^cand) (T115 candidate reduction)
  -> S_orient^cand (R1 target slot class)
```

This composition is still candidate-level and must not be promoted into an
actual strict-core theta export or into a discharge of `T2`.

## Hard limits

`T117` must not claim:

1. strict derivation or uniqueness of the sign mask `b_{i,k}`,
2. strict derivation or uniqueness of the parameter `eps`,
3. strict derivation of the resulting `E_pair`,
4. that the composition discharges `N302` / `T2` / `QW-2191`,
5. actual `theta_1`, `theta_2`,
6. actual pair population,
7. ToE closure.
