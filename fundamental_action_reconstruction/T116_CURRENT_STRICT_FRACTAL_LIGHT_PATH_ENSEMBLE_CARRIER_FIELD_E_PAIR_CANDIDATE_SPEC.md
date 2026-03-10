# T116 Current Strict Fractal-Light Path-Ensemble Carrier Field `E_pair` Candidate Spec

Status: `T116_CURRENT_STRICT_FRACTAL_LIGHT_PATH_ENSEMBLE_CARRIER_FIELD_E_PAIR_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

After `T115/F268/N380`, the repo exports one strict-side **pair-map-rule
candidate form** that reduces a finite path ensemble carrier into
`(theta_1^cand, theta_2^cand)` via the strict kernel coupling channel.

However, that candidate form is only meaningful if the lane also supplies one
explicit **pair-indexed path-ensemble carrier field**:

```text
E_pair = (E_pair1, E_pair2)
```

where each `E_pairi` is a finite set of nonnegative distances with
nonnegative weights.

`T116` specifies what an admissible **candidate** carrier field `E_pair` is,
and what a minimal template instance looks like, while remaining:

1. noncyclic (`N18`: no `theta` input, no populated-instance input),
2. observer-free (`AX9` discipline: no `K_obs` primary source),
3. below: actual theta export, population, `E_orient`, `S_sel_int`,
   selector closure, and ToE closure.

## Carrier-field type

On the designated minimal pair family:

```text
pair_family := [pair1, pair2]
```

an admissible carrier field is a pair-indexed object:

```text
E_pair := {
  pair1: {paths: [ {d,w}, ... ]},
  pair2: {paths: [ {d,w}, ... ]}
}
```

with constraints (candidate-level contracts):

1. `d >= 0`,
2. `w >= 0`,
3. `sum_k w_k = 1` for each pair slot,
4. supplied without using `theta_1, theta_2` as input,
5. supplied without using a populated basis-pair instance as input,
6. no observer-indexed selection as primary source.

## Minimal template instance (nad12 uniform)

One admissible **template** (not strict-derived) carrier instance is:

```text
E_pair_nad12_uniform_template_v1 :=
  for each pair slot in [pair1,pair2]:
    paths = [ (d=k, w=1/12) for k in {0,1,...,11} ].
```

This uses the existing 12-slot recurrence motif as a finite carrier scaffold
only; it does not claim any strict derivation of the octave labels or of the
weights.

## Hard limits

`T116` must not claim:

1. strict derivation of `E_pair`,
2. that any one template is physically/ontologically canonical,
3. that this resolves `N302`,
4. that this discharges `QW-2191`,
5. that this upgrades `T115` into actual theta export or pair population,
6. ToE closure.

