# F270 First Actual Strict Sigma-Int Driven `E_pair` Generator Candidate Packet

Status: `F270_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_DRIVEN_E_PAIR_GENERATOR_CANDIDATE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `N381`, the strict-side fractal-light lane has:

1. an explicit reduction form `E_pair -> (theta_1^cand, theta_2^cand)` (`T115`),
2. one explicit template carrier instance `E_pair_nad12_uniform_template_v1`
   (`T116/F269`),

but it still lacks any explicit **source-side** (internal-datum driven)
generator for `E_pair`.

`F270` packages one explicit **candidate** generator:

```text
G_sigma_int_to_E_pair_generator_candidate_v1
```

that maps an internal `Z2` sigma-int datum value `sigma_int_input ∈ {+1,-1}`
into a finite, normalized, pair-indexed carrier field `E_pair`,
without using `theta` inputs and without using populated-instance inputs.

## Inputs reused

1. `B4`
   - `sigma_int_candidate := chi_FR(gamma_pi1) ∈ {+1,-1}` (candidate object; hybrid-supported).
2. `F307/N418`
   - `sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1) ∈ {+1,-1}` (strict-side source upgrade; explicit premise provenance).
3. `T117`
   - explicit admissible generator form for `E_pair` typed over `sigma_int_input ∈ {+1,-1}`.

## Packet result

`F270` exports one actual packaged generator candidate:

```text
G_sigma_int_to_E_pair_generator_candidate_v1
```

defined only as:

```text
G_sigma_int_to_E_pair_generator_candidate_v1 :=
(
  input = sigma_int_input ∈ {+1,-1},
  parameter_eps = eps ∈ [0,1],
  sign_masks = {
    b_{1,k} = (-1)^k,
    b_{2,k} = (-1)^(k+1)
  } for k=0..11,
  output = E_pair with:
    d_{i,k}=k,
    w_{i,k}=(1 + sigma_int_input*eps*b_{i,k})/12,
  contracts = [
    noncyclic_no_theta_input,
    noncyclic_no_populated_instance_input,
    observer_free_no_K_obs_primary_selection
  ],
  strict_derivation_present = false,
  strict_uniqueness_present = false
)
```

## Meaning

This packet means only:

1. the repo now contains one explicit internal-datum-driven candidate rule for
   generating a finite `E_pair` carrier field from a `Z2` sigma-int input value,
2. this is stronger than having only a fixed uniform template instance,
3. but it remains strictly below:
   - strict derivation of the generator,
   - strict bridge/export-map discharge,
   - strict-core theta export,
   - ToE closure.

No identification theorem between `sigma_int_candidate` (`B4`) and
`sigma_int_strict_derived_v1` (`F307/N418`) is used or implied here; both are
simply admissible instantiations of the abstract input
`sigma_int_input ∈ {+1,-1}`.

## What F270 does not claim

`F270` does not claim:

1. discharge of `T2`,
2. resolution of `N302`,
3. any strict-core selector closure or `QW-2191` discharge,
4. actual `theta_1`, `theta_2`,
5. actual pair population,
6. ToE closure.
