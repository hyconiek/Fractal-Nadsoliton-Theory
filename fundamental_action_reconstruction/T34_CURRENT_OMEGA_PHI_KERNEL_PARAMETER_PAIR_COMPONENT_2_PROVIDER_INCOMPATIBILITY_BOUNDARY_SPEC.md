# T34 Current Omega-Phi Kernel Parameter Pair Component 2 Provider Incompatibility Boundary Spec

Status: `T34_CURRENT_OMEGA_PHI_KERNEL_PARAMETER_PAIR_COMPONENT_2_PROVIDER_INCOMPATIBILITY_BOUNDARY_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N297`, one concrete third-provider-class target already exists on the
residual-datum / `sigma_int_candidate` lane.

An alternative idea now proposes one genuinely new blocker-cut:

```text
treat the internal strict-kernel parameter pair (omega, phi)
as the pair-indexed minimal carrier for component 2
```

The question is not whether `omega` and `phi` exist.

They do exist.

The narrower question is:

```text
can the current repo honestly treat (omega, phi)
as a pair-indexed population anchor
for T26 component 2
```

without false pass?

## Scope

`T34` is scoped only to the proposal:

```text
(omega, phi)
  -> pair-indexed population anchor for component 2 ?
```

using only current repo material.

It does **not** decide:

1. whether `omega` and `phi` are useful at all,
2. whether a future typed map from `(omega,phi)` to `theta_1, theta_2` could
   ever be exported,
3. impossibility in principle of all kernel-parameter-based routes,
4. impossibility in principle of component 2.

## Reused support

`T34` may reuse only already exported material:

1. `F163/N285`
   - `(omega,phi)` already support one local source-side derivative datum for
     component 1,
2. `N289`
   - current component-1 blocker still includes missing object carrier and
     missing pair-indexed seed carrier,
3. `T26/N284`
   - component 2 requires one pair-indexed noncyclic population anchor for at
     least `[pair1,pair2]`,
4. `R1/C48/C49`
   - the codomain-side type for component 2 remains
     `theta_1/theta_2 -> u_1/u_2 -> S_orient_cand`,
5. `F2`
   - `K_strict_gate` remains a later-pipeline operational strict working
     kernel and may not silently inherit stronger roles without explicit
     bridging or typed export.

## Exact decision question

Can the current repo honestly export anything stronger than:

```text
(omega, phi) are one real internal parameter pair of K_strict_gate,
but they remain only source-side scalar kernel parameters
and not an actual pair-indexed population anchor for component 2
```

for the present proposal?

## Boundary target

If the answer is negative, freeze:

```text
OmegaPhi_component_2_pair_indexed_population_anchor_incompatibility_boundary_v1
```

with the intended meaning:

```text
the pair (omega, phi) may already support local source-side seed calculations,
but the current repo exports no typed reduction
(omega, phi) -> (theta_1, theta_2),
no pair-indexed carrier over [pair1,pair2],
and no basis-pair population rule derived from that pair,
so (omega, phi) cannot yet serve as the component-2 pair-indexed population
anchor on the current repo state
```

## Hard limits

`T34` must not claim:

1. that `(omega,phi)` are useless,
2. that future typed transport from `(omega,phi)` is impossible in principle,
3. actual support for component 2,
4. actual `theta_1`, `theta_2`,
5. actual populated basis-pair instance,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. ToE closure.
