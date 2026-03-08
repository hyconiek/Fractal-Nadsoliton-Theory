# N256 Current First Actual Source Topology Basis Independent Promotion Witness Theorem

Status: `N256_DISCHARGED_CURRENT_FIRST_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`N256` packages the current `F148/P236` result into one theorem-level
current-state statement, without upgrading it into:

1. a quotient-safe `QW-2191` resolution theorem,
2. a current selector closure theorem,
3. a current global `QW-2191` discharge theorem.

## Fixed theorem statement

```text
N256_Current_First_Actual_SourceTopology_BasisIndependent_Promotion_Witness_Theorem

On the current repo state, one actual source-side basis-independent
selector-promotion witness is exported:

  Upsilon_sel_basis_actual_witness_v1 :
    tau_src_candidate_v1 -> Sigma_sel_basis_free_target_v1

supported by the current-repo-state packet

  W_src_basis_promotion_support_packet_v1
    = (
        tau_src_candidate_v1,
        Theta_src_nontriv_actual_discharge_witness_v1,
        Pi_sel_src_actual_witness_v1,
        Q_basis_sel_v1,
        selector_axis_basis_free_class_v1,
        selector_signed_split_basis_free_class_v1,
        preobserver_basis_free_scope_tag_v1,
        observer_downstream_only
      )

This witness remains:
  - observer-free in the witness domain,
  - source-side only,
  - basis-independent at the class level,
  - below quotient-safe QW-2191 resolution,
  - below current selector closure,
  - below current global QW-2191 discharge.
```

## Why this is the honest theorem

Because on the current repo state:

1. `F136/P224/N244` export only a future basis-independent promotion target,
2. `F146/P234/N254` add one actual full source-topology nontriviality
   witness,
3. `F147/P235/N255` add one actual source-side selector witness,
4. `F148/P236` add only the class-level basis reduction from that actual
   selector witness into `Sigma_sel_basis_free_target_v1`,
5. but quotient-safe `QW-2191`, selector closure, and global discharge remain
   open,
6. observer-side asymmetry remains downstream witness only.

Therefore the strongest honest theorem is:

```text
the current repo exports one actual source-side basis-independent
selector-promotion witness for tau_src_candidate_v1
```

and nothing stronger.

## Consequence

After `N256`, the `T14` frontier is narrower:

1. the basis-independent selector-promotion layer no longer remains only a
   future target,
2. it now contains one actual source-side basis-independent promotion witness
   for `tau_src_candidate_v1`,
3. but quotient-safe `QW-2191` resolution and selector closure remain
   downstream of that still-limited witness.

## Hard limits

`N256` does not discharge:

1. quotient-safe `QW-2191` resolution,
2. current selector closure,
3. current global `QW-2191` discharge,
4. ToE closure.
