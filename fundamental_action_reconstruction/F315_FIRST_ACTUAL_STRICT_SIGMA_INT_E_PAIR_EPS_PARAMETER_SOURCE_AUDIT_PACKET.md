# F315 First Actual Strict Sigma-Int `E_pair` Eps-Parameter Source Audit Packet

Status: `F315_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_E_PAIR_EPS_PARAMETER_SOURCE_AUDIT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `F312/N423` and `F314/N425`, the repo exports strict-input instantiation
artifacts of the sigma-int → residual-datum projection-candidate pipeline,
each fixing one explicit amplitude parameter choice:

```text
eps = 1/2
```

This packet performs the next honest anti-false-pass move:

```text
audit whether the current declared strict core exports any internal strict-derived
or strict-source-upgraded value object supplying the amplitude parameter eps
used in the sigma-int-driven E_pair generator weight law.
```

If no such internal eps source exists, then:

1. all E_pair generator uses of `eps` remain explicitly parameter choices,
2. all downstream theta outputs remain candidate-only and noncanonical with
   respect to `eps`,
3. no strict-core theta-source upgrade may be implied from any one chosen `eps`.

## Inputs reused (strict-admissible)

1. `T117/F270/N382`
   - sigma-int driven `E_pair` generator candidate:
     `w_{i,k} := (1 + sigma_int * eps * b_{i,k})/12`,
     with `parameter_eps = eps ∈ [0,1]`,
     and explicit non-claims of strict derivation/uniqueness.
2. `F312/N423`
   - strict-input instantiation fixing `eps = 1/2` as an explicit parameter choice
     (candidate-only output).
3. `F314/N425`
   - strict-input positive-window instantiation deriving `delta_d` from the strict kernel
     tuple, while still fixing `eps = 1/2` as an explicit parameter choice.
4. `A10`
   - anti-overclaim boundary.

## Audit scope (repo-wide)

The audit asks one narrow question:

```text
does the current strict core export a dedicated eps value object
with strict provenance (strict-derived or strict-source-upgraded),
distinct from unrelated epsilon/radius symbols used elsewhere in the repo?
```

Repo-wide grep was performed for explicit `eps` value-object exports
(`eps_strict_*`, `eps_*_derived`, `eps_source`, etc.) inside:

- `fundamental_action_reconstruction/`,
- `fundamental_action_reconstruction/generated/`.

### Semantic separation (critical)

The repo does export multiple `epsilon_*` radius objects (e.g. local barrier
radii such as `epsilon_src_local_barrier_radius_v1`), but these belong to
different semantics (angle perturbation/stability radii) and are **not**
exported as the amplitude parameter `eps` of the `E_pair` weight law.

No strict bridge identifying those radius epsilons with the generator amplitude
parameter `eps` is exported.

## Finding

On the current repo state:

1. no dedicated strict-core eps value object is exported,
2. `eps` appears only as:
   - a declared candidate parameter `eps ∈ [0,1]` in the generator definition, and
   - an explicit numeric choice in candidate instantiation artifacts.

Therefore `eps` remains a free candidate parameter, and must not be treated as
strict-derived or canonically fixed by the current strict core.

## Persisted summary

This packet persists the audit summary artifact:

```text
fundamental_action_reconstruction/generated/f315_strict_sigma_int_e_pair_eps_parameter_source_audit_summary.json
```

## Status discipline

This packet does **not** claim:

1. strict-core derivation or uniqueness of `eps`,
2. actual strict-core `theta_1`, `theta_2` export,
3. actual populated basis-pair instance,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

