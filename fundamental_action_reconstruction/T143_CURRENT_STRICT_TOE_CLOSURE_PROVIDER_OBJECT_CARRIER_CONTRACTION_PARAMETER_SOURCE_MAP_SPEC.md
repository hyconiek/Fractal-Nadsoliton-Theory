# T143 Current Strict ToE Closure Provider-Object Carrier Contraction-Parameter Source Map Spec

Status: `T143_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_CONTRACTION_PARAMETER_SOURCE_MAP_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T142/N409` isolate one concrete upstream gap on the strict provider-object
carrier lane:

```text
the contraction parameters (a,b) used by the provider-object carrier candidate
lane are not yet exported as strict-derived source-side outputs of tau_src_candidate_v1
```

The next honest move is therefore not a promotion of the provider-object lane
into realization, but a narrow upstream discharge attempt:

```text
export one explicit source-side map tau_src_candidate_v1 -> (a,b),
derived from already exported source topology control data,
so that (a,b) are no longer free external knobs on the provider-object carrier lane.
```

`T143` specifies exactly such a source-side contraction-parameter map, using
only:

1. the exported strict source packet `tau_src_candidate_v1` (`F127/N235`),
2. the exported strict-kernel core flow amplitude `cos(phi)` already carried by
   that packet (`F74/F127`),
3. the exported barrier-protected sign and margin witnesses ensuring
   `0 < cos(phi) < 1` on the declared core branch (`F141/N249`).

This is kernel-split-safe: it uses `K_strict_gate` only through the already
exported core-limit control datum `K_strict(0)=cos(phi)` and does not claim any
legacy-to-strict bridge.

## Map object

Export one explicit source-side map:

```text
A_strict_provider_object_contraction_parameter_source_map_v1 :
  tau_src_candidate_v1 -> (a,b)
```

defined by:

```text
(a,b) := (cos(phi), cos(phi))
```

where `cos(phi)` is the strict-kernel core flow amplitude carried by
`tau_src_candidate_v1` via:

```text
T_flow^(0) := cos(phi) * e_topo.
```

## Contraction acceptance facts (declared; relies only on exported witnesses)

From `F141/N249` the repo already exports:

1. `sign(cos(phi)) = +1`,
2. a positive barrier margin
   `delta_src_barrier_sign_margin_v1 := (pi/2) - |phi| > 0`,
3. `phi` is fixed and nonzero on the declared core branch.

Therefore on that declared strict branch:

```text
0 < cos(phi) < 1
```

so the exported parameters satisfy:

```text
0 < |a| < 1,
0 < |b| < 1.
```

## Noncyclic / observer-free discipline

`T143` is admissible only if:

1. no `theta_{1,2}` are used as inputs,
2. no populated basis-pair instance is used as input,
3. no `K_obs`-indexed selection is used as a primary source of uniqueness,
4. no discharge of `N302` is implied.

## Hard limits

`T143` must not claim:

1. discharge of `Epsilon_strict_provider_object_carrier_layer_target_v1`,
2. discharge of `N302` (actual residual bridge/export-map object support),
3. any export-map object export (`N300`),
4. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
5. ToE closure.

