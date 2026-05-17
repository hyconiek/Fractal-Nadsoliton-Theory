# P1955 Strict Minimal hAA Vertex Export Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__MINIMAL_TREE_hAA_VERTEX_EXPORTED_NO_FALSE_PASS`
As of: `2026-05-17`

## Pre-Execution Grep

Before execution, the repository was searched in English and Polish for:

```text
hAA, graviton gauge, metric perturbation, sqrt(-g), F_{mu nu},
physical-state projector, wierzcholek, amplituda ubrana, metryka,
zaburzenie metryczne, projektor stanu fizycznego
```

The search found an important prior source:

```text
P1657: sqrt(-g)*[-1/4 g^{mu alpha}g^{nu beta}F_{mu nu}F_{alpha beta}]
```

Therefore the honest next move is not another nonavailability packet, but a
minimal tree-level `hAA` vertex export.

## Result

`P1955` derives the metric linearization:

```text
g_mu_nu = eta_mu_nu + kappa*h_mu_nu
g^{mu nu} = eta^{mu nu} - kappa*h^{mu nu} + O(kappa^2)
sqrt(-g) = 1 + (kappa/2)*h + O(kappa^2)
```

from the minimal gauge density:

```text
L_gauge = -sqrt(-g)/4 * g^{mu alpha}g^{nu beta}F_{mu nu}F_{alpha beta}.
```

The exported vertex density is:

```text
L_hAA_minimal_tree = (kappa/2) h^{mu nu} T^gauge_{mu nu}
```

where

```text
T^gauge_{mu nu}
  = Z_gauge * (F_{mu rho}F_nu^rho
      - (1/4) eta_{mu nu}F_{rho sigma}F^{rho sigma}).
```

## Machine Check

The script verifies:

```text
linearized_L_gauge_coefficient - (1/2)h^{mu nu}T^gauge_{mu nu} = 0
```

and also checks symmetry and tracelessness of the Maxwell stress tensor in 4D.

## Scope

This discharges the `P1954` missing `M1` item only in the minimal tree-level
field-strength sense.

It does not provide:

1. `R F^2` or higher-curvature contact terms,
2. BRST physical-state projection,
3. dressed graviton propagator residues,
4. same-scheme `DiscM=CutSum`,
5. global `UR_link`.

## Outputs

- `p1955_s905_strict_minimal_hAA_vertex_export.py`
- `generated/p1955_s905_strict_minimal_hAA_vertex_export.json`

## Next Honest Step

Build `P1956`: attach a same-scheme BRST physical-state projector and
polarization-sum certificate for the `gauge_gauge` final state, or export the
exact missing ghost/projector data obstruction.
