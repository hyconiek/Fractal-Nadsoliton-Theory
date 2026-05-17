# P1958 Strict Local Abelian Gauge-Fixing Ghost-Action Seed Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__LOCAL_ABELIAN_GHOST_SEED_EXPORTED_NO_BRST_PASS`
As of: `2026-05-17`

## Pre-Execution Grep

Before execution, the repository was searched in English and Polish for:

```text
gauge sector, L_gauge, F_{mu nu}, Yang/Mills, structure constants,
gauge fixing, Faddeev/Popov, ghost action, antighost,
BRST differential, gauge transformation,
sektor cechowania, cechowanie, ustalenie cechowania,
duchy, antyduch, akcja duchow, transformacja cechowania,
pochodna kowariantna
```

Relevant existing sources:

```text
P1657: minimal covariant gauge-metric density, with gauge fixing not discharged.
P1907: SM gauge-sector registry, with ghost/BRST constraints OPEN_SYMBOLIC.
P1764: covariant operator-level E_A^mu, with BRST/Cutkosky still open.
P1957: full BV/BRST ghost-sector witness pack not available.
```

No theorem-grade structure-constant table, `BV_BRST_operator_map`, or nonproxy
ghost/antighost action was found. `TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf`
was ignored.

## Result

`P1958` exports the smallest gauge-fixing/ghost object that can be honestly
derived from the currently exported `F^2` gauge sector without importing
unexported non-Abelian data.

Scope:

```text
single local Abelianized gauge generator
flat B1 tangent patch
eta_mu_nu = diag(-1,1,1,1)
A_mu lower-index convention
delta A_mu = partial_mu alpha
```

The Lorenz gauge functional is:

```text
G[A] = partial_mu A^mu
     = -partial_t A_0 + partial_x A_1 + partial_y A_2 + partial_z A_3
```

The gauge-fixing Lagrangian is:

```text
L_GF = -(1/(2*xi))*G[A]^2
```

The Faddeev-Popov operator is computed by:

```text
M_FP alpha = d/deps G[A_mu + eps*partial_mu alpha] | eps=0
```

The machine result is:

```text
M_FP alpha - Box alpha = 0
```

The ghost action seed is:

```text
L_ghost = - cbar * Box c
```

and the by-parts identity is machine-checked:

```text
L_ghost - (partial_mu cbar partial^mu c + boundary_divergence) = 0
```

## Machine Checks

The script verifies:

```text
delta F_{mu nu} = 0
M_FP = Box
operator ghost form = by-parts ghost form + boundary term
all_local_checks_zero = true
```

## Scope Guard

This updates `P1957` only as:

```text
B1_ghost_sector_nonproxy_export =
  PARTIAL_SEED_ONLY_FOR_LOCAL_FREE_ABELIANIZED_GAUGE_FIELD
```

Still open:

```text
B2_BV_BRST_operator_map
B3_explicit_BRST_charge_Q
B4_Q_squared_simplified_zero
B5_physical_state_cohomology
B6_gauge_gauge_ghost_cancellation_trace
```

Therefore:

```text
TG2_BRST_GLOBAL_NILPOTENCY = NOT_PROMOTED
TG3_CUTKOSKY_GLOBAL_UNITARITY = NOT_PROMOTED
```

## Outputs

- `p1958_s908_strict_local_abelian_gauge_fixing_ghost_action_seed.py`
- `generated/p1958_s908_strict_local_abelian_gauge_fixing_ghost_action_seed.json`

## Next Honest Step

Build `P1959` with high reasoning: decide whether the strict lane has enough
exported gauge-group data to extend this seed to `SU(3)xSU(2)xU(1)` structure
constants and BRST differential, or prove a non-Abelian structure-data
obstruction.
