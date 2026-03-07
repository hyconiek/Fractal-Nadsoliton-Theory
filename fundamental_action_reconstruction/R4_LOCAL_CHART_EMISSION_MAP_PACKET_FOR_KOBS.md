# R4 Local-Chart Emission Map Packet For K_obs

Status: `R4_EXECUTED_LOCAL_CHART_EMISSION_MAP_PACKET_FOR_KOBS_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

`P7` localized one remaining operator-chain object:

```text
explicit_emission_map_E_from_M_pair_to_L_int
```

`R4` does not try to build the full route

```text
existing kernel feedback -> H3 operator chain -> selector-facing K_obs
```

and it does not claim:

- strict-core selector reduction,
- factorization from current kernel feedback,
- pair-target privilege,
- or theorem-level discharge.

`R4` does something narrower:

- materialize one explicit current-pair representative
  `E_1 : V_1 -> L_1`,
- identify it as the current `M_pair -> L_int` representative only for the
  actual pair-level test on `pair1`,
- keep explicit that the map is only a local-chart preoriented packet built from
  the already exported `psi0` coordinate embedding.

## Inputs reused

1. `H30`
   - `psi0` is a deterministic kernel-invariant anchor candidate
2. `H31`
   - legal coordinate embedding
     `u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1`
3. `H33`
   - `pair1` is only a deterministic local chart
4. `H34`
   - no basis-covariance / target-independence discharge
5. `H35`
   - no strict physical axis selection inside `pair1`
6. `H36`
   - no strict directed-axis selection
7. `H37`
   - no sign-sensitive selector-state object on `pair1`
8. `R3`
   - explicit `G_light^(1)` packet on `L_1 = span{ell_+, ell_-}`
9. `P1`
   - explicit extension-lane matrix on `pair1` for the same
     `psi0 + anisotropic retardation` geometry

## Current-pair carrier identification

For the present actual pair-level test only, define:

```text
M_pair,current := V_1 = span{c_1, s_1}
L_int,current  := L_1 = span{ell_+, ell_-}
```

This packet does **not** claim:

- that `pair1` is the unique physical selector target,
- that `M_pair` is globally reduced to `pair1`,
- or that the resulting map is chart-independent.

## Local-chart emission map

Let

```text
u_psi0_pair1 := cos(psi0)c_1 + sin(psi0)s_1
v_psi0_pair1 := -sin(psi0)c_1 + cos(psi0)s_1
```

where `u_psi0_pair1` is exactly the legal local representative already exported
in `H31`, and `v_psi0_pair1` is its orthogonal chart complement inside
`V_1 = span{c_1, s_1}`.

Define `E_1 : V_1 -> L_1` by:

```text
E_1(u_psi0_pair1) = ell_+
E_1(v_psi0_pair1) = ell_-
```

In the ordered bases

```text
(c_1, s_1)  for V_1
(ell_+, ell_-) for L_1
```

this gives the explicit matrix

```text
E_1 = R(-psi0)
    = [[ cos(psi0),  sin(psi0)],
       [-sin(psi0),  cos(psi0)]]
```

## Consistency check with `R3` and `P1`

Using the explicit `R3` packet

```text
G_light^(1) = diag(lambda_+, lambda_-)
```

the partial pullback

```text
E_1^* G_light^(1) E_1
```

is a real symmetric `2x2` matrix on `V_1`.

`R4` verifies that this partial pullback exactly matches the already executed
extension-lane matrix from `P1/Test C configured_path_split`.

This is only a coherence check between:

- the local-chart emission packet,
- the light propagator packet,
- and the old `psi0 + anisotropic retardation` extension geometry.

It is **not** a proof that the full `H3` chain is instantiated.

## Result of `R4`

`R4` establishes:

1. one explicit current-pair emission map packet now exists,
2. it is serialized as a real orthogonal `2x2` matrix,
3. it is explicitly preoriented by `psi0`,
4. it is explicitly current-pair/local-chart scoped,
5. it coheres numerically with the already computed `P1` extension matrix when
   combined only with `R3`,
6. but it is still not a factorization of current kernel feedback and still not
   a selector-source discharge.

## Honest frontier after `R4`

`R4` does **not** establish:

- `R_mat : L_1 -> Q_mat`,
- `O_obs : Q_mat -> Q_mat`,
- an equivalence/factorization map from current kernel feedback to the `H3`
  chain,
- a full `H3` selector-facing projected block,
- strict-core selector closure,
- theorem-level closure,
- full ToE closure.

So the honest residual frontier becomes:

- `R4_B1 := an explicit current-pair local-chart emission map packet E_1 now exists, but it remains preoriented by psi0, current-pair scoped, and not yet a factorization of current kernel feedback into a full selector-facing H3 chain`
- `H33_B1 := pair1 is available as a deterministic local mode chart for the primary psi0 lane, but no strict-core justification elevates it to a uniquely selector-relevant reduction target`
- `H34_B1 := strict core contains local chart embeddings for psi0 but no basis-covariance or target-independence argument elevating those embeddings beyond chart dependence`
- `H35_B1 := strict core supports only a coordinate-level direction u_psi0_pair1 inside pair1 and contains no strict argument that psi0 selects a physically privileged axis there`
- `H36_B1 := strict core supports only a coordinate-level undirected axis representative u_psi0_pair1 inside pair1 and contains no strict argument selecting a directed orientation on that axis`
- `H37_B1 := strict core contains no sign-sensitive state object or observable on pair1 and therefore does not distinguish u from -u as physically different selector states`

## What `R4` does not claim

`R4` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that `pair1` is physically privileged,
- that `psi0 = theta_1`,
- that `E_1` is chart-independent,
- that `E_1` is already factorized from current kernel feedback,
- that `QW-2191` is discharged,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the exact current route with explicit `E_1` and `G_light^(1)` in scope,
2. accept only:
   - a blocker-set reduced by exactly one object,
   - or a review result if the route assumptions changed,
3. keep the route negative unless `R_mat`, `O_obs`, or the factorization map
   are actually instantiated.
