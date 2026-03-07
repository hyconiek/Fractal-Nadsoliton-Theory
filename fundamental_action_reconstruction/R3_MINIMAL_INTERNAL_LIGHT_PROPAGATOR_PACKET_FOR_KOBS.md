# R3 Minimal Internal Light Propagator Packet For K_obs

Status: `R3_EXECUTED_MINIMAL_INTERNAL_LIGHT_PROPAGATOR_PACKET_FOR_KOBS_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

`P6` isolated one missing operator-chain object:

```text
explicit_internal_light_propagator_G_light_on_L_int
```

`R3` does not try to build the full chain

```text
existing kernel feedback -> H3 operator chain -> selector-facing K_obs
```

and it does not claim any factorization from existing kernel feedback.

`R3` does something narrower:

- materialize one explicit finite internal-light propagator packet
  `G_light^(1) : L_1 -> L_1`,
- derive it only from already persisted retardation/anisotropy proxies,
- keep explicit that:
  - no emission map `E` exists yet,
  - no `R_mat` exists yet,
  - no `O_obs` exists yet,
  - no factorization/equivalence map to the current kernel feedback exists yet,
  - no selector-sector projected block is exported.

## Inputs reused

1. `QW-1952`
   - `omega`
   - `retard_phase`
   - `anisotropy_strength`
2. `QW-1955`
   - persistence of `retard_phase` across the minimal-repair layer
3. `QW-1956`
   - persistence of `retard_phase` across the repaired-operator layer
4. `H29`
   - old retardation proxies modulate a pre-oriented channel and do not
     themselves generate the selector anchor
5. `H42`
   - minimal retardation law
     `phi_ret(L;c)=omega*L/c`
   - anisotropic light-channel form appears only after anisotropic path data

## Minimal carrier and operator

Define the minimal finite internal-light carrier

```text
L_1 := span{ell_+, ell_-}
```

with basis order:

```text
(ell_+, ell_-)
```

where `ell_+` and `ell_-` are only ordered light eigenchannels associated with
the scalar anisotropy split `1 +/- anisotropy_strength`.

This packet does **not** identify these channels with a spatial orientation in
the selector sector.

Using the already persisted retardation data, define:

```text
L0    := c * retard_phase / omega
L_+   := L0 * (1 + anisotropy_strength)
L_-   := L0 * (1 - anisotropy_strength)
lambda_+ := cos(omega * L_+ / c)
lambda_- := cos(omega * L_- / c)
```

in the same unit convention already used by the executable pair-level probe:

```text
c = 1
```

Then the minimal explicit light propagator packet is

```text
G_light^(1) = diag(lambda_+, lambda_-)
```

on `L_1`.

## Result of `R3`

`R3` establishes:

1. one explicit finite internal-light propagator packet now exists,
2. it is serialized as a real `2x2` matrix on `L_1`,
3. the matrix uses only `omega`, `retard_phase`, `anisotropy_strength`,
4. the matrix itself does **not** use `psi0`,
5. the packet is therefore explicit at operator level,
   but still not selector-facing and still not factorized from current kernel
   feedback.

## Honest frontier after `R3`

`R3` does **not** establish:

- `E : M_pair -> L_1`,
- `R_mat : L_1 -> Q_mat`,
- `O_obs : Q_mat -> Q_mat`,
- an equivalence/factorization map from current kernel feedback to the `H3`
  chain,
- a selector-sector projected `2x2` block,
- strict-core closure,
- theorem-level closure,
- full ToE closure.

So the honest residual frontier becomes:

- `R3_B1 := an explicit finite internal-light propagator packet G_light^(1) now exists on L_1, but it remains only an eigenchannel-level operator packet and not an identified factorization of current kernel feedback into a selector-facing H3 chain`
- `H14_B1 := existing kernel feedback is real but no explicit equivalence map or selector-sector reduction identifies it with the H-lane operator K_obs`
- `H15_B1 := existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback`
- `H29_B1 := old wave/retardation/memory proxies only modulate a pre-oriented anisotropic channel and do not by themselves generate an internal strict-core orientation anchor for theta_i`

## What `R3` does not claim

`R3` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that current kernel feedback already equals `K_obs`,
- that `G_light^(1)` is already strict-core selector data,
- that `psi0` is internalized,
- that `QW-2191` is discharged,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the exact `P6` route with `G_light^(1)` explicitly in scope,
2. accept only:
   - a blocker-set reduced by exactly one object,
   - or a review result if the route assumptions changed,
3. keep the route negative unless the remaining objects or the factorization map
   are actually instantiated.
