# P2019 S969 Strict Cutkosky First Transverse Tree-Amplitude Component Theorem

Status: `P2019_FIRST_LOCAL_TRANSVERSE_TREE_COMPONENT_EXPORTED_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

Move beyond governance-only audits by computing one concrete strict-side local
amplitude component for the channel

```text
graviton -> gauge_gauge.
```

P2019 uses the already exported P1955 minimal `hAA` vertex and the P1956 local
transverse two-gauge polarization projector.  It does **not** add a new
quadrature ansatz.

## Construction

In the massless two-body cut frame with `E=1`, P2019 takes

```text
k1 = (1,0,0, 1),
k2 = (1,0,0,-1),
eps_x = (0,1,0,0),
eps_y = (0,0,1,0),
eta = diag(-1,1,1,1).
```

From the P1955 minimal vertex

```text
L_hAA = (kappa/2) h^{mu nu} T^gauge_{mu nu},
```

P2019 computes the cross stress tensor for two gauge waves:

```text
T_cross_{mu nu}
  = F1_{mu rho} F2_nu^rho + F2_{mu rho} F1_nu^rho
    - (1/2) eta_{mu nu} F1_{rho sigma} F2^{rho sigma}.
```

It then contracts this with two local graviton transverse polarizations:

```text
h_plus  = diag(0,1,-1,0),
h_cross has h_xy=h_yx=1.
```

The output amplitudes are reported in units of `kappa*Z_gauge`.

## Result

The local transverse tree-level amplitude matrices are:

```text
M_plus  = [[-2, 0], [0, 2]],
M_cross = [[ 0,-2], [-2,0]].
```

The local polarization-summed square is:

```text
sum |M_tree_transverse|^2 / (kappa^2 Z_gauge^2) = 16.
```

## No-false-pass boundary

P2019 is only a local tree-level transverse component.  It is not:

1. a dressed loop amplitude,
2. a BRST cohomology theorem,
3. a same-scheme `MSbar_B1_seed` dressed object,
4. a `DiscM = CutSum` proof,
5. all-state unitarity,
6. `QW-2191` discharge,
7. ToE closure.

Therefore P2019 updates the P1953 contract only partially:

```text
M_dressed_common_basis = PARTIAL_TREE_TRANSVERSE_COMPONENT_ONLY_NOT_DRESSED.
```

## Next honest step

Promote the local component toward the P1953 contract by adding either:

1. a same-scheme dressing/residue layer, or
2. a genuine BRST cohomology projector and ghost-cancellation trace,

then rerun P2018 to see which contract gates can be truthfully upgraded.
