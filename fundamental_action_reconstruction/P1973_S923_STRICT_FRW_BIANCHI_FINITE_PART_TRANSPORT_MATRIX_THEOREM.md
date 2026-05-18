# P1973 / S923 — Strict FRW/Bianchi-I Finite-Part Transport Matrix Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Local transport residual: `PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

After `P1972` closed the local B1 counterterm-cancellation residual, the next
honest background-independence move is not to claim global covariance.  The
honest move is to build the first explicit finite-part transport matrix for the
already exported `P1897` FRW/Bianchi-I `nu` branch and test whether the B1 zero
renormalization residual remains locked under that local transport.

## Strict source

`P1897` exports the preclosure bridge form:

```text
G_BI_b1 = (1 + nu*sigma2) G_BI_frw.
```

Define:

```text
q := nu*sigma2,
m(lambda) := 1 + lambda*q,
lambda in [0,1].
```

The local admissible branch is:

```text
1 + q > 0.
```

This excludes the singular transport point and keeps the local map invertible.

## Transport equation

For a finite part `F(lambda)` transported from FRW to Bianchi-I, define:

```text
F(lambda) = m(lambda) F_FRW.
```

Then:

```text
dF/dlambda = q F_FRW
A(lambda)  = q/(1 + lambda*q)
A(lambda) F(lambda) = q F_FRW
```

so the residual is exactly:

```text
dF/dlambda - A(lambda) F(lambda) = 0.
```

For the four B1 curvature coefficients, the local matrix transport is:

```text
T(lambda) = m(lambda) I_4,
A_4(lambda) = q/(1 + lambda*q) I_4.
```

The matrix residual is:

```text
dT/dlambda - A_4(lambda) T(lambda) = 0_4.
```

At `lambda=1`:

```text
T(1) = (1 + q) I_4,
det T(1) = (1 + q)^4.
```

Thus the local branch is invertible whenever `1 + q != 0`, and orientation
preserving on the stricter branch `1 + q > 0`.

## Link to P1972

`P1972` exports the B1 renormalization residual vector:

```text
[0,0,0,0].
```

Because linear transport maps the zero vector to the zero vector,

```text
T(lambda) [0,0,0,0]^T = [0,0,0,0]^T.
```

Therefore the B1 divergent counterterm cancellation remains locked under this
local scalar finite-part transport branch.

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.json
```

The witness includes:

1. symbolic scalar transport residual,
2. symbolic four-by-four matrix residual,
3. determinant condition for `T(1)`,
4. replay via `scipy.linalg.expm` for representative positive local `q` values,
5. explicit false-pass boundaries.

## What is now honestly improved

The repo now has a concrete local transport equation connecting the B1
renormalization finite-part lock to the P1897 FRW/Bianchi-I preclosure bridge.
This advances the background-independence task from a bridge statement to a
machine-checkable local transport residual.

## False-pass boundary

This packet does **not** prove:

1. global background-independence,
2. full variational covariance of every EOM tensor,
3. transport through the singular branch `1 + nu*sigma2 = 0`,
4. BRST closure,
5. Cutkosky/unitarity closure,
6. `QW-2191` selector discharge,
7. ToE closure.

The status therefore remains `OPEN_OBSTRUCTION_WITH_TRACE` globally.

## Next honest step

Lift this scalar finite-part transport to the full FRW/Bianchi-I variational
operator bundle.  Concretely, compute componentwise EOM residual tensors for the
`nu`-dependent anisotropic terms and verify whether each tensor component maps
covariantly, not merely by a scalar multiplier.

## Lay explanation

Po skasowaniu jednej grupy nieskończoności sprawdzamy, czy ten wynik przenosi
się między dwoma sposobami opisu geometrii kosmosu.  Dla najprostszego lokalnego
mnożnika most działa dokładnie.  To jest ważny test techniczny, ale jeszcze nie
dowód pełnej niezależności od tła, bo pełne równania mają więcej składników niż
jeden mnożnik.
