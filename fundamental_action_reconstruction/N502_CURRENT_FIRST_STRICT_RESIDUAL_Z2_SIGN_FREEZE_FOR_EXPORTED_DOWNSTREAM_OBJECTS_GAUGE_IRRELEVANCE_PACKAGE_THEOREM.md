# N502 Current First Strict Residual `Z2` Sign Freeze for Exported Downstream Objects — Gauge‑Irrelevance Package Theorem (No False‑PASS)

Status: `N502_DISCHARGED_CURRENT_FIRST_STRICT_RESIDUAL_Z2_SIGN_FREEZE_FOR_EXPORTED_DOWNSTREAM_OBJECTS_GAUGE_IRRELEVANCE_PACKAGE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

After the strict lane‑scoped `O(2) -> Z2` cuts (`N487`, `N496`) and the corresponding internal orientation datum exports (`N492`, `N500`),
the remaining ambiguity on each degenerate pair plane is residual sign:

```text
u  ->  -u.
```

This theorem packages the narrow strict clarification required to avoid a false next blocker:

```text
for the currently exported strict downstream objects actually used in the strict closure stack,
residual Z2 sign is a tracked gauge/convention layer and does not change the downstream object/audit result.
```

This does **not** derive a sign‑sensitive physical orientation datum, does **not** discharge `QW-2191` globally, and does **not** claim ToE closure.

## Strict-admissible inputs reused

1. `N493` / `N495`
   - residual `Z2` sign flips (and even full `O(2)` rotations) are conjugation/basis gauge for the `QW-2190` embedding audits.
2. `N501`
   - residual sign flips are irrelevant for the `R1` target slot because the slot semantics is span/projector‑based.
3. `F456`
   - strict downstream operator object `A_1(pair1) := |u_1><u_1|` constructed from the exported `u_1`, explicitly sign‑gauge invariant.
4. `F457`
   - strict lane-scoped transition angle export `alpha_12 := (theta_2-theta_1) mod 2π` and its axis-only reduction `alpha_12 mod π`.
5. `A10`
   - anti-overclaim boundary.

## Theorem (residual sign is frozen as gauge for the exported downstream objects)

### Claim 1. Residual sign does not change the `R1` target-slot inhabitant object.

By `R1`, the target slot is the subspace object:

```text
S(theta_1,theta_2) := span{u_1(theta_1), u_2(theta_2)}.
```

By `N501`, for any `s_1,s_2 ∈ {+1,-1}`:

```text
span{s_1 u_1, s_2 u_2} = span{u_1, u_2},
```

so residual sign does not change the `R1` target-slot object (equivalently, its span projector is unchanged). ∎

### Claim 2. Residual sign does not change the exported strict operator `A_1(pair1)`.

By `F456`, the exported operator is the rank‑one projector:

```text
A_1(pair1) := |u_1><u_1|.
```

For the residual sign flip `u_1' := -u_1`:

```text
|u_1'><u_1'| = |-u_1><-u_1| = |u_1><u_1|.
```

Therefore residual sign does not change the exported operator object `A_1(pair1)`. ∎

### Claim 3. Residual sign (and even full `O(2)` rotations) does not change the `QW-2190` embedding audits.

By `N493`, residual sign flips of an orthonormal basis inside a declared `QW-2190` embedding subspace act only by conjugation
on the embedded generators and preserve all invariance and Lie‑closure audits.

By `N495`, the same holds for the full continuous `O(2)` basis freedom: it is likewise conjugation‑only gauge for the `QW-2190` embedding audits.

Therefore, within the `QW-2190` audit scope, residual sign does not change the audit results. ∎

### Claim 4. Residual sign does not change the exported axis-only transition angle `alpha_12 mod π`.

By `F457`, the exported transition data includes:

```text
alpha_12_mod_pi := (theta_2 - theta_1) mod π.
```

Under residual sign flips `u_i -> -u_i`, the local phase representatives may shift by `theta_i -> theta_i + π` (mod `2π`),
but the axis-only difference is invariant:

```text
( (theta_2 + k2·π) - (theta_1 + k1·π) ) mod π  =  (theta_2 - theta_1) mod π
```

for any `k1,k2 ∈ {0,1}`. Therefore residual sign does not change the exported `alpha_12 mod π`. ∎

### Conclusion (B3 continuation, option A)

In the strict scope of the currently exported downstream objects actually used by the strict closure stack:

- `R1` target-slot inhabitant is a span/projector object,
- `A_1(pair1)` is a projector,
- `alpha_12 mod π` is an axis-only transition angle,
- `QW-2190` embedding audits are conjugation‑invariant,

so residual `Z2` sign can be **frozen** as a tracked gauge/convention layer without affecting these downstream objects/audits.

This closes the need to *lift* residual sign **for these specific objects**; any future claim that depends on a sign‑sensitive physical orientation
must still either:

1. prove sign gauge‑irrelevance for that specific observable, or
2. derive a genuinely sign‑sensitive strict datum without hidden marked‑direction inputs,

before promotion. ∎

## What `N502` does not claim

`N502` does not claim:

1. a sign‑sensitive physical orientation datum derived (residual sign remains a convention layer outside the stated objects),
2. strict-core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.
