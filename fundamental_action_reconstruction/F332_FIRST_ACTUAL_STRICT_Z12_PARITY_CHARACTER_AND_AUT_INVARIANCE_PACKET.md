# F332 First Actual Strict `Z_12` Parity Character + `Aut(Z_12)` Invariance Packet

Status: `F332_EXECUTED_FIRST_ACTUAL_STRICT_Z12_PARITY_CHARACTER_AND_AUT_INVARIANCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `F329/N450` and `F331/N453` the repo exports:

1. a typed cyclic group carrier `Z_12_v1`,
2. the typed gauge symmetry `Aut_Z12_v1 ≅ (Z/12Z)^× = {1,5,7,11}` acting on the `Z_12` / phase-embedding layer.

One honest next question for the `AX20/T162` typed lane is:

```text
which ingredients are already invariant under Aut(Z_12),
so they do not introduce any hidden generator/orientation slot?
```

This packet executes one narrow, strict, audit-safe move:

```text
export the unique nontrivial order-2 character χ_parity : Z_12_v1 -> {+1,-1},
and export the theorem-level fact that χ_parity is invariant under Aut_Z12_v1.
```

This does **not** attempt any density-operator rigidity, any Berry/holonomy ingredient, any theta export,
any `O(2)` cut, or any `QW-2191` discharge.

## Inputs reused (strict-admissible)

1. `F329/N450`
   - typed cyclic group object `Z_12_v1`.
2. `F331/N453`
   - typed automorphism group `Aut_Z12_v1` and its action on `Z_12` / `Phase_12`.
3. `N435`
   - strict sigma-int lane already uses the sign-mask form `b_{1,k} := s^k` with `s=-1`.

## Packet result (typed parity character + invariance)

### 1) Typed parity character on `Z_12_v1`

Export the unique nontrivial order-2 character:

```text
chi_parity_Z12_v1 : Z_12_v1 -> {+1,-1}
chi_parity_Z12_v1(k) := (-1)^(k mod 2).
```

Equivalently (factorization through the unique quotient `Z_12 -> Z_2`):

```text
Z_12_v1  --(mod 2)-->  Z_2  --(sign)-->  {+1,-1}.
```

Persisted artifact:

`fundamental_action_reconstruction/generated/chi_parity_z12_v1_character.json`

### 2) Induced parity character on `Phase_12_v1`

Under the presentation `Phase_12_v1 = {ζ^k}`, export the induced character:

```text
chi_parity_Phase12_v1(ζ^k) := (-1)^(k mod 2).
```

Persisted artifact:

`fundamental_action_reconstruction/generated/chi_parity_phase_12_v1_character.json`

### 3) `Aut_Z12_v1` invariance (gauge-robustness)

Export the strict invariance statement:

```text
for all u ∈ Aut_Z12_v1 and all k ∈ Z_12_v1:
  chi_parity_Z12_v1(u·k) = chi_parity_Z12_v1(k).
```

Reason (finite arithmetic): every unit `u ∈ (Z/12Z)^×` is odd, so `(u·k mod 2) = (k mod 2)`.

## Meaning

This packet means only:

1. the sign pattern `(-1)^k` is now an explicitly exported **typed character** of `Z_12_v1`,
2. this character is **gauge-invariant** under the exported `Aut_Z12_v1` symmetry,
3. therefore any construction that depends only on the parity class `k mod 2` is already
   `Aut(Z_12)`-quotient safe at the phase-embedding layer (no hidden generator/orientation slot).

## Hard limits (no false pass)

`F332` does not claim:

1. any canonical choice of `Z_12_v1 -> Phase_12_v1` embedding (it only exports an invariant),
2. any density-operator `1/2` rigidity ingredient,
3. any Berry/holonomy construction,
4. slot-free sigma-int → theta construction (`T162`) or strict-core upgrade (`T159`),
5. strict theta export, strict-core selector closure, or `QW-2191` discharge,
6. ToE closure.

