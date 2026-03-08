# N186 Current Additive Preobserver Source Object Nonreduction Theorem

Status: `N186_DISCHARGED_CURRENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_NONREDUCTION_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Theorem statement

Within the current repo state, the first additive preobserver source-object
construction attempt

```text
S_preLM_additive_candidate_v1 := exp(A_up) u_T
```

is not equal to the same-basis packaged realization of the older `F75` target:

```text
S_preLM_target_packaged_realization_v1(d_*=1)
```

under the fixed basis realization:

```text
P_NL^(0) u_T := u_L
P_LM(d_*) u_L := I_mat(d_*) u_M
```

## Proof basis

1. `F76` exports the additive attempt with closed form
   `u_T + cos(phi) u_L + (cos(phi)/4) u_M`.
2. `F79` exports the same-basis packaged target realization with closed form
   `cos(phi) u_T + cos(phi) u_L + (cos(phi)/2) u_M`.
3. `F79` computes the explicit difference
   `Delta_preLM_nonreduction_v1 = (1 - cos(phi)) u_T - (cos(phi)/4) u_M`.
4. `P167` confirms that this difference is explicitly nonzero.

Therefore, on the current repo state, one explicit nonreduction witness exists.

## What this theorem does mean

`N186` supports the limited conclusion that:

```text
S_preLM_additive_candidate_v1
!=
S_preLM_target_packaged_realization_v1(d_*=1)
```

on the fixed same-basis realization.

## What this theorem does not mean

`N186` does not claim:

- that the first admissibility clause is already satisfied,
- that `S_preLM_additive_candidate_v1` is already a genuinely new strict-core
  source object,
- that the object is already exported as `constructed_source_object`,
- admissible `S_sel_int`,
- admissible `E_orient`,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
