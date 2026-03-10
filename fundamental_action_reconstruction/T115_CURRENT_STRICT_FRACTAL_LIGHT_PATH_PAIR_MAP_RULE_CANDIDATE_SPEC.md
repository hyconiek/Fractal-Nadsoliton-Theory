# T115 Current Strict Fractal-Light Path Pair-Map Rule Candidate Spec

Status: `T115_CURRENT_STRICT_FRACTAL_LIGHT_PATH_PAIR_MAP_RULE_CANDIDATE_SPEC_NO_FALSE_PASS`
As of: `2026-03-10`

## Goal

After `N292/N293/N294` and `F194`, the fractal branch is in a sharply stated
state:

```text
actual carrier object present
actual interface support present
actual fractal-to-pair map rule absent
```

At the same time, the strict closure lane (`N327`) isolates the dominant
missing ingredient as one genuine source-side, observer-free, pair-indexed,
noncyclic selector/provider object-carrier layer.

`T115` proposes one **strict-side candidate** map-rule form intended to:

1. keep the ontology discipline `nadsoliton -> light -> matter -> observer`,
2. treat the oscillatory ("light") motif as part of the coupling,
3. treat the fractal substrate as the multiplicity/geometry carrier,
4. remain **noncyclic** (no `theta` input, no populated-instance input),
5. remain **observer-free** (no `K_obs` as primary source),
6. remain below: actual theta export, actual pair population, and closure.

`T115` does not claim an actual map rule exists or is discharged.

## Admissible move

The admissible move here is not:

1. actual `theta_1`, `theta_2`,
2. actual populated basis-pair instance,
3. actual component-2 support witness,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. selector closure,
7. ToE closure.

The admissible move is weaker:

```text
export one explicit strict-side pair-map-rule candidate form
```

that is compatible with:

- noncyclic guardrails (`N18`),
- observer-free ordering (`AX9` discipline),
- and the strict kernel as an operational coupling channel (`K_strict_gate`).

## Candidate law-form (fractal-light phasor reduction)

Reuse strict kernel notation:

```math
\Theta(d):=\omega d+\phi,
\qquad
D(d):=1+\beta d^{\eta},
\qquad
K_{strict\_gate}(d)=\frac{\cos(\Theta(d))}{D(d)}.
```

Introduce, for each designated pair slot `i ∈ {1,2}`, a finite path-ensemble
carrier:

```text
E_i(pair) := { (d_{i,k}, w_{i,k}) }_{k=1..n_i}
```

with:

1. distances `d_{i,k} >= 0`,
2. weights `w_{i,k} >= 0`,
3. a normalization convention (candidate): `sum_k w_{i,k} = 1`.

Define the phasor components:

```math
X_i^{cand}
:=
\sum_{k=1}^{n_i} w_{i,k}\frac{\cos(\Theta(d_{i,k}))}{D(d_{i,k})},
\qquad
Y_i^{cand}
:=
\sum_{k=1}^{n_i} w_{i,k}\frac{\sin(\Theta(d_{i,k}))}{D(d_{i,k})}.
```

Then define the candidate strict-side phase output:

```math
\theta_i^{cand}:=\operatorname{atan2}\!\left(Y_i^{cand},X_i^{cand}\right).
```

Finally define:

```math
u_i^{cand}:=\cos(\theta_i^{cand})c_i+\sin(\theta_i^{cand})s_i,
\qquad
S_{orient}^{cand}:=\operatorname{span}\{u_1^{cand},u_2^{cand}\}.
```

Interpretation discipline:

1. the oscillatory numerator (`cos/sin`) is the light-coupling motif,
2. the denominator `D(d)` is the fractal-geometry damping motif,
3. the ensemble `{(d_{i,k},w_{i,k})}` is the finite carrier of fractal path
   multiplicity,
4. nothing here presumes an observer; everything is pair-indexed and
   source-side by construction intent.

## Noncyclic / observer-free constraints

`T115` is guardrail-safe only if:

1. `E_i(pair)` is supplied without using `theta_1, theta_2` as inputs,
2. `E_i(pair)` is supplied without using a populated basis-pair instance as
   input,
3. no `K_obs`-indexed selection is used as a primary source,
4. no selector-closure claim is made from the presence of a candidate rule.

## Degeneracy frontier (honest failure mode)

If for some `i`:

```math
X_i^{cand}=0 \ \text{and}\ Y_i^{cand}=0,
```

then `atan2` is undefined and the candidate does not produce `theta_i^{cand}`.

This is not a defect to hide. It is a **frontier**:

1. it signals missing non-degeneracy structure in the carrier `E_i(pair)`,
2. it is compatible with the known strict-core uniqueness/selector obstruction
   (`QW-2191`) and must not be silently bypassed.

## Candidate object

If the repo admits the above form as a strict-side candidate rule, export:

```text
M_fractal_light_path_pair_map_rule_candidate_v1
```

with intended meaning:

```text
a noncyclic, observer-free, pair-indexed candidate rule form
for reducing a finite fractal path-ensemble carrier E_pair
through the strict kernel coupling channel
into (theta_1^cand, theta_2^cand) and S_orient^cand,
without claiming actual theta export or population
```

## Hard limits

`T115` must not claim:

1. actual fractal-to-pair map rule,
2. actual `theta_1`, `theta_2`,
3. actual populated basis-pair instance,
4. actual component-2 support,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. `QW-2191` discharge,
9. ToE closure.

