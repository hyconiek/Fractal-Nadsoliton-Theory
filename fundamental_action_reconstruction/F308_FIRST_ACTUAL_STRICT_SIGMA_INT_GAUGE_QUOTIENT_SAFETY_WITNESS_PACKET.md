# F308 First Actual Strict Sigma-Int Gauge-Quotient Safety Witness Packet

Status: `F308_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_GAUGE_QUOTIENT_SAFETY_WITNESS_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T148/P388` keep one strict prerequisite explicit for any honest strict-core
bridge/export-map claim from sigma-int into the residual orientation target
slot:

```text
theorem-level gauge-quotient safety for the strict sigma-int datum (T123/N388).
```

After `F307/N418`, the repo now exports sigma-int as a strict-core source datum
on a declared strict domain, with explicit premise-based provenance and no
hybrid FR reuse.

This packet executes one narrow, audit-safe move:

```text
export one declared gauge action G ⟲ C_v1
and export one explicit invariance/quotient-level witness showing
sigma_int is gauge-quotient-safe on that declared strict domain,
without gauge fixing and without implied selector closure.
```

## Inputs reused (strict-admissible)

1. `F306/N417`
   - declared strict configuration space `C_v1`,
   - strict witness `pi_1(C_v1) ≅ Z_2` with generator label `gamma_pi1_v1`.
2. `F307/N418`
   - strict FR-sign map object `chi_FR_strict_v1`,
   - strict sigma-int source-upgrade object `sigma_int_strict_derived_v1`,
   - explicit premise-based provenance (no hybrid reuse).
3. `T123`
   - acceptance tests for gauge-quotient safety.

## Declared gauge action (no gauge fixing)

Declare:

```text
C_v1 := Map_*^deg=1(S^3, SU(2)).
G    := Map_*(S^3, SU(2)).
```

Define a (strict, observer-free) gauge action by pointwise conjugation:

```text
(g ⋅ U)(x) := g(x) U(x) g(x)^(-1),
```

for all `g ∈ G`, `U ∈ C_v1`, `x ∈ S^3`.

This action preserves the based condition `U(∞)=1` and preserves the degree
component `deg=1` (no silent exit from the declared strict domain).

No gauge fixing is used in this packet.

## Gauge-quotient semantics for sigma-int

On the declared strict domain, `sigma_int_strict_derived_v1` is exported as:

```text
sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1) = -1,
```

and therefore it is a domain-level constant datum on `C_v1` (not a coordinate-
dependent gauge choice).

Define the induced class function:

```text
sigma_int_strict_value : C_v1 -> {+1,-1}
sigma_int_strict_value(U) := sigma_int_strict_derived_v1.
```

Then `sigma_int_strict_value` descends to a well-defined function on the orbit
space `C_v1/G`.

## Invariance witness (theorem-level, topological)

For any `g ∈ G`, the map `a_g : C_v1 -> C_v1`, `a_g(U)=g⋅U`, is a
homeomorphism. It therefore induces an automorphism on the fundamental group:

```text
(a_g)_* : pi_1(C_v1) -> pi_1(C_v1).
```

From `N417`, `pi_1(C_v1) ≅ Z_2`, so `pi_1(C_v1)` has a unique nontrivial
element. Therefore every automorphism fixes the nontrivial class, hence:

```text
(a_g)_*(gamma_pi1_v1) = gamma_pi1_v1.
```

Because `chi_FR_strict_v1` is the unique nontrivial character
`pi_1(C_v1)->{+1,-1}` exported in `F307`, it is invariant under such
automorphisms, and thus:

```text
chi_FR_strict_v1((a_g)_*(gamma_pi1_v1)) = chi_FR_strict_v1(gamma_pi1_v1).
```

Equivalently, the strict sigma-int value is gauge-invariant:

```text
sigma_int_strict_value(U) = sigma_int_strict_value(g⋅U)   for all g ∈ G, U ∈ C_v1.
```

This discharges the `T123` gauge-quotient safety acceptance requirement on the
declared strict domain without any gauge fixing and without relying on axiom-
lane promotion.

## Persisted artifact

`fundamental_action_reconstruction/generated/sigma_int_gauge_quotient_safety_witness_v1.json`

## Status discipline

This packet does **not** claim:

1. selector-track identification beyond overlay-only (`T147/N414`),
2. export of any residual-datum bridge/export-map object (`N301`) or discharge
   of `N300`,
3. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
4. ToE closure.

