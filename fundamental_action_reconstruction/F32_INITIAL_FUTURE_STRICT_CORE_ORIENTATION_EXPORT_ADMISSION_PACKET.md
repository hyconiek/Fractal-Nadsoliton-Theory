# F32 Initial Future Strict-Core Orientation Export Admission Packet

Status: `F32_EXECUTED_INITIAL_FUTURE_STRICT_CORE_ORIENTATION_EXPORT_ADMISSION_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N128`, the last positive branch is already reduced to one initial seed:

```text
S_sel_int -> E_orient
```

`F32` asks the next honest constructive question:

```text
what admissible contract must E_orient satisfy if the future seed is ever to count?
```

## Orientation export contract

A future `E_orient` may count as the required strict-core orientation export
only if all of the following are satisfied:

1. **derived from future source object**
   - `E_orient` is exported from the future strict-core source object
     `S_sel_int`, not imported from outside the seed chain,
2. **strict-core only**
   - `E_orient` is exported in strict core and not only on
     control / extension / axiom lanes,
3. **internal orientation datum or theorem-level equivalent**
   - it supplies the missing internal orientation datum, residual selector
     datum, or a theorem-level equivalent selector-bearing source,
4. **selector-bearing without external anchor**
   - it does not rely on imported `psi0`, explicit selector control, or theory-
     level acceptance outside strict core as if those were the export itself,
5. **quotient / gauge safe**
   - it is compatible with quotient/gauge safety strongly enough to remain a
     genuine internal export rather than an overlay-only witness,
6. **bridge-ready for B_sel**
   - it is typed/export-stable strongly enough that a future strict-core bridge
     `B_sel` could honestly start from it,
7. **no silent kernel substitution**
   - it does not silently identify `K_legacy_ont` with `K_strict_gate`.

## Why this contract is forced

The contract is forced by the current repo state:

1. `F29/N126` already freeze the admission gate for a future source object,
2. `B2` still keeps the internal orientation datum absent in current strict
   core,
3. `N4/N5` keep `psi0` from counting as a strict-core selector-source export,
4. `N7/N8` keep the current residual-datum bridge absent,
5. `N123` forbids silent legacy-to-strict substitution,
6. `N125` keeps selector acceptance outside strict core,
7. `N128` reduces the first future seed to `S_sel_int -> E_orient`.

Therefore `E_orient` cannot be left as a vague placeholder; it needs a
strictly admissible export contract.

## What F32 does not claim

`F32` does not claim:

- that `E_orient` already exists,
- that `S_sel_int` already exists,
- that the seed is already constructible,
- that `B_sel`, `R_sel`, or `O_sel` are solved,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. test whether the current repo already reduces the last positive branch to
   one explicit package
   `admissible S_sel_int + E_orient export contract`,
2. and if so, freeze that package as the only honest initial future
   construction package.
