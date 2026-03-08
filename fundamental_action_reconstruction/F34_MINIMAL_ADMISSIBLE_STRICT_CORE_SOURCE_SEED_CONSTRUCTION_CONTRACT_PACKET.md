# F34 Minimal Admissible Strict-Core Source Seed Construction Contract Packet

Status: `F34_EXECUTED_MINIMAL_ADMISSIBLE_STRICT_CORE_SOURCE_SEED_CONSTRUCTION_CONTRACT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N130`, the next honest constructive question is:

```text
what is the minimal admissible construction contract for S_sel_int itself?
```

## Minimal source-seed construction contract

A future object may count as an admissible constructed `S_sel_int` only if all
of the following are satisfied:

1. **genuinely new strict-core source object**
   - it is added as a new strict-core exported object rather than reusing
     current `psi0`, FR, `sigma_int`, overlay-fit, or axiom-lane artifacts,
2. **carrier-typed enough for later export**
   - it has enough strict-core carrier/type structure that a later
     `E_orient` export can be meaningfully defined from it,
3. **source-seed only**
   - it is counted only as construction of `S_sel_int` itself and not as a
     hidden discharge of `E_orient`, `B_sel`, `R_sel`, or `O_sel`,
4. **strict-core only**
   - it is exported in strict core rather than on control / extension / axiom
     lanes,
5. **non-substitutive**
   - it does not silently identify `K_legacy_ont` with `K_strict_gate`,
6. **selector-acceptance independent**
   - theory-level selector acceptance outside strict core may not count as the
     construction of `S_sel_int`,
7. **future-bridge compatible**
   - it is stable enough that later `E_orient` and then `B_sel` could honestly
     start from it.

## Why this contract is forced

The contract is forced by the current repo state:

1. `N124/N126` keep current strict core without any admissible selector source,
2. `N130` reduces the first remaining branch to construction of `S_sel_int`,
3. `E_orient` remains later and may not be laundered into the source step,
4. therefore the source step itself needs a minimal admissible contract.

## What F34 does not claim

`F34` does not claim:

- that `S_sel_int` already exists,
- that `S_sel_int` is already constructible,
- that `E_orient` already exists,
- that downstream `B_sel -> R_sel -> O_sel` is solved,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. test whether the current repo already reduces the last positive branch to
   one explicit minimal construction contract for `S_sel_int`,
2. and if so, freeze that contract as the only honest next constructive move.
