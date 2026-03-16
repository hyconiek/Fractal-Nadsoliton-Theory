# F473 Current Strict `T164` `Z_12` Generator/Orientation Canonical‑Fixing Datum — Export Packet (No False‑PASS)

Status: `F473_EXECUTED_CURRENT_STRICT_T164_Z12_GENERATOR_ORIENTATION_CANONICAL_FIXING_DATUM_EXPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`T164` names the missing strict-core ingredient required to make any **directed / sign-sensitive** continuation honest:

- export an explicit fixing datum that reduces the `Aut(Z_12)` ambiguity in a tracked way,
- so later work cannot silently treat “choose generator 1” / “choose an orientation” as already canonical (`N462`).

This packet exports one **premise-based strict provenance** datum object implementing the minimal content of the `T164` acceptance tests:

```text
Kappa_Z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1
```

This is not a strict derivation of canonicity; it is an explicit premise recorded as such (no false pass).

## Inputs (strict)

1. `Z_12_v1` (typed group object)  
   `fundamental_action_reconstruction/generated/z_12_v1_group.json`
2. `Aut_Z12_v1` (typed automorphism group object)  
   `fundamental_action_reconstruction/generated/aut_z12_v1_units_group.json`
3. Boundary theorem:
   - `N462` (no `Aut(Z_12)`‑invariant canonical generator/orientation choice from typed structure alone)
4. Target spec:
   - `T164`

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/f473_current_strict_t164_z12_generator_orientation_canonical_fixing_datum_export_packet.py
```

Exports:

1. fixing datum (premise‑based strict provenance):
   - `fundamental_action_reconstruction/generated/kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json`
2. summary:
   - `fundamental_action_reconstruction/generated/f473_current_strict_t164_z12_generator_orientation_canonical_fixing_datum_export_packet_summary.json`

## Meaning (no false‑PASS)

This packet means only:

1. an explicit fixing datum is now exported as a strict provenance object (premise‑based),
2. the datum reduces the full `Aut(Z_12)` ambiguity by declaration to the subgroup preserving the fix (here: the identity),
3. downstream directed/sign‑sensitive constructions must cite this datum (or remain non‑strict).

## Hard limits

This packet does **not** claim:

1. `Aut(Z_12)`‑invariant canonicity (ruled out by `N462`),
2. discharge of `T163`/`T162`,
3. strict theta export or selector closure,
4. global discharge of `QW-2191`,
5. ToE closure.
