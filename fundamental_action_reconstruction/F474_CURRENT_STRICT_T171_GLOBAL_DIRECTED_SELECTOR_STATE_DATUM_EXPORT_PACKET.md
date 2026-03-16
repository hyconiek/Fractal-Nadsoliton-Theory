# F474 Current Strict `T171` Global Directed Selector State Datum — Export Packet (No False‑PASS)

Status: `F474_EXECUTED_CURRENT_STRICT_T171_GLOBAL_DIRECTED_SELECTOR_STATE_DATUM_EXPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After:

1. `T170` discharge (`F469/N515`): global selector atlas + transition/gluing objects on `C_v1`,
2. `H39` discharge (`F470/N516`): global **projective** selector state object on `C_v1`,
3. `T164` discharge (`F473/N523`): an explicit premise‑based strict provenance fixing datum reducing `Aut(Z_12)` ambiguity,

the post‑projective frontier `T171` can be addressed honestly:

```text
export a sign-sensitive / directed selector state datum (vector-level lift) compatible with the already exported projective state.
```

This packet exports:

1. a strict directed sign‑distinction observable on `pair1`,
2. a strict global **directed** selector state object on `C_v1` (vector‑level lift),

while keeping all hard limits explicit (no implied selector closure; no `QW-2191` discharge; no ToE closure).

## Inputs (strict)

1. `alpha_geo_strict_derived_v1 := 4 ln 2` (`F309/N420`)
2. `Kappa_Z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1` (`F473/N523`)
3. global projective selector state object:
   - `SelectorState_global_C_v1_projective_strict_v1` (`F470/N516`)
4. convention-layer oriented vector section on `{pair1..pair5}`:
   - `U_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1` (`F467/N511`)

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/f474_current_strict_t171_global_directed_selector_state_datum_export_packet.py
```

Exports:

1. directed sign-distinction observable on `pair1`:
   - `fundamental_action_reconstruction/generated/s_dir_pair1_strict_v1.json`
2. global directed selector state object on `C_v1` (vector-level lift):
   - `fundamental_action_reconstruction/generated/selector_state_global_c_v1_directed_strict_v1.json`
3. summary:
   - `fundamental_action_reconstruction/generated/f474_current_strict_t171_global_directed_selector_state_datum_export_packet_summary.json`

## Meaning (no false‑PASS)

This packet means only:

1. residual `Z2` sign on the exported vector representative is now treated as *directed/physical* in the declared scope by citing an explicit fixing datum (`T164`),
2. the directed state descends to the already exported projective state (projectors/spans match),
3. no global selector closure or `QW-2191` discharge is implied.

## Hard limits

This packet does **not** claim:

1. `Aut(Z_12)`‑invariant canonicity (ruled out by `N462`; the fix is premise‑based),
2. strict‑core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.
