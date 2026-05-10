# Experimental Unfreeze Direction Sensitivity Boundary

Status: `NONSTRICT_PARAMETER_UNFREEZE_DIRECTION_SENSITIVITY_BOUNDARY`
As of: `2026-03-27`

## Goal

Test the narrowest honest question opened by the experimental unfreeze lane:

```text
can parameter motion in (omega, phi, beta, eta) inside the current strict
distance-only kernel family, by itself, generate the inversion-sensitive
source-side rule required by F960/T183?
```

## Input State

This note is anchored to four already-frozen repo facts:

1. `K1` keeps the strict working family as

```text
K_strict_gate(d) := cos(omega*d + phi) / (1 + beta*d^eta)
```

that is, a family whose explicit argument is still just `d`.
2. `P729/N725` already localize the surviving strict frontier split as the
   opposite residual-datum branches `delta_k` versus `delta_-k`.
3. `P730/N726` already prove that the current direction-free lane fails
   exactly because inversion partners stay score-equal under the present
   source-side inputs.
4. `F960` already freezes the active missing object as
   `ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1`,
   and its own recommended next move asks for a genuinely new
   inversion-sensitive source-side provider class.

## Honest Test

The experimental lane may move `(omega, phi, beta, eta)`, but that motion is
still confined to the same family:

```text
distance-only input -> K_strict_gate(d; omega, phi, beta, eta)
```

Changing the tuple can retune amplitudes over `d`, but it does not by itself:

1. add a signed orbit-direction variable,
2. add a branch tag distinguishing `k` from `-k`,
3. add a non-distance source-side datum,
4. or change the provider class from direction-erasing to
   inversion-sensitive.

## Conclusion

The strongest honest current reading is:

```text
parameter unfreeze inside the present K_strict_gate(d) family can be used as
an exploratory pressure test, but it cannot by itself export the exact F960/T183
residual-datum orbit-direction selection bridge
```

Therefore this lane remains admissible only for:

1. `bridge-pressure only`,
2. `fit-dependence audit`,
3. `operational-only retune`,
4. or other explicitly non-strict heuristic uses.

It is not admissible to treat parameter motion alone as:

1. `T183` discharge,
2. `T176` discharge,
3. `QW-2191` closure,
4. or a strict physical orientation selector.

## Next Consequence

If this experimental lane continues, the honest question is no longer:

```text
can one more retune of omega/phi/beta/eta solve the strict selector problem?
```

The honest question becomes:

```text
does parameter motion expose bridge-pressure, fit-dependence, or a clue toward
a genuinely new inversion-sensitive provider class beyond the present
distance-only kernel family?
```
