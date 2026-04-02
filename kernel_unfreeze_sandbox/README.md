# Kernel Unfreeze Sandbox

This directory contains a deliberately non-strict experimental lane for
unfreezing the current strict working tuple

```text
(omega, phi, beta, eta) = (0.18575, 0.16250, 1.0, 1.8)
```

used by the later operational strict kernel.

Why this sandbox exists:

1. `QW-2038` is an explicit refreeze scan.
2. `QW-2039` is an explicit refrozen gate.
3. `QW-2041` confirms canonical/refrozen semantic drift.
4. `QW-2049` selects the current working gate tuple.
5. `ANALIZA_FITTINGU_I_TRIKOW_KOMPENSACYJNYCH.md` records a broader project
   history in which fitting and compensation tricks must be tracked openly.

Therefore, unfreezing these parameters is admissible only as:

```text
fit-derived exploratory lane
```

and not as a new strict theorem.

This sandbox does not claim:

1. `T173` discharge,
2. `T176` export,
3. `QW-2191` selector closure,
4. a bridge between `K_legacy_ont` and `K_strict_gate`,
5. inheritance of legacy physical-role claims onto the strict kernel.

Files:

1. `EXPERIMENTAL_UNFREEZE_OMEGA_PHI_BETA_ETA_LANE.md`
   - human-readable experimental scope and stop rules,
2. `unfreeze_omega_phi_beta_eta_probe.py`
   - parses provenance reports and emits a machine-readable manifest,
3. `EXPERIMENTAL_UNFREEZE_DIRECTION_SENSITIVITY_BOUNDARY.md`
   - explains why parameter motion inside the same distance-only kernel family
     cannot by itself export the active `F960/T183` inversion-sensitive bridge,
4. `unfreeze_direction_sensitivity_boundary_probe.py`
   - checks the structural boundary against `K1`, `P729`, `P730/N726`, and
     `F960`,
5. `generated/unfreeze_omega_phi_beta_eta_manifest.json`
   - provenance probe output,
6. `generated/unfreeze_direction_sensitivity_boundary_summary.json`
   - boundary probe output.
