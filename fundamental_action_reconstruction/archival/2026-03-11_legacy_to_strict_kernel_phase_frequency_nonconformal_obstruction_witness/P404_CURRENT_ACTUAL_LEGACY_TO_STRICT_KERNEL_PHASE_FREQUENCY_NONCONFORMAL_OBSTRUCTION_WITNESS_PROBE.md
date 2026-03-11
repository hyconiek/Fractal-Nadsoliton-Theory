# (ARCHIVAL) P404 Legacy→Strict Kernel Phase/Frequency Nonconformal Obstruction Witness Probe (Draft)

Status: `P404_ARCHIVAL_SCRATCH_DRAFT_NO_FALSE_PASS`  
As of: `2026-03-11`

## Archival warning

This probe draft is preserved only for provenance in the archival lane.
It does **not** represent an active current-repo-state probe, because `S2`
freezes the legacy→strict bridge/non-bridge route after the author’s retirement
decree for `K_legacy_ont`.

## Intended goal (draft)

Check whether the archival draft packet computation (`F326_ARCHIVAL`) produces a
consistent “phase/frequency bridge absent” obstruction verdict on the current
export set, **without** promoting it into any active-lane discharge.

## Intended probe checks (draft)

| Check | Verdict (draft) |
|---|---|
| `F326_ARCHIVAL` output summary exists in this archival directory | required for a “YES” |
| `explicit_phase_frequency_bridge_present == false` (from `P47` summary) | required |
| `phase_frequency_layer_obstructed_on_current_exports == true` | required |

## Implementation note (archival-only)

The accompanying script `p404_current_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_probe.py`
evaluates the local archival `generated/` outputs and writes its probe report
into the same archival directory.

## Hard limits

This archival probe draft must not be cited as:

1. an active executed probe,
2. a strict-core discharge,
3. a revival of the frozen legacy→strict bridge/non-bridge frontier,
4. a selector closure / `QW-2191` discharge / ToE closure claim.

