# P2623/S1573 Wilson-loop flux orientation-source boundary

Status: `P2623_WILSON_FLUX_ORIENTATION_SOURCE_BOUNDARY_NO_SELECTOR_SOURCE_NO_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE`

## Content-first anti-duplication grep audit

Mode: `content-first semantic patterns, not ticket/name lookup`.
- `closed_loop_flux_holonomy_content`: 1015 hits; samples retained in JSON certificate.
- `gauge_link_transport_content`: 793 hits; samples retained in JSON certificate.
- `orientation_reversal_content`: 137 hits; samples retained in JSON certificate.
- `selector_source_content`: 555 hits; samples retained in JSON certificate.
- `guardrail_content`: 12568 hits; samples retained in JSON certificate.

## Theorem export

**Claim.** A gauge-invariant Wilson loop can supply an orientation-odd sign only after three typed inputs are present: a gauge-safe connection on a closed cycle, nonzero flux with Im W != 0, and an independently sourced orientation of that cycle.  The gauge-invariant but unoriented datum {W,W^{-1}} is sign-blind, so Wilson/holonomy content alone does not repair the P2620 orientation atom.

Positive retained content:
- Closed-cycle Wilson products are gauge-invariant by telescoping link gauge shifts.
- Reversing the oriented cycle conjugates the Wilson loop and flips sign(Im W) when the flux is not 0 or pi.
- If a future strict theorem supplies gauge-safe connection, nonzero flux, and cycle orientation, sigma=sign(Im W) is a mathematically valid orientation-odd selector candidate.

Obstructions:
- Without a cycle-orientation source the physical datum is the conjugacy pair {W,W^{-1}}, which contains both signs.
- A Wilson loop can be a gauge-invariant flux diagnostic, but it cannot create the missing orientation by itself.
- No nonlinear damping completion source is supplied by this selector-side analysis.

## Computational certificate

- Gauge-invariance defect: `4.965068306494546e-16`.
- Reversal conjugacy defect: `0.0`.
- Oriented sign flips under reversal: `True`.
- Unoriented conjugacy class contains both signs: `True`.
- Orientation-source lattice accepting rows: `1` of 8.
- P2620 accepting rows now: `0`.

## Next admissible targets
- search for an internal strict source of oriented closed cycles rather than another holonomy value
- derive a gauge-safe connection/parallel transport from strict nadsoliton field content
- keep nonlinear damping completion as a separate bridge atom

Not licensed:
- unconditional orientation_odd_selector_source
- promotion of unoriented Wilson-loop data to strict selector source
- P2620 bridge-source cut repair
- role-bearing L_total
- QW-2191 discharge
- ToE closure

Certificate hash: `ded5f83438c06618ca14dc093c59e6c2582816abf613196043b091ea12690f42`
