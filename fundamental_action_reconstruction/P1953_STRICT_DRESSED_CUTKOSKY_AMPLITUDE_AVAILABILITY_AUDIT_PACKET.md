# P1953 Strict Dressed Cutkosky Amplitude Availability Audit Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__DRESSED_AMPLITUDE_INTERFACE_MISSING_NO_FALSE_PASS`
As of: `2026-05-17`

## Goal

After `P1951/P1952`, the honest next question is:

```text
Can K_cut_seed be replaced by a full dressed graviton->gauge_gauge amplitude?
```

`P1953` audits the current repository for the minimum theorem-grade Cutkosky
inputs.

## Result

The local verdict is:

```text
DRESSED_CUTKOSKY_AMPLITUDE_NOT_AVAILABLE__INTERFACE_CONTRACT_EXPORTED
```

The current repo contains:

1. seed/proxy Cutkosky data,
2. seed physical-pole projection,
3. common-basis placeholders,
4. open BRST/Cutkosky contracts.

It does not contain the full dressed amplitude object required for theorem
closure.

## Blocking Missing Objects

The audit blocks theorem promotion because the following objects are missing
or only seed-level:

1. `|M_dressed(graviton->gauge_gauge)|^2` in a common basis,
2. same-scheme `DiscM_common_basis` and `CutSum_common_basis`,
3. computed dressed propagator poles and residues,
4. BRST physical-state projector for the channel,
5. a locked common scheme between `MSbar_B1_seed` and later common-basis rows.

## Exported Interface Contract

`P1953` exports the required interface:

```text
DressedCutkoskyAmplitude_graviton_to_gauge_gauge_strict_B1_v1
```

with required fields including:

```text
M_dressed_common_basis
AbsM_dressed_squared_common_basis
DiscM_common_basis
CutSum_common_basis
DiscM_minus_CutSum_simplified
dressed_graviton_propagator_pole_list
residue_values_per_pole
ghost_sector_exclusion_trace
BRST physical-state projector
```

## No False Pass

`P1953` does not weaken `P1950/P1951/P1952`. Those remain valid in their local
scope. It only prevents their seed-level outputs from being promoted into a
global Cutkosky theorem.

## Outputs

- `p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.py`
- `generated/p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json`

## Next Honest Step

Build `P1954` with high-reasoning support: derive the same-scheme dressed
amplitude from `L_total` and exported vertices/projectors, or export a formal
nonavailability theorem identifying the exact missing vertex/projector data.
