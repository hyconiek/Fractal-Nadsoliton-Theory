# P663 Current Actual Legacy-to-Strict Kernel Bifurcated Frontier Probe (v2, post-`N554`)

Status: `P663_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PROBE_V2_NO_FALSE_PASS`  
As of: `2026-03-17`

## Purpose

Audit that the repo exports the current `F663` frontier status packet summary
and that it reflects the post-`N554` state:

- bridge branch explicit + future-only,
- nonbridge target explicit + **not** future-only (actual discharge present),
- no branch-selection claim.

## Inputs

- `generated/f663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_packet_v2_summary.json`

## Output

- `generated/p663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_probe_v2_summary.json`

