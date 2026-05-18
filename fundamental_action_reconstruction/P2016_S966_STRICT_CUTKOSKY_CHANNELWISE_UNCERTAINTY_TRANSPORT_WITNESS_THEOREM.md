# P2016 S966 Strict Cutkosky Channelwise Uncertainty Transport Witness Theorem

Status: `P2016_EXECUTED_CHANNELWISE_UNCERTAINTY_TRANSPORT_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

Strengthen strict-lane unitarity diagnostics by moving from whole-cut multiplicative
uncertainty (P2015) to channel-resolved uncertainty transport on `gg`, `gh`, `hh`.

## Construction

For each `s` sample from `P1997`, keep `ImM(s)` fixed and scan independent channel
perturbations:

- `eps_gg ∈ [-0.05, 0.05]`,
- `eps_gh ∈ [-0.06, 0.06]`,
- `eps_hh ∈ [-0.08, 0.08]`.

For every triple, compute transported cut:

`Cut_transported = gg*(1+eps_gg) + gh*(1+eps_gh) + hh*(1+eps_hh)`

and record `Delta_opt = ImM - Cut_transported`, interval statistics, and channelwise
positivity flags.

## Result

`P2016` exports a channel-resolved uncertainty table with gate checks for:

1. upstream witness availability,
2. table completeness,
3. residue positivity under the scan,
4. bounded max-abs delta in scanned window.

No global closure claim is made.

## Next honest step

Replace proxy channel perturbations by explicit loop-derived channel amplitudes
from backend integrals and rerun transport on a wider energy/RG atlas.
