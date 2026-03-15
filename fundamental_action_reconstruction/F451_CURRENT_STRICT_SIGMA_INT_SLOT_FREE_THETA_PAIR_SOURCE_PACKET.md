# F451 Current Strict Sigma‑Int Slot‑Free Theta‑Pair Source Packet (No False‑PASS)

Status: `F451_EXECUTED_CURRENT_STRICT_SIGMA_INT_SLOT_FREE_THETA_PAIR_SOURCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

The strict sigma‑int → theta lane currently has two exposed selector slots in its **old** candidate class:

1. `eps ∈ [0,1]` (generator amplitude; `T117`),
2. `delta_d ∈ (0, delta_max]` (positive‑window corridor step; `T119`),

and `N441/N437` prove the theta output depends on those choices. Therefore strict‑core upgrade (`T159`) cannot be
obtained by “invariance under slots” on the current class (`N443`), and strict‑derived slot selection (`T160/T161`)
is not exported.

This packet executes the **construction‑class change** route (`T162`):

```text
export one slot-free strict theta-pair source:
  sigma_int_strict_derived_v1 -> (theta_1, theta_2),
using the already discharged strict Shannon selector ingredient class (T165 via F446/N480),
extended to pair2 by N488,
with no eps/delta_d parameters and no theta inputs.
```

This packet is about **theta supply only**. It does not claim:

1. admissible `S_sel_int` or strict‑core selector closure,
2. global `QW-2191` discharge,
3. ToE closure.

## Strict‑admissible inputs reused

1. `F307/N418` + `F308/N419`
   - `sigma_int_strict_derived_v1 ∈ {+1,-1}` with gauge‑quotient safety,
2. `F309/N420`
   - `alpha_geo_strict_derived_v1 := 4 ln 2` (strict‑derived Shannon amplitude),
3. `F329` / `QW-2190`
   - typed `Z_12` scaffold + real Fourier pairs `pair1, pair2`,
4. `F446/N480`
   - direction‑free element‑order reference `r_ord` + theorem‑level `pair1` minimizer set `θ = π/2 (mod π)`,
5. `N488`
   - theorem‑level `pair2` minimizer set `θ = 0 (mod π)` for the same reference/objective family.

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/f451_current_strict_sigma_int_slot_free_theta_pair_source_packet.py
```

Exports:

1. strict theta‑pair source artifact:
   - `fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json`
2. summary:
   - `fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1_summary.json`

## Meaning (no false‑PASS)

This packet means only:

1. the repo exports one **slot‑free** sigma‑int → theta‑pair source object (no `eps`, no `delta_d`, no theta inputs),
2. the `pair1` and `pair2` angles are fixed by theorem‑level Shannon objectives built from strict objects:
   `alpha_geo_strict_derived_v1` + `ord_Z12` (direction‑free),
3. the residual `Z2` sign on `pair1` is explicitly tracked via the existing sigma‑int sign convention
   (`F311` maps `sigma_int=-1` to `theta_1 -> theta_1 + π`), while the target‑slot object remains a span (sign‑free).

This packet does **not** claim:

1. discharge of `T160/T161` (strict‑derived eps/delta selection),
2. discharge of `T168/T169` (scalar→per‑site lift),
3. global selector closure, global `QW-2191` discharge, or ToE closure.

