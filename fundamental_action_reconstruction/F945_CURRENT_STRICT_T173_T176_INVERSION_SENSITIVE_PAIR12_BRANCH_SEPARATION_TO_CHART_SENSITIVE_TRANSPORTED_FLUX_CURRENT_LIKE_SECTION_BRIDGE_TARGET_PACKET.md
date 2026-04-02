# F945 Current Strict `T173/T176` Inversion-Sensitive Pair1/Pair2 Branch-Separation to Chart-Sensitive Transported Flux/Current-Like Section Bridge Target Packet

Status: `F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS`
As of: `2026-03-21`

## Goal

Freeze the exact missing bridge object now exposed on the active `QW-2191`
frontier:

```text
the current repo already exports inversion-sensitive pair1/pair2 branch
separation,
but still does not export the bridge from that separation
to a chart-sensitive transported flux/current-like section on full C_v1
```

## Why This Packet Exists

1. `P945` shows that the current repo already exports a real
   inversion-sensitive `pair1/pair2` branch separation.
2. `P721/P722` keep the `T176/T177` chart-sensitive provider bridge
   explicitly unexported.
3. The honest next move is therefore to freeze the exact bridge target,
   not to pretend that `T176`, `T177`, or `T185` are already discharged.

## Exported Target

```text
inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1
```

## Hard Limits

`F945` does not claim:

1. that `T176` is discharged,
2. that `T177` is discharged,
3. that `T185` is discharged,
4. that `F647` is admissible `S_sel_int`,
5. that the rooted `w_break` sign lift is a strict physical orientation datum,
6. kernel-alone/global `QW-2191` discharge,
7. ToE closure.
