# P646 Current Strict Witness‑Provider Scan Probe (Seed‑v1)

Status: `P646_EXECUTABLE_CURRENT_STRICT_WITNESS_PROVIDER_SCAN_PROBE_FOR_SEED_V1_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Mechanically scan `generated/*.json` for any export matching the `F646`
witness‑provider signature (seed‑v1).

## Allowed conclusions (no false pass)

This probe supports exactly the following current‑repo‑state conclusions:

```text
either:  no strict witness provider matching the F646 signature is exported (seed v1)
or:      one-or-more strict witness providers matching the F646 signature are exported (seed v1)
```
