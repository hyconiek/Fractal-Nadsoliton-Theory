# P702 Current Non‑Strict Standard‑Model Host Matching From P696 Channel Proxy Probe (No False‑PASS)

Status: `P702_CURRENT_NONSTRICT_STANDARD_MODEL_HOST_MATCHING_FROM_P696_CHANNEL_PROXY_PROBE_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Provide a **non‑strict** (externally‑seeded) host‑matching harness that compares the exported strict projective channel proxy spectrum (`P696`)
to an explicitly provided Standard‑Model mass target dataset.

This probe is intentionally **not strict** because it depends on external target data and on an explicit identification policy
(mapping channels to particle masses and choosing a global unit scale).

The probe computes:

1. a best-fit **bijection** between the selected P696 channels and the provided target list (assignment on centered log‑mass²),
2. the resulting best-fit **global scale** from proxy units to the external target units,
3. diagnostic metrics (log‑RMSE, max relative error),
4. an optional pass/fail under **user-supplied explicit thresholds** (if provided in the dataset file).

If the external target dataset file is missing, the probe returns `NOT_COMPUTABLE` and points to the template file.

## Hard limits

This probe does **not** claim:

1. any strict Standard Model identification or host matching discharge,
2. any directed/sign-sensitive physical orientation datum in strict core,
3. kernel-alone/global `QW-2191` discharge,
4. strict-core ToE closure or global ToE closure.

