# P717 Current Strict `T176` `pair4` Exact Branch-Split Release-7 OS Gauge-Irrelevance Bridge Audit Probe (No False-PASS)

Status: `P717_CURRENT_STRICT_T176_PAIR4_EXACT_BRANCH_SPLIT_RELEASE_7_OS_GAUGE_IRRELEVANCE_BRIDGE_AUDIT_PROBE_SPEC_NO_FALSE_PASS`
As of: `2026-03-18`

## Goal

After `P716`, the next honest continuation is not to pretend that the strict-core
branch problem is solved.

The sharper question is:

```text
does the currently localized exact branch split still matter for the concrete
Release-7 OS observables actually used downstream?
```

This probe audits that bridge, and nothing more.

## Inputs

1. `generated/p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json`
   - exact/projective rooted-branch structure after parity completion.
2. `generated/p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe_summary.json`
   - localization of the remaining exact split to `pair4`.
3. `generated/p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json`
   - Release-7 OS gauge-irrelevance audit over all residual sign patterns.

## Audit rule

Test whether the exact two-branch situation isolated by `P715/P716`

```text
reference exact branch
vs
pair4-root exact branch = global negation of the reference branch
```

is already covered by the stronger existing Release-7 OS gauge-irrelevance
result `P709`, which audits all residual sign patterns for:

- `P694` mass proxy,
- `P696` selector-aligned channel spectrum proxy invariants,
- `F704` basis-invariant mass observable.

## Acceptance

The probe passes only if it can show all of:

1. `P715` still exhibits the exact branch split,
2. `P716` still localizes that split to the `pair4` negative cosine-axis role,
3. `P709` still passes sign-gauge irrelevance for the full sign-pattern family,
4. therefore the current exact `pair4` branch split is gauge-irrelevant for
   the concrete Release-7 OS observables already exported downstream.

## Hard limits

This probe does **not**:

- claim a strict-core directed/sign-sensitive physical orientation datum,
- claim `T176` discharge,
- claim kernel-alone/global `QW-2191` discharge,
- claim gauge-irrelevance for arbitrary future observables,
- claim ToE closure.
