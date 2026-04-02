# P802 Current Strict Alpha_s Provider-Class Reorganization Audit Probe

Status: `P802_CURRENT_STRICT_ALPHA_S_PASSIVE_PROVIDER_SKELETON_SUPPORTED_ACTIVE_PROVIDER_ACTION_RULE_BLOCKED`
As of: `2026-03-19`

## Goal

After `F801`, the next honest question is:

```text
can the current same-domain provider skeleton
already be reorganized
into a real alpha_s provider class
without importing foreign-domain reference structures
or new host semantics?
```

## Scope

`P802` does not build that provider class.
It only separates two claims:

1. passive same-domain provider skeleton,
2. active provider action rule.

## Main Checks

1. confirm `F801` already exports the passive same-domain skeleton,
2. confirm the same-domain lane already supplies the passive roles needed under a future provider class,
3. test whether the same-domain lane already exports any active provider action rule or supply schema,
4. keep foreign-domain reference structures and non-strict calibration excluded from the reorganization.

## Result

`P802` finds another asymmetric split:

1. the passive same-domain provider skeleton is **supported/exported**,
2. the active provider action rule remains **blocked / nonexport**.

So the blocker narrows again:

```text
not "do we have a same-domain provider base",
but "what active rule makes that base act as a provider class"
```

## Hard Limit

`P802` does not claim that passive support already implies
an active provider class.
It only blocks silent promotion from the first to the second.
