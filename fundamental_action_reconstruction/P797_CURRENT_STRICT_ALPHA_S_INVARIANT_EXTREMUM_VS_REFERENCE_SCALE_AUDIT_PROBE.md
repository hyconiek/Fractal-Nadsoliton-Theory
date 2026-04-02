# P797 Current Strict Alpha_s Invariant Extremum vs Reference-Scale Audit Probe

Status: `P797_CURRENT_STRICT_ALPHA_S_INVARIANT_EXTREMUM_CANDIDATE_SUPPORTED_REFERENCE_SCALE_RULE_BLOCKED`
As of: `2026-03-19`

## Goal

After `P796/F796`, the next honest question is:

```text
does the current strict repo already support
the F704 maximum as an invariant extremum candidate,
even though it still does not support it
as a reference-scale point?
```

## Scope

`P797` does not promote the F704 maximum.
It only separates two claims:

1. invariantly defined numeric extremum candidate,
2. reference-scale point.

## Main Checks

1. confirm the F704 maximum is well-defined on the exported basis-invariant
   spectrum,
2. confirm current Release-7 invariance theorems keep that extremum stable
   under the admitted gauge/convention layer,
3. confirm current meaning theorems still stop before reference-scale
   semantics,
4. identify the sharper remaining blocker.

## Result

`P797` finds another asymmetric split:

1. the F704 maximum is **candidate-supported** as an invariant numeric
   extremum,
2. the upgrade to **reference-scale point** remains **blocked / nonexport**.

So the blocker narrows again:

```text
not "is the F704 maximum stable",
but "what rule upgrades that stable extremum into a reference-scale point"
```

## Hard Limit

`P797` does not claim that invariant extremum support already implies
reference-scale semantics.
It only blocks silent promotion from the first to the second.
