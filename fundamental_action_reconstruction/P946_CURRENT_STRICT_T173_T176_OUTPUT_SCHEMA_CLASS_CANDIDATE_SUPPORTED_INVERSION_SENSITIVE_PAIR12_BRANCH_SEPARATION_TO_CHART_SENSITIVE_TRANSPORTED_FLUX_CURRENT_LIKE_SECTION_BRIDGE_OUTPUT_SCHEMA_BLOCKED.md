# P946 Current Strict `T173/T176` Output-Schema Class Candidate Supported Inversion-Sensitive Pair1/Pair2 Branch-Separation to Chart-Sensitive Transported Flux/Current-Like Section Bridge Output Schema Blocked

Status: `P946_CURRENT_STRICT_T173_T176_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_BLOCKED`
As of: `2026-03-21`

## Goal

After `F945`, the next honest question is:

```text
does the repo already export
one exact bridge_output_schema
for the new inversion-sensitive pair1/pair2 to chart-sensitive transported
flux/current-like section bridge,
or does it only preserve neighboring output-schema class
at candidate grade around that lane?
```

## Scope

`P946` does not export the bridge output.
It only audits whether the frozen `F945` bridge target already has one exact
output-schema object, or whether the repo still preserves only neighboring
partial bridge classes, local chart-sensitive atlas classes, and pair1/pair2
branch-sensitive support around that lane.

## Main Checks

1. confirm `F945` already names `bridge_output_schema` as one exact missing
   field,
2. confirm `P945` keeps the bridge itself unexported while preserving the
   already-exported inversion-sensitive `pair1/pair2` branch separation,
3. confirm `P721/P722` preserve a physics-facing but chart-blind source-side
   provider class and the exact chart-sensitive transported section target
   class,
4. confirm `P742/P747` preserve neighboring partial bridge classes on the
   pair1/pair2 typed-carrier side and on the local chart-sensitive atlas side,
5. confirm `P729/P731` preserve exact `pair1/pair2 = delta_k / delta_-k`
   branch sensitivity together with antisymmetric witness-side separation,
6. confirm `F647/P684` keep the lane below admissible `S_sel_int` and below
   any strict physical orientation datum,
7. confirm no exact new `bridge_output_schema` object is exported on the
   current repo state.

## Result

`P946` shows:

```text
bridge_output_schema class is candidate-supported,
but the exact bridge_output_schema itself remains blocked
```

So the sharp blocker is now the exact output-schema object required by `F945`.

## Hard Limit

`P946` does not claim:

1. that the exact `bridge_output_schema` already exists,
2. that neighboring partial bridge classes silently discharge the frozen
   `F945` bridge target,
3. that `F647` becomes admissible `S_sel_int`,
4. that the rooted `w_break` sign lift becomes a strict physical orientation
   datum,
5. that `T176`, `T177`, or `T185` are discharged,
6. kernel-alone/global `QW-2191` discharge,
7. ToE closure.
