# P250 Current Actual Nonstrict Declared-Scope Selector Closure Probe

Status: `P250_EXECUTED_CURRENT_ACTUAL_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Test whether the current repo now exports one actual non-strict
declared-scope selector-closure witness without pretending strict or global
closure.

## Input

`P250` reads:

1. `N125`,
2. `N258`,
3. `N269`,
4. `F160`.

## Probe question

Does the current repo export:

```text
C_sel_declared_scope_nonstrict_actual_witness_v1 :
  tau_src_candidate_v1
    -> axiom_augmented_declared_scope_selector_closure_target_v1
```

such that:

1. declared-scope source-topology selector theorem is real,
2. selector requirement is explicitly accepted in `axiom_augmented_only`
   scope,
3. bridge/nonbridge is not treated as a mandatory gate,
4. strict-core closure is still not claimed,
5. global closure is still not claimed?

## Expected outcome

If the packet is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_WITNESS_AFTER_P250
```

## Hard limits

Passing `P250` does not mean:

1. strict-core selector closure is proved,
2. global selector closure is proved,
3. global `QW-2191` discharge is proved,
4. ToE closure is proved.
