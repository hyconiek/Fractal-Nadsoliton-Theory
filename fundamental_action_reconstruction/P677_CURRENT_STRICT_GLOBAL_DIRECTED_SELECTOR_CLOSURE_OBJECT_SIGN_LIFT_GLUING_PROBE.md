# P677 Current Strict Global Directed Selector Closure Object Sign-Lift Gluing Probe

Status: `P677_CURRENT_STRICT_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OBJECT_SIGN_LIFT_GLUING_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit that the repo exports one explicit **directed** global selector closure object on `C_v1`:

```text
SelectorClosure_global_C_v1_directed_strict_v1
```

and that it is well-defined **only** after applying the explicitly exported per-chart sign-lift/section choice (no hidden sign fixing).

## Acceptance (no false pass)

The probe must:

1. validate the object exists and is marked `no_false_pass`,
2. validate the closure scope is `directed_vector_state`,
3. validate the sign-lift payload is present and covers `{pair1..pair5}` with values in `{±1}`,
4. validate the object contains an explicit well-definedness certificate and that it passes,
5. keep hard limits explicit (no strict-core selector closure claim, no kernel-alone/global `QW-2191` discharge claim, no operator-level groupoid promotion).

