# P673 Current Strict Global Projective Selector Closure Object Gluing Probe

Status: `P673_EXECUTABLE_CURRENT_STRICT_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OBJECT_GLUING_PROBE_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit that the exported object:

```text
SelectorClosure_global_C_v1_projective_strict_v1
```

is **well-defined at projector/section level** on the declared global cover `{pair1..pair5}` of `C_v1`, meaning its
chartwise output projectors `B_out(pair_m)` agree within tolerance.

## Positive rule

The probe may return a positive result only if:

1. the closure object file exists and is typed as projective closure on `C_v1`,
2. the embedded well-definedness certificate passes:
   `max_abs_diff_across_chartwise_output_projectors <= tolerance`,
3. all referenced inputs exist (global atlas, transition, state, output operator, and chartwise `A_m` operators).

## Hard limits

Even a positive result does **not** imply:

- strict-core selector closure,
- global `QW-2191` discharge,
- ToE closure,
- any operator-level transition groupoid promotion beyond projector/section level (`N512` boundary).

