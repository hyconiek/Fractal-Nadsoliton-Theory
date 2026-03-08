# P207 — Current Next Actual Emergent Observer Closure Readout Operator Probe

## Goal

Check whether the repo exports one admissible next downstream readout operator
from `AI_obs_actual_closure_support_preLM_v1`.

## Expected positive condition

The probe passes if:

1. `AI_obs_actual_closure_support_preLM_v1` is already admissible,
2. `AJ_obs_actual_closure_readout_preLM_v1` is derived only from that map,
3. the new readout operator remains strict-core only,
4. it stays downstream-only,
5. the commit amplitude is positive,
6. the gap channel stays zero.

## Non-claim

Passing this probe still does **not** mean:

- actual emergent-observer closure,
- `QW-2191` discharge,
- strict-core selector closure,
- final ToE closure.
