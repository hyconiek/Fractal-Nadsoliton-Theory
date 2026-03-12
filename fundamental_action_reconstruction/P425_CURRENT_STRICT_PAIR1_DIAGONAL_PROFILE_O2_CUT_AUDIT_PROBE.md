# P425 Current Strict Pair1 Diagonal Profile O(2)-Cut Audit Probe

Status: `P425_EXECUTED_CURRENT_STRICT_PAIR1_DIAGONAL_PROFILE_O2_CUT_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `N465`, the strict frozen host `A = K_total + m0^2 I` cannot break `O(2)` on `pair1`.

The next honest mathematical question is:

```text
what exactly must a diagonal/local sector contribute in order to break O(2) on pair1?
```

`P425` is a pure audit probe that numerically verifies the closed-form criterion stated in `N466`
on a few explicit diagonal profiles on the strict 12-site scaffold.

This probe does **not** claim any of those profiles is physically realized.

## Strict-admissible evidence reused

1. `QW-2118`
   - fixes the strict `n=12` scaffold size for the current kernel-mode lane.

No other physics input is used; the probe is purely linear-algebraic.

## Computation artifact

This probe is executed by:

- `fundamental_action_reconstruction/p425_current_strict_pair1_diagonal_profile_o2_cut_audit_probe.py`

and it writes:

- `fundamental_action_reconstruction/generated/p425_current_strict_pair1_diagonal_profile_o2_cut_audit_probe.json`
- `fundamental_action_reconstruction/generated/p425_current_strict_pair1_diagonal_profile_o2_cut_audit_probe_summary.json`

## Verdict

On the current repo state, the numerical audit matches the closed-form formula from `N466` within
tolerance on the tested profiles:

1. constant diagonal profiles are isotropic on `pair1`,
2. profiles with a nonzero `mode-2` Fourier component produce `Δ1 ≠ 0` and therefore break `O(2)` on `pair1`.

## Hard limits (no false pass)

`P425` does not:

1. export any physically derived diagonal profile,
2. discharge any direct-family witness (`g4/g6/gY/m2`),
3. discharge `QW-2191`,
4. export a strict-core selector ingredient,
5. claim ToE closure.

