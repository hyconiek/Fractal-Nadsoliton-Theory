# P424 Current Strict QW‑2118 Host Pair1 Block Isotropy Audit Probe

Status: `P424_EXECUTED_CURRENT_STRICT_QW2118_HOST_PAIR1_BLOCK_ISOTROPY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Take one concrete, strict-admissible step toward the “physical accelerator / symmetry breaking” discussion
without false pass:

```text
compute the pair1 (c1,s1) 2×2 restriction of the certified host operator A = K_total + m0^2 I,
and check whether that restriction can itself supply an O(2)-cut (anisotropy) on pair1.
```

This probe does **not** claim any `K_obs`-lane factorization or selector closure. It audits only the
host-side symmetry constraint.

## Strict-admissible evidence reused

1. `QW-2118`
   - frozen `12`-octave cyclic-distance `K_total` distance profile.
2. `QW-2124`
   - branch-resolved scalar vacuum floor `m0^2` (broken-floor).
3. `R8`
   - host scope schema: `A = K_total + m0^2 I` on the declared `12`-slot carrier.

For comparison only (not used as a strict source), the probe also reads the already exported computed
current-pair block from:

- `P10` (`generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain.json`).

## Computation artifact

This probe is executed by:

- `fundamental_action_reconstruction/p424_current_strict_qw2118_host_pair1_block_isotropy_audit_probe.py`

and it writes:

- `fundamental_action_reconstruction/generated/p424_current_strict_qw2118_host_pair1_block_isotropy_audit_probe.json`
- `fundamental_action_reconstruction/generated/p424_current_strict_qw2118_host_pair1_block_isotropy_audit_probe_summary.json`

## Exact verdict (current repo state)

On the current frozen host (`A = K_total + m0^2 I`), the restriction to `pair1 = span{c1,s1}` is
numerically isotropic (scalar on pair1) within the declared tolerance.

Therefore the host operator alone does **not** supply an `O(2)`-cut on `pair1`.

In particular, the anisotropic current-pair `A_1^(current)` block exported in `P10` cannot equal the
host `pair1` restriction without importing additional symmetry-breaking structure.

## Hard limits (no false pass)

`P424` does not:

1. identify `P10`’s current-pair chain as a strict reduction of existing kernel feedback,
2. discharge `QW-2191`,
3. export a strict-core selector ingredient,
4. claim strict-core theta export,
5. claim ToE closure.

