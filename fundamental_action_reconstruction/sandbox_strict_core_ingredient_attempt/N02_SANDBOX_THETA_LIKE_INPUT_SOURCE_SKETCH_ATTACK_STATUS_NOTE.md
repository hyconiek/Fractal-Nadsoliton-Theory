# N02 Sandbox Theta-Like Input Source Sketch Attack Status Note

Status: `N02_SANDBOX_THETA_LIKE_INPUT_SOURCE_SKETCH_ATTACK_STATUS_NOTE_NO_FALSE_PASS`
As of: `2026-03-09`

## Strongest honest conclusion

The sandbox attack on theta-like inputs does succeed in one narrow,
repo-consistent sense:

```text
rho_int_orientation_request_slot_v1 now becomes
rho_int_orientation_request_slot_v2,
which is not only target-slot-aligned
but also theta-source-sketch-aware
```

More precisely, the refined slot is now bound to:

1. the required-input pair `theta_1`, `theta_2` from `R1`,
2. the class-level orientation-slice candidate from `C47`,
3. the minimal basis-pair export skeleton from `C48`,
4. the exact `C35` split:
   strict-core actual theta absent,
   axiom-augmented theta branch present but non-strict.

## What did not happen

This attack did **not** derive:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure.

## Why the attack is still useful

It improves the sandbox in one practical way:

1. the orientation slot is now connected to the real phase-input frontier
   instead of only to an abstract target slot,
2. future sandbox work can now attack the exact place where strict core still
   fails:
   no strict-core actual theta source.

## Honest next move

If the sandbox is continued, the next clean move is:

1. either attack whether a non-placeholder strict-core theta-source skeleton
   can be written without importing the axiom-augmented branch,
2. or stop here and switch to an incompatibility boundary saying that the
   theta-source frontier itself does not presently enter strict core.
