# N711 Current Strict `T176` Parity-Completed Dual-Anchor All-Root Support but Global `Z2` Orbit-Split Boundary Theorem

Status: `N711_CURRENT_STRICT_T176_PARITY_COMPLETED_DUAL_ANCHOR_ALL_ROOT_SUPPORT_BUT_GLOBAL_Z2_ORBIT_SPLIT_BOUNDARY_THEOREM_NO_FALSE_PASS`
As of: `2026-03-17`

## Goal

After `P714/N710`, the next honest narrowing step is to test the smallest
parity-completed extension already present in the exported strict-core payload:

```text
w_break_by_x  +  w_ref_unnormalized_by_x
```

`N711` packages the strongest honest reading of that test.

## Theorem-level conclusion

From `P715`, the current repo supports the following joint statement:

1. the parity-completed dual-anchor rule does remove the old root-support gap:
   every chart root on the exported `pair1..pair5` lane admits a nonzero root
   anchor using the fixed slot-free rule
   `w_break first, w_ref_unnormalized fallback`;
2. the resulting rooted sections agree **projectively** across all roots:
   every supported root recovers the same section up to one global sign;
3. nevertheless, the result still does **not** yield an exact directed
   section:
   relative to the reference section, the `pair4` root produces the globally
   negated section, while `pair1`, `pair2`, `pair3`, and `pair5` recover the
   same exact branch.

Therefore the strongest honest current statement is:

```text
parity completion is enough to globalize root support and recover one
projective root-orbit,
but it is not enough to discharge the final directed-sign branch:
one residual global Z2 orbit split remains.
```

## Consequence

This narrows the next strict provider requirement further:

1. adding an even anchor component is **not** by itself sufficient,
2. the remaining missing ingredient must fix the last global `Z2` orbit split
   without smuggling a marked generator, orientation slot, or `T164`-style
   premise,
3. otherwise the honest continuation remains to freeze this residual orbit as
   gauge/convention and keep physics claims in projective/basis-invariant
   scope.

## What this theorem does not prove

`N711` does **not** prove:

1. `T176` discharge,
2. kernel-alone/global `QW-2191` discharge,
3. any directed/sign-sensitive physical orientation datum in strict core,
4. ToE closure.
