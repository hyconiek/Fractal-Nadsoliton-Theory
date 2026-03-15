# P439 Current Strict QW‑2191 Weighted KL‑to‑Reference Objective O(2)‑Cut Audit Probe

Status: `P439_EXECUTED_CURRENT_STRICT_QW2191_WEIGHTED_KL_REFERENCE_OBJECTIVE_O2_CUT_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

After `N463/N464`, a strict Shannon/KL objective built only from squared site‑amplitude distributions
and **permutation/translation invariant** structure cannot yield a unique `O(2)` cut.

`P439` audits one concrete *escape hatch* explicitly permitted by `N464`:

```text
introduce an explicit non‑translation‑invariant reference distribution r(x)
on the 12-site ring (an anchored or group-structure-derived datum), then test KL(p_theta || r)
for uniqueness on the QW-2191 O(2) family.
```

This probe answers only:

1. does such a **non‑permutation‑invariant** KL-to-reference objective vary with `theta` (nontriviality), and
2. does it yield a unique global minimizer cluster on a dense audit grid (candidate `O(2)` cut),

for several **explicit reference distributions** `r(x)` constructed from strict-side artifacts (e.g. `R14` kernel matrix row)
and from intrinsic `Z_12` group structure (e.g. element-order weights).

`P439` does **not** claim that any such reference is already an admissible strict-core selector ingredient.
It is an audit probe to decide whether the *objective-shape class* is even viable before attempting a theorem-level
`T165` discharge.

## Strict-admissible evidence reused

1. `QW-2190`
   - fixed 12‑octave ring scaffold and real Fourier mode basis.
2. `QW-2191`
   - continuous `O(2)` rotation family in degenerate two‑mode subspaces.
3. `N463/N464`
   - theorem-level boundaries for permutation-invariant objectives; and explicit permission that a marked-site/non-translation-invariant datum
     would be required for uniqueness.
4. `R14`
   - frozen strict `K_total` specialization packet (used only as a strict-sourced **reference-shape** provider for `r(x)`).

## What is tested (objective family)

Let `u_1,theta := cos(theta)c1 + sin(theta)s1` and define the squared site distribution `p_theta(x)=u_1,theta(x)^2`.

For each declared reference distribution `r(x)` on sites `x=0..11`, test objectives such as:

```text
J(theta) := KL( p_theta || r )
J_avg(theta) := 0.5 * ( KL(p_theta||r) + KL(q_theta||r) )   where q_theta is the orthogonal rotated partner v_1,theta squared
```

The tested references include non-uniform vectors derived from strict artifacts (e.g. the absolute `K_total[0,x]` profile from `R14`),
normalized into a probability distribution.

## Computation artifact

This probe is executed by:

- `fundamental_action_reconstruction/p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.py`

and it writes:

- `fundamental_action_reconstruction/generated/p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.json`
- `fundamental_action_reconstruction/generated/p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe_summary.json`

## Hard limits (no false pass)

`P439` does **not**:

1. export a strict-core selector ingredient,
2. discharge `T165`,
3. discharge `QW-2191`,
4. claim any marked-site/reference datum is already strict-admissible as physics,
5. claim strict-core theta export or ToE closure.
