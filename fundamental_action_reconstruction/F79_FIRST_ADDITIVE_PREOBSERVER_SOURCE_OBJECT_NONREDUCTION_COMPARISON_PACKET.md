# F79 First Additive Preobserver Source Object Nonreduction Comparison Packet

Status: `F79_EXECUTED_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_NONREDUCTION_COMPARISON_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N185`, the next honest question is:

```text
is S_preLM_additive_candidate_v1 merely a direct repackaging
of the older F75 target profile, or is it already nontrivially different?
```

This packet answers only that comparison question. It does not claim that the
first admissibility clause is already satisfied.

## Fixed same-basis realization

Work in the same ordered basis used by `F76`:

```text
u_T, u_L, u_M
```

with:

```text
u_T := e_topo
```

Freeze the first matter slice:

```text
d_* := 1
I_mat(d_*) = 1/2
```

Now realize the `F75` packaged target in this same basis by the canonical
reading:

```text
P_NL^(0) u_T := u_L
P_LM(d_*) u_L := I_mat(d_*) u_M
```

so that:

```text
S_preLM_target_packaged_realization_v1(d_*)
  := cos(phi) u_T + cos(phi) u_L + (cos(phi)/2) u_M
```

## Additive attempt reused

Reuse the additive attempt from `F76`:

```text
S_preLM_additive_candidate_v1
  = u_T + cos(phi) u_L + (cos(phi)/4) u_M
```

## Nonreduction witness

Define the same-basis difference:

```text
Delta_preLM_nonreduction_v1
  := S_preLM_additive_candidate_v1
   - S_preLM_target_packaged_realization_v1(d_*)
```

Hence:

```text
Delta_preLM_nonreduction_v1
  = (1 - cos(phi)) u_T - (cos(phi)/4) u_M
```

This is the first explicit nonreduction witness on the preobserver lane.

## Meaning

If `Delta_preLM_nonreduction_v1 != 0`, then `S_preLM_additive_candidate_v1`
is not merely the direct tuple packaging from `F75` under the same basis
realization and same fixed slice `d_* = 1`.

That still does **not** imply:

- genuinely-new strict-core source object,
- constructed source object export,
- admissible `S_sel_int`,
- admissible `E_orient`,
- closure.

It only removes the simplest reduction-to-packaging reading.

## Hard limits

`F79` does not discharge:

- the first admissibility clause,
- a constructed source object,
- admissible `S_sel_int`,
- admissible `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191`,
- ToE closure.
