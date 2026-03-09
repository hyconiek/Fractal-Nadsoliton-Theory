# F164 First Actual Nonstrict Declared-Scope ToE Local Derivative Candidate Support Packet

Status: `F164_EXECUTED_FIRST_ACTUAL_NONSTRICT_DECLARED_SCOPE_TOE_LOCAL_DERIVATIVE_CANDIDATE_SUPPORT_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N272`, the repo already exports:

1. one actual non-strict declared-scope ToE preclosure support packet,
2. one explicit future-only non-strict declared-scope ToE target,
3. one local source-side derivative calculation packet from `F163`,
4. one theorem-level boundary theorem from `N273` saying that this local
   derivative is still not actual ToE closure.

So the strongest honest next move is still not:

1. actual non-strict declared-scope ToE closure,
2. actual strict-core ToE closure,
3. actual global ToE closure.

It is narrower:

```text
freeze one actual candidate-support packet
integrating the local derivative datum
into the already existing non-strict declared-scope ToE lane
```

while keeping the result explicitly below actual closure.

## Fixed support packet

Reuse:

1. `N271`
   one actual non-strict declared-scope ToE preclosure support packet,
2. `N272`
   one future-only non-strict declared-scope ToE closure target,
3. `N273`
   one theorem-level boundary theorem saying that the new derivative datum is
   only a candidate ingredient and not an actual closure discharge,
4. `F163`
   one local source-side derivative calculation witness.

Freeze one actual candidate-support packet:

```text
Sigma_nonstrict_declared_scope_toe_local_derivative_candidate_support_v1 :=
(
  Lambda_nonstrict_declared_scope_toe_preclosure_support_v1,
  C_toe_declared_scope_nonstrict_future_target_v1,
  chi_src_local_derivative_calculation_witness_v1,
  local_derivative_candidate_only_boundary
)
```

## Result

`F164` exports one actual candidate-support packet:

```text
Sigma_nonstrict_declared_scope_toe_local_derivative_candidate_support_v1
```

meaning only:

1. the non-strict declared-scope ToE lane now contains one explicit local
   derivative candidate datum in addition to the earlier preclosure support
   packet,
2. the candidate datum is now packaged at theorem-route level rather than left
   as an isolated calculation,
3. the result still remains below any actual non-strict declared-scope ToE
   closure theorem.

## Hard limits

`F164` does not discharge:

1. actual non-strict declared-scope ToE closure,
2. actual strict-core selector closure,
3. actual global selector closure,
4. actual global `QW-2191` discharge,
5. actual strict-core ToE closure,
6. actual global ToE closure.

## Recommended next move

The correct next move is:

1. test whether the current repo really exports this new candidate-support
   packet,
2. keep the derivative datum explicitly local and source-side only,
3. do not relabel the packet as already discharged ToE closure.
