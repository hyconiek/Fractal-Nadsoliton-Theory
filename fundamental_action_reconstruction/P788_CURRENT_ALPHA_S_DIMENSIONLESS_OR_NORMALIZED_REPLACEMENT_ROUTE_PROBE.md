# P788 Current Alpha_s Dimensionless Or Normalized Replacement Route Probe

Status: `P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE`
As of: `2026-03-19`

## Goal

Test one narrow question after `F787`:

```text
does the current repo already export a real dimensionless or explicitly
normalized replacement route for the alpha_s boundary, so that mu0_gev can be
removed from the minimal strict bridge?
```

## Scope

`P788` does not derive a new alpha_s boundary.
It only checks whether the replacement route already exists in the current repo
state.

## Main Checks

1. confirm the strict mass proxy layer exists and remains dimensionless,
2. confirm the current alpha_s input lane is still GeV-level,
3. confirm the current QW-2093 alpha_s formulas still route through
   `m_bottom` and `m_top/m_bottom`,
4. confirm the current QCD running baseline still accepts GeV-scale inputs,
5. confirm no exported dimensionless/normalized alpha_s adapter object is
   presently available.

## Hard Limit

If the answer is negative, `P788` must reduce the gap to explicit missing
interfaces rather than silently promote any proxy or ansatz into a replacement
route.
