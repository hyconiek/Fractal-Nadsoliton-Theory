# P462 Current Strict Shannon Element‑Order Reference `Z_24` Mode‑Index Assignment Candidate Probe

Status: `P462_EXECUTED_Z24_MODE_INDEX_ASSIGNMENT_CANDIDATE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`P461` performed a probe-only scan of the element-order reference defect

```text
F_{2m}(ord) := Σ_{x=0}^{n-1} ord_{Z_n}(x) * exp(i * 4π m x / n),
```

and found `F_{2m}(ord_{Z_n}) ≠ 0` for all Fourier-degenerate pairs on each scanned
`n ∈ {6,8,10,12,14,16,18,20,24}`.

This probe takes one cautious next step:

- pick one carrier with `n ≠ 12` (here: `n=24`),
- build the real Fourier basis on `Z_24`,
- export a **probe-level** mode-index assignment basis candidate induced by the same diagonal-defect angle rule
  `θ_* = (1/2) atan2(Im F_{2m}, Re F_{2m})`,
- without promoting any `n≠12` claim to theorem level.

## Inputs

- Pure group-structure computation on `Z_24` (no physics role-transfer claim).
- The reference profile is `ord_{Z_24}(x)` with `ord_{Z_24}(0)=1` and `ord_{Z_24}(x)=24/gcd(x,24)` for `x≠0`.

## Output artifacts

- `fundamental_action_reconstruction/generated/p462_current_strict_shannon_element_order_reference_z24_mode_index_assignment_candidate_probe.json`
- `fundamental_action_reconstruction/generated/p462_current_strict_shannon_element_order_reference_z24_mode_index_assignment_candidate_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. any theorem-level scope-extension of the `n=12` uniqueness/assignment packages to `n=24` (note: the general direction-freeness lemma
   `ord_{Z_n}` is `Aut(Z_n)`‑invariant is recorded as `N503`, but this does not by itself promote `n=24` into strict-core closure),
2. any strict-core promotion of `n=24` into the `QW-2190` physical mode scaffold,
3. any discharge of `QW-2191` beyond the declared `n=12` lanes,
4. any strict-core selector closure / admissible `S_sel_int`,
5. any ToE closure.
