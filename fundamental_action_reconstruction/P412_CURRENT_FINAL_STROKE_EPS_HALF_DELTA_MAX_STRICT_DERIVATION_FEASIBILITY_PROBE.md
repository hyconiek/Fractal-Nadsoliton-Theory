# P412 Current “Final Stroke” (eps=1/2 + delta_d=delta_max) Strict-Derivation Feasibility Probe

Status: `P412_EXECUTED_CURRENT_FINAL_STROKE_EPS_HALF_DELTA_MAX_STRICT_DERIVATION_FEASIBILITY_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Audit one proposed “Final Stroke” strategy in **strict** mode:

1. derive `eps = 1/2` from a “charge parity balance / zero-charge point” principle, and
2. derive `delta_d = delta_max` from a “maximum information packing / void saturation” principle,
   optionally citing the strict Shannon-derived `alpha_geo = 4 ln 2` (`N420`),

and decide whether any such derivation already exists in current strict exports, without false pass.

## Strict-admissible evidence reused

1. `T117`
   - the strict sigma-int driven `E_pair` generator candidate uses free `eps ∈ [0,1]`.
2. `T119`
   - the strict positive-window corridor admits free `delta_d ∈ (0, delta_max]`, with `delta_max := d_local/11`.
3. `F317/N428`
   - exports `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2` as `strict_source_upgraded` (premise).
4. `F328/N440`
   - exports `delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max` as `strict_source_upgraded` (premise).
5. `N441` and `N437`
   - prove theta-pair dependence on admissible `eps` and `delta_d` choices (both are real selector slots).
6. `P408`
   - strict-admissibility audit: no strict theorem derives `eps = 1/2` nor uniquely selects `delta_d = delta_max` today.
7. `N446` and `N447`
   - theorem-level nonderivation packaging for the two proposed strict moves.
8. `T160/T161` and `P410/P411` and `N444/N445`
   - the strict-derived slot-selection targets are named and remain not discharged.

## Probe table

| Proposed “Final Stroke” move | Strict verdict | Evidence |
|---|---|---|
| derive `eps = 1/2` from “charge parity balance / zero-charge point” | **NO** | no typed strict charge-parity-balance law exported (`P408`); literal `Z2` balance constraints on `T117` do not derive `1/2` (they do not constrain or force `eps=0`) (`N446`); eps remains a selector slot (`N441`) |
| derive `delta_d = delta_max` from “maximum information packing / void saturation” (+ cite `alpha_geo`) | **NO** | `delta_max` is a corridor bound, not an objective optimizer (`T119`); no typed strict information-packing objective exported (`P408`); `alpha_geo` derivation does not supply `delta_d` (`N420`, packaged in `N447`); delta_d remains a selector slot (`N437`) |

## Exact verdict (no false pass)

On the current strict sigma-int lane:

```text
eps = 1/2: premise-based strict provenance export exists (F317/N428),
           but strict-derived derivation from charge parity balance is NOT exported (N446, P408, T160/P410/N444).

delta_d = delta_max: premise-based strict provenance export exists (F328/N440),
                     but strict-derived derivation from information packing is NOT exported (N447, P408, T161/P411/N445).
```

Therefore the “Final Stroke” as a **strict-derived** discharge is not admissible today and cannot be used
to clear strict non-claims (strict-core theta export, strict-core canonical `O(2)`-cut, or `QW-2191` discharge).

## Next honest step

If the goal is strict-core upgrade (`T159`) rather than premise-based representative work:

1. export genuinely `strict_derived` (not premise-only) slot-selection ingredients for `eps` and `delta_d`
   (`T160` and `T161`), **or**
2. export a new strict sigma-int → theta construction class that does not contain the `eps` / `delta_d` slots,
   **or**
3. proceed in a separated non-strict scope by declaring explicit selector premises (no promotion to strict-derived).

