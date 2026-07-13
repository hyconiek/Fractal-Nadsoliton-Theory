# P3138/S2088 Non-Fourier D_HL extremum joint-source audit

Status: `P3138_NONFOURIER_E_DHL_EXTREMUM_JOINT_SOURCE_BOUNDED_NO_GO`

## Constructed object
- `E_DHL zero-crossing/curvature-extremum joint source candidate`
- Formula: `E_DHL(D) = (argmax_x Delta D(x), zero-crossing set, argmax_x |Delta^2 D(x)|) with attempted polarity from slope sign`
- Non-Fourier discipline: No Fourier coefficients, characters, projectors, or phase gauges are used; only local cyclic samples and finite differences are audited.

## Finite certificate
- profiles tested: `120`
- translation-covariant receiver rows under t=1: `120`
- inversion-paired rows: `36`
- tie-free positive-slope rows: `0`
- multi-tie positive-slope rows: `120`
- accepted import-free joint sources: `0`

## Gate table
- `non_fourier_formula`: `True` — uses local values, forward differences, second differences, and zero crossings only
- `computes_nonempty_origin_receiver`: `True` — every tested profile has local slope/zero/extremum receiver data
- `tie_free_unique_cell`: `False` — fails on nonprimitive/higher-mode profiles where several equal extrema occur
- `translation_quotient_invariant`: `False` — the receiver set shifts with diagonal Z12 translation instead of descending to an absolute support representative
- `inversion_polarity_unpaired`: `False` — inversion exchanges positive and negative slope/curvature polarity classes
- `import_free_cell_order`: `False` — single-cell extraction requires a labelled cyclic order/tie convention

## Decision
P3138 constructs one non-Fourier joint source candidate E_DHL from local zero-crossings, slopes, and curvature extrema of the P3134 D_HL profiles. It is a genuine receiver and can locate local structure, but it remains translation-covariant rather than translation-quotient invariant; higher modes have equal-extremum ties; and inversion pairs the polarity classes. Thus extracting one absolute (r,lambda) still imports a labelled support order/orientation convention. No import-free joint source is exported.

## Recommendation
After P3133-P3138, the D_HL lane has now tested legacy/torsion shape, joint origin-polarity matrices, Fourier frame sourcing, and one non-Fourier extremum source. Unless a genuinely new strict source law supplies absolute support origin and polarity without imported labels, the next proof-grade move should be a no-new-live-frontier reconciliation for the D_HL selector/symmetry-breaking lane.
