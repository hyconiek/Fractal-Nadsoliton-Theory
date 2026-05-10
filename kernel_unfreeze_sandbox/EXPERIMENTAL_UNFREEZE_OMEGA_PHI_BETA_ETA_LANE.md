# Experimental Unfreeze Omega Phi Beta Eta Lane

Status: `NONSTRICT_FIT_DERIVED_PARAMETER_UNFREEZE_SANDBOX`
Branch: `experimental/unfreeze-omega-phi-beta-eta`

## Goal

Open one honest exploratory lane in which the current strict working tuple

```text
(omega, phi, beta, eta) = (0.18575, 0.16250, 1.0, 1.8)
```

is treated as unfrozen for testing purposes.

The lane exists because the repo itself shows the following chain:

```text
QW-2038 refreeze scan
-> QW-2039 refrozen gate
-> QW-2041 semantic drift audit
-> QW-2049 working gate selection
```

So the present tuple is not treated here as a timeless first-principles
necessity, but as a later pipeline working choice with scan/gate provenance.

## Why This Is Allowed

The lane is allowed only because it stays explicitly non-strict.

It does not:

1. identify `K_strict_gate` with `K_legacy_ont`,
2. import legacy electroweak or gravity claims onto the strict kernel,
3. claim strict selector closure,
4. claim that parameter motion itself solves `QW-2191`,
5. replace the active strict frontier with a fit-only story.

## Source Evidence

### QW-2038

- explicit `DERIVATION_COMPATIBLE_KERNEL_REFREEZE_PASS`,
- best row:

```text
omega/phi/beta/eta = 0.185750 / 0.162500 / 1.000000 / 1.800000
```

### QW-2039

- explicit `DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_PASS`,
- selected refrozen tuple:

```text
omega/phi/beta/eta = 0.185750 / 0.162500 / 0.920000 / 1.800000
```

### QW-2041

- explicit `CANONICAL_REFROZEN_REPARAMETERIZATION_FAIL`,
- canonical/refrozen semantic drift is already documented internally.

### QW-2049

- explicit `SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS`,
- selected working tuple:

```text
omega/phi/beta/eta = 0.185750 / 0.162500 / 1.000000 / 1.800000
```

## Sandbox Objectives

This lane is allowed to test only the following questions:

1. whether local motion of `(omega, phi, beta, eta)` yields a cleaner
   inversion-sensitive source-side provider hypothesis than the frozen tuple,
2. whether any such motion reduces the `pair1/pair2` residual orbit-direction
   ambiguity heuristically,
3. whether the working tuple sits near a merely operational gate optimum or
   near a broader structurally stable region,
4. whether unfreezing exposes that the strict tuple is too fit-dependent to
   support stronger selector claims.

## Search Windows

### Local continuity window

```text
omega in [0.14575, 0.22575]
phi   in [0.12250, 0.20250]
beta  in [0.85000, 1.15000]
eta   in [1.60000, 2.00000]
```

### Canonical-tension window

```text
omega in [0.18575, 0.785398]
phi   in [0.16250, 0.523599]
beta  in [0.01000, 1.00000]
eta   in [1.00000, 1.80000]
```

### Hard admissibility hints from QW-2049 micro support

```text
beta CI95 overlap: [0.000054, 106.682605]
eta  CI95 overlap: [1.000000, 3.000000]
```

These are recorded as provenance, not as recommended search ranges.

## Score Axes

Every exploratory result should be judged on four axes:

1. operational continuity to `QW-2049`,
2. semantic tension to the canonical tuple from `QW-2041`,
3. selector relevance,
4. anti-fitting discipline.

## Stop Rules

Stop the lane if any one of the following happens:

1. better observable fit appears without any new selector-relevant structure,
2. the lane starts recycling old compensation tricks,
3. parameter motion is used to imply strict closure language,
4. the lane silently reimports legacy semantics into the strict kernel.

## Honest Outcome Labels

Allowed labels:

1. `non-strict exploratory improvement`,
2. `fit-recycling warning`,
3. `selector-relevant heuristic`,
4. `operational-only retune`,
5. `bridge-pressure only`.

Forbidden labels:

1. `strict discharge`,
2. `T176 exported`,
3. `QW-2191 solved`,
4. `legacy/strict bridge proven`.
