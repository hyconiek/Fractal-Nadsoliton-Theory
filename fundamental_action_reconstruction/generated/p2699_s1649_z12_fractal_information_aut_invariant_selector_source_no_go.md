# P2699/S1649 Z12 fractal-information Aut-invariant selector-source no-go

Status: `P2699_Z12_FRACTAL_INFORMATION_AUT_INVARIANT_SELECTOR_SOURCE_NO_GO`

## Finite Aut(Z12) calculation
- units: `[1, 5, 7, 11]`
- orbits: `[[0], [1, 5, 7, 11], [2, 10], [3, 9], [4, 8], [6]]`
- fixed points / invariant translation strides: `[0, 6]`
- directed generator orbit of +1: `[1, 5, 7, 11]`

## No-go matrix
- `unique_directed_generator_from_aut_invariant_fractal_information`: passes=`False`; +1 has a four-element Aut orbit, so invariant information cannot select it as a unique directed generator.
- `distinguish_plus_one_from_minus_one_without_premise`: passes=`False`; +1 and -1 lie in the same Aut orbit, so a scalar Aut-invariant fractal-information source cannot provide directed sign choice.
- `derive_existing_premise_based_direction_as_nonpremise_source`: passes=`False`; The existing directed lane is real but premise-based; P2699 does not add a new non-premise source theorem.
- `strict_core_selector_or_qw2191_closure`: passes=`False`; The finite Aut calculation supplies no missing strict-core selector source and does not change P739/P740 nonexport boundaries.

## Decision
A pure internal fractal-information candidate constrained to Aut(Z12)-invariant data cannot choose a unique directed generator or distinguish +1 from -1; the real directed lane remains premise-based rather than a new strict source theorem.

## Next honest step
No proof-grade closure move is unlocked by P2699.  The next admissible move must introduce a genuinely new non-premise strict selector/orientation source beyond Aut-invariant Z12/fractal-information data, or pivot to a new typed object outside the closed selector/direct/bridge/Lagrangian lanes.  Otherwise preserve the P2697-P2699 no-new-live-frontier certificate.
