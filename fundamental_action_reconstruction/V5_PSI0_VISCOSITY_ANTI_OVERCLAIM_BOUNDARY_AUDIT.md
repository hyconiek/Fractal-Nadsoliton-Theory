# V5 psi0+viscosity anti-overclaim boundary audit

Data: 2026-03-06
Status: `PASS_PARTIAL_BOUNDARY_CERTIFIED_NO_FALSE_PROMOTION`

## Cel

Jawnie odgraniczyc lane `psi0 + anizotropowa viscosity` od:
- strict core,
- niezaleznego generatora selektora,
- theorem-level/full-closure claims.

## Wejscia

- `H30`: `psi0` jako deterministic kernel-invariant anchor candidate.
- `H31`: formalny embedding `psi0 -> pair1`.
- `V3`: pair-level informational viscosity operator.
- `V4`: sprzezenie `psi0 + anizotropowa viscosity` daje pair-level efekt anchor-amplifying.

## Wynik

Lane `psi0 + viscosity` jest metodologicznie dopuszczalny tylko pod nastepujacymi warunkami:

1. `psi0` pozostaje importowanym anchorem.
2. anizotropowa viscosity dziala jedynie jako `anchor-amplifying / anchor-refining`.
3. lane nie moze byc reinterpretowany jako niezalezne zrodlo `theta_i`.
4. lane nie moze byc reinterpretowany jako strict-core closure.
5. lane nie moze byc promowany do `theorem-level PASS` ani `full-closure PASS`.

## Frontier

- `V5_B1 := psi0_coupled_viscosity_is_now_boundary_certified_as_a_secondary_anchor_amplifying_lane_only_and_cannot_be_promoted_to_a_primary_selector_source_or_strict_core_closure`

## Forbidden claims

- `psi0_plus_viscosity_generates_selector_anchor_from_scratch`
- `psi0_plus_viscosity_discharges_QW_2191`
- `psi0_plus_viscosity_is_part_of_base_strict_core`
- `psi0_plus_viscosity_achieves_theorem_level_or_full_closure`

## Najuczciwszy wniosek

Lane `psi0 + viscosity` pozostaje pomocniczym mechanizmem wzmacniajacym juz dostarczony anchor i nie stanowi samodzielnego domkniecia selektora.
