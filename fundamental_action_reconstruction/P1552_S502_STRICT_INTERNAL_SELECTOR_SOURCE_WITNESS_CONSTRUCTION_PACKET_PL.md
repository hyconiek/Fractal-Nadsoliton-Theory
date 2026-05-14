# P1552 S502 Strict Internal Selector Source Witness Construction Packet (No Legacy Bridge)

Status: `P1552_PROPOSED_STRICT_INTERNAL_SELECTOR_SOURCE_WITNESS_CONSTRUCTION_PACKET`
As of: `2026-05-14`

## Cel

Wykonać pierwszy konstrukcyjny krok bezpośrednio pod rozwiązanie `QW-2191` na trasie:

`F_Nadsoliton => L_SM + L_GR`

w trybie strict-only i z pełną fizyczną interpretowalnością.

## Decyzja profesorska

Następny uczciwy krok to eksport **kandydata witness**:

`S_sel_int_strict_source_witness_v1_candidate`

z jawnym bilansem dowodowym i testami, ale bez claimu closure.

## Zakres fizyczny (minimalny)

Kandydat musi naraz zawierać:

1. komponent źródła selektora wewnętrznego (`selector_source_component`),
2. mapę zgodności orientacji dla obu kanałów (`sm_channel`, `gr_channel`),
3. test stabilności przy perturbacji sygnału wejściowego,
4. wynik testu rozdzielczości uniqueness (czyli brak nieusuwalnej degeneracji),
5. jawny status: `QW-2191` pozostaje open, jeśli uniqueness nie jest theorem-level.

## PASS/FAIL

PASS tylko jako:

`STRICT_SOURCE_WITNESS_CANDIDATE_EXPORTED_WITH_PHYSICAL_TESTS`

FAIL jeśli:

- kandydat nie ma osi SM/GR,
- brak testu perturbacyjnego,
- pojawia się closure claim bez theorem-level uniqueness.

## Co to znaczy dla ToE

ToE nadal otwarte. `S502` ma dostarczyć brakujący, fizycznie sprawdzalny filar
pod późniejsze theorem-level domknięcie `QW-2191`.

## Omówienie dla laika

To etap budowy „prototypu klucza” do zamka `QW-2191`.
Sprawdzamy, czy klucz faktycznie pasuje do dwóch zamków naraz (SM i GR)
i czy nie psuje się przy drobnych zmianach warunków.
Dopiero wtedy można iść do finalnego dowodu.
