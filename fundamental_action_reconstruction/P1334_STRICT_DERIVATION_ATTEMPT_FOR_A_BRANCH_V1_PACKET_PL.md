# P1334 — Strict derivation attempt for `A_branch_v1` (PL)

Status: `STRICT_DERIVATION_NOT_ESTABLISHED`
As of: `2026-05-12`
Depends on: `P1333`

## Cel
Wykonać kolejny krok po `P1333`: spróbować podnieść aksjomat
`A_branch_v1` z poziomu non-strict do strict poprzez sprawdzenie, czy istnieje
już wewnętrzne źródło strict dla tej reguły.

## Artefakt wykonawczy
- skrypt: `p1334_strict_derivation_attempt_for_a_branch_v1.py`
- raport: `generated/p1334_strict_derivation_attempt_for_a_branch_v1_report_v1.json`

## Wynik
- `v4_replay_adversarial_ready = true`,
- `l1_exported = true`,
- `axiom_non_strict_closure_exists = true`,
- `strict_internal_source_for_A_branch_v1_exported = false`,
- `strict_derivation_succeeded = false`.

## Decyzja profesorska
To jest klarowny wynik: mamy skuteczną wersję non-strict i mocne wsparcie
numeryczne, ale brak strict-side derivation aksjomatu `A_branch_v1`.

Dlatego status pozostaje:
- `QW-2191 strict = NOT_CLOSED`,
- `QW-2191 non-strict = CLOSED_NON_STRICT`.

## Konsekwencja
Nie wolno ogłaszać strict closure. Program wymaga eksportu wewnętrznego źródła
strict dla reguły gałęzi near-boundary.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że non-strict closure jest bezwartościowe.

## Rekomendowany następny uczciwy krok
Uruchomić **P1335 internal-source construction for A_branch_v1**:
1. jawnie skonstruować strict-side mechanizm wyboru gałęzi,
2. pokazać jego kompatybilność z `R_common_v1` i `v4`,
3. zaktualizować checker O3 po tym eksporcie.

## Dla laika
Mamy działającą regułę graniczną, ale wciąż jako "dodatkowe założenie".
Następny krok to pokazać, że ta reguła wynika naturalnie z samej teorii,
a nie jest dopisana z zewnątrz.
