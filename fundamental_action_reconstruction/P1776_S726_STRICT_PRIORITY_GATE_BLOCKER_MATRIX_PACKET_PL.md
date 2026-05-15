# P1776 — S726
## STRICT PRIORITY GATE BLOCKER MATRIX PACKET (PL)

## Cel

Zsynchronizować aktualny status krytycznych bramek priorytetu (E_Aμ, E_H, EL_g,
H1 4D, Bianchi/Ward, BRST/Cutkosky) w jednym strict-only artefakcie bez fałszywego PASS.

## Wynik

- `E_A^μ` i `E_H`: jawne nonproxy eksporty kowariantne istnieją, ale nadal OPEN na poziomie komponentowym.
- `EL_g^{μν}`: jawny nonproxy eksport kowariantny istnieje, ale residual `EL_g-E_{μν}` nadal OPEN komponentowo.
- `boundary-term control`: kontrakt istnieje i jest reużyty.
- `Bianchi/Ward`: blokada praktyczna przez niedomknięty odbiór `W1` (`OBSTRUCTION_W1_NOT_YET_FULL_EXPORT`).
- `BRST/Cutkosky`: pozostają zablokowane przez warunki wejścia z BW i brak theorem-level witness.

## Co jest dowiedzione

1. Priorytetowe eksporty operatorowe (E_Aμ, E_H, EL_g) są obecne jako nonproxy na poziomie kowariantnym.
2. W1 nadal nie spełnia warunków `FULL_EXPORT`, więc GBW rerun i dalsze bramki nie są dopuszczalne.

## Co pozostaje OPEN

1. komponentowe świadectwo H1,
2. komponentowy residual `EL_g-E_{μν}` na bazie `B1/B2/B3/C1/C2`,
3. theorem-level Bianchi/Ward,
4. theorem-level BRST nilpotency i Cutkosky unitarity.

## Ryzyko false-pass

- Najwyższe ryzyko: traktowanie „operator exported” jako równoważne theorem-level closure.
- Drugie ryzyko: przeskoczenie blokady W1 i uruchamianie BRST/Cutkosky „na skróty”.

## Następny uczciwy krok

Dowieźć W1 do `FULL_EXPORT` (brakujące współczynniki B2/B3 + normalizacja dywergencji H(R2)),
a potem uruchomić komponentowe testy H1 i `EL_g-E_{μν}` wraz z raportem Bianchi/Ward.

## Wyjaśnienie dla laika

To etap porządkowania sygnalizacji: wiemy, które elementy już działają, a które są tylko częściowo gotowe.
Dopiero po domknięciu jednego brakującego modułu (W1) można bezpiecznie wejść w końcowe testy spójności.
