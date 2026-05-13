# P1509 — S4.59 Strict Contradiction-Branch Stress For F⇒LSM+LGR Packet (PL)

Status: `P1509_EXECUTED_STRICT_CONTRADICTION_BRANCH_STRESS`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1508`: uruchomić gałąź obaleniową
(contradiction branch) dla draftu twierdzenia `F (Nadsoliton) ⇒ L_SM + L_GR`
na dopuszczonych perturbacjach strict-only.

## Decyzja profesorska

Nie podnosimy claimu. Zamiast tego testujemy, czy draft twierdzenia przeżywa
kontrolowane odchylenia parametrów bez łamania guardrail.

## Protokół stresowy (strict-only)

Testujemy trzy rodziny perturbacji:

1. `selector_scale` (zmiana siły `S_internal_v1`),
2. `weight_shift` (przesunięcie `w_SM`, `w_GR` przy zachowaniu sumy bliskiej 1),
3. `orientation_flip_probe` (próba odwrócenia wspólnej orientacji selektora).

Każda próba sprawdza A1..A5 z `P1508`.

## Kryterium wyniku

- jeśli pojawi się kontrprzykład A1..A5 w dopuszczonej rodzinie, draft odpada,
- jeśli nie pojawi się kontrprzykład, draft utrzymuje status roboczy,
- w obu przypadkach: `qw2191_closed = false`.

## Wynik P1509

Publikujemy raport stresu gałęzi obaleniowej jako uczciwy krok strict-only
na torze `F⇒LSM+LGR`, bez bridge do legacy.
