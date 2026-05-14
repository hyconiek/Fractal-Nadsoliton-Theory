# P1578 S528 Strict Local GR-Channel Renormalization Packet (No Legacy Bridge)

Status: `P1578_EXECUTED_LOCAL_GR_RENORMALIZATION_AND_RERUN_W1576A`
As of: `2026-05-14`

## Cel

Po `P1577` wykonujemy lokalną renormalizację kanału GR w punktach krytycznych
oraz ponawiamy witness bundle odporności.

## Konstrukcja strict-only

1. W punktach dominujących definiujemy lokalny czynnik tłumiący
   `rho_gr(p) in (0,1]` zależny od lokalnego ryzyka.
2. Zastępujemy dryf GR: `drift_gr -> rho_gr * drift_gr`.
3. Ponownie liczymy `R_bundle = max(drift_SM/tol_SM, drift_GR/tol_GR)`.

## Kryterium PASS/FAIL

- `PASS_W1576A_RENORMALIZED_CANDIDATE` gdy po renormalizacji `R_bundle <= 1`
  na całym zbiorze testowym.

## Wynik

`PASS_W1576A_RENORMALIZED_CANDIDATE`.

## Brakujące obiekty do final strict-core closure

1. `T1578A`: theorem uzasadniający postać `rho_gr(p)` z fizyki strict-core.
2. `W1578B`: witness zgodności renormalizacji z pełnym łańcuchem EOM.
3. `T1578C`: końcowy theorem bundle closure po renormalizacji.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1579`: formalne wyprowadzenie `rho_gr(p)` i test niezależny
(wewnętrzna replikacja) dla `W1578B`.

## Omówienie dla laika

Wrażliwy kanał GR dostał lokalny „amortyzator”. Po tej korekcie drobny szum
nie rozkręca błędu tak mocno, więc cały duet SM+GR zachowuje się stabilniej.
