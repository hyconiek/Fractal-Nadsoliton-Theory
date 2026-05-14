# P1572 S522 Strict Semiglobal Chart Patching And Regional Bound Packet (No Legacy Bridge)

Status: `P1572_EXECUTED_CHART_PATCHING_AND_REGIONAL_ERROR_BOUND_CANDIDATE`
As of: `2026-05-14`

## Cel

Po `P1571` domykamy następny uczciwy krok:

1. podział domeny na karty (chart patching),
2. odseparowanie obszaru złego uwarunkowania,
3. eksport regionalnych boundów błędu inwersji.

## Konstrukcja strict-only

Dla mapy 4-obserwowalnej wyznaczamy na siatce:
- `good_chart`: punkty z `||J^{-1}||_inf <= cond_cut`,
- `buffer_chart`: punkty przejściowe,
- `bad_chart`: punkty przekraczające próg.

Następnie raportujemy regionalne boundy `||delta p|| <= B_region * ||delta F||`.

## Kryterium PASS/FAIL

- `PASS_T1571A_CANDIDATE`:
  1. istnieją niepuste `good_chart` i `buffer_chart`,
  2. `bad_chart` jest jawnie odseparowany,
  3. regionalne boundy są skończone i monotoniczne (`B_good < B_buffer < B_bad`).

## Wynik

`PASS_T1571A_CANDIDATE`.

## Brakujące obiekty do final strict-core closure

1. `T1572B`: formalny dowód ciągłości map przejścia między kartami.
2. `W1572C`: witness zgodności patchingu z pełnym bundlem `L_SM + L_GR`.
3. `T1572D`: globalny theorem błędu dla klasy eksperymentów EOM.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1573`: formalizacja map przejścia chart->chart i dowód ciągłości
(`T1572B`) plus integracja z pełnym torem Lagrangian/EOM.

## Omówienie dla laika

To jak mapa miasta podzielona na strefy: w jednej jedzie się pewnie, w innej
trzeba większej ostrożności. Dzięki temu nie udajemy, że wszędzie jest tak samo
łatwo, tylko uczciwie podajemy gdzie model jest najbardziej wrażliwy.
