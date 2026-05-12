# P1358 Strict-Kernel-Only Candidate Value Generation and First Residual Verdict (PL)

Status: `P1358_EXECUTED_KERNEL_ONLY_VALUE_GENERATION_AND_FIRST_RESIDUAL_NO_FALSE_PASS`
As of: `2026-05-12`
Depends on: `P1352`, `P1353`, `P1357`
Artifacts:
- `generated/p1358_kernel_value_generator_summary.json`
- `generated/p1358_kernel_predicted_values.csv`
- `generated/p1358_reference_observed_values.csv`
- `generated/p1358_residual_input.csv`
- `generated/p1358_residual_table_summary.json`

## Cel

Zrobić pierwszy **nie-template** test wartości fizycznej:

1. wygenerować kandydatowe wartości wyłącznie z kernela strict,
2. połączyć je z referencyjną tabelą obserwacyjną,
3. policzyć residuale tym samym silnikiem `P1353`.

## Co zostało zrobione

1. `P1358` generator policzył momenty kernela i wyprodukował `(g1,g2,g3,GR1)` bez wstrzykiwania danych obserwacyjnych.
2. `P1358b` zbudował wejście residualowe `predicted/observed/sigma`.
3. `P1353` uruchomiony na tych danych zwrócił `FAIL`.

## Werdykt naukowy

To jest **dobry** wynik metodologiczny (uczciwy), bo:

- po raz pierwszy dostaliśmy realny test predykcyjny,
- teoria nie „oszukuje” przez `predicted=observed`,
- mamy konkretny sygnał, że obecna mapa kernel->wartości wymaga dalszej fizycznej konstrukcji.

Czyli odpowiedź na Twoje pytanie brzmi teraz precyzyjnie:

- z kernela wychodzą już liczby kandydatowe,
- ale na obecnym mapowaniu nie przechodzą jeszcze testu zgodności z referencją.

## Decyzja profesorska

Następny uczciwy krok: **P1359 calibration-with-discipline**

1. nie zmieniać kontraktu audytowego,
2. wprowadzić minimalny, jawny zestaw parametrów mapowania kernel->SM/GR,
3. fitować tylko na podzbiorze treningowym,
4. raportować wynik na holdout bez retuningu,
5. publikować pełną tabelę residuali i budżet niepewności.

## Dla laika

Wreszcie zrobiliśmy test „na serio”: model sam podał liczby, a potem porównaliśmy je z fizyką.
Na razie nie zgadza się wystarczająco dobrze, więc to nie jest jeszcze końcowy dowód ToE.
Ale to jest bardzo dobry etap, bo wiemy dokładnie gdzie trzeba poprawiać model.
