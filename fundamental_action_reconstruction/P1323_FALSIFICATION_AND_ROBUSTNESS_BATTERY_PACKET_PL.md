# P1323 — Falsification and robustness battery packet for `S_sel_strict_v2` (PL)

Status: `CANDIDATE_FALSIFIED_ON_EXTENDED_WINDOW`
As of: `2026-05-12`
Depends on: `P1322`

## Cel
Wykonać rekomendowany krok z `P1322`: rozszerzony test falsyfikacyjny
`S_sel_strict_v2` poza lokalny zestaw C1–C4.

## Artefakt wykonawczy
- skrypt: `p1323_falsification_robustness_battery.py`
- raport: `generated/p1323_falsification_robustness_battery_report_v1.json`

## Wynik
- `sample_count = 200`
- `negative_count = 15`
- status: `COUNTEREXAMPLES_FOUND`
- `robust_positive_on_test_window = false`

## Decyzja profesorska
`S_sel_strict_v2` **nie przechodzi** rozszerzonej baterii robustności.
To znaczy: kandydat był lokalnie obiecujący, ale nie jest stabilnym prawem
selektora w szerszym oknie testowym.

## Konsekwencja
- `S_sel_strict_v2` nie może zostać awansowany na
  `strict_core_selector_source_exported = true`.
- `QW-2191` strict pozostaje `NOT_CLOSED`.
- Trzeba wrócić do etapu konstrukcji nowego kandydata (`v3`) lub
  przeformułować przestrzeń cech.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że żadna wersja `S_sel_strict` nie może działać.

## Rekomendowany następny uczciwy krok
Uruchomić **P1324 candidate-v3 construction**:
1. dodać nowy składnik stabilizujący (np. termin nieliniowy lub warunek
   monotoniczny),
2. jawnie opisać klasę dopuszczalności,
3. powtórzyć pełną baterię P1323 jako test wejścia.

## Dla laika
Nowy kompas działał na małej próbce, ale gdy sprawdziliśmy go szerzej,
pojawiły się przypadki, gdzie wskazuje zły kierunek. To normalny etap nauki:
teraz budujemy lepszą wersję kompasu i testujemy od nowa.
