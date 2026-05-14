# P1540 S490 TB_THM_MAIN Composition And Anti-Loop Progress Packet (No Legacy Bridge)

Status: `P1540_EXECUTED_TB_THM_MAIN_COMPOSITION_AND_ANTI_LOOP_PROGRESS_PROVISIONAL`
As of: `2026-05-14`

## Cel

Odpowiedzieć formalnie na pytanie: "czy idziemy do przodu, czy kręcimy się w
kółko?" oraz wykonać następny krok po `P1539`.

## Diagnoza postępu

W strict-only pipeline postęp uznajemy za realny, jeśli jednocześnie:

1. rośnie liczba lematów podniesionych do `theorem_level_candidate`,
2. maleje liczba `critical_open_gaps`,
3. pojawia się nowy obiekt kompozycji (`TB_THM_MAIN` composition record),
4. brak restartu do legacy bridge.

## Zakres S490

`S490`:

1. tworzy rekord kompozycji `TB_THM_MAIN` z aktualnych lematów,
2. liczy metrykę `progress_delta` vs poprzedni checkpoint,
3. klasyfikuje stan jako `forward_progress` albo `loop_risk`.

## Kontrakt wyjścia

- `tb_thm_main_composition_record`,
- `progress_metrics`,
- `progress_classification` in `{forward_progress, loop_risk}`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli `progress_classification=forward_progress` oraz brak legacy transfer.

FAIL jeśli `progress_delta <= 0` i nie ma nowego obiektu kompozycji.
