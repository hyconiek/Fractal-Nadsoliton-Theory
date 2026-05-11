# P1214 Provisional Uncertified W4 Path Packet

Status: `P1214_EXECUTED_PROVISIONAL_UNCERTIFIED_PATH_WITH_HARD_NONCLOSURE`
As of: `2026-05-11`

## Goal

Na żądanie operacyjne dopuszczamy artefakt do ścieżki W4 bez pełnej certyfikacji
`P1213`, ale wyłącznie jako ścieżkę prowizoryczną/eksploracyjną.

## Professor-level decision

Modyfikujemy `P1209` o jawny tryb:

- `--allow-uncertified-artifact`.

W tym trybie gate może się otworzyć przy `base_ready=true` nawet jeśli
`certificate_ready=false`, ale status musi być wyraźnie oznaczony jako:

- `gate_mode=PROVISIONAL_UNCERTIFIED`,
- `strict_closure_claim_allowed=false`,
- `theory_closure_status=OPEN`.

## Honest boundary

To NIE jest domknięcie teorii. To jedynie kontrolowane dopuszczenie wejścia do
próby W4 przy zachowaniu twardego zakazu false-pass i closure claim.
