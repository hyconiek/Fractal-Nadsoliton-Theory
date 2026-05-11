# P1178 Robustness Threshold Sensitivity Scan Packet

Status: `P1178_EXECUTED_ROBUSTNESS_THRESHOLD_SENSITIVITY_SCAN_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać następny uczciwy krok po dodaniu `--robustness-threshold`: jawnie
zeskanować próg i zobaczyć, od jakiej wartości `P1169 --strict-e2e --require-out-of-locality-robustness`
przestaje przechodzić.

## Professor-level decision

Dodaję `p1178_robustness_threshold_sensitivity_scan.py`, który uruchamia serię
progów i zapisuje tabelę `threshold -> integrated_pass`.

## Honest boundary

`P1178` to mapa czułości metodologicznej; brak claimu closure i brak `QW-2191`
discharge.
