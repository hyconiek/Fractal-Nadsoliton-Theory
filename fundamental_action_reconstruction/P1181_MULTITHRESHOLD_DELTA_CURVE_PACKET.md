# P1181 Multithreshold Delta Curve Packet

Status: `P1181_EXECUTED_MULTITHRESHOLD_DELTA_CURVE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wyznaczyć zakres progów robustness, gdzie kandydat uplift ma przewagę nad
baseline pod identycznym strict-E2E gate.

## Professor-level decision

Uruchamiam krzywą progową (`p1181_multithreshold_delta_curve.py`) i oznaczam
progi, gdzie `uplift_pass=true` oraz `baseline_pass=false`.

## Honest boundary

`P1181` pozostaje analizą porównawczą; brak claimu closure i brak `QW-2191`
discharge.
