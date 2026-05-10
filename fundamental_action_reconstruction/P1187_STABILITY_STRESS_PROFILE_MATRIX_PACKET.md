# P1187 Stability Stress Profile Matrix Packet

Status: `P1187_EXECUTED_STABILITY_STRESS_PROFILE_MATRIX_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Zbudować macierz stabilności `stability_rate(threshold, candidate)` dla
kandydatów baseline i uplift, aby sprawdzić odporność procesu poza jednym
ustawieniem.

## Professor-level decision

Dodaję `p1187_stability_stress_profile_matrix.py`, który dla każdego punktu
(candidate, threshold) wykonuje kilka powtórzeń i raportuje `stability_rate`.

## Honest boundary

`P1187` to profil stresowy procesu; brak claimu closure i brak `QW-2191`
discharge.
