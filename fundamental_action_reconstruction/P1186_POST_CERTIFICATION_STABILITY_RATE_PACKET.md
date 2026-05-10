# P1186 Post Certification Stability Rate Packet

Status: `P1186_EXECUTED_POST_CERTIFICATION_STABILITY_RATE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Ocenić stabilność po certyfikacji nie punktowo, ale jako częstość sukcesu na
powtarzanych rerunach drift-checka.

## Professor-level decision

Dodaję `p1186_post_certification_stability_rate.py`, który powtarza `P1185`
(w tej wersji: 5 razy) i raportuje `stability_rate`.

## Honest boundary

`P1186` to statystyka procesu; brak claimu closure i brak `QW-2191` discharge.
