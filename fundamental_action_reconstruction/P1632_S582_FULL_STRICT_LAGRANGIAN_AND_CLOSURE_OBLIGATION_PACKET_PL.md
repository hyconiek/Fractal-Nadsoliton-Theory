# P1632 / S582 — Full strict Lagrangian export and closure obligations map

## Cel
Dowieźć pełny (nie-szkieletowy) zapis toru strict-only:
`K_strict -> współczynniki -> L_total -> EOM`,
oraz jawnie wyprowadzić listę brakujących eksportów/witnessów/theoremów do strict-core closure.

## Wejścia
- `generated/p1631_s581_cover_wise_jacobian_invertibility_summary.json`

## Wyjście
- `generated/p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json`

## Zakres rygoru
- strict-only; zero legacy-bridge.
- pełny łańcuch fizyczny: Nadsoliton => L_SM + L_GR + sektor strict-skalar + sprzężenia mieszane.
- domknięcie ToE pozostaje `OPEN`, dopóki braki theorem-level nie są dowiedzione.

## Następny uczciwy krok
Udowodnić globalną kompatybilność overlapów na coverze (noncyclic) i z niej wyprowadzić twierdzenie o globalnej jednoznaczności selektora (QW-2191 strict-core).
