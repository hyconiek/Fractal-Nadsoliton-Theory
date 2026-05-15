# P1696 S646 Strict Gauge+Higgs Covariant EOM Export Packet (PL)

Status: `P1696_EXECUTED_STRICT_GAUGE_HIGGS_COVARIANT_EOM_EXPORT_NO_FALSE_PASS`  
As of: `2026-05-14`

## Cel

Po `P1695` wykonać pierwszy realny eksport kowariantnego bloku EOM
(zamiast samego planu), zachowując strict-only tor:

`kernel strict -> współczynniki -> pełny lagranżian -> EOM`.

## Co wyeksportowano

1. Blok `gauge + Higgs` wyprowadzony symbolicznie z kotwicy współczynników strict.
2. Jawny blok lagranżjanu (`L_gauge_higgs_reduced`) z markerem krzywizny tła `R0`.
3. Równania EOM dla pola gauge `A` oraz pola Higgs `h`.
4. Utrzymany globalny status:

`KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Rygor

- strict-only discipline utrzymane,
- brak legacy bridge,
- brak fałszywego final-pass: to nadal częściowy blok, nie pełny theorem-level closure.

## Dla laika

To pierwszy „prawdziwy fragment” pełnych równań ruchu policzony z tej wersji
teorii. Pokazujemy, że da się przejść od parametrów strict do konkretnych
równań dla pól. Ale żeby zamknąć całość, trzeba jeszcze dodać pozostałe części
(np. metrykę i fermiony) oraz globalne dowody kwantowe.

## Następny uczciwy krok (rekomendacja)

Dołączyć do tego eksportu pełny blok metryczny i spinorowy oraz zapiąć wszystko
z theorem-level witnessami: counterterm-flow, BRST/Cutkosky, background-independence.
