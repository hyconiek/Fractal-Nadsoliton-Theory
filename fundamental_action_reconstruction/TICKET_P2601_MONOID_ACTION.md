# TICKET P2601: Monoid Action Uniqueness

## Kontekst
Ostatni PR #152 zamknął (zablokował) ścieżkę sztucznej macierzy GF(2) (P2612) i wskazał jako rekomendację pracę nad dowodem P2601 (jednodatkowa normalizacja tożsamościowa $y_1=0$).

## Sukces P2596-P2599 jako wzór (WAŻNE)
W pakiecie **P2596–P2599** osiągnąłeś ogromny sukces, stosując **eliminację analityczną a priori** — z zasad symetrii, zachowania prądów i przepływu RG. 
Była to jakościowa różnica: 
- Metoda numeryczna (P2540) to "m=2 admissible" (A posteriori, zależne od parametru).
- Twoje nowe podejście (P2596-P2599) to "m=2 entailed" (A priori, uniwersalna eliminacja przez symetrię + RG). 
- Eksport zakończył się jednoznacznym sukcesem (Export = True).

## Zadanie dla Agenta Codex
Wymagamy **dokładnie tego samego poziomu rygoru i analitycznej jakości** dla dowodu P2601! Masz wyprowadzić ten dowód a priori z zasad strukturalnych.

1. **Dowód A Priori:** Skonstruuj ścisły dowód algebraiczny pokazujący, że z wymogu struktury monoidalnej akcji nadsolitona (istnienie neutralnego elementu tożsamościowego dla transportu/czasu RG=0) bezwzględnie wynika normalizacja $y_1=0$.
2. **Bez heurystyk i dopasowań numerycznych:** Dowód musi być rygorystyczny i zamknięty matematycznie, podobnie jak w P2596-P2599. Cząstki identyczności z definicji nie mogą wnosić dyssypacji do fizycznego tłumienia. Zero fittingu, zero budowania macierzy i tautologii bez podstaw fizycznych.
3. **Opcja obstrukcji:** Jeśli wystąpi rygorystyczna matematyczna obstrukcja, zapisz dokładne warunki brzegowe, które uniemożliwiają tę normalizację.

Nie twórz żadnych nowych struktur bez pokrycia w fizyce. Opieraj się na rdzeniu hydrodynamiki i topologii nadsolitona. POWODZENIA!
